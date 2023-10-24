#' Aggregate eBird data over cells within BCRs
#' @param data filtered erd
#' @param bcrs sf representation of BCRs
#' @param .year single year for which index is desired
#' @param tgrid_min minimum week of focal period
#' @param tgrid_max maximum week of focal period
#' @param time_window one of "weekly" or "full"; controlls whether the bootstrap is stratified by week
#' @param min_lat minimum latitude
#' @param max_lat maximum latitude
#' @param min_lon minimum longitude
#' @param small_grid resolution of dggrid for small cells to use in bootstrapping
#' @param .cores number of cores for parallel processing; used only if time_window = "weekly"
#' @return abundance data
#' @export
get_bcr_data <- function (data, bcrs, .year, tgrid_min, tgrid_max, time_window = "full",
                          min_lat, max_lat, min_lon, max_lon, small_grid = 11,
                          .cores = 4)
{
  bcrs <- sf::st_transform(bcrs, 4326)
  cell_data <- list()
  zfd <- data |>
    dplyr::filter(year == .year) |>
    dplyr::filter(latitude > min_lat & latitude < max_lat) |>
    dplyr::filter(longitude > min_lon & longitude < max_lon) |>
    merge(get_tgrid(), by = "day_of_year", all = F) |>
    filter(tgrid >= tgrid_min & tgrid <= tgrid_max)

  if (time_window == "full") {
    zfd$tgrid <- 1
  }

  dg_small <- dggridR::dgconstruct(res = small_grid)
  zfd2 <- sf::st_as_sf(zfd, coords = c("longitude", "latitude"), crs = 4326)
  tt <- sf::st_intersects(zfd2, bcrs)
  lengths <- lapply(tt, length)
  assertthat::assert_that(max(unlist(lengths)) == 1)
  tt[lengths == 0] <- NA
  ttt <- unlist(tt)
  zfd$bcr <- bcrs$BCR[ttt]
  zfd$cells_small <- dggridR::dgGEO_to_SEQNUM(dg_small, zfd$longitude,
                                              zfd$latitude)$seqnum
  bcrs_used <- bcrs[unique(ttt), ]

  conus <- spData::us_states
  conus_raster <- fasterize::fasterize(conus, raster::raster(ncol = 1000,
                                                             nrow = 1000, xmn = -125, xmx = -66, ymn = 24, ymx = 50))
  conus_coords <- raster::as.data.frame(conus_raster, xy = T)
  conus_coords <- conus_coords[!is.na(conus_coords$layer),
  ]
  conus_cells_small <- data.frame(cell = unique(dggridR::dgGEO_to_SEQNUM(dg_small,
                                                                         conus_coords$x, conus_coords$y)[[1]]))
  conus_cells_small$lat <- dggridR::dgSEQNUM_to_GEO(dg_small,
                                                    conus_cells_small$cell)$lat_deg
  conus_cells_small$lon <- dggridR::dgSEQNUM_to_GEO(dg_small,
                                                    conus_cells_small$cell)$lon_deg

  ccs <- sf::st_as_sf(conus_cells_small, coords = c("lon", "lat"), crs = 4326)
  uu <- sf::st_intersects(ccs, bcrs_used)
  uu2 <- lapply(uu, length) %>% unlist()

  all_cells_small <- conus_cells_small$cell[uu2 > 0]
  if (time_window == "weekly") {
    cl <- parallel::makeCluster(.cores, "FORK")
    doParallel::registerDoParallel(cl)
    tgrid <- tgrid_min:tgrid_max
    cell_data <- foreach(t = 1:length(tgrid)) %dopar% {
      get_cell_data(zfd, tgrid[t], all_cells_small)
    }
    parallel::stopCluster(cl = cl)
  }
  else if (time_window == "full") {
    tgrid <- 1
    cell_data <- list(get_cell_data(zfd, 1, all_cells_small))
  }
  names(cell_data) <- paste0("tgrid_", tgrid)
  attr(cell_data, "species") <- attributes(data)$species
  attr(cell_data, "year") <- .year
  attr(cell_data, "small_grid") <- small_grid
  attr(cell_data, "min_lat") <- min_lat
  attr(cell_data, "max_lat") <- max_lat
  attr(cell_data, "min_lon") <- min_lon
  cell_data
}

#' Get the small pixels corresponding to each BCR, each of which represents a large pixel
#' @param bcrs sf representation of BCRs
#' @param dg_small a dggrid representing the small pixels
#' @return dataframe giving the bcr for each small cell
get_pixels_bcr <- function (bcrs, dg_small)
{
  conus <- spData::us_states
  conus_raster <- fasterize::fasterize(conus, raster::raster(ncol = 1000,
                                                             nrow = 1000, xmn = -125, xmx = -66, ymn = 24, ymx = 50))
  conus_coords <- raster::as.data.frame(conus_raster, xy = T)
  conus_coords <- conus_coords[!is.na(conus_coords$layer), ]

  conus_cells_small <- data.frame(cell = unique(dggridR::dgGEO_to_SEQNUM(dg_small,
                                                                         conus_coords$x, conus_coords$y)[[1]]))
  conus_cells_small$lat <- dggridR::dgSEQNUM_to_GEO(dg_small,
                                                    conus_cells_small$cell)$lat_deg
  conus_cells_small$lon <- dggridR::dgSEQNUM_to_GEO(dg_small,
                                                    conus_cells_small$cell)$lon_deg
  ccs2 <- sf::st_as_sf(conus_cells_small, coords = c("lon", "lat"), crs = 4326)
  tt <- sf::st_intersects(ccs2, bcrs)
  lengths <- lapply(tt, length)
  assertthat::assert_that(max(unlist(lengths)) == 1)
  tt[lengths == 0] <- NA
  ttt <- unlist(tt)
  conus_cells_small$cell_large <- bcrs$BCR[ttt]
  conus_cells_small
}

#' get bootstrapped abundance estimates for bcrs
#' @param sp_data species abundance data outputted by get_bcr_data
#' @param bcrs sf object holding the bcrs
#' @param n_rep number of bootstrap replicates
#' @param .cores number of cores for parallel processing
#' @return matrix of bootstrapped abundances
#' @export
get_abun_bcr <- function (sp_data, bcrs, n_rep, .cores = 4)
{
  dg_small <- dggridR::dgconstruct(res = attributes(sp_data)$small_grid)

  conus <- spData::us_states
  conus_raster <- fasterize::fasterize(conus, raster::raster(ncol = 1000,
                                                             nrow = 1000, xmn = -125, xmx = -66, ymn = 24, ymx = 50))
  conus_coords <- raster::as.data.frame(conus_raster, xy = T)
  conus_coords <- conus_coords[!is.na(conus_coords$layer), ]
  conus_cells_small <- data.frame(cell = unique(dggridR::dgGEO_to_SEQNUM(dg_small,
                                                                         conus_coords$x, conus_coords$y)[[1]]))
  cells <- get_pixels_bcr(bcrs = bcrs, dg_small = dg_small)
  cells <- cells[cells$lat > attributes(sp_data)$min_lat &
                   cells$lat < attributes(sp_data)$max_lat & cells$lon >
                   attributes(sp_data)$min_lon, ]
  big_cells <- unique(cells$cell_large)
  cl <- parallel::makeCluster(.cores, "FORK")
  doParallel::registerDoParallel(cl)
  tgrid_list <- foreach(t = seq_along(sp_data)) %dopar% {
    get_abun_tgrid_slice(sp_data, t, cells, big_cells, n_rep)
  }
  parallel::stopCluster(cl = cl)
  do.call(rbind, tgrid_list)
}
