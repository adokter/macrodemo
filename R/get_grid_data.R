#' Get eBird summaries on cells of a small grid, tagged with membership in a large grid.
#' @param data zero-filled data (output from import_from_erd())
#' @param year the year
#' @param tgrid_min the minimum time grid index (week number if `time_grid=7`).
#' @param tgrid_max the maximum time grid index (week number if `time_grid=7`).
#' @param time_window "gridded" or "full"
#' @param min_lat the minimum latitude (decimal degrees) where the prediction is desired
#' @param max_lat the maximum latitude
#' @param min_lon the minimum longitude
#' @param max_lon the maximum longitude
#' @param large_grid the dggridR resolution of the large grid
#' @param small_grid the dggridR resolution of the small grid
#' @param time_grid the temporal grid size in days (default 7, one week).
#' @param .cores number of parallel cores for `get_cell_data()` calls
#' @param roi region of interest
#' @export
#' @return a list with grid-sampled data
#' @details
#' The returned list object either contains a single list of lists (when time_window equals full)
#' or a multiple lists of lists (when time_window equals weekly), one for each week. The deepest
#' list elements contains count information for the small grid, with list elements labled as
#' `cell_X` with `X` the small cell index. These small cell data elements contain the following
#' list elements returned by [get_cell_data()]:
#' \describe{
#' \item{`n`}{number of observations (number of zero-filled checklist available in the small cell)}
#' \item{`n_z`}{number of zero counts}
#' \item{`n_x`}{number of presence only observations (X's)}
#' \item{`n_value`}{number of observations with a count value (zero or a count, i.e. non-X observations)}
#' \item{`mean_positive`}{average count of observations with a count > 0}
#' \item{`stixel_mean_small`}{average count for the small cell, as calculated with [stixel_mean_small()]}
#' }
#' The average count calculated with [stixel_mean_small()] currently assumes for presence only observations
#' (X's) a count equal to `mean_positive`.
get_grid_data <- function(data, year,
                     tgrid_min, tgrid_max, time_window = "gridded",
                     min_lat, max_lat, min_lon, max_lon = Inf,
                     large_grid = 6, small_grid = 11, time_grid = 7, roi,
                     .cores = 4) {
  cell_data <- list()

  zfd <- data
  my_year=year # dummy to use in next line
  zfd <- zfd[year == my_year]
  zfd <- zfd[latitude > min_lat & latitude < max_lat & longitude > min_lon & longitude < max_lon]
  zfd <- merge(zfd, get_tgrid(time_grid), by = "day_of_year", all = FALSE)
  zfd <- zfd[tgrid >= tgrid_min & tgrid <= tgrid_max]

  if (time_window == "full") {
    zfd$tgrid <- 1
  }

  dg_large <- dggridR::dgconstruct(res=large_grid)
  dg_small <- dggridR::dgconstruct(res=small_grid)

  pixels <- get_pixels(dg_large = dg_large, dg_small = dg_small, roi = roi)
  all_cells_small <- pixels$cells_small$cell[pixels$cells_small$cell_large %in% zfd$seqnum_large]

  if (time_window == "gridded") {
    cl <- parallel::makeCluster(.cores, "FORK")
    doParallel::registerDoParallel(cl)
    tgrid <- tgrid_min:tgrid_max

    cell_data <- foreach (t = 1:length(tgrid)
    ) %dopar% {
      get_cell_data(zfd, tgrid[t], all_cells_small)
    }

    parallel::stopCluster(cl = cl)
  } else if (time_window == "full") {
    tgrid <- 1
    cell_data <- list(get_cell_data(zfd, 1, all_cells_small))
  }

  names(cell_data) <- paste0("tgrid_", tgrid)

  attr(cell_data, "species") <- attributes(data)$species
  attr(cell_data, "year") <- year
  attr(cell_data, "large_grid") <- large_grid
  attr(cell_data, "small_grid") <- small_grid
  attr(cell_data, "min_lat") <- min_lat
  attr(cell_data, "max_lat") <- max_lat
  attr(cell_data, "min_lon") <- min_lon
  attr(cell_data, "max_lon") <- max_lon

  return(cell_data)
}





#' get a 7-day grid over ordinal days
#' @param days time grid resolution in days
#' @return a two-column data-frame that tags each day-of-year to a tgrid cell
get_tgrid <- function(days=7) {
  out <- data.frame(day_of_year = 1:366)
  out$tgrid <- floor((out$day_of_year - 1)/days) + 1
  return(out)
}


#' get a hexagonal grid over the contiguous US
#' @param dg_large large dggridR grid
#' @param dg_small small dggridR grid
#' @param roi a raster object defining the Region Of Interest (ROI)
get_pixels <- function(dg_large, dg_small, roi) {
  roi_coords <- raster::as.data.frame(roi, xy = T)
  roi_coords <- roi_coords[!is.na(roi_coords$layer), ]
  roi_cells <- data.frame(cell = unique(dggridR::dgGEO_to_SEQNUM(dg_large, roi_coords$x, roi_coords$y)[[1]]))
  roi_cells$lat <- dggridR::dgSEQNUM_to_GEO(dg_large, roi_cells$cell)$lat_deg
  roi_cells$lon <- dggridR::dgSEQNUM_to_GEO(dg_large, roi_cells$cell)$lon_deg

  roi_cells_small <- data.frame(cell = unique(dggridR::dgGEO_to_SEQNUM(dg_small, roi_coords$x, roi_coords$y)[[1]]))
  roi_cells_small$lat <- dggridR::dgSEQNUM_to_GEO(dg_small, roi_cells_small$cell)$lat_deg
  roi_cells_small$lon <- dggridR::dgSEQNUM_to_GEO(dg_small, roi_cells_small$cell)$lon_deg
  roi_cells_small$cell_large <- dggridR::dgGEO_to_SEQNUM(dg_large, roi_cells_small$lon, roi_cells_small$lat)[[1]]

  out <- list(cells_large = roi_cells, cells_small = roi_cells_small)
  return(out)
}


#' Get mean value for a stixel
#' @param x data for a single cell
#' @return average abundance for the cell
stixel_mean_small <- function(x){
  if (x$n == 0) {
    return(NA)
  } else if (x$n_x + x$n_value == 0) {
    return(0)
  } else if (x$n_value == 0) {
#    return(weighted.mean(x = c(0,2), w = c(x$n_z, x$n_x)))
    return(NA)
  } else {
    # Rescale the numbers of zeros by the proportion non-X non-zero observations
    # This has the effect of assuming the mean positive value for the presence only (X) observations
    n_z_rescaled <- x$n_z * x$n_value / (x$n_value + x$n_x)
    return(weighted.mean(x = c(0, x$mean_positive), w = c(n_z_rescaled, x$n_value)))
  }
}


#' Get data for a stixel
#' @param zfd zero-filled data
#' @param t the tgrid element
#' @param all_cells_small the list of all cells
#' @return row vector
get_cell_data <- function(zfd, t, all_cells_small) {
  out <- list()
  zfd_t <- zfd[zfd$tgrid == t, ]
  for (g in seq_along(all_cells_small)) {
    cd <- list()
    obs <- zfd_t[zfd_t$seqnum_small == all_cells_small[g], ]
    cd$n <- nrow(obs)
    cd$n_z <- sum(obs$obs_count == "0" & obs$only_presence_reported == 0)
    cd$n_x <- sum(obs$only_presence_reported)
    not_x <- as.integer(obs$obs_count[obs$only_presence_reported == 0])
    cd$n_value <- sum(not_x > 0)
    cd$mean_positive <- mean(not_x[not_x > 0])
    cd$stixel_mean_small <- stixel_mean_small(cd)
    out[[g]] <- cd
    names(out)[g] <- paste0("cell_", all_cells_small[g])
  }
  return(out)
}

