#' sample a species given thresholds on effort, space and time
#' @param species_code six letter species code
#' @param erd_path path to erd
#' @param checklists filtered checklists
#' @param bcrs sf object holding bcrs
#' @param effort_thresholds effort thresholds
#' @param extent_space spatial extent
#' @param extent_time temporal extent
#' @param time_window stratify by time grid or not
#' @param small_grid resolution of small grid
#' @param large_grid resolution of large grid
#' @param time_grid resolution of time
#' @param .cores cores for parallel computation
#' @param quiet if TRUE, suppress informational messages (warnings/errors still returned)
#' @return large list object
#' @export
sample_bcr_abun <- function(
    species_code, erd_path, checklists,
    effort_thresholds, extent_space, extent_time,
    time_window="full", small_grid=11, bcrs,
    time_grid=7, .cores=4, quiet = TRUE){
  # verify input arguments
  assertthat::assert_that(is.character(species_code))
  assertthat::assert_that(file.exists(erd_path))
  assertthat::assert_that(is.data.frame(checklists))
  assertthat::assert_that(is.data.frame(effort_thresholds))
  assertthat::assert_that(all(c("dist_max","time_min","time_max","cci_min") %in% colnames(effort_thresholds)), msg="missing threshold value(s) for dist_max, time_min, time_max, cci_min")
  assertthat::assert_that(effort_thresholds$time_min < effort_thresholds$time_max)
  assertthat::assert_that(effort_thresholds$dist_max > 0)
  assertthat::assert_that(is.data.frame(extent_space))
  assertthat::assert_that(all(c("min_lon","max_lon","min_lat","max_lat") %in% colnames(extent_space)), msg="missing threshold value(s) for min_lon, max_lon, min_lat, max_lat")
  assertthat::assert_that(extent_space$min_lon < extent_space$max_lon)
  assertthat::assert_that(extent_space$min_lat < extent_space$max_lat)
  assertthat::assert_that(is.data.frame(extent_time))
  assertthat::assert_that(all(extent_time$year_min <= extent_time$year_max))
  assertthat::assert_that(all(extent_time$tgrid_min < extent_time$tgrid_max))
  assertthat::assert_that(all(c("period","tgrid_min","tgrid_max","year_min", "year_max") %in% colnames(extent_time)), msg="missing threshold value(s) for period, tgrid_min, tgrid_max, year_min, year_max")
  assertthat::assert_that(time_window %in% c("gridded", "full"))

  # load species data from ERD
  # takes ~ 2-3 mins for carwre (Carolina Wren)
  sp_data <- import_from_erd(
    species_code,
    erd_path = erd_path,
    checklists = checklists
  )

  # loop over time periods
  data_grid <- data_abun <- list()
  for (i in 1:nrow(extent_time)) {
    # loop over years
    years=extent_time$year_min[i]:extent_time$year_max[i]
    # initialize year lists
    data_grid[[i]] <- data_abun[[i]] <- list()
    for (y in seq_along(years)) {
      if(!quiet){
        print(paste("grid sampling",species_code,"data for",extent_time$period[i],years[y],"..."))
      }
      data_grid[[i]][[y]] <-
        get_bcr_data(
          data = sp_data, .year = years[y],
          bcrs=bcrs,
          tgrid_min = extent_time$tgrid_min[i],
          tgrid_max = extent_time$tgrid_max[i],
          time_window = time_window, min_lat = extent_space$min_lat,
          max_lat = extent_space$max_lat, min_lon = extent_space$min_lon,
          max_lon = extent_space$max_lon,
          small_grid=small_grid,
          .cores = .cores
        )

      data_abun[[i]][[y]] <- get_abun_bcr(data_grid[[i]][[y]], bcrs2, n_rep=100)
    }
    names(data_grid[[i]]) <- years
    names(data_abun[[i]]) <- years
  }
  names(data_grid) <- extent_time$period
  names(data_abun) <- extent_time$period

  output <- list(grid=data_grid, abun=data_abun)
  output
}


