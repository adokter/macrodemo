#' Get log-ratios
#' @param data abundance index data
#' @param cells_all cells to compute
#' @param period period to compute
#' @param n_small_min minimum number of small cells required to calculate an abundance ratio for a large cell
#' @export
get_ratios <- function(data, cells_all, period=c("spring", "fall"), n_small_min=10){
  # verify we have all periods (spring and fall) available
  assert_that(all(period %in% names(data)))
  assert_that(length(period)==2)
  assert_that(is.character(period))
  # verify that loaded data contains no unknown grid cells
  # QUESTION: in what situation could this evaluate to FALSE?
  cells_present <- sapply(period, function(x) assert_that(all(data[[x]][[1]]$cell %in% cells_all)))
  assert_that(all(cells_present), msg="loaded data contains unknown grid cells")
  
  message(paste("calculating ratios with period1 =",period[1],"and period2 =",period[2]))
  
  # Extract cell-specific abundance data
  spring_abun_summary <- get_abun_summary(data[[period[1]]], n_small_min)
  fall_abun_summary <- get_abun_summary(data[[period[2]]], n_small_min)
  cell_timeseries <- get_cell_timeseries(cells_all,
                                         spring_abun_summary,
                                         fall_abun_summary)
  # Take the log-ratios along the timeseries to get a ratio series
  # cell_ratio_series is a summary of the bootstrap uncertainty, and
  # cell_ratio_series_full gives all bootstrap replicates.
  cells <- unique(spring_abun_summary[[1]]$cell)
  cell_ratio_series <- cell_ratio_series_full <- list()
  ratio_series_length <- 2*length(spring_abun_summary)-1
  
  for (i in 1:length(cell_timeseries)) {
    if (cells_all[i] %in% cells) {
      # calculate logarithmic ratios (lrats)
      lrats <- apply(cell_timeseries[[i]][ ,2:ncol(cell_timeseries[[i]])], 2, function(x){log(stocks::ratios(x))})
      cell_ratio_series[[i]] <- list()
      
      # cell / year / reference period (named after the latest period on which the ratio is based)
      cell_ratio_series[[i]]$cell <- rep(cells_all[i],ratio_series_length)
      cell_ratio_series[[i]]$year <- as.numeric(sapply(names(spring_abun_summary), function(x) rep(x,times=2)))[-1]
      cell_ratio_series[[i]]$period <- rep(rev(period), length.out=ratio_series_length)
      
      cell_ratio_series[[i]]$median <- apply(lrats, 1, median)
      cell_ratio_series[[i]]$avg <- apply(lrats, 1, mean)
      cell_ratio_series[[i]]$sd <- apply(lrats, 1, sd)
      cell_ratio_series[[i]]$q10 <- apply(lrats, 1, function(x){quantile(x, .1, na.rm = T)})
      cell_ratio_series[[i]]$q90 <- apply(lrats, 1, function(x){quantile(x, .9, na.rm = T)})
      cell_ratio_series[[i]]$skewness <- apply(lrats, 1, moments::skewness)
      cell_ratio_series[[i]]$kurtosis <- apply(lrats, 1, moments::kurtosis)
  
      cell_ratio_series_full[[i]] <- lrats
      # set rownames to year_period
      # colnames already set in format "rep_i"
      rownames(cell_ratio_series_full[[i]]) <- paste0(cell_ratio_series[[i]]$year,"_",cell_ratio_series[[i]]$period)
    } else {
      cell_ratio_series[[i]] <- cell_ratio_series_full[[i]] <- NA
    }
  }
  
  # make output a named list:
  names(cell_ratio_series) <- names(cell_timeseries)
  names(cell_ratio_series_full) <- names(cell_timeseries)
  
  return(list(summary=cell_ratio_series, replicates=cell_ratio_series_full))
}
