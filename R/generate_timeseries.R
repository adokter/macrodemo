#' Function to extract cell-specific abundance data
#' @param abun_data abundance data
#' @param n_small_min minimum number of small cells to compute abundance index for large cell
#' @return data.frame of abundances
#' @export
#' @details
#' The output data.frame has columns `cell` for the cell index, `average` for the average abundance, and columns `rep_i`
#' containing the replicates of sampled abundance on which `average` is based.
abun_data_bycell <- function(abun_data, n_small_min = 10) {
  cells <- unique(abun_data$cell)
  # number of columns in abun_data:
  ncol <- ncol(abun_data)
  ad2 <- abun_data[abun_data$n_small >= n_small_min, ]
  abun_cols <- as.data.frame(matrix(nrow = length(cells), ncol = ncol - 5))
  names(abun_cols) <- c("average", paste0("rep_", c(1:(ncol - 6))))
  out <- cbind(data.frame(cell = cells), abun_cols)
  for(i in seq_along(cells)) {
    ad3 <- ad2[ad2$cell == cells[i], ]
    if (nrow(ad3) == 1) {
      out[i, 2:(ncol-4)] <- colMeans(ad3[6:ncol])
    }
  }
  return(out)
}

#' Abundance summary across years
#' @inheritParams abun_data_bycell
#' @return a list of data.frames for each year calculated by [abun_data_bycell()].
#' @export
#' @details Each list element for this function contains a data.frame for a year, provided by `abun_data`
get_abun_summary <- function(abun_data, n_small_min) {
  lapply(abun_data, abun_data_bycell, n_small_min = n_small_min)
}

#' Merge spring and fall indices into a timeseries
#' @param cells_all numeric vector of cell ids
#' @param spring_abun_summary spring abundance summary
#' @param fall_abun_summary fall abundance summary
#' @param quiet if TRUE suppress informational messages (warnings & errors still shown)
#' @return a list of data.frames, with list elements corresponding to cells
#' @export
#' @details This functions pastes the spring and full abundance summaries together
#' in one data.frame for each cell, with the spring on the uneven column indices and fall on the
#' even column indices, with in total 2x (number of years) columns
get_cell_timeseries <- function (cells_all,
                                 spring_abun_summary, fall_abun_summary,
                                 quiet = TRUE) {
  cells <- unique(spring_abun_summary[[1]]$cell)
  cell_timeseries <- list()
  # column indices containin sampled abundance data (rep_i) and their average (average)
  # dropping first column containing cell name
  col_abun = 2:dim(spring_abun_summary[[1]])[2]
  for (i in seq_along(cells_all)) {
    if (cells_all[i] %in% cells) {
      if(!quiet) print(paste("calculating cell",i,"..."))

      cell_data <- spring_abun_summary[[1]][spring_abun_summary[[1]]$cell == cells_all[i], col_abun]
      cell_data <- rbind(cell_data, fall_abun_summary[[1]][fall_abun_summary[[1]]$cell == cells_all[i], col_abun])
      years <- as.numeric(names(spring_abun_summary))
      for(j in 2:length(years)) {
        cell_data <- rbind(cell_data, spring_abun_summary[[j]][spring_abun_summary[[j]]$cell == cells_all[i], col_abun])
        cell_data <- rbind(cell_data, fall_abun_summary[[j]][fall_abun_summary[[j]]$cell == cells_all[i], col_abun])
      }
      cell_timeseries[[i]] <- cell_data
    } else {
      cell_timeseries[[i]] <- NA
    }

  }
  names(cell_timeseries) = cells_all
  return(cell_timeseries)
}

#' function to determine which cell-years are analyzeable
#' @param ratio_series a ratio series
#' @param uncertainty_high_grade maximum allowable uncertainty in cell-year
#' @param inf_exclude exclude all cell years when the ratio is infinite in any year
#' @param element_1_exclude exclude the first year of the series
#' @param n_min_prod minimimum number of productivity indices necessary to include
#' any years/seasons at all
#' @param n_min_surv minimum number of survicial indices necessary to include any
#' years/seasons at all
#' @param n_min_full minimum number of individual years with both a productivity and
#' a survical index to include any years/seasons at all
#' @return a logical vector of which elements to include
#' @export
use_cell_years <- function (ratio_series, uncertainty_high_grade = Inf,
                            inf_exclude = F, element_1_exclude = T,
                            n_min_prod = 5, n_min_surv = 5,
                            n_min_full = 5) {
  if(!identical(ratio_series, NA)){
    if(any(is.infinite(ratio_series$median))) {
      contains_inf <- TRUE
      # set infinite values to NA
      for(i in 1:length(ratio_series)){
        ratio_series[[i]][is.infinite(ratio_series[[i]])] <- NA
      }
    } else {contains_inf <- FALSE}
    
    insufficient_return <- rep(FALSE, length(ratio_series$median))
    if(inf_exclude & contains_inf) {
      return(insufficient_return)
    }
    
    idx_prod <- seq(from=1,to=length(ratio_series$median),by=2)
    idx_surv <- seq(from=2,to=length(ratio_series$median),by=2)
    
    if(sum(!is.na(ratio_series$median[idx_prod])) < n_min_prod) {
      return(insufficient_return)
    }
    
    if(sum(!is.na(ratio_series$median[idx_surv])) < n_min_surv) {
      return(insufficient_return)
    }
    
    if(any(!is.na(ratio_series$median[idx_prod])) &
       any(!is.na(ratio_series$median[idx_surv]))){
      # QUESTION: why are we dropping the first productivity index?
      prod_dif <- max(ratio_series$median[idx_prod[-1]], na.rm = T) -
        min(ratio_series$median[idx_prod[-1]], na.rm = T)
      surv_dif <- max(ratio_series$median[idx_surv], na.rm = T) -
        min(ratio_series$median[idx_surv], na.rm = T)
      use <- rep(1, length(ratio_series$median))
      if(element_1_exclude){
        use[1] <- 0
      }else if(is.na(ratio_series$median[1])){
        use[1] <- 0
      }else{
        uncertainty <- ratio_series$q90[1] - ratio_series$q10[1]
        if(uncertainty > (uncertainty_high_grade * prod_dif)) {
          use[1] <- 0
        }
      }
      
      for (i in 1:length(idx_surv)) {
        p_i <- 1 + 2*i
        s_i <- 2*i

        if(is.na(ratio_series$median[p_i])){
          use[p_i] <- 0
        } else {
          uncertainty <- ratio_series$q90[p_i] - ratio_series$q10[p_i]
          if(uncertainty > (uncertainty_high_grade * prod_dif)) {
            use[p_i] <- 0
          }
        }

        if(is.na(ratio_series$median[s_i])){
          use[s_i] <- 0
        } else {
          uncertainty <- ratio_series$q90[s_i] - ratio_series$q10[s_i]
          #QUESTION: is this type of thresholding still in use?
          if(uncertainty > (uncertainty_high_grade * surv_dif)) {
            use[s_i] <- 0
          }
        }
      }
      use <- as.logical(use)
    } else {
      use <- rep(FALSE, length(ratio_series$median))
    }
    
    if(sum(use[idx_prod]) < n_min_prod) {
      use <- rep(FALSE, length(ratio_series$median))
    } else if(sum(use[idx_surv]) < n_min_surv) {
      use <- rep(FALSE, length(ratio_series$median))
    } else if ((sum(use[idx_prod] * use[idx_surv])) < n_min_full) {
      use <- rep(FALSE, length(ratio_series$median))
    }

    return(use)
  }
}
