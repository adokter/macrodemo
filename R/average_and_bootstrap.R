#' Get abundance data
#' @inheritParams get_abun_tgrid_slice
#' @inheritParams get_pixels
#' @param .cores number of parallel cores to use
#' @return a dataframe
#' @export
get_abun <- function(sp_data, n_rep, roi, .cores = 4){
  dg_large <- dggridR::dgconstruct(res=attributes(sp_data)$large_grid)
  dg_small <- dggridR::dgconstruct(res=attributes(sp_data)$small_grid)
  cells <- get_pixels(dg_large = dg_large, dg_small = dg_small, roi)$cells_small
  cells <- cells[cells$lat > attributes(sp_data)$min_lat & 
                   cells$lat < attributes(sp_data)$max_lat &
                   cells$lon > attributes(sp_data)$min_lon, ]
  
  big_cells <- unique(cells$cell_large)
  
  cl <- parallel::makeCluster(.cores, "FORK")
  doParallel::registerDoParallel(cl)
  
  tgrid_list <- foreach (t = seq_along(sp_data)
  ) %dopar% {
    get_abun_tgrid_slice(sp_data, t, cells, big_cells, n_rep)
  }
  parallel::stopCluster(cl = cl)
  
  out <- do.call(rbind, tgrid_list)
  return(out)
}


#' Get abundance for one timeslice
#' @param sp_data Data for a particular species and year from get_grid_data()
#' @param t index for the timeslice
#' @param cells output of get_pixels() subset to ROI
#' @param big_cells list of the big-cell identities
#' @param n_rep number of bootstrap replicates
get_abun_tgrid_slice <- function(sp_data, t, cells, big_cells, n_rep){
  tgrid_data <- as.data.frame(do.call(rbind, lapply(sp_data[[t]], data.frame)))
  out1 <- data.frame(tgrid = names(sp_data)[t], cell = big_cells, tot_small = NA, n_small = NA, n_list = NA,  mean = NA)
  out2 <- as.data.frame(matrix(data = 0, nrow = length(big_cells), ncol = n_rep))
  names(out2) <- paste0("rep_", 1:n_rep)
  out3 <- cbind(out1, out2)
  for (i in seq_along(big_cells)) { # loop over the big cells
    cells_small <- paste0("cell_", cells$cell[cells$cell_large == big_cells[i]])
    tgi <- tgrid_data[rownames(tgrid_data) %in% cells_small,]
    out3$tot_small[i] <- nrow(tgi)
    out3$n_small[i] <- sum(!is.na(tgi$stixel_mean_small))
    out3$n_list[i] <- sum(tgi$n)
    
    if(sum(!is.na(tgi$stixel_mean_small)) == 0) {
      out3$mean[i] <- NA
      out3[i, 6+(1:n_rep)] <- NA
    } else {
      tgi2 <- tgi[!is.na(tgi$stixel_mean_small), ]
      
      
      out3$mean[i] <- mean(tgi2$stixel_mean_small)
      
      if(out3$mean[i] != 0) {
        for(j in 1:n_rep){
          out3[i, j + 6] <- weighted.mean(tgi2$stixel_mean_small,
                                          w = gtools::rdirichlet(1, rep(1, length(tgi2$stixel_mean_small))))
          
        }
      }
    }
  }
  return(out3)
}





