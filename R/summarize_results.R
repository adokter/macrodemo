#' Summarize computed results from get_abun
#' @param spring output from get_abun for a spring
#' @param fall output from get_abun for a fall of the same year
#' @param year_remove remove all stixels that are NA in any year after
#'                      year_remove.  Remove these stixels from ALL years
#' @param remove if T, remove all years after year_remove.  Useful scatterplots
#'                      of north versus south or index versus duration, where 
#'                      the year information isn't visually apparent
#' @param min_small the minimum number of small stixels with at least one checklist
#'                      require to include a big stixel in the analysis
#' @param min_list  the minimum number of checklists required to include a big
#'                      stixel in the analysis.  Guaranteed to have no effect if 
#'                      less than or equal to min_small
#' @param ci Credible interval width for quantities reported with credible intervals.
#' @return a list of means and CIs for spring and fall indices and ratios
#' @export
summarize_avg_and_rep <- function(spring, fall, year_remove=4, remove = T, 
                                  min_small = 1, min_list = 1, ci = .8) {
  if(length(spring) != length(fall)) {stop("spring and fall have different lengths")}
  
  lq <- (1-ci)/2
  uq <- 1-lq
  
  spring_remove <- fall_remove <- vector()
  for (i in year_remove:length(spring)) {
    spring_remove <- c(spring_remove, which(spring[[i]]$n_small < min_small | spring[[i]]$n_list < min_list))
    fall_remove <- c(fall_remove, which(fall[[i]]$n_small < min_small | fall[[i]]$n_list < min_list))
  }
  spring_remove <- unique(spring_remove)
  fall_remove <- unique(fall_remove)
  
  spring_keep <- fall_keep <- 1:nrow(spring[[1]])
  if(length(spring_remove) > 0) {
    spring_keep <- spring_keep[-spring_remove]
  }
  if(length(fall_remove) > 0) {
    fall_keep <- fall_keep[-fall_remove]
  }
  
  if(remove){
    spring <- spring[(year_remove + 1):length(spring)]
    fall <- fall[(year_remove + 1):length(fall)]
  }
  
  
  year_vals_spring <- year_vals_fall <- list()
  mean_spring <- lci_spring <- uci_spring <-
    mean_fall <- lci_fall <- uci_fall <- vector()
  for (i in 1:length(spring)) {
    year_vals_spring[[i]] <- colMeans(spring[[i]][spring_keep, 6:ncol(spring[[i]])], na.rm = T)
    year_vals_fall[[i]] <- colMeans(fall[[i]][fall_keep, 6:ncol(spring[[i]])], na.rm = T)
    mean_spring[i] <- year_vals_spring[[i]]["mean"]
    mean_fall[i] <- year_vals_fall[[i]]["mean"]
    lci_spring[i] <- quantile(year_vals_spring[[i]][2:(ncol(spring[[i]]) - 5)], lq)
    uci_spring[i] <- quantile(year_vals_spring[[i]][2:(ncol(spring[[i]]) - 5)], uq)
    lci_fall[i] <- quantile(year_vals_fall[[i]][2:(ncol(spring[[i]]) - 5)], lq)
    uci_fall[i] <- quantile(year_vals_fall[[i]][2:(ncol(spring[[i]]) - 5)], uq)
  }
  mean_vec <- rep(NA, 2*length(spring))
  mean_vec[2*c(1:length(spring)) - 1] <- mean_spring
  mean_vec[2*c(1:length(spring))] <- mean_fall
  mean_lrat <- log(stocks::ratios(mean_vec))
  
  rep_lrat <- matrix(NA, nrow = 2*length(spring) - 1, ncol = ncol(spring[[1]]) - 6)
  for (i in 1:(ncol(spring[[1]]) - 6)){
    rep_vec <- rep(NA, 2*length(spring))
    for (j in 1:length(spring)) {
      rep_vec[2*j - 1] <- year_vals_spring[[j]][1+i]
      rep_vec[2*j] <- year_vals_fall[[j]][1+i]
    }
    rep_lrat[ , i] <- log(stocks::ratios(rep_vec))
  }
  
  lrat_lci <- apply(rep_lrat, 1, function(x){quantile(x, lq)})
  lrat_uci <- apply(rep_lrat, 1, function(x){quantile(x, uq)})
  
  return(list(mean_spring = mean_spring, lci_spring = lci_spring, uci_spring = uci_spring,
              mean_fall = mean_fall, lci_fall = lci_fall, uci_fall = uci_fall,
              mean_lrat = mean_lrat, lci_lrat = lrat_lci, uci_lrat = lrat_uci))
}


#' plotting function for summary returned by summarize_avg_and_rep
#' @param results_summary summary returned by summarize_avg_and_rep
#' @param speed_data_spring spring speed data from median_passage
#' @param speed_data_fall fall speed data from median_passage
#' @export
plot_one_summary <- function (results_summary, speed_data_spring = NULL, speed_data_fall = NULL) {
  plot(results_summary$mean_spring, xlab = "year", ylab = "index",
       ylim = c(min(results_summary$lci_spring), max(results_summary$uci_spring)),
       main = "spring index")
  for(i in 1:length(results_summary$mean_spring)) {
    lines(c(i,i), c(results_summary$lci_spring[i], results_summary$uci_spring[i]))
  }

  plot(results_summary$mean_fall, xlab = "year", ylab = "index",
       ylim = c(min(results_summary$lci_fall), max(results_summary$uci_fall)),
       main = "fall index")
  for(i in 1:length(results_summary$mean_spring)) {
    lines(c(i,i), c(results_summary$lci_fall[i], results_summary$uci_fall[i]))
  }

  if(!is.null(speed_data_spring)) {
    duration <- speed_data_spring$duration[(1+nrow(speed_data_spring)-length(results_summary$mean_spring)):nrow(speed_data_spring)]
    plot(results_summary$mean_spring ~ duration, xlab = "duration", ylab = "index",
         ylim = c(min(results_summary$lci_spring), max(results_summary$uci_spring)),
         main = "spring index vs duration")
    for(i in 1:length(results_summary$mean_spring)) {
      lines(c(duration[i],duration[i]), c(results_summary$lci_spring[i], results_summary$uci_spring[i]))
    }
  }

  if(!is.null(speed_data_fall)) {
    duration <- speed_data_fall$duration[(1+nrow(speed_data_fall)-length(results_summary$mean_fall)):nrow(speed_data_fall)]
    plot(results_summary$mean_fall ~ duration, xlab = "duration", ylab = "index",
         ylim = c(min(results_summary$lci_fall), max(results_summary$uci_fall)),
         main = "fall index vs duration")
    for(i in 1:length(results_summary$mean_fall)) {
      lines(c(duration[i],duration[i]), c(results_summary$lci_fall[i], results_summary$uci_fall[i]))
    }
  }
  
}

#' plotting function for two summaries returned by summarize_avg_and_rep,
#' one for a northern region and one for a southern region.
#' @param results_summary_n summary for the northern region
#' @param results_summary_s summary for the southern region
#' @export
plot_ns_summary <- function (results_summary_n, results_summary_s) {
  if(length(results_summary_n$mean_spring) != length(results_summary_s$mean_spring)){
    stop("input summaries have different numbers of years")
  }
  plot(results_summary_n$mean_spring, xlab = "year", ylab = "index", 
       ylim = c(min(c(results_summary_n$lci_spring, results_summary_s$lci_spring)), 
                max(c(results_summary_s$uci_spring, results_summary_s$uci_spring))),
       main = "spring index")
  for(i in 1:length(results_summary_n$mean_spring)) {
    points(i+.2, results_summary_s$mean_spring[i], pch = 16, col = "darkgoldenrod2")
    lines(c(i,i), c(results_summary_n$lci_spring[i], results_summary_n$uci_spring[i]), col = "darkmagenta")
    lines(c(i+.2,i+.2), c(results_summary_s$lci_spring[i], results_summary_s$uci_spring[i]), col = "darkgoldenrod2")
  }
  
  plot(results_summary_n$mean_fall, xlab = "year", ylab = "index", 
       ylim = c(min(c(results_summary_n$lci_fall, results_summary_s$lci_fall)), 
                max(c(results_summary_s$uci_fall, results_summary_s$uci_fall))),
       main = "fall index")
  for(i in 1:length(results_summary_n$mean_fall)) {
    points(i+.2, results_summary_s$mean_fall[i], pch = 16, col = "darkgoldenrod2")
    lines(c(i,i), c(results_summary_n$lci_fall[i], results_summary_n$uci_fall[i]), col = "darkmagenta")
    lines(c(i+.2,i+.2), c(results_summary_s$lci_fall[i], results_summary_s$uci_fall[i]), col = "darkgoldenrod2")
  }
  
  plot(results_summary_n$mean_lrat, xlab = "", ylab = "log-ratio",
       ylim = c(min(c(results_summary_n$lci_lrat, results_summary_s$lci_lrat)), 
                max(c(results_summary_n$uci_lrat, results_summary_s$uci_lrat))),
       main = "log_ratios")
  for (i in 1:length(results_summary_n$mean_lrat)) {
    points(i+.2, results_summary_s$mean_lrat[i], pch = 16, col = "darkgoldenrod2")
    lines(c(i,i), c(results_summary_n$lci_lrat[i], results_summary_n$uci_lrat[i]), col = "darkmagenta")
    lines(c(i+.2,i+.2), c(results_summary_s$lci_lrat[i], results_summary_s$uci_lrat[i]), col = "darkgoldenrod2")
    
  }
  
  plot(results_summary_n$mean_lrat ~ results_summary_s$mean_lrat, xlab = "south", ylab = "north", main = "log-ratios",
       col = rep(c("blue", "darkorange"), length(results_summary_n$mean_spring))[1:length(results_summary_n$mean_lrat)], pch = 16)
  for(i in 1:length(results_summary_n$mean_lrat)) {
    lines(rep(results_summary_s$mean_lrat[i], 2), c(results_summary_n$lci_lrat[i], results_summary_n$uci_lrat[i]))
    lines(c(results_summary_s$lci_lrat[i], results_summary_s$uci_lrat[i]), rep(results_summary_n$mean_lrat[i], 2))
  }
  points(results_summary_n$mean_lrat ~ results_summary_s$mean_lrat,
       col = rep(c("blue", "darkorange"), length(results_summary_n$mean_spring))[1:length(results_summary_n$mean_lrat)], pch = 16)
}


