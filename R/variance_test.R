#' perform formal test of whether productivity or survival variance is larger
#' @param cell_ratios cell ratios object
#' @param cells_all all cells in cell_ratios
#' @param cells the cells over which output is desired
#' @param years years underlying the cell ratios object
#' @param min_n minimum sample size for which to estimate a variance
#' @return dataframe. prob_surv gives posterior probability that survival variance
#'  is higher than productivity variance
#'  coef is the average log-scale effect.
#' @export
variance_test <- function(cell_ratios, cells_all, cells, years, min_n = 5){
  cell_ratio_series <- cell_ratios$summary
  cell_ratio_series_full <- cell_ratios$replicates
  years_len <- length(years) - 1

  sd_holder_prod <- sd_holder_surv <- matrix(nrow = length(cells_all), ncol = 100)
  colnames(sd_holder_prod) <- paste0("prod_sd_rep_", 1:100)
  colnames(sd_holder_surv) <- paste0("surv_sd_rep_", 1:100)

  cell_lrat_sd <- cbind(data.frame(cell = cells_all, n_prod = NA, n_surv = NA),
                        as.data.frame(sd_holder_prod), as.data.frame(sd_holder_surv),
                        var_p = NA, var_d = NA, diagnostic_pass = NA)

  pb <- txtProgressBar(min = 1, max = length(cells_all), initial = 1)
  for(i in seq_along(cells_all)){
    setTxtProgressBar(pb, i)
    if(!identical(cell_ratio_series[[i]], NA)){
      lrats_avg <- cell_ratio_series[[i]]$avg
      lrats_avg[!use_cell_years(cell_ratio_series[[i]], inf_exclude=T)] <- NA
      if (cells_all[i] %in% cells & !all(is.na(lrats_avg))) {
        # This should be guaranteed by upstream data munging
        assertthat::assert_that(sum(!is.na(lrats_avg)) >= 5)
        prod_rats <- lrats_avg[1 + 2*seq_len(years_len)]
        surv_rats <- lrats_avg[2*seq_len(years_len)]
        cell_lrat_sd$n_prod[i] <- sum(!is.na(prod_rats))
        cell_lrat_sd$n_surv[i] <- sum(!is.na(surv_rats))
        lrat_means <- rowMeans(cell_ratio_series_full[[i]])
        lrat_sds <- apply(cell_ratio_series_full[[i]], 1, sd)

        prod_df <- data.frame(ratio = lrat_means[1 + 2*seq_len(years_len)],
                              sd = lrat_sds[1 + 2*seq_len(years_len)],
                              season = "prod")
        prod_df$ratio[is.infinite(prod_df$ratio) | is.nan(prod_df$ratio)] <- NA

        surv_df <- data.frame(ratio = lrat_means[2*c(1:13)],
                              sd = lrat_sds[2*c(1:13)],
                              season = "surv")
        surv_df$ratio[is.infinite(surv_df$ratio) | is.nan(surv_df$ratio)] <- NA


        model_df <- rbind(prod_df, surv_df)

        ##### model formula #####
        # modeling the ratio with known Gaussian measurement error `| resp_se...`
        # plus an additional modeled residual term `sigma = TRUE`
        # response varies by season `~ season`
        # residual sd also varies by season `sigma ~ season`
        mod_formula <- bf(ratio | resp_se(sd, sigma = TRUE) ~ season,
                          sigma ~ season)

        if(sum(!is.na(prod_df$ratio)) >= min_n & sum(!is.na(surv_df$ratio)) >= min_n){
          capture.output(suppressWarnings(suppressMessages(
            mod <- brm(mod_formula, data = model_df, family = gaussian(),
                       iter = 2000, warmup = 1000, chains = 3, refresh = 0,
                       backend = "cmdstanr")
          )))

          # check for divergences and poor r-hat (modern version of Gelman-Rubin diagnostic)
          # that checks split, folded, rank-normalized r-hat
          # if divergences occur, increase adapt_delta to adapt smaller step-sizes in leapfrog integrator (default .8)
          # if rhat is bad, run chains for longer (default 2000)
          diagnostics <- check_brmsfit_diagnostics(mod)
          if (!all(diagnostics)) {
            adapt_delta <- .8
            iter <- 2000
            if (!diagnostics[1]) {
              adapt_delta <- .99
            }
            if (!diagnostics[2]) {
              iter <- 4000
            }
            capture.output(suppressWarnings(suppressMessages(
              mod <- brm(mod_formula, data = model_df, family = gaussian(),
                         iter = iter, warmup = 1000, chains = 3, refresh = 0,
                         adapt_delta = adapt_delta,
                         backend = "cmdstanr")
            )))
            diagnostics <- check_brmsfit_diagnostics(mod)
          }
          if (all(diagnostics)) {
            cell_lrat_sd$diagnostic_pass[i] <- TRUE

            d <- as_draws_df(mod)
            # the coefficient of the sd (log link) is b_sigma_seasonsurv
            cell_lrat_sd$prob_surv[i] <- mean(d$b_sigma_seasonsurv > 0)
            cell_lrat_sd$coef[i] <- mean(d$b_sigma_seasonsurv)
          } else {
            cell_lrat_sd$diagnostic_pass[i] <- FALSE
          }
        }
      }
    }
  }
  cell_lrat_sd
}
