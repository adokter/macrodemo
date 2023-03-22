#' perform weather regressions. currently expects to see data (and perform regressions)
#' for winter temp, winter swe, and summer temp, but this would be easily updated by
#' setting NULL defaults for these inputs plus a bit of control flow.
#' @param cell_ratios cell ratios object
#' @param cells_all all cells in cell_ratios
#' @param years years underlying the cell ratios object
#' @param min_n minimum sample size for which to estimate a slope
#' @param janfeb_temp jan/feb temps
#' @param julaug_temp july/august temps
#' @param decmar_swe dec-mar swe
#' @return 2-element list.  prob_surv gives posterior probability that survival variance
#'  is higher than productivity variance
#'  coef is the average log-scale effect.
#' @export
weather_regressions <- function(
  cell_ratios, cells_all, years, min_n,
  janfeb_temp, julaug_temp, decmar_swe
){
  plotting_data2 <- data.frame(
    cell = cells_all,
    surv_length = NA,
    prod_length = NA,
    full_length = NA,
    janfeb_mean_bayes = NA,
    janfeb_median_bayes = NA,
    janfeb_sd_bayes = NA,
    janfeb_skew_bayes = NA,
    janfeb_kurt_bayes = NA,
    janfeb_p_bayes = NA,
    julaug_mean_bayes = NA,
    julaug_median_bayes = NA,
    julaug_sd_bayes = NA,
    julaug_skew_bayes = NA,
    julaug_kurt_bayes = NA,
    julaug_p_bayes = NA,
    swe_mean_bayes = NA,
    swe_median_bayes = NA,
    swe_sd_bayes = NA,
    swe_skew_bayes = NA,
    swe_kurt_bayes = NA,
    swe_p_bayes = NA,
    janfeb_flag = NA, julaug_flag = NA, swe_flag = NA
  )

  cell_ratio_series <- cell_ratios$summary
  cell_ratio_series_full <- cell_ratios$replicates
  years_len <- length(years) - 1

  pb <- txtProgressBar(min = 1, max = length(cells_all), initial = 1)
  for (i in seq_along(cells_all)) {
    setTxtProgressBar(pb, i)
    if(i == 233){next} # This is a cell out over open ocean
    if (!identical(cell_ratio_series[[i]], NA)) {
      crsi <- cell_ratio_series[[i]]
      plotting_data2$surv_length[i] <- sum(is.finite(crsi$median[2*seq_len(years_len)]))
      plotting_data2$prod_length[i] <- sum(is.finite(crsi$median[2*seq_len(years_len + 1) - 1]))
      plotting_data2$full_length[i] <- sum(is.finite(crsi$median[1 + 2*seq_len(years_len)] + crsi$median[2*seq_len(years_len)]))

      janfeb_temp_slopes <- julaug_temp_slopes <- janfeb_swe_slopes <- vector()

      lrat_means <- rowMeans(cell_ratio_series_full[[i]])
      lrat_sds <- apply(cell_ratio_series_full[[i]], 1, sd)

      surv_means <- prod_means <- lrat_means
      surv_sds <- prod_sds <- lrat_sds
      surv_means[!use_cell_years(cell_ratio_series[[i]],
                                 inf_exclude = F,
                                 n_min_prod = 0,
                                 n_min_surv = min_n,
                                 n_min_full = 0)] <- NA
      surv_means <- surv_means[2*c(1:13)]
      surv_sds[!use_cell_years(cell_ratio_series[[i]],
                               inf_exclude = F,
                               n_min_prod = 0,
                               n_min_surv = min_n,
                               n_min_full = 0)] <- NA
      surv_sds <- surv_sds[2*seq_len(years_len)]

      prod_means[!use_cell_years(cell_ratio_series[[i]],
                                 inf_exclude = F,
                                 n_min_prod = min_n,
                                 n_min_surv = 0,
                                 n_min_full = 0)] <- NA
      prod_means <- prod_means[1 + 2*c(0:13)]
      prod_sds[!use_cell_years(cell_ratio_series[[i]],
                               inf_exclude = F,
                               n_min_prod = min_n,
                               n_min_surv = 0,
                               n_min_full = 0)] <- NA
      prod_sds <- prod_sds[2*seq_len(years_len + 1) - 1]

      if (sum(!is.na(surv_means)) > min_n) {
        surv_df <- data.frame(mean = surv_means, sd = surv_sds,
                              janfeb_temp = unlist(janfeb_temp[[i]][1:13]),
                              decmar_swe = sqrt(unlist(decmar_swe[[i]][1:13])))
        capture.output(suppressWarnings(suppressMessages(
          janfeb_mod <- brm(bf(mean | resp_se(sd, sigma = TRUE) ~ janfeb_temp),
                          data = surv_df, family = gaussian(),
                          prior = prior(std_normal(), class = "b"),
                          iter = 2000, warmup = 1000, chains = 3, refresh = 0,
                          backend = "cmdstanr")
        )))
        diagnostics <- check_brmsfit_diagnostics(janfeb_mod)
        if (!all(diagnostics)) {
          adapt_delta <- .8
          iter <- 2000
          if (!diagnostics[1]) {
            adapt_delta <- .99
          }
          if (!diagnostics[2]) {
            iter <- janfeb_iter[i] <- 4000
          }
          capture.output(suppressWarnings(suppressMessages(
            janfeb_mod <- brm(bf(mean | resp_se(sd, sigma = TRUE) ~ janfeb_temp),
                            data = surv_df, family = gaussian(),
                            prior = prior(std_normal(), class = "b"),
                            iter = iter, warmup = 1000, adapt_delta = adapt_delta,
                            chains = 3, refresh = 0,
                            backend = "cmdstanr")
          )))
          diagnostics <- check_brmsfit_diagnostics(janfeb_mod)
          if(!all(diagnostics)){
            plotting_data2$janfeb_flag[i] <- TRUE
          }
        }

        janfeb_temp_slopes <- as_draws_df(janfeb_mod)$b_janfeb_temp

        plotting_data2$janfeb_mean_bayes[i] <- mean(janfeb_temp_slopes)
        plotting_data2$janfeb_median_bayes[i] <- median(janfeb_temp_slopes)
        plotting_data2$janfeb_sd_bayes[i] <- sd(janfeb_temp_slopes)
        plotting_data2$janfeb_skew_bayes[i] <- moments::skewness(janfeb_temp_slopes)
        plotting_data2$janfeb_kurt_bayes[i] <- moments::kurtosis(janfeb_temp_slopes)
        plotting_data2$janfeb_p_bayes[i] <- mean(janfeb_temp_slopes > 0)

        # Here we hard-code that in addition to meeting the thresholds for inclusion,
        # we also need to see at least two years with nonzero snow.
        if (sum((surv_df$decmar_swe > 0) & (!is.na(surv_df$mean))) > 2) {
          capture.output(suppressWarnings(suppressMessages(
            swe_mod <- brms::brm(bf(mean | resp_se(sd, sigma = TRUE) ~ decmar_swe),
                               prior = prior(std_normal(), class = "b"),
                               data = surv_df, family = gaussian(),
                               iter = 2000, warmup = 1000, chains = 3, refresh = 0,
                               backend = "cmdstanr")
          )))
          diagnostics <- check_brmsfit_diagnostics(swe_mod)
          if (!all(diagnostics)) {
            adapt_delta <- .8
            iter <- 2000
            if (!diagnostics[1]) {
              adapt_delta <- .99
            }
            if (!diagnostics[2]) {
              iter <- swe_iter[i] <- 4000
            }
            capture.output(suppressWarnings(suppressMessages(
              swe_mod <- brm(bf(mean | resp_se(sd, sigma = TRUE) ~ decmar_swe),
                           data = surv_df, family = gaussian(),
                           prior = prior(std_normal(), class = "b"),
                           iter = iter, warmup = 1000, adapt_delta = adapt_delta,
                           chains = 3, refresh = 0,
                           backend = "cmdstanr")
            )))
            diagnostics <- check_brmsfit_diagnostics(swe_mod)
            if(!all(diagnostics)){
              plotting_data2$swe_flag[i] <- TRUE
            }
          }

          swe_slopes <- as_draws_df(swe_mod)$b_decmar_swe

          plotting_data2$swe_mean_bayes[i] <- mean(swe_slopes)
          plotting_data2$swe_median_bayes[i] <- median(swe_slopes)
          plotting_data2$swe_sd_bayes[i] <- sd(swe_slopes)
          plotting_data2$swe_skew_bayes[i] <- moments::skewness(swe_slopes)
          plotting_data2$swe_kurt_bayes[i] <- moments::kurtosis(swe_slopes)
          plotting_data2$swe_p_bayes[i] <- mean(swe_slopes > 0)
        }
      }

      if (sum(!is.na(prod_means)) > min_n) {
        prod_df <- data.frame(mean = prod_means, sd = prod_sds,
                              julaug_temp = unlist(julaug_temp[[i]]))
        capture.output(suppressWarnings(suppressMessages(
          julaug_mod <- brms::brm(bf(mean | resp_se(sd, sigma = TRUE) ~ julaug_temp),
                                data = prod_df, family = gaussian(),
                                prior = prior(std_normal(), class = "b"),
                                iter = 2000, warmup = 1000, chains = 3, refresh = 0,
                                backend = "cmdstanr")
        )))
        diagnostics <- check_brmsfit_diagnostics(julaug_mod)
        if (!all(diagnostics)) {
          adapt_delta <- .8
          iter <- 2000
          if (!diagnostics[1]) {
            adapt_delta <- .99
          }
          if (!diagnostics[2]) {
            iter <- julaug_iter[i] <- 4000
          }
          capture.output(suppressWarnings(suppressMessages(
            julaug_mod <- brm(bf(mean | resp_se(sd, sigma = TRUE) ~ julaug_temp),
                            data = surv_df, family = gaussian(),
                            prior = prior(std_normal(), class = "b"),
                            iter = iter, warmup = 1000, adapt_delta = adapt_delta,
                            chains = 3, refresh = 0,
                            backend = "cmdstanr")
          )))
          diagnostics <- check_brmsfit_diagnostics(julaug_mod)
          if(!all(diagnostics)){
            plotting_data2$julaug_flag[i] <- TRUE
          }
        }
        julaug_temp_slopes <- as_draws_df(julaug_mod)$b_julaug_temp

        plotting_data2$julaug_mean_bayes[i] <- mean(julaug_temp_slopes)
        plotting_data2$julaug_median_bayes[i] <- median(julaug_temp_slopes)
        plotting_data2$julaug_sd_bayes[i] <- sd(julaug_temp_slopes)
        plotting_data2$julaug_skew_bayes[i] <- moments::skewness(julaug_temp_slopes)
        plotting_data2$julaug_kurt_bayes[i] <- moments::kurtosis(julaug_temp_slopes)
        plotting_data2$julaug_p_bayes[i] <- mean(julaug_temp_slopes > 0)
      }
    }
  }
  plotting_data2
}
