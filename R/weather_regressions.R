#' Perform weather regressions.
#'
#' @param tidy_ratios ratios object, typically output of \link{make_ratios_tidy}
#' @param data_daymet weather data, typically output of \link{daymet_set_extract}
#' @param params_daymet data.frame specifying Daymet data extraction, see Details section of \link{daymet_set_extract}
#' @param min_n minimum sample size for which to estimate a slope
#' @return a data.frame with regression results
#' @export
#' @inheritParams compare_ratio_variances
weather_regressions <- function(tidy_ratios, data_daymet, params_daymet, min_n, warmup=1000, iter=2000, chains=3, quiet=FALSE){

  #initialize output data.frame
  df_out <- data.frame()

  tidy_ratio_series <- tidy_ratios$summary %>%
    filter(cell != 233) %>% # cell 233 is out over open ocean
    left_join(data_daymet, by=c("cell","year")) # add weather data

  years=unique(tidy_ratios$summary$year)
  years_len <- length(years) - 1
  cells_all <- sort(unique(tidy_ratios$summary$cell))

  # initiate a progress bar when quiet = T
  if(quiet) pb <- txtProgressBar(min = 1, max = length(cells_all), initial = 1)

  for (i in seq_along(cells_all)){
    # update the progress bar
    if(quiet) setTxtProgressBar(pb, i)
    # select data for a single cell
    # crsi stands for 'cell ratio series i'

    adapt_delta <- 0.8
    converged <- FALSE
    tries <- 0

    for(j in seq_along(nrow(params_daymet))){
      # the daymet parameter to use in regression:
      par = params_daymet$label[j]
      # the period (demographic index) to use in regression:
      period_demographic = params_daymet$period[j]
      # select the data for the regression:
      data_regression <- tidy_ratio_series %>%
        filter(cell==cells_all[i]) %>%
        filter(period==period_demographic) %>%
        filter(!is.na(avg))
      # initialize output
      output=data.frame(period=period_demographic,label=par,n=nrow(data_regression),
                        mean=NA,median=NA,sd=NA,skewness=NA,kurtosis=NA,p_value=NA,converged=NA)

      if(nrow(data_regression)>min_n){
        # add column `response` containing a copy of the data we want to regress:
        data_regression$predictor=data_regression[[par]]

        # construct brms regression formula:
        brms_formula <- bf(avg | resp_se(sd, sigma = TRUE) ~ predictor)

        while(!converged & tries < 2){

          if(!quiet) print(paste("starting estimation for cell",cells_all[i],"and parameter", par))

          # Here we hard-code that in addition to meeting the thresholds for inclusion,
          # we also need to see at least two years with nonzero snow.
          if(par=="swe" & (sum((data_regression$predictor > 0) & (!is.na(data_regression$avg))) > 2)) break

          # run brms model
          capture.output(suppressWarnings(suppressMessages(
            mod <- brm(brms_formula, data = data_regression, family = gaussian(),
                       prior = prior(std_normal(), class = "b"),
                       ter = iter, warmup = warmup, chains = chains, refresh = 0,
                       backend = "cmdstanr", silent=ifelse(quiet,2,1))
          )))

          diagnostics <- check_brmsfit_diagnostics(mod)
          if(all(diagnostics)) converged = TRUE

          if (!diagnostics[1]) {
            adapt_delta <- .99
          }
          if (!diagnostics[2]) {
            iter <- iter*2
          }

          tries = tries + 1

        } # while

        if(converged){
          # sample the posterior distribution:
          d <- as_draws_df(mod)

          slopes <- d$b_predictor

          output=data.frame(period=period_demographic,
                            label=par,
                            n=nrow(data_regression),
                            mean=mean(slopes),
                            median=median(slopes),
                            sd=sd(slopes),
                            skewness=moments::skewness(slopes),
                            kurtosis=moments::kurtosis(slopes),
                            p_value=mean(slopes > 0),
                            converged=converged
          )
        } else{
          output$converged=FALSE
          if(!quiet) print(paste("diagnostic failure for cell", cell_index))
        }
      } # if
    df_out <- rbind(df_out,output)
    } # for loop across daymet parameters
  } # for loop across cells
  df_out
}
