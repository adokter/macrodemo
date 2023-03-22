#' reformat the output of get_ratios() into tidy format
#' @param cell_ratios data output from [get_ratios()]
#' @return a list with two tidy data frames
#' @export
make_ratios_tidy <- function(cell_ratios){
  ############################
  # tidy up replicate data 
  ############################
  tidy_rep <- function(data_rep){
    assert_that(is.list(data_rep))
    # function expects a list of length one (in order to keep the list element name available)
    assert_that(length(data_rep)==1)
    output <- as_tibble(data_rep[[1]])
    replicate_cols <- colnames(output)
    output$cell=names(data_rep)[1]
    output$year <- as.numeric(substr(rownames(data_rep[[1]]),1,4))
    output$period <- substr(rownames(data_rep[[1]]),6,10^6)
    output <- tidyr::pivot_longer(output,all_of(replicate_cols),names_to = "replicate")
    output$replicate <- as.numeric(gsub("rep_","",output$replicate))
    output
  }
  data_rep <- cell_ratios$replicates[!is.na(cell_ratios$replicates)]
  output_rep <- do.call(rbind,lapply(seq_along(data_rep), function(i) tidy_rep(data_rep[i])))
  ############################
  # tidy up summary data 
  ############################ 
  # convert lists to tibble dataframe
  output_sum <- do.call(rbind,lapply(cell_ratios$summary[!is.na(cell_ratios$summary)], as_tibble))
  
  # add information on whether there are Inf values in a given cell time series
  output_sum <- left_join(output_sum, output_sum %>% 
    group_by(cell) %>%
    summarise(has_inf=sum(is.infinite(avg))>0), by="cell")
  
  # add information on number of productivity, survival and total indices per year
  # replace Inf values with NA
  
  period1=unique(output_sum$period)[1]
  period2=unique(output_sum$period)[2]
  
  left_join(output_sum, output_sum %>%
              mutate(across(everything(), ~ replace(., is.infinite(.),NA))) %>% # change infinite values to NA
              group_by(cell) %>%
              summarise(n_prod=sum(is.finite(avg) & period==period1), # typically evaluates to period=="fall
                        n_surv=sum(is.finite(avg) & period==period2), # typically evaluates to period=="fall
                        n=sum(is.finite(avg))),
            by="cell") -> output_sum
  
  # return both as a list:
  list(summary=output_sum, replicates=output_rep)
}
