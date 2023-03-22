#' plot the demographic ratio data for a large cell
#' @param cell_index the index of a large cell
#' @param data a tidy dataframe of ratio data, typically output from [make_ratios_tidy]
#' @param log whether to plot ratios on a logarithmic scale
#' @export
plot_ratios <- function(cell_index, data, log=TRUE){
  if(!log){
    data$median=exp(data$median)
    data$q10=exp(data$q10)
    data$q90=exp(data$q90)
    ylabel="ratio"
  } else{
    ylabel="log-ratio"
  }
  data %>% filter(cell==cell_index) %>%
    group_by(season) %>%
    ggplot(aes(x=year,y=median, col=season)) + 
    geom_errorbar(aes(ymin=q10, ymax=q90), width=.1,
                  position=position_dodge(0.2)) + 
    geom_point() + 
    geom_line() + 
    ggtitle(paste("cell index =",cell_index)) + 
    ylab(ylabel) 
}
