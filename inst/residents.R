
####
#### Parameters -------------------------------------------------------------
####
params <- list()
if(Sys.info()[["user"]]=="amd427"){
  params$erd_path <- "~/Dropbox/macrodemography/erd/erd.db"
  params$output_path <- "~/Dropbox/macrodemography_refactor/data/residents"
} else if(Sys.info()[["user"]]=="bg423"){
  params$erd_path <- "~/Documents/macrodemography/data/erd.db"
  params$output_path <- "~/Documents/macrodemography/data"
}
params$years <- c(2006:2019)
params$extent_space <-  data.frame( min_lon=-125, max_lon=-66, min_lat=24, max_lat=50 )
params$period <- c("spring", "fall")
params$time_grid <- 7
params$tgrid_min <- c(13, 40)
params$tgrid_max <- c(16, 43)
params$max_altitude <- 2000
params$max_altitude_above_lat42 <- 1500
params$effort_thresholds <- data.frame( dist_max=3, time_min=5/60, time_max=1, cci_min=0 )
params$hexagon_area_large <- 70000
params$hexagon_area_small <- 300
params$n_small_min <- 10
params$n_year_min <- 5
params$daymet <- data.frame( label=c("tmax_winter","tmax_summer","swe"),
                             variable=c("tmax","tmax","swe"),
                             date_min=c("01-01","07-01","12-01"),
                             date_max=c("02-28","08-31","03-15"),
                              period=c("spring","fall","spring") )
params$species_to_process <- c("carwre", "norcar")
params$always_import_checklists  <- FALSE
params$always_filter_checklists  <- FALSE
params$always_resample_bootstrap <- FALSE
params$always_run_variance_test  <- FALSE
params$always_download_weather   <- FALSE
params$always_run_regressions    <- FALSE
params$quiet                     <- TRUE
params$region <- "eastern_us"
params$plotting_xlim <- c(-107, -65)
params$extent_time <-
  data.frame(
    period = params$period,
    tgrid_min = params$tgrid_min,
    tgrid_max = params$tgrid_max,
    year_min = min(params$years),
    year_max = max(params$years)
  )

####
#### Load packages -------------------------------------------------------------
####
# Install an older version of dggridR from source. See:
# https://github.com/r-barnes/dggridR/issues/63#issuecomment-1454929653
# remotes::install_github("r-barnes/dggridR", ref = "ec2a040")

# we make sure we have the same version of the erdPackage as this .rmd
package_path <- strsplit(rstudioapi::getSourceEditorContext()$path, "inst")[[1]][1]
install.packages(package_path, repos = NULL, type = 'source')

# package load
library(macrodemography)
library(data.table) # don't worry about openMP warnings on Mac.
library(brms)
library(ggplot2)
library(dplyr)
library(magrittr)
library(sf)
library(dtplyr) # enable dplyr for data.table
library(assertthat)
library(ebirdst)
library(fasterize)
library(dggridR)
library(rgee)
library(tidyr)
library(gridExtra)

# this package also needs a working version of rcmdstan:
# see https://mc-stan.org/cmdstanr/
cmdstanr::check_cmdstan_toolchain()
assert_that(cmdstanr::cmdstan_version() >= "2.31")

####
#### Color scales / themes -------------------------------------------------------------
####
# set color scales
cols_bd <- c(hsv(seq(0,.17,length.out = 100),1,seq(.9,.6,length.out = 100)), hsv(seq(.45,.65, length.out = 100),1,seq(.6,1,length.out = 100)))
cols_bd2 <- c(hsv(seq(0,.17,length.out = 100),seq(1, .2, length.out = 100),.9), hsv(seq(.45,.65, length.out = 100),seq(.2, 1, length.out = 100),.9))

# blank ggplot2 theme
blank_theme <-
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )

####
#### Regions of interest -------------------------------------------------------------
####
western_states <- c("arizona", "california", "colorado", "idaho", "montana",
                    "nevada", "new mexico", "oregon", "utah", "washington",
                    "wyoming")

if(params$region == "eastern_us"){
  region_of_interest <- spData::us_states %>% filter(!(tolower(NAME) %in% western_states))
} else if(params$region == "western_us") {
  region_of_interest <- spData::us_states %>% filter(toupper(NAME) %in% western_states)
} else if(params$region == "conus") {
  region_of_interest <- spData::us_states
} else if(params$region == "north_america") {
  region_of_interest <- ne_states(country =  c("United States of America","Canada"), returnclass = "sf") %>%
  filter(region != "Northern Canada") %>%
  filter(!name %in% c("Hawaii", "Alaska"))
} else{
  region_of_interest <- spData::us_states %>% filter(tolower(NAME) %in% params$region)
}

# NOTE: this adds additional cropping of the region of interest
# adjust if necessary
raster_of_interest <- fasterize::fasterize(region_of_interest,
                                   raster::raster(ncol=1000, nrow = 1000,
                                                  xmn = -125, xmx = -66,
                                                  ymn =24, ymx = 50))
####
#### Define hexagon grids -------------------------------------------------------------
####
# for an area of ~ 70000 km^2, we get a resolution of 6:
grid_large <- dggridR::dgconstruct(area = params$hexagon_area_large)
# for an area of ~ 300 km^2, we get a resolution of 11:
grid_small <- dggridR::dgconstruct(area = params$hexagon_area_small)

####
#### Import checklists -------------------------------------------------------------
####
# print some species info
ebirdst::ebirdst_runs %>% filter(substr(species_code,1,6) %in% species_codes$six)

checklists_path <- paste0(params$output_path, "/checklists.RDS")
filtered_checklists_path <- paste0(params$output_path, "/filtered_checklists.RDS")

if(params$always_import_checklists | !file.exists(checklists_path)){
  checklists <- import_checklists(params$erd_path)
  saveRDS(checklists, checklists_path)
}

if(params$always_filter_checklists | !file.exists(filtered_checklists_path)){
  checklists <- readRDS(checklists_path) %>%
    filter(latitude > params$extent_space$min_lat &
             latitude < params$extent_space$max_lat &
             longitude > params$extent_space$min_lon &
             longitude < params$extent_space$max_lon &
             year >= min(params$years) &
             year <= max(params$years) &
             ELEV_30M_MEDIAN < params$max_altitude &
             ((ELEV_30M_MEDIAN < params$max_altitude_above_lat42) | (latitude < 42)) &
             cci > params$effort_thresholds$cci_min &
             effort_distance_km <= params$effort_thresholds$dist_max &
             effort_hrs >= params$effort_thresholds$time_min &
             effort_hrs <= params$effort_thresholds$time_max
    ) %>%
    mutate(
      seqnum_large = dggridR::dgGEO_to_SEQNUM(
        grid_large,
        .$longitude,
        .$latitude
        )[[1]],
      seqnum_small = dggridR::dgGEO_to_SEQNUM(
        grid_small,
        .$longitude,
        .$latitude
        )[[1]]
      )
  saveRDS(checklists, filtered_checklists_path)
} else {
  checklists <- readRDS(filtered_checklists_path)
}

# extract unique large hexagons
cells_all <- unique(checklists$seqnum_large)

####
#### Bootstrap abundances  -------------------------------------------------------------
####
for(species_code in params$species_to_process){
  file_out <- paste0(params$output_path, "/abun_data/", species_code , ".rds")

  if(params$always_resample_bootstrap | !file.exists(file_out)){
    data <-
      sample_grid_abun(
        species_code, params$erd_path, checklists, params$effort_thresholds,
        params$extent_space, params$extent_time, time_window="full",
        small_grid=grid_small$res, large_grid=grid_large$res, time_grid=7,
        roi = raster_of_interest, quiet = FALSE
      )
    if(!dir.exists(dirname(file_out))) dir.create(dirname(file_out), recursive = TRUE)
    saveRDS(data, file_out)
  }
}

####
#### Calculate demographic indices (spring/fall log-ratios)  ------------------------------------------------------
####

for(species_code in params$species_to_process){
  ##### load abundance data
  file_species <- paste0(params$output_path, "/abun_data/", species_code , ".rds")

  print(paste("loading data from file", file_species,"..."))
  data <- readRDS(file_species)

# Get demographic indices
  cell_ratios <- get_ratios(data$abun, cells_all, n_small_min = params$n_small_min, quiet=params$quiet)
  tidy_ratios <- make_ratios_tidy(cell_ratios)

  # rename fall/spring ratios to productivity/recruitment:
  tidy_ratios$summary %>%
    mutate(season=ifelse(period=="fall", "prod","surv")) -> tidy_ratios$summary

  # save the ratio data
  saveRDS(tidy_ratios, paste0(params$output_path, "/abun_data/", species_code, "_ratios.rds"))
}

####
#### Plot demographic indices ------------------------------------------------------
####

# carwre
species_code <- c("norcar")
species_code <- c("carwre")
tidy_ratios <- readRDS(paste0(params$output_path, "/abun_data/", species_code, "_ratios.rds"))

# select cells with at least 20 valid ratios
tidy_ratios$summary %>%
  group_by(cell) %>%
  summarize(sufficient_data = sum(is.na(median))<20) %>%
  filter(sufficient_data) %>% pull(cell) -> cells_select

# plot selected cell 20 (example)
plot_ratios(cells_select[10], data=tidy_ratios$summary)

# number of analyzable timeseries
df <- tidy_ratios$summary
table(df$cell, df$season)

# number of cells
df2 <- df %>% group_by(cell, season) %>% summarise(n=n(), non_na=sum(!is.na(median))) %>%
  filter(season=="surv", non_na>=5)
df2 <- df %>% group_by(cell, season) %>% summarise(n=n(), non_na=sum(!is.na(median))) %>%
  filter(season=="prod", non_na>=5)

# Plotting demographic indices for all cells

# Initialize a list to hold the plots
plots <- list()
# Loop through each cell
for (i in seq_along(cells_select)) {
  # Generate a plot for the current cell
  p <- plot_ratios(cells_select[i], data=tidy_ratios$summary)+theme(legend.position = "none") +ylim(-1.5, 1.5)
  # Add the plot to the list of plots
  plots[[length(plots) + 1]] <- p
}

# Arrange the plots in a grid and display
p <- grid.arrange(grobs = plots, ncol = ceiling(sqrt(length(plots))), nrow=round(sqrt(length(plots))))

# save plot panel
ggsave(plot=p, paste0(params$output_path, species_code, "/plots/", "surv_prod_perYearCell.png"), width = 80, height= 40, units = "cm")


####
#### Verify normality assumptions ------------------------------------------------------
####
# Check the higher moments of the abundance log-ratios to verify adequacy of Gaussian approximations.
for(species_code in params$species_to_process){
  tidy_ratios <- readRDS(paste0(params$output_path, "/abun_data/", species_code, "_ratios.rds"))
  # excess kurtosis = kurtosis - 3
  ggplot(tidy_ratios$summary, aes(x=kurtosis-3)) +
    geom_histogram(binwidth=.1) +
    xlab("excess kurtosis") +
    ggtitle(species_code)
  # skewness
  ggplot(tidy_ratios$summary, aes(x=skewness)) +
    geom_histogram(binwidth=.1) +
    xlab("skewness") +
    ggtitle(species_code)
}

####
#### Compare recruitment vs mortality variances  ------------------------------------------------------
####

for(species_code in params$species_to_process){
  tidy_ratios <- readRDS(paste0(params$output_path, "/abun_data/", species_code, "_ratios.rds"))

  # compare variances in productivity and survival for each cell across years:
  var_save_path <- paste0(params$output_path, "/variance_results/", species_code,"/variance_test.rds")
  dir.create(dirname(var_save_path), recursive = TRUE, showWarnings = FALSE)
  if(params$always_run_variance_test | !file.exists(var_save_path)){
      print(paste0("processing_species ", species_code))
    data_compare_ratios <- lapply(sort(unique(tidy_ratios$summary$cell)),
                                  compare_ratio_variances, data=tidy_ratios$summary,
                                  n_ratio_min=params$n_year_min,
                                  quiet=params$quiet) %>% do.call(rbind, .)
    saveRDS(data_compare_ratios,file=var_save_path)
  } else{
    data_compare_ratios <- readRDS(file=var_save_path)
  }
}

####
#### Plot variances  ------------------------------------------------------
####

# helper function to plot cells on a map
plot_cells_on_map <- function(data, param, color_scale){
  ggplot() + blank_theme + coord_fixed() +
    geom_sf(data=region_of_interest, fill=NA, color="black") +
    geom_polygon(data=data, aes(x=long, y=lat, group=group, fill=eval(parse(text=param))), alpha=0.7) +
    geom_path(data=data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
    color_scale +
    labs(fill=param) +
    xlim(params$plotting_xlim)
}

plot1 <- list()
plot2 <- list()
plot3 <- list()
plot4 <- list()

for(species_code in params$species_to_process){
  var_save_path <- paste0(params$output_path, "/variance_results/", species_code,"/variance_test.rds")
  cell_lrat_sd <- readRDS(var_save_path)

  tidy_ratios$summary %>%
    select(cell, n_prod,n_surv,n) %>%
    group_by(cell) %>%
    filter(row_number()==1) %>%
    right_join(cell_lrat_sd, by="cell") %>%
    filter(n_prod >= params$n_year_min & n_surv >= params$n_year_min) %>%
    filter(!is.na(p_survival_variance_higher)) -> plotting_data

  dggridR::dgcellstogrid(
      grid_large,
      plotting_data$cell,
      frame=TRUE,
      wrapcells=TRUE
      ) %>%
    mutate(cell=as.numeric(cell)) %>%
    left_join(plotting_data, by="cell") -> grid2

  # define color scales
  scale_viridis <- viridis::scale_fill_viridis(limits = c(params$n_year_min, length(params$years)))
  scale_blue_red <- scale_fill_gradientn(colours = cols_bd2, na.value=NA, limits = c(0, 1))
#  fl <- max(abs(grid2$effect_size_log), na.rm = T) + .1
  fl <- 6.7

  # plot number of years with a productivity index
  p1<- plot_cells_on_map(grid2, "n_prod", color_scale=scale_viridis)+theme(legend.position = "none")+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  print(p1)

  # plot number of years with a survival index
  p2=plot_cells_on_map(grid2, "n_surv", color_scale=scale_viridis)+theme(legend.position = "none")+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  print(p2)

  # plot probability that survival variance is higher than productivity
  p3=plot_cells_on_map(grid2, "p_survival_variance_higher", color_scale=scale_blue_red)+
    theme(legend.position = "none")+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  print(p3)

  # plot the difference in survival and recruitment variance (log-effect size)
  # scale opacity by the probability that the survival variance is higher.
  p4=ggplot() + coord_fixed() + blank_theme +
    geom_sf(data=region_of_interest, fill=NA, color="black") +
    geom_polygon(data=grid2, aes(x=long, y=lat, group=group, fill = effect_size_log), alpha = 2*abs(grid2$p_survival_variance_higher - 0.5))   +
    geom_path(data=grid2, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
    scale_fill_gradientn(colours = cols_bd, na.value=NA, limits = c(-fl, fl), oob=scales::squish) +
    xlim(params$plotting_xlim)+theme(legend.position = "none")
  print(p4)

plot1[[species_code]] <- gridExtra::grid.arrange(p1)
plot2[[species_code]] <- gridExtra::grid.arrange(p2)
plot3[[species_code]] <- gridExtra::grid.arrange(p3)
plot4[[species_code]] <- gridExtra::grid.arrange(p4)

}

do.call(gridExtra::grid.arrange, c(plot1, ncol=2))
do.call(gridExtra::grid.arrange, c(plot2, ncol=2))
do.call(gridExtra::grid.arrange, c(plot3, ncol=2))
do.call(gridExtra::grid.arrange, c(plot4, ncol=2))


# re-producing figure 3 of the manuscript
# ---------------------------------------
library(cowplot)
grid1 <- plot_grid(plotlist = plot4, ncol=2, labels = "auto")

# Get legends from one of the plots
p4 <-
  ggplot() + coord_fixed() + blank_theme +
  geom_sf(data=region_of_interest, fill=NA, color="black") +
  geom_polygon(data=grid2, aes(x=long, y=lat, group=group, fill = effect_size_log), alpha = 2*abs(grid2$p_survival_variance_higher - 0.5))   +
  geom_path(data=grid2, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradientn(colours = cols_bd, na.value=NA, limits = c(-fl, fl), oob=scales::squish) +
  xlim(params$plotting_xlim)+labs(fill="log-sd difference \n (survival-productivity)")+
  guides(fill=guide_colorbar(title.position = "right", title.hjust = .5))

legend <- get_legend(p4+theme(legend.key.size = unit(.4, "cm")))

# Add common legend to the plots
fig3 <- grid.arrange(grid1, right=legend)
ggsave("~/Documents/macrodemography/data/carwre_norcar_plots/fig3.png", plot=fig3, width = 22, height= 7, units = "cm")


# re-producing figure S3 of the manuscript
# ---------------------------------------
library(cowplot)
grid1 <- plot_grid(plotlist = plot3, ncol=2, labels = "auto")

# Get legends from one of the plots
legend <- get_legend(plot_cells_on_map(grid2, "p_survival_variance_higher", color_scale=scale_blue_red)+
                       labs(fill="p(survival)")+theme(legend.key.size = unit(.4, "cm")))

# Add common legend to the plots
figS3 <- grid.arrange(grid1, right=legend)
ggsave("~/Documents/macrodemography/data/carwre_norcar_plots/figS3_3.png", plot=figS3, width = 22, height= 8, units = "cm")

# # Get legends from one of the plots
# legend <- get_legend(plot_cells_on_map(grid2, "p_survival_variance_higher", color_scale=scale_blue_red)+
#                        labs(fill="p(survival)")+theme(legend.position = c(.5, 1.5), legend.direction = "horizontal",
#                                                       legend.key.width = unit(1, "cm"))+
#                        guides(fill = guide_colorbar(title.position = "top", title.hjust=.5)))
#
# # Add common legend to the plots
# figS3 <- grid.arrange(grid1, legend, nrow = 2, heights = c(15, 1))
# ggsave("~/Documents/macrodemography/data/carwre_norcar_plots/figS3_2.png", plot=figS3, width = 22, height= 10, units = "cm")

# re-producing figure S2 of the manuscript
# ---------------------------------------
library(cowplot)
grid1 <- plot_grid(plotlist = c(plot1, plot2), nrow=2, labels = "auto")

# Get legends from one of the plots
legend <- get_legend(plot_cells_on_map(grid2, "n_prod", color_scale=scale_viridis)+labs(fill="timeseries\nlength (years)")+
                       scale_fill_viridis_c(limits=c(5,13)))

# Add common legend to the plots
figS2 <- grid.arrange(grid1, legend, ncol = 2, widths = c(3, .5))
ggsave("~/Documents/macrodemography/data/carwre_norcar_plots/figS2.png", plot=figS2, width = 27, height= 17, units = "cm")


####
#### Compare bayesian/frequentist averages  ------------------------------------------------------
####

# Quick comparisons of inverse-sd weighted averages and the Bayesian estimates.
# Using the last generated `tidy_ratios` object for now.

# average indices over cells across years (using inverse weighting by sd)
# (note that n_prod,n_surv and n are identical within the grouping, so
# the weighted.mean just repeats the group value)
tidy_ratios$summary %>%
  group_by(cell,season) %>%
  filter(is.finite(avg)) %>%
  summarise(across(c("median","avg","q10","q90","sd","n_prod","n_surv","n", "has_inf"), \(x) weighted.mean(x, 1/sd, na.rm=TRUE)), .groups="drop_last") %>%
  mutate(has_inf=as.logical(has_inf)) -> data_cell

# join variance test results to year-averaged cell data:
data_cell <- left_join(data_cell,data_compare_ratios,by="cell")


####
#### Load weather data  ------------------------------------------------------
####
dir.create(paste0(params$output_path, "/weather"), showWarnings = FALSE)
weather_file <- paste0(params$output_path, "/weather/",paste0(params$daymet$label, collapse = "-"),".rds")

if(!params$quiet) print(paste("loading/processing weather file",basename(weather_file)))

# Warning messages below are from `erdPackage::daymet_extract` which calls
# `sp::getSpPPolygonsIDSlots`, which raises the deprecation warning
# if(params$always_download_weather | !file.exists(weather_file)){
#   # initialize google earth engine
#   ee_Initialize()
#
#   # loop over cells and years
#   data_daymet=data.frame()
#   for (cell in cells_all) {
#     for (year in params$years){
#       print(paste("year = ",year,"cell = ",cell))
#       data_daymet <- rbind(data_daymet, daymet_set_extract(year, cell, grid_large, params$daymet))
#
#     }
#   }
#   saveRDS(data_daymet, weather_file)
# } else{
#   data_daymet <- readRDS(weather_file)
# }


flag_file <- paste0(params$output_path, "/weather/flag.rds")

if (params$always_download_weather | !file.exists(weather_file) | !file.exists(flag_file)) {
  # initialize google earth engine
  ee_Initialize()

  # load or create flag variable to indicate whether data extraction was successful or not
  flag <- if (file.exists(flag_file)) readRDS(flag_file) else data.frame(cell = character(), year = integer(), success = logical())

  # check if weather file exists and load data
  data_daymet <- if (file.exists(weather_file)) readRDS(weather_file) else data.frame()

  # loop over cells and years
  for (cell in cells_all) {
    for (year in params$years) {
      # check if data extraction was successful for this cell and year
      if (nrow(flag) == 0 || !any(flag$success[flag$cell == cell & flag$year == year])) {
        # attempt data extraction
        tryCatch({
          print(paste("year =", year, "cell =", cell))
          data_daymet <- rbind(data_daymet, daymet_set_extract(year, cell, grid_large, params$daymet))

          # update flag variable to indicate success
          flag <- rbind(flag, data.frame(cell = cell, year = year, success = TRUE))
          saveRDS(flag, flag_file)

          # save weather data after each iteration of the inner loop
          saveRDS(data_daymet, weather_file)
        }, error = function(e) {
          # if an error occurred during data extraction, print an error message
          message(paste("Error occurred for year =", year, "cell =", cell))
          # update flag variable to indicate failure
          flag <- rbind(flag, data.frame(cell = cell, year = year, success = FALSE))
          saveRDS(flag, flag_file)
        })
      }
    }
  }
} else {
  # load weather data if always_download_weather is false and the weather file exists
  data_daymet <- readRDS(weather_file)
}







# Perform the regressions of fluctuations (from a single cell, single season,
# across years) against weather variables.

####
#### Weather regressions  ------------------------------------------------------
####
for(species_code in params$species_to_process){
  regression_save_path <- paste0(params$output_path, "/regression_results/", species_code)
  dir.create(regression_save_path, recursive = TRUE, showWarnings = FALSE)

  if(params$always_run_regressions | !file.exists(paste0(regression_save_path, "/regressions.rds"))){
    print(paste0("processing_species ", species_code))

    file_path <- paste0(params$output_path, "/abun_data/", species_code, "_ratios.rds")
    tidy_ratios <- readRDS(file_path)

    data_regression <- weather_regressions(tidy_ratios, data_daymet,params$daymet,params$n_year_min, quiet=FALSE)

    saveRDS(data_regression, paste0(regression_save_path, "/regressions.rds"))
  }
}

####
#### Plot weather regressions------------------------------------------------------
####

# select a species:
species_code="carwre"
# load the data for that species:
regression_save_path <- paste0(params$output_path, "/regression_results/", species_code)
data_regression <- readRDS(paste0(regression_save_path, "/regressions.rds"))

# plot function for regression data:
plot_regression <- function(data, label_daymet, moment, x_lim, fill_lim="auto", alpha="auto", labels=FALSE, colors=cols_bd) {
  assert_that(moment %in% names(data),msg=paste("column",variable,"not found in data"))
  assert_that(label_daymet %in% unique(data_regression$label),msg=paste("daymet variable",variable,"not found in data"))
  data <- filter(data, label==label_daymet)
  # get the spatial polygons for the cells
  grid <- dggridR::dgcellstogrid(grid_large,data$cell,frame=TRUE,wrapcells=TRUE)
  # merge grid and regression data
  grid_data <- merge(grid,data,by="cell")
  # calculate average lat,long for cells:
  grid %>%
    group_by(cell) %>%
    summarize(long=mean(long), lat=mean(lat)) -> data_label
  # define opacity based on value of alpha parameter:
  if(is.number(alpha)) grid_data$alpha=alpha
  # when alpha is NULL, color according to the certainty that the coefficient is either positively or
  # negatively different from zero
  if(alpha=="auto") grid_data$alpha=2*abs(grid_data$p_value-0.5)

  if(identical(fill_lim,"auto")){
    max_moment = max(abs(grid_data[[moment]]), na.rm = T) + 0.1
    fill_lim = c(-max_moment,max_moment)
  }
  p <- ggplot() + coord_fixed() + blank_theme +
    geom_sf(data=region_of_interest, fill=NA, color="black")   +
    geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = !!sym(moment)), alpha=grid_data$alpha) +
    geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
    scale_fill_gradientn(colours = colors, na.value=NA, limits=fill_lim, oob=scales::squish) + xlim(x_lim)
  if(labels) p = p + geom_text(aes(x=long,y=lat,label=cell), data=data_label)
  print(p)

}

# plot excess kurtosis and skewness of the regression coefficient posterior:
data_regression %>%
  mutate(excess_kurtosis=kurtosis-3) %>%
  tidyr::pivot_longer(c(skewness, excess_kurtosis)) %>%
  group_by(label,name) %>%
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(c("name","label"))


# plot number of years of data on which the regression was based
# also add the cell number labels
plot_regression(data_regression, "tmax_winter", "n", params$plotting_xlim, fill_lim=c(0,15), alpha=.8, labels=TRUE, colors=viridisLite::viridis(100))
plot_regression(data_regression, "swe", "n", params$plotting_xlim, fill_lim=c(0,15), alpha=.8, labels=TRUE, colors=viridisLite::viridis(100))
plot_regression(data_regression, "tmax_summer", "n", params$plotting_xlim, fill_lim=c(0,15), alpha=.8, labels=TRUE, colors=viridisLite::viridis(100))

# plot the mean regression slopes
plot_regression(data_regression, "tmax_winter", "mean", params$plotting_xlim, alpha=.8)
plot_regression(data_regression, "swe", "mean", params$plotting_xlim, alpha=.8)
plot_regression(data_regression, "tmax_summer", "mean", params$plotting_xlim, alpha=.8)

# plot the posterior probability that the regression coefficient is larger than zero (p_value)
plot_regression(data_regression, "tmax_winter", "p_value", params$plotting_xlim, alpha=.8, fill_lim=c(0,1))
plot_regression(data_regression, "swe", "p_value", params$plotting_xlim, alpha=.8, fill_lim=c(0,1))
plot_regression(data_regression, "tmax_summer", "p_value", params$plotting_xlim, alpha=.8, fill_lim=c(0,1))

# plot the mean regression slopes, with transparency denoting the probability of being non-zero
# essentially a combining the previous two types of plots in one figure
plot_regression(data_regression, "tmax_winter", "mean", params$plotting_xlim)
plot_regression(data_regression, "swe", "mean", params$plotting_xlim)
plot_regression(data_regression, "tmax_summer", "mean", params$plotting_xlim)

# plot the kurtosis (note for Guassian distribution the kurtosis = 3)
plot_regression(data_regression, "tmax_winter", "kurtosis", params$plotting_xlim, fill_lim=c(3,9))
plot_regression(data_regression, "swe", "kurtosis", params$plotting_xlim, fill_lim=c(3,9))
plot_regression(data_regression, "tmax_summer", "kurtosis", params$plotting_xlim, fill_lim=c(3,9))

# plot the skewness
plot_regression(data_regression, "tmax_winter", "skewness", params$plotting_xlim)
plot_regression(data_regression, "swe", "skewness", params$plotting_xlim)
plot_regression(data_regression, "tmax_summer", "skewness", params$plotting_xlim)



# Re-producing weather regression figure(s) for ms (Unsmoothed slope + probability)
# ------------------------------------------------------------------------------------

plot_list1 <- list()
plot_list2 <- list()
plot_list3 <- list()

for (species_code1 in params$species_to_process) {

  # load the data for that species:
  regression_save_path <- paste0(params$output_path, "/regression_results/", species_code1)
  data_regression <- readRDS(paste0(regression_save_path, "/regressions.rds"))

  # plot the mean regression slopes, with transparency denoting the probability of being non-zero
  # essentially a combining the previous two types of plots in one figure
  p1 <- plot_regression(data_regression, "tmax_winter", "mean", params$plotting_xlim, fill_lim = c(-1,1))+
    theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, -2, 0), "cm"))

  p2 <- plot_regression(data_regression, "swe", "mean", params$plotting_xlim, fill_lim = c(-1,1))+
    theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

  p3 <- plot_regression(data_regression, "tmax_summer", "mean", params$plotting_xlim, fill_lim = c(-1,1))+
    theme(legend.position = "none")+theme(plot.margin = unit(c(-2, 0, 0, 0), "cm"))

  plot_list1[[species_code1]] <- gridExtra::grid.arrange(p1)
  plot_list2[[species_code1]] <- gridExtra::grid.arrange(p2)
  plot_list3[[species_code1]] <- gridExtra::grid.arrange(p3)
}

library(gridExtra)
library(cowplot)

# Get legends from one of the plots
legend <- get_legend(plot_regression(data_regression, "swe", "mean", params$plotting_xlim, fill_lim = c(-1,1))+
                       labs(fill="estimated regression\nslope (unsmoothed)"))

# Create the grid with titles
plots <- plot_grid(plotlist = c(plot_list1,plot_list2, plot_list3),
                              ncol = 2, nrow = 3, labels = "auto",
                   label_x = c(.05,.05,.05,.05,.05,.05),
                   label_y = c(.85,.85,.95,.95,1.05,1.05))

# Add common legend to the plots
figSX <- grid.arrange(plots,
                      right=legend,
                      top=grid::textGrob(c("Carolina Wren", "Northern Cardinal"), x=c(.20, .70), y=c(-.1,-.1)),
                      left=grid::textGrob(c("productivity ~ summer temperature", "survival ~ snow cover",
                                            "survival ~ winter temperature"), y=c(.22,.50,.80), x=.7, rot=90, gp=grid::gpar(fontface=3, fontsize=10)))


ggsave("~/Documents/macrodemography/data/carwre_norcar_plots/WeatherReg_slope_prob_plotArrange_unsmoothed2.png", plot=figSX, width = 24, height= 25, units = "cm")


# Re-producing weather regression figure(s) for ms (Unmoothed probability)
# ------------------------------------------------------------------------------------

plot_list1 <- list()
plot_list2 <- list()
plot_list3 <- list()

for (species_code1 in params$species_to_process) {

  # load the data for that species:
  regression_save_path <- paste0(params$output_path, "/regression_results/", species_code1)
  data_regression <- readRDS(paste0(regression_save_path, "/regressions.rds"))

  # plot the posterior probability that the regression coefficient is larger than zero (p_value)
p1 <- plot_regression(data_regression, "tmax_winter", "p_value", params$plotting_xlim, alpha=.8, fill_lim=c(0,1), colors = cols_bd2)+
  theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, -2, 0), "cm"))

p2 <- plot_regression(data_regression, "swe", "p_value", params$plotting_xlim, alpha=.8, fill_lim=c(0,1), colors = cols_bd2)+
  theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

p3 <- plot_regression(data_regression, "tmax_summer", "p_value", params$plotting_xlim, alpha=.8, fill_lim=c(0,1), colors = cols_bd2)+
  theme(legend.position = "none")+theme(plot.margin = unit(c(-2, 0, 0, 0), "cm"))

  plot_list1[[species_code1]] <- gridExtra::grid.arrange(p1)
  plot_list2[[species_code1]] <- gridExtra::grid.arrange(p2)
  plot_list3[[species_code1]] <- gridExtra::grid.arrange(p3)
}

library(gridExtra)
library(cowplot)

# Get legends from one of the plots
legend <- get_legend(plot_regression(data_regression, "tmax_winter", "p_value", params$plotting_xlim, alpha=.8, fill_lim=c(0,1), colors = cols_bd2)+
                       labs(fill="posterior probability\nof positive slope\n(unsmoothed model)"))

# Create the grid with titles
plots <- plot_grid(plotlist = c(plot_list1,plot_list2, plot_list3),
                   ncol = 2, nrow = 3, labels = "auto",
                   label_x = c(.05,.05,.05,.05,.05,.05),
                   label_y = c(.85,.85,.95,.95,1.05,1.05))

# Add common legend to the plots
figSX <- grid.arrange(plots,
                      right=legend,
                      top=grid::textGrob(c("Carolina Wren", "Northern Cardinal"), x=c(.20, .70), y=c(-.1,-.1)),
                      left=grid::textGrob(c("productivity ~ summer temperature", "survival ~ snow cover",
                                            "survival ~ winter temperature"), y=c(.22,.50,.80), x=.7, rot=90, gp=grid::gpar(fontface=3, fontsize=10)))


ggsave("~/Documents/macrodemography/data/carwre_norcar_plots/WeatherReg_prob_plotArrange_unsmoothed2.png", plot=figSX, width = 24, height= 25, units = "cm")


####
#### Smoothing using ICAR/CAR models: latitude regression----------------------------------------------
####

# prepare CAR model data from weather regression data
data_regression %>%
  mutate(lat = dggridR::dgSEQNUM_to_GEO(grid_large, cell)$lat_deg,
         lon = dggridR::dgSEQNUM_to_GEO(grid_large, cell)$lon_deg) %>%
  filter(lon > params$plotting_xlim[1]) %>%
  # remove a region where we have sparse data
  filter(!(lat > 38.7 & lon < -99.5)) %>%
  filter(!is.na(mean)) %>%
  group_by(label) %>%
  # rescale the data (to help with model convergence)
  mutate(slope_scaled = scale(`mean`),
         lat_scaled = scale(lat),
         known_se = `sd`/sd(`mean`)) -> car_data

# verify kurtosis and skewness
car_data %>%
  mutate(excess_kurtosis=kurtosis-3) %>%
  pivot_longer(c(skewness, excess_kurtosis)) %>%
  group_by(label,name) %>%
  ggplot(aes(value)) +
  geom_density() +
  facet_wrap(c("name","label"))

# filter data, add upper lower bounds to the slope
car_data %>%
  filter(label=="tmax_winter") %>%
  arrange(lat) %>%
  mutate(lower=slope_scaled-known_se) %>%
  mutate(upper=slope_scaled+known_se) %>%
  mutate(lat_jit = lat + rnorm(n(), 0, .5)) -> car_data_select

macrodemography:::get_adjacency_matrix(car_data_select$cell, grid_large) -> adjacency_mat

icar_fit_lat <-
  brm(
    slope_scaled | resp_se(known_se, sigma = TRUE) ~
      lat_scaled + car(M, gr = cell, type = "icar"),
    data = car_data_select,
    data2 = list(M = adjacency_mat),
    backend = 'cmdstanr',
    iter = 12000, warmup = 2000,
    cores = 4, adapt_delta = .8, max_treedepth = 11, refresh = 0)

summary(icar_fit_lat)

# backtransform the model output using the inverse of the scale() operation
post_summ <- data.frame(posterior_summary(icar_fit_lat, variable = c("b_Intercept", "b_lat_scaled", "sdcar")))
post_summ2 <- post_summ * sd(car_data_select$mean) +
  mean(car_data_select$mean)

# back-transform according to Jacob's suggestion

# The relevant ICAR model regresses slope_scaled ~ lat_scaled.
# slope_scaled is in units of scaled(logarithms per degree C).
# Let's call the standard deviation associated with that scaling operation sd1.
# lat_scaled is in units of scaled(latitude). Call the associated standard deviation sd2.
# To go from the model-reported slope, in units of scaled(logarithms per degree C) per scaled(latitude)
# to units of logarithms per degree C per latitude, we need to multiply the numerator by sd1
# and multiply the denominator by sd2, or in other words multiply the slope by sd1/sd2.
# We don't need to add any means back in because translations don't affect the slope.

post_summ <- data.frame(posterior_summary(icar_fit_lat, variable = c("b_Intercept", "b_lat_scaled", "sdcar")))
post_summ3 <- post_summ*(sd(car_data_select$mean)/sd(car_data_select$lat))

# verify BFMI statistic is sufficiently high
rstan::get_bfmi(icar_fit_lat$fit)

# QUESTION: why re_formula = NA, i.e. not including group-level effects
# QUESTION: why incl_autor = F, i.e. do not include correlation structures in the predictions
# sample the posterior distribution:
draws <-  posterior_epred(icar_fit_lat, re_formula = NA, incl_autocor = F)
# extract quantiles for each cell
draws_quantiles <- as_tibble(t(sapply(1:ncol(draws),function(idx) quantile(draws[,idx],c(.5,seq(.05,1,.1))))))
# bind with original data
car_data_select %>%
  select(cell, lat, slope_scaled, known_se, lower,upper, lat_jit) %>%
  rename(slope=slope_scaled) %>%
  cbind(draws_quantiles) -> q_frame

# scaled certainty, inverse of the standard error:
certainty <- min(q_frame$known_se)/q_frame$known_se

# backtransform the data using the inverse of the scale() operation
q_frame[, ! (names(q_frame) %in% c("label","cell","lat", "lat_jit"))] <-
  q_frame[, ! (names(q_frame) %in% c("label","cell","lat", "lat_jit"))] * sd(car_data_select$mean) +
  mean(car_data_select$mean)
# bring up points outside of the plotting range
q_frame$lower[q_frame$lower < -.5] <- -.5
# color settings
fill_color <- "salmon2"
alpha <- .5
point_color <- "gray25"
# plot the data and the scatter points, opacity by their certainty:
ggplot(q_frame) + theme_classic() +
  geom_ribbon(aes(x = lat, ymin = `5%`, ymax = `95%`), fill = fill_color, alpha = alpha) +
  geom_ribbon(aes(x = lat, ymin = `15%`, ymax = `85%`), fill = fill_color, alpha = alpha) +
  geom_ribbon(aes(x = lat, ymin = `25%`, ymax = `75%`), fill = fill_color, alpha = alpha) +
  geom_ribbon(aes(x = lat, ymin = `35%`, ymax = `65%`), fill = fill_color, alpha = alpha) +
  geom_ribbon(aes(x = lat, ymin = `45%`, ymax = `55%`), fill = fill_color, alpha = alpha) +
  geom_point(aes(x = lat_jit, y = slope, alpha = certainty), color = point_color) +
  geom_segment(aes(x = lat_jit, y = lower, xend = lat_jit, yend = upper, alpha = certainty), color = point_color) +
  geom_line(aes(x = lat, y = `50%`)) +
  ylim(c(-.5, .5))
# plot only the model fit:
ggplot(q_frame) + theme_classic() +
  geom_ribbon(aes(x = lat, ymin = `5%`, ymax = `95%`), fill = fill_color, alpha = alpha) +
  geom_line(aes(x = lat, y = `50%`)) +
  ylim(c(-.5, .5))


####
#### Smoothing using ICAR/CAR models: smoothed maps----------------------------------------------
####

icar_regression <- function(data, adjancency_mat, warmup=2000, cores=4, iter=12000){
  assert_that(nrow(adjacency_mat) == ncol(adjacency_mat))
  assert_that(nrow(data) == nrow(adjacency_mat))
  # ICAR model without latitude
  print("Estimating icar model ...")
  icar_fit <-
    brm(slope_scaled | resp_se(known_se, sigma = TRUE) ~
          car(M, gr = cell, type = "icar"),
        data = data,
        data2 = list(M = adjacency_mat),
        backend = 'cmdstanr',
        iter = iter, warmup = warmup, cores = cores, refresh = 0)
  # number of cells:
  npt <- nrow(data)
  # number of posterior draws:
  n_samples = cores*(iter-warmup)
  # initialize matrix to contain posterior draws for each cell:
  true_values <- matrix(nrow = n_samples, ncol = npt)

  print("Drawing from posterior ...")
  for(i in 1:npt) {
    print(paste("cell",i,"/",npt,"..."))
    # the (unsmoothed) posterior slope:
    M <- data$slope_scaled[i]
    # the posterior uncertainty from the weather regression:
    se <- data$known_se[i]
    # get posterior draws of the linear predictor
    # re.form=NULL (default), i.e. include all group-level effects
    # incl_autor=TRUE (default), i.e. include correlation structures in prediction
    LP <-  posterior_linpred(icar_fit, re.form = NULL, incl_autocor = T)[ , i]
    # get draws of the unobserved latent uncertainty:
    sigma <- as_draws_df(icar_fit)$sigma
    # error propagation, combining sigma form CAR model and se from weather regression
    # QUESTION: why the linear predictor here
    mu = (M*sigma^2 + LP * se^2)/(sigma^2 + se^2)
    scale = 1/sqrt(1/sigma^2 + 1/se^2)
    true_values[ , i] <- rnorm(n_samples, mu, scale)
  }
  unscale <- function(x){x + mean(data$mean)/sd(data$mean)}
  smooth_prob=apply(true_values, 2,function(x){mean(unscale(x) > 0)})
  smooth_mean=apply(true_values, 2,function(x){mean(unscale(x))})
  smooth_data=data.frame(cell=data$cell, smooth_prob=smooth_prob, smooth_mean=smooth_mean)
  output_data <- merge(data,smooth_data, by = "cell")
  output_data
}

# plots the mean probability for a nonzero effect
plot_smooth_prob <- function(grid_data, region_of_interest){
  ggplot() + coord_fixed() + blank_theme +
    geom_sf(data=region_of_interest, fill=NA, color="black")   +
    geom_polygon(data=grid_data, aes(x=long, y=lat.x, group=group, fill = smooth_prob), alpha = .8)   +
    geom_path(data=grid_data, aes(x=long, y=lat.x, group=group), alpha=0.4, color="white") +
    scale_fill_gradientn(colours = cols_bd2, na.value=NA, limits = c(0,1), oob=scales::squish) +
    xlim(params$plotting_xlim)
}

# plots the mean smoothed effect size, with opacity according to probability
plot_smooth_mean <- function(grid_data, region_of_interest){
  fl <- max(abs(grid_data$smooth_mean), na.rm = T) + .1
  fl <- 1. # Baagi has changed this for consistency for manuscript!
  ggplot() + coord_fixed() + blank_theme +
    geom_sf(data=region_of_interest, fill=NA, color="black")   +
    geom_polygon(data=grid_data, aes(x=long, y=lat.x, group=group, fill = smooth_mean), alpha = 2*abs(grid_data$smooth_prob - 0.5))   +
    geom_path(data=grid_data, aes(x=long, y=lat.x, group=group), alpha=0.4, color="white") +
    scale_fill_gradientn(colours = cols_bd, na.value=NA, limits = c(-fl, fl), oob=scales::squish) +
    xlim(params$plotting_xlim)
}

merge_grid_to_data <- function(data, grid_large){
  grid <- dggridR::dgcellstogrid(grid_large,data$cell,frame=TRUE,wrapcells=TRUE)
  grid <- merge(grid,data,by=c("cell"))
  # remove smoothed data equal to NA
  grid[!is.na(grid$smooth_prob), ]
}


# First, we need to run ICAR regressions and save both "output data" into respective folders

for(species_code in params$species_to_process) {

  # load the data for that species:
  regression_save_path <- paste0(params$output_path, "/regression_results/", species_code)
  data_regression <- readRDS(paste0(regression_save_path, "/regressions.rds"))

  # prepare CAR model data from weather regression data
  data_regression %>%
    mutate(lat = dggridR::dgSEQNUM_to_GEO(grid_large, cell)$lat_deg,
           lon = dggridR::dgSEQNUM_to_GEO(grid_large, cell)$lon_deg) %>%
    filter(lon > params$plotting_xlim[1]) %>%
    # remove a region where we have sparse data
    #  filter(!(lat > 38.7 & lon < -99.5)) %>%
    filter(!is.na(mean)) %>%
    group_by(label) %>%
    # rescale the data (to help with model convergence)
    mutate(slope_scaled = scale(`mean`),
           lat_scaled = scale(lat),
           known_se = `sd`/sd(`mean`)) -> car_data

  # tmax_winter
  car_data %>% filter(label=="tmax_winter") -> car_data_select
  adjacency_mat <- macrodemography:::get_adjacency_matrix(car_data_select$cell, grid_large)
  icar_smooth_tmax_winter <- icar_regression(car_data_select, adjacency_mat)
  grid_data_surv_tmax <- merge_grid_to_data(icar_smooth_tmax_winter, grid_large)

  wintT_pval_sm    <- plot_smooth_prob(grid_data = grid_data_surv_tmax, region_of_interest = region_of_interest)+labs(subtitle = species_code)
  wintT_sl_pval_sm <- plot_smooth_mean(grid_data = grid_data_surv_tmax, region_of_interest = region_of_interest)+labs(subtitle = species_code)

  # tmax_summer
  car_data %>% filter(label=="tmax_summer") -> car_data_select
  adjacency_mat <- macrodemography:::get_adjacency_matrix(car_data_select$cell, grid_large)
  icar_smooth_tmax_summer <- icar_regression(car_data_select, adjacency_mat)
  grid_data_prod_tmax <- merge_grid_to_data(icar_smooth_tmax_summer, grid_large)

  summerT_pval_sm    <- plot_smooth_prob(grid_data = grid_data_prod_tmax, region_of_interest = region_of_interest)+labs(subtitle = species_code)
  summerT_sl_pval_sm <- plot_smooth_mean(grid_data = grid_data_prod_tmax, region_of_interest = region_of_interest)+labs(subtitle = species_code)

  # swe
  car_data %>% filter(label=="swe") -> car_data_select
  adjacency_mat <- macrodemography:::get_adjacency_matrix(car_data_select$cell, grid_large)
  icar_smooth_swe <- icar_regression(car_data_select, adjacency_mat)
  grid_data_surv_snow <- merge_grid_to_data(icar_smooth_swe, grid_large)

  snow_pval_sm    <- plot_smooth_prob(grid_data = grid_data_surv_snow, region_of_interest = region_of_interest)+labs(subtitle = species_code)
  snow_sl_pval_sm <- plot_smooth_mean(grid_data = grid_data_surv_snow, region_of_interest = region_of_interest)+labs(subtitle = species_code)


  # merge and save all output data
  weatherReg_icar_data <- rbind(grid_data_surv_tmax, grid_data_prod_tmax, grid_data_surv_snow)
  saveRDS(weatherReg_icar_data, paste0(params$output_path, "/regression_results/", species_code, "/weather_regressions_smoothed.rds"))

}


# Re-producing weather regression figure(s) for ms (Smoothed slope + probability)
# ------------------------------------------------------------------------------------

plot_list1 <- list()
plot_list2 <- list()
plot_list3 <- list()

for (species_code in params$species_to_process) {

# load the data for that species:
regression_save_path <- paste0(params$output_path, "/regression_results/", species_code)
data_regression <- readRDS(paste0(regression_save_path, "/weather_regressions_smoothed.rds"))

# tmax_winter
grid_data <- data_regression %>% filter(label=="tmax_winter")
p1 <- plot_smooth_mean(grid_data = grid_data, region_of_interest = region_of_interest)+
  theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, -2, 0), "cm"))

# swe
grid_data <- data_regression %>% filter(label=="swe")
p2 <- plot_smooth_mean(grid_data = grid_data, region_of_interest = region_of_interest)+
  theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

# tmax summer
grid_data <- data_regression %>% filter(label=="tmax_summer")
p3 <- plot_smooth_mean(grid_data = grid_data, region_of_interest = region_of_interest)+
  theme(legend.position = "none")+theme(plot.margin = unit(c(-2, 0, 0, 0), "cm"))

  plot_list1[[species_code]] <- gridExtra::grid.arrange(p1)
  plot_list2[[species_code]] <- gridExtra::grid.arrange(p2)
  plot_list3[[species_code]] <- gridExtra::grid.arrange(p3)
}

library(gridExtra)
library(cowplot)

# Get legends from one of the plots
legend <- get_legend(plot_smooth_mean(grid_data = grid_data, region_of_interest = region_of_interest)+
                       labs(fill="estimated regression\nslope (smoothed)"))

# Create the grid with titles
plots <- plot_grid(plotlist = c(plot_list1,plot_list2, plot_list3),
                   ncol = 2, nrow = 3, labels = "auto", label_x = c(.05,.05,.05,.05,.05,.05),
                                                        label_y = c(.85,.85,.95,.95,1.05,1.05))

# Add common legend to the plots
figSX <- grid.arrange(plots,
                      right=legend,
                      top=grid::textGrob(c("Carolina Wren", "Northern Cardinal"), x=c(.20, .70), y=c(-.1,-.1)),
                      left=grid::textGrob(c("productivity ~ summer temperature", "survival ~ snow cover",
                                            "survival ~ winter temperature"), y=c(.22,.50,.80), x=.7, rot=90, gp=grid::gpar(fontface=3, fontsize=10)))

ggsave("~/Documents/macrodemography/data/carwre_norcar_plots/WeatherReg_slope_prob_plotArrange_Smoothed.png", plot=figSX, width = 24, height= 25, units = "cm")


# Re-producing weather regression figure(s) for ms (probability)
# ------------------------------------------------------------------------------------

plot_list1 <- list()
plot_list2 <- list()
plot_list3 <- list()

for (species_code in params$species_to_process) {

  # load the data for that species:
  regression_save_path <- paste0(params$output_path, "/regression_results/", species_code)
  data_regression <- readRDS(paste0(regression_save_path, "/weather_regressions_smoothed.rds"))

  # tmax_winter
  grid_data <- data_regression %>% filter(label=="tmax_winter")
  p1 <- plot_smooth_prob(grid_data = grid_data, region_of_interest = region_of_interest)+
    theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, -2, 0), "cm"))

  # swe
  grid_data <- data_regression %>% filter(label=="swe")
  p2 <- plot_smooth_prob(grid_data = grid_data, region_of_interest = region_of_interest)+
    theme(legend.position = "none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

  # tmax summer
  grid_data <- data_regression %>% filter(label=="tmax_summer")
  p3 <- plot_smooth_prob(grid_data = grid_data, region_of_interest = region_of_interest)+
    theme(legend.position = "none")+theme(plot.margin = unit(c(-2, 0, 0, 0), "cm"))

  plot_list1[[species_code]] <- gridExtra::grid.arrange(p1)
  plot_list2[[species_code]] <- gridExtra::grid.arrange(p2)
  plot_list3[[species_code]] <- gridExtra::grid.arrange(p3)
}

library(gridExtra)
library(cowplot)

# Get legends from one of the plots
legend <- get_legend(plot_smooth_prob(grid_data = grid_data, region_of_interest = region_of_interest)+
                       labs(fill="posterior probability\nof positive slope\n(smoothed model)"))

# Create the grid with titles
plots <- plot_grid(plotlist = c(plot_list1,plot_list2, plot_list3),
                   ncol = 2, nrow = 3, labels = "auto", label_x = c(.05,.05,.05,.05,.05,.05),
                   label_y = c(.85,.85,.95,.95,1.05,1.05))

# Add common legend to the plots
figSX <- grid.arrange(plots,
                      right=legend,
                      top=grid::textGrob(c("Carolina Wren", "Northern Cardinal"), x=c(.20, .70), y=c(-.1,-.1)),
                      left=grid::textGrob(c("productivity ~ summer temperature", "survival ~ snow cover",
                                            "survival ~ winter temperature"), y=c(.22,.50,.80), x=.7, rot=90, gp=grid::gpar(fontface=3, fontsize=10)))

ggsave("~/Documents/macrodemography/data/carwre_norcar_plots/WeatherReg_prob_plotArrange_Smoothed.png", plot=figSX, width = 24, height= 25, units = "cm")


