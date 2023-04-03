

# Parameters ----

params <- list()

params$erd_path <- "~/Documents/erd/erd.db"
params$output_path <- "~/Documents/macrodemography/data/"
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
params$daymet <- data.frame( label=c("tmax_winter","tmax_summer","swe"), variable=c("tmax","tmax","swe"), date_min=c("01-01","07-01","12-01"), date_max=c("02-28","08-31","03-15"), period=c("spring","fall","spring") )
params$species_to_process <- c("carwre")
params$always_import_checklists  <- FALSE
params$always_filter_checklists  <- FALSE
params$always_resample_bootstrap <- FALSE
params$always_run_variance_test  <- FALSE
params$always_download_weather   <- FALSE
params$always_run_regressions    <- FALSE
params$quiet                     <- TRUE
params$region <- "eastern_us"
params$plotting_xlim <- c(-107, -65)

# Packages ----
# Install an older version of dggridR from source. See:
# https://github.com/r-barnes/dggridR/issues/63#issuecomment-1454929653
remotes::install_github("r-barnes/dggridR", ref = "ec2a040")

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

# this package also needs a working version of rcmdstan:
# see https://mc-stan.org/cmdstanr/
cmdstanr::check_cmdstan_toolchain()
assert_that(cmdstanr::cmdstan_version() >= "2.31")


# misc-setup ----
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

# define time extent ----
# we can't define this in the run params because years isn't yet defined
extent_time <-
  data.frame(
    period = params$period,
    tgrid_min = params$tgrid_min, 
    tgrid_max = params$tgrid_max, 
    year_min = min(params$years), 
    year_max = max(params$years)
  )

# create dggrids----

# for an area of ~ 70000 km^2, we get a resolution of 6:
grid_large <- dggridR::dgconstruct(area = params$hexagon_area_large)
# for an area of ~ 300 km^2, we get a resolution of 11:
grid_small <- dggridR::dgconstruct(area = params$hexagon_area_small)

# species info ----
ebirdst::ebirdst_runs %>% filter(substr(species_code,1,6) %in% species_codes$six)

# import checklists----
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

# bootstrap-abundance ----

for(species_code in params$species_to_process){
  file_out <- paste0(params$output_path, "/abun_data/", species_code , ".rds")

  if(params$always_resample_bootstrap | !file.exists(file_out)){
    data <- 
      sample_grid_abun(
        species_code, params$erd_path, checklists, params$effort_thresholds, 
        params$extent_space, extent_time, time_window="full", 
        small_grid=grid_small$res, large_grid=grid_large$res, time_grid=7, 
        roi = raster_of_interest, quiet = FALSE
      )
    if(!dir.exists(dirname(file_out))) dir.create(dirname(file_out), recursive = TRUE)
    saveRDS(data, file_out)
  }
}

# Calculate spring/fall log-ratios ----

for(species_code in params$species_to_process){
  ##### load abundance data ----
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

# plot demographic indices for specific cells ----

species_code="carwre"
tidy_ratios <- readRDS(paste0(params$output_path, "/abun_data/", species_code, "_ratios.rds"))

# select cells with at least 20 valid ratios
tidy_ratios$summary %>%
  group_by(cell) %>%
  summarize(sufficient_data = sum(is.na(median))<20) %>%  
  filter(sufficient_data) %>% pull(cell) -> cells_select

# plot selected cell 20 (example)
plot_ratios(cells_select[20], data=tidy_ratios$summary)

# check normality assumptions ----
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

# compare recruitment/mortality variances ----

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


# plot variances ----

# I don't think there's anything random in this chunk; setting a seed just in case
set.seed(5)

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

  plot_cells_on_map <- function(x, color_scale){
    ggplot() + blank_theme + coord_fixed() +
    geom_sf(data=region_of_interest, fill=NA, color="black") +
    geom_polygon(data=grid2, aes(x=long, y=lat, group=group, fill=eval(parse(text=x))), alpha=0.7)    +
    geom_path(data=grid2, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
    color_scale +
    labs(fill=x) +
    xlim(params$plotting_xlim)
    }
  
  # define color scales
  scale_viridis <- viridis::scale_fill_viridis(limits = c(params$n_year_min, length(params$years) - 1))
  scale_blue_red <- scale_fill_gradientn(colours = cols_bd2, na.value=NA, limits = c(0, 1))
  fl <- max(abs(grid2$effect_size_log), na.rm = T) + .1
  
  # plot number of years with a productivity index
  plot_cells_on_map("n_prod", color_scale=scale_viridis)
  # plot number of years with a survival index
  plot_cells_on_map("n_surv", color_scale=scale_viridis)
  # plot probability that survival variance is higher than productivity
  plot_cells_on_map("p_survival_variance_higher", color_scale=scale_blue_red)
  # plot the difference in survival and recruitment variance (log-effect size)
  # scale opacity by the probability that the survival variance is higher.
  ggplot() + coord_fixed() + blank_theme +
    geom_sf(data=region_of_interest, fill=NA, color="black") +
    geom_polygon(data=grid2, aes(x=long, y=lat, group=group, fill = effect_size_log), alpha = 2*abs(grid2$p_survival_variance_higher - 0.5))   +
    geom_path(data=grid2, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
    scale_fill_gradientn(colours = cols_bd, na.value=NA, limits = c(-fl, fl)) +
    xlim(params$plotting_xlim)
}



# compare bayesian/frequentist averages ----
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

# Download or load weather data ----
# Download weather data for all relevant cells, if necessary.

dir.create(paste0(params$output_path, "/weather"), showWarnings = FALSE)
weather_file <- paste0(params$output_path, "/weather/",paste0(params$daymet$label, collapse = "-"),".rds")

if(!params$quiet) print(paste("loading/processing weather file",basename(weather_file)))

# Warning messages below are from `erdPackage::daymet_extract` which calls
# `sp::getSpPPolygonsIDSlots`, which raises the deprecation warning
if(params$always_download_weather | !file.exists(weather_file)){
  # initialize google earth engine
  ee_Initialize()
  
  # loop over cells and years
  data_daymet=data.frame()
  for (cell in cells_all) {
    for (year in params$years){
      print(paste("year = ",year,"cell = ",cell))
      data_daymet <- rbind(data_daymet, daymet_set_extract(year, cell, grid_large, params$daymet))
    }
  }
  saveRDS(data_daymet, weather_file)
} else{
  data_daymet <- readRDS(weather_file)
}

# Perform the regressions of fluctuations (from a single cell, single season, 
# across years) against weather variables.

# do regressions ----
for(species_code in params$species_to_process){
  regression_save_path <- paste0(params$output_path, "/regression_results/", species_code)
  dir.create(regression_save_path, recursive = TRUE, showWarnings = FALSE)
  
  if(params$always_run_regressions | !file.exists(paste0(regression_save_path, "/regressions.rds"))){
    print(paste0("processing_species ", species_code))
    
    file_path <- paste0(params$output_path, "/abun_data/", species_code, "_ratios.rds")
    tidy_ratios <- readRDS(file_path)
    
    data_regression <- weather_regressions(tidy_ratios, data_daymet,params$daymet,params$n_year_min, quiet=TRUE)
    
    saveRDS(data_regression, paste0(regression_save_path, "/regressions.rds"))
  }
}

