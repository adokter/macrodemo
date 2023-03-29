#' Extract Daymet summaries over a spatiotemporal cell of interest
#' @param cell the cell id from a dggrid
#' @param res the resolution of the dggrid
#' @param variable the name of the daymet variable
#' @param mindate the minimum date over which to summarize
#' @param maxdate the maximum date over which to summarize
#' @return the median (over pixels in the spatial cell) of the mean (over days
#'  within the pixels) of the variable of interest.
#' @export
daymet_extract <- function(cell = 621, res = 6, variable = "swe", mindate = "2014-01-01", maxdate = "2014-03-01"){
  dg <- dggridR::dgconstruct(res = res)
  global <- dggridR::dgearthgrid(dg, frame=FALSE)
  global.cell <- data.frame(cell=sp::getSpPPolygonsIDSlots(global), row.names=sp::getSpPPolygonsIDSlots(global))
  global <- sf::st_as_sf(SpatialPolygonsDataFrame(global, global.cell))
  polygon = sf_as_ee(global[cell,]) # rgee already imported
  daymet <- ee$
    ImageCollection("NASA/ORNL/DAYMET_V4")$
    filterBounds(polygon)$
    select(variable)$
    filterDate(mindate, maxdate)$
    mean()$
    reduceRegion(
      reducer = ee$Reducer$median(),
      geometry = polygon,
      scale = 1000,
      maxPixels = 1e10
    )
  return(daymet$getInfo())
}

#' Extract a set of Daymet parameters for a cell and year
#'
#' Extracts a set of Daymet data according to a specification of variables, minimum and maximum dates
#' using \link{daymet_extract}
#' @param year the year of interest
#' @param cell the large grid index
#' @param grid a grid specification obtained via \link[dggridR]{dgconstruct}
#' @param params_daymet a data.frame specifying daymet parameters and averaging periods
#' @return a data.frame with cell-averaged daymet parameters
#' @export
#' @details
#' Parameter `params_daymet` should be a `data.frame` with the following columns:
#' \describe{
#' \item{`label`}{a user-defined label to define the entry}
#' \item{`variable`}{a valid Daymet variable}
#' \item{`date_min`}{a minimum date formatted as mm-dd (without year)}
#' \item{`date_max`}{a maximum date formatted as mm-dd (without year)}
#' \item{`period`}{the corresponding demographic period (to associate a regression using \link{weather_regressions})}
#' }
daymet_set_extract <- function(year, cell, grid, params_daymet){
  # initialize the data.frame for daymet data:
  daymet_data=data.frame(matrix(nrow=1,ncol=nrow(params_daymet)+2,dimnames=list(NULL,c("cell","year",params_daymet$label))))
  daymet_data$cell=cell
  daymet_data$year=year
  # fill data.frame with averaged daymet data:
  for(lab in params_daymet$label){
    # filter parameters to use for daymet extraction:
    par_daymet <- params_daymet %>% filter(label==lab)
    # construct min and max dates
    maxdate=paste0(year,"-",par_daymet$date_max)
    mindate=paste0(year,"-",par_daymet$date_min)
    # address cases in which the min to max period includes a year change:
    if(as.Date(maxdate)<as.Date(mindate)) mindate=paste0(year-1,"-",par_daymet$date_min)
    # average over large cell, and store in data.frame
    daymet_data[par_daymet$label] <- daymet_extract(cell=cell, res=grid$res, variable=par_daymet$variable, mindate=mindate, maxdate=maxdate)
  }
  return(daymet_data)
}

