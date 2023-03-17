#' extract Daymet summaries over a spatiotemporal cell of interest
#' @import rgee 
#' @param cell the cell id from a dggrid
#' @param res the resolution of the dggrid
#' @param variable the name of the daymet variable
#' @param mindate the minimum date over which to summarize
#' @param maxdate the maximum date over which to summarize
#' @return the median (over pixels in the spatial cell) of the mean (over days
#'  within the pixels) of the variable of interest.
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
