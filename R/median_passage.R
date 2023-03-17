#' Compute median passage dates for a species
#' @param sp_data species data returned from import_from_erd
#' @param season the season
#' @param min_day minimum day of year
#' @param max_day maximum day of year
#' @param lat1 southern boundary of ROI
#' @param lat2 northern boundary of ROI
#' @param bandwidth width of latitudinal bands to get median, in degrees
#' @return dataframe
#' @export
median_passage <- function(sp_data, season, min_day, max_day, lat1, lat2, bandwidth = 1) {
  years <- as.numeric(unique(sp_data$year))
  years <- years[order(years)]
  
  out <- data.frame(year = years, 
                    south = NA, north = NA,
                    south_effort = NA, north_effort = NA
                    )
  for (i in seq_along(years)) {
    ydata <- sp_data[sp_data$year == years[i] & day_of_year >= min_day & day_of_year <= max_day]
    ydata$reported <- as.numeric((ydata$obs_count + ydata$only_presence_reported) > 0)
    out$south[i] <- median(ydata$day_of_year[ydata$latitude > lat1 &
                                                ydata$latitude < (lat1 + bandwidth) &
                                                ydata$reported == 1])
    out$south_effort[i] <- median(ydata$day_of_year[ydata$latitude > lat1 &
                                                ydata$latitude < (lat1 + bandwidth)])
    
    out$north[i] <- median(ydata$day_of_year[ydata$latitude < lat2 &
                                                ydata$latitude > (lat2 - bandwidth) &
                                                ydata$reported == 1])
    out$north_effort[i] <- median(ydata$day_of_year[ydata$latitude < lat2 &
                                                       ydata$latitude > (lat2 - bandwidth)])
    if (season == "spring") {
      out$duration <- out$north - out$south
    } else if (season == "fall") {
      out$duration <- out$south - out$north
    }
  }
  
  return(out)
}



