% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/effort_weather_correction.R
\name{effort_weather_correction}
\alias{effort_weather_correction}
\title{A helper function to generate effort (detectability) corrected counts for each checklist}
\usage{
effort_weather_correction(sp_data)
}
\arguments{
\item{sp_data}{zero-filled 'sp_data' data.table object, output from \code{import_from_erd} function}
}
\value{
A \code{data.table} object with zero-filled abundance (effort corrected) and presence-only reported colums
as well as additional columns for latitude, longitude, year, month, day_of_year, hours_of_day,
protocol_id, is_stationary, is_traveling, effort_hrs, effort_distance_km, cci, local weather
}
\description{
This function runs on zero-filled data and runs XGboost to predict count for
original conditions (lat,long, day of year, checlist duration, distance, hour of day, number of observer, mode of checklist,
cci, temperature, precipitation, wind), as well for standarndized condition.
}
