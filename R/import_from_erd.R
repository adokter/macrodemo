#' Import and zero-fill data from the ebird reference dataset
#' @param sp_code A six-letter species code from eBird (e.g. "tenwar")
#' @param erd_path The path to the eBird reference dataset .db file
#' @param checklists An optional path to a previously extracted checklists file
#' @param obs An optional path to a previously extracted observation file
#' @return A \code{data.table} object with zero-filled abundance and presence-only reported colums
#'    as well as additional columns for latitude, longitude, year, month, day_of_year, hours_of_day,
#'    protocol_id, is_stationary, is_traveling, effort_hrs, effort_distance_km, cci
#' @export

import_from_erd <- function(sp_code, erd_path = "/Users/jacobsocolar/Dropbox/Work/macrodemography/erd/erd.db", 
                            checklists = NULL, obs = NULL) {
  db <- DBI::dbConnect(RSQLite::SQLite(), erd_path)
  if (is.null(checklists)) {
    message("querying checklists")
    checklists <- import_checklists(erd_path = erd_path)
  }
  
  if (is.null(obs)) {
    message("querying observations")
    obs_query <- DBI::dbSendQuery(db, 
                             paste0("SELECT * FROM obs
      WHERE species_code='", sp_code, "'"))
    obs <- DBI::dbFetch(obs_query)
  }
  DBI::dbClearResult(obs_query)
  DBI::dbDisconnect(db)
  
  message("converting to data table and joining")
  checklists_dt <- data.table::data.table(checklists)
  obs_dt <- data.table::data.table(obs)
  zf <- data.table::merge.data.table(checklists_dt, obs_dt, by = "sampling_event_id", all.x = T)
  
  zf$obs_count[is.na(zf$obs_count)] <- 0
  zf$only_presence_reported[is.na(zf$only_presence_reported)] <- 0
  
  attr(zf, "species") <- sp_code
  
  return(zf)
}

#' Import checklists data from the erd
#' @param erd_path The path to the eBird reference dataset .db file
#' @return a data.frame of the checklists
#' @export
import_checklists <- function(erd_path = "/Users/jacobsocolar/Dropbox/Work/macrodemography/erd/erd.db") {
  db <- DBI::dbConnect(RSQLite::SQLite(), erd_path)
  checklists_query <- DBI::dbSendQuery(db, 
                                       "SELECT sampling_event_id, latitude, longitude, year, month, day_of_year,
                                        hours_of_day, protocol_id, is_stationary, is_traveling, 
                                        effort_hrs, effort_distance_km, cci, ntl_mean, ELEV_30M_MEDIAN
                                     FROM checklists")
  checklists <- DBI::dbFetch(checklists_query)
  DBI::dbClearResult(checklists_query)
  DBI::dbDisconnect(db)
  return(checklists)
}
