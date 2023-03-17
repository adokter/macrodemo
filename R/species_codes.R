#' four and six letter bird species codes for focal species
#'
#' @format A data.frame with 22 rows and 2 columns
#' @source this script

species_codes <- 
  data.frame(
    four = c("bcch", "bhnu", "cach", "carw", "noca", "piwo", "tuti",
             "canw", "cacw", "wren", "calt", "cant", "cath", "cbth", "lbth",
             "btgn", "verd", "pyrr", "oati", "juti", "casj", "wosj"
            ),
    six = c("bkcchi", "bnhnut", "carchi", "carwre", "norcar", "pilwoo", 
            "tuftit", "canwre", "cacwre", "wrenti", "caltow", "cantow", 
            "calthr", "cubthr", "lobthr", "bktgna", "verdin", "pyrrhu", 
            "oaktit", "juntit", "cowscj", "wooscj"
            )
    )

usethis::use_data(species_codes, overwrite = TRUE)