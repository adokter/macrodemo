% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_bcr_abun.R
\name{sample_bcr_abun}
\alias{sample_bcr_abun}
\title{sample a species given thresholds on effort, space and time}
\usage{
sample_bcr_abun(
  species_code,
  erd_path,
  checklists,
  effort_thresholds,
  extent_space,
  extent_time,
  time_window = "full",
  small_grid = 11,
  bcrs,
  time_grid = 7,
  .cores = 4,
  quiet = TRUE
)
}
\arguments{
\item{species_code}{six letter species code}

\item{erd_path}{path to erd}

\item{checklists}{filtered checklists}

\item{effort_thresholds}{effort thresholds}

\item{extent_space}{spatial extent}

\item{extent_time}{temporal extent}

\item{time_window}{stratify by time grid or not}

\item{small_grid}{resolution of small grid}

\item{bcrs}{sf object holding bcrs}

\item{time_grid}{resolution of time}

\item{.cores}{cores for parallel computation}

\item{quiet}{if TRUE, suppress informational messages (warnings/errors still returned)}

\item{large_grid}{resolution of large grid}
}
\value{
large list object
}
\description{
sample a species given thresholds on effort, space and time
}
