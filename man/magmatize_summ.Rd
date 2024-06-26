% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/magma_summ.R
\name{magmatize_summ}
\alias{magmatize_summ}
\title{Summarize model output}
\usage{
magmatize_summ(
  ma_out = NULL,
  ma_dat,
  summ_level,
  which_dist = NULL,
  fst_files = NULL,
  save_trace = "in_memory"
)
}
\arguments{
\item{ma_out}{MAGMA output}

\item{ma_dat}{MAGMA input data}

\item{summ_level}{Summarize at district or subdistrict level}

\item{which_dist}{Function format raw magma output one district at a time.
Identify district as 1, 2, ... Default = NULL will summarize all districts.}

\item{fst_files}{Fst files location if MAGMA model output was saved.}

\item{save_trace}{default = "in_memory" to have trace history as a part of summary. Or specify the path of a directory to save trace history as fst files.}
}
\value{
Summary tables for reporting groups and/or age classes.
}
\description{
Summarize model output
}
\examples{
\dontrun{
# format data
wd <- getwd() # path to data folder
magma_data <- magmatize_data(wd = wd, save_data = FALSE)

# model run
magma_out <- magmatize_mdl(magma_data, nreps = 50, nburn = 25, thin = 1, nchains = 2)

# summary using output as an object
magma_summ <- magmatize_summ(ma_out = magma_out, ma_dat = magma_data, summ_level = "district", which_dist = 1, save_trace = wd)

# summary using output saved as Fst files
magma_summ <- magmatize_summ(ma_dat = magma_data, summ_level = "district", which_dist = 1, fst_files = wd)
}

}
