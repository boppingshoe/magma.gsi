% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/magma_summ_dt.R
\name{magmatize_summ}
\alias{magmatize_summ}
\title{Summarize model output}
\usage{
magmatize_summ(ma_out, ma_dat, summ_level, which_dist = NULL, type = NULL)
}
\arguments{
\item{ma_out}{MAGMA output}

\item{ma_dat}{MAGMA input data}

\item{summ_level}{Summarize at district or subdistrict level}

\item{which_dist}{Function format raw magma output one district at a time.
Identify district as 1, 2, ... Default = NULL will summarize all districts.}

\item{type}{Identify "pop" or "age" to summarize only populations or age class.
if you don't specify a "type", it will summarize both pop and age at the same time.}
}
\value{
Summary tables for reporting groups and/or age classes.
}
\description{
Use this function to summarize "smaller" output (i.e., not TBR).
}
\examples{
\dontrun{
# format data
wd <- getwd() # path to data folder
magma_data <- magmatize_data(wd = wd, save_data = FALSE)

# model run
magma_out <- magmatize_mdl(magma_data, nreps = 50, nburn = 25, thin = 1, nchains = 2)

# summary
magma_summ <- magmatize_summ(ma_out = magma_out, ma_dat = magma_data, summ_level = "district", which_dist = 1, type = "pop")
}

}
