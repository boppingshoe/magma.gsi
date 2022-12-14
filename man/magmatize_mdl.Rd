% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/magma_mdl.R
\name{magmatize_mdl}
\alias{magmatize_mdl}
\title{Run MAGMA model}
\usage{
magmatize_mdl(
  dat_in,
  nreps,
  nburn,
  thin,
  nchains,
  nadapt = 50,
  keep_burn = FALSE,
  age_priors = "weak",
  cond_gsi = TRUE,
  out_path = NULL,
  seed = NULL
)
}
\arguments{
\item{dat_in}{Input data list.}

\item{nreps}{Amount of simulations to run.}

\item{nburn}{Number of burn-in simulations to discard.}

\item{thin}{Interval to keep the simulations.}

\item{nchains}{Number of independent MCMC chains run in the simulation.}

\item{nadapt}{Amount of warm-up/adapt runs before the simulation.}

\item{keep_burn}{Keep the burn-ins in the output or not. The default is \code{FALSE}.}

\item{cond_gsi}{Option to use conditional GSI model. The default is \code{TRUE}.}

\item{out_path}{File path for saving the output.
The default is \code{NULL} for not saving the output.}

\item{seed}{Option to initialize a pseudo-random number generator (set random seed)
so the output can be reproduced exactly.
Just pick a seed number and make note of it for future reference.
The default is \code{NULL}.}

\item{age_prior}{Specify the level of prior influence for age.
The default is "weak", which has sum of age priors = 1.
If specify \code{age_prior = "strong"}, it will have sum of age priors
= number of specified age classes.}
}
\value{
The raw output of MAGMA is a multi-way array that need to be summarized
using summary functions.
}
\description{
Run MAGMA model
}
\examples{
# format data
wd <- "D:/bobby_adfg/backup_013122/projects/magma/test_TBR" # path to data folder
magma_data <- magmatize_data(wd = wd, save_data = FALSE)

# model run
magma_out <- magmatize_mdl(magma_data, nreps = 50, nburn = 25, thin = 1, nchains = 3)

}
