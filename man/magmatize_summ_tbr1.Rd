% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/magma_summ_tbr_1.R
\name{magmatize_summ_tbr1}
\alias{magmatize_summ_tbr1}
\title{TBR output format wrapper step 1}
\usage{
magmatize_summ_tbr1(
  which_dist = NULL,
  outraw,
  ma_dat,
  nreps,
  nburn,
  thin,
  nchains,
  keep_burn = FALSE,
  malia = FALSE
)
}
\arguments{
\item{which_dist}{Function format raw magma output one district at a time.
Identify district as 1, 2, ... Default = NULL will summarize all districts.}

\item{outraw}{MAGMA output.}

\item{ma_dat}{The same MAGMA input data for model run.}

\item{nreps}{The same as \emph{nreps} in MAGMA model run.}

\item{nburn}{The same as \emph{nburn} in MAGMA model run.}

\item{thin}{The same as \emph{thin} in MAGMA model run.}

\item{nchains}{The same as \emph{nchains} in MAGMA model run.}

\item{keep_burn}{The same as \emph{keep_burn} in MAGMA model run.}

\item{malia}{= TRUE if the output is from program malia.}

\item{summ_level}{Summarize at district or subdistrict level.}

\item{type}{Identify "pop" or "age" to summarize populations or age class.}
}
\value{
Model output in multiway array and subset of metadata as a list object.
}
\description{
TBR output format wrapper step 1
}
\examples{
# format data
wd <- "D:/bobby_adfg/backup_013122/projects/magma/test_TBR" # path to data folder
magma_data <- magmatize_data(wd = wd, save_data = FALSE)

# model run
magma_out <- magmatize_mdl(magma_data,
  nreps = 50, nburn = 25, thin = 1, nchains = 3)

# summary step 1
tbr1 <- magmatize_summ_tbr1(which_dist = 3,
  outraw = magma_out,
  ma_dat = magma_data,
  nreps = 50, nburn = 25, thin = 1, nchains = 3)

}
