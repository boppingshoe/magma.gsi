# Run MAGMA model

Run MAGMA model

## Usage

``` r
magmatize_mdl(
  dat_in,
  nreps,
  nburn,
  thin,
  nchains,
  nadapt = 50,
  keep_burn = FALSE,
  age_prior = "flat",
  cond_gsi = TRUE,
  file = NULL,
  seed = NULL,
  iden_output = FALSE
)
```

## Arguments

- dat_in:

  Input data list.

- nreps:

  Amount of simulations to run.

- nburn:

  Number of burn-in simulations to discard.

- thin:

  At what interval to keep the simulations.

- nchains:

  Number of independent MCMC chains run in the simulation.

- nadapt:

  Amount of warm-up/adapt runs before the simulation (only for fully
  Bayesian mode).

- keep_burn:

  Logical (default = `FALSE`). To keep the burn-ins in the output or
  not.

- age_prior:

  Option to adjust prior weight on the age class proportions.

  - `flat` (default): conventional setup that puts a flat prior on the
    age proportions.

  - `zero_out`: prior weight will concentrate on the major age groups
    that are observed in metadata (i.e., "zero out" the unobserved age
    classes).

  - `weak_flat`: a less influential flat prior than the conventional one

- cond_gsi:

  Logical (default = `TRUE`). Option to use conditional GSI model. See
  vignette for details.

  - `TRUE`: run MAGMA with a hybrid algorithm of conditional GSI and
    fully Bayesian.

  - `FALSE`: run MAGMA with fully Bayesian algorithm.

- file:

  File path for saving the output. The default is `NULL` for not saving
  the output.

- seed:

  Option to initialize a pseudo-random number generator (set random
  seed) so the output can be reproduced exactly. Just pick a seed number
  and make note of it for future reference. The default is `NULL`.

- iden_output:

  Option to have trace history for individual group membership
  assignments included in the final output. Default is FALSE.

## Value

A list object contains:

- the raw output of MAGMA as a list/multi-way array that need to be
  further summarized using summary functions,

- specifications for the model run (information needed for summary),

- and individual group membership assignment history (optional).

## Examples

``` r
if (FALSE) { # \dontrun{
# format data
wd <- getwd() # path to data folder
magma_data <- magmatize_data(wd = wd, save_data = FALSE)

# model run
magma_out <- magmatize_mdl(magma_data, nreps = 50, nburn = 25, thin = 1, nchains = 2)
} # }
```
