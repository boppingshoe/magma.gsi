# Plot MCMC trace

Plot MCMC trace

## Usage

``` r
magmatize_tr_plot(obj, nburn = 0, thin = 1, name_order = NULL)
```

## Arguments

- obj:

  Trace from the model output.

- nburn:

  Number of burn-in you set up when you ran the model. Default is 0 if
  you didn't save the burn-ins (keep_burn = FALSE).

- thin:

  Number of thinning you set up when you ran the model. Default is 1 (no
  thinning).

- name_order:

  Arrange the reporting groups as you wish. Leave it empty if you want
  to accept the default.

## Value

Trace plot in ggplot

## Examples

``` r
if (FALSE) { # \dontrun{
# format data
wd <- getwd() # path to data folder
magma_data <- magmatize_data(wd = wd, save_data = FALSE)

# model run
magma_out <- magmatize_mdl(magma_data,
  nreps = 50, nburn = 25, thin = 1, nchains = 2)

# summary
magma_summ <- magmatize_summ(which_dist = 1,
  ma_out = magma_out,
  ma_dat = magma_data,
  summ_level = "district",
  type = "pop")

# trace plot
magmatize_tr_plot(obj = magma_summ$pop_prop[[1]])
} # }
```
