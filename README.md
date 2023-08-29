
<!-- README.md is generated from README.Rmd. Please edit that file -->

# magma.gsi <img src="man/figures/logo.png" align="right" height="97" alt="" />

<!-- badges: start -->
<!-- badges: end -->

*MAGMA* is a slowing moving but powerful model that combine mark and age
information with genetic mixture analysis. Mainly, *MAGMA* estimates two
sets of parameters: stock and age. A composition of the two, information
that is often required for run reconstruction models, is simply the
product of these two sets of parameters. After all age/stock
compositions in other strata are calculated in the same fashion, they
are multiplied by the harvest proportions of their corresponding strata
and summed up to get a weighted average age/stock composition.

## Installation

You can install the development version of *magma.gsi* from
[GitHub](https://github.com/boppingshoe/magma.gsi) with:

``` r
# install.packages("devtools")
devtools::install_github("boppingshoe/magma.gsi", build_vignettes = TRUE)
```

## Example

Once you have the data sets up at a designated directory, you can
compile the input object:

``` r
library(magma.gsi)

wd <- getwd() # path to data folder
magma_data <- magmatize_data(wd = paste0(wd, "/vignettes"), save_data = FALSE)
#> Compiling input data, may take a minute or two...
#> No missing hatcheries
#> Time difference of 8.80534 secs
```

Run the model:

``` r
magma_out <- magmatize_mdl(magma_data, nreps = 50, nburn = 25, thin = 1, nchains = 3)
#> Running model (and the category is... Live, Werk, Pose!)
#> Time difference of 1.923433 secs
#> 2023-08-28 17:07:52.779157
```

Summarize the results:

``` r
magma_summ <- magmatize_summ(ma_out = magma_out, ma_dat = magma_data, summ_level = "district", type = "age")
#> Preparing output (patience grasshopper...)
#> Time difference of 0.5067289 secs
#> 2023-08-28 17:07:53.325211

magma_summ$age_summ[1]
#> $`(1)_D1_Koyukuk`
#> # A tibble: 9 × 10
#>   group   age     mean  median     sd      ci.05  ci.95    p0    GR n_eff
#>   <chr>   <chr>  <dbl>   <dbl>  <dbl>      <dbl>  <dbl> <dbl> <dbl> <dbl>
#> 1 Koyukuk 11    0.0212 0.00672 0.0428 0.0000656  0.0757 0.08   1.03  75  
#> 2 Koyukuk 12    0.0158 0.00458 0.0261 0.00000451 0.0547 0.187  1.08  68.4
#> 3 Koyukuk 13    0.0200 0.00722 0.0361 0.0000141  0.0683 0.133  1.02  73.6
#> 4 Koyukuk 21    0.0147 0.00508 0.0235 0.0000249  0.0480 0.12   1.05  75  
#> 5 Koyukuk 22    0.0657 0.0542  0.0550 0.00803    0.186  0      1.11  75  
#> 6 Koyukuk 23    0.512  0.545   0.118  0.294      0.666  0      1.15  97.0
#> 7 Koyukuk 31    0.323  0.309   0.115  0.162      0.530  0      1.12  75  
#> 8 Koyukuk 32    0.0141 0.00363 0.0210 0.0000108  0.0573 0.173  1.07  75  
#> 9 Koyukuk 33    0.0136 0.00337 0.0197 0.0000391  0.0572 0.133  1.06  75
```

There’s a function in the package to make trace plots and inspect mixing
of chains.

``` r
tr_plot(magma_summ$age_prop[[1]])
```

<img src="man/figures/README-example_trace_plot-1.png" width="100%" />

To see more examples on using *magma.gsi*, you can call the manual using
`vignette("magma-vignette", package = "magma.gsi")` after you installed
*magma.gsi*.
