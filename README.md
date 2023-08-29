
<!-- README.md is generated from README.Rmd. Please edit that file -->

# magma.gsi <img src="man/figures/logo.png" align="right" height="200" alt="" />

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
#> Time difference of 8.675583 secs
```

Run the model:

``` r
magma_out <- magmatize_mdl(magma_data, nreps = 50, nburn = 25, thin = 1, nchains = 3)
#> Running model (and the category is... Live, Werk, Pose!)
#> Time difference of 2.01657 secs
#> 2023-08-28 17:15:41.888492
```

Summarize the results:

``` r
magma_summ <- magmatize_summ(ma_out = magma_out, ma_dat = magma_data, summ_level = "district", type = "age")
#> Preparing output (patience grasshopper...)
#> Time difference of 0.5365539 secs
#> 2023-08-28 17:15:42.465052

magma_summ$age_summ[1]
#> $`(1)_D1_Koyukuk`
#> # A tibble: 9 × 10
#>   group   age     mean  median     sd     ci.05  ci.95     p0    GR n_eff
#>   <chr>   <chr>  <dbl>   <dbl>  <dbl>     <dbl>  <dbl>  <dbl> <dbl> <dbl>
#> 1 Koyukuk 11    0.0190 0.00882 0.0257 0.000163  0.0780 0.0533 1.07   75  
#> 2 Koyukuk 12    0.0216 0.00609 0.0376 0.0000110 0.0873 0.12   0.999  75  
#> 3 Koyukuk 13    0.0158 0.00369 0.0329 0.0000149 0.0797 0.147  1.02   75  
#> 4 Koyukuk 21    0.0223 0.00773 0.0364 0.0000600 0.0995 0.0667 1.03   70.8
#> 5 Koyukuk 22    0.0541 0.0350  0.0469 0.00467   0.139  0      0.989  75  
#> 6 Koyukuk 23    0.516  0.518   0.115  0.343     0.718  0      1.02   62.8
#> 7 Koyukuk 31    0.318  0.327   0.100  0.179     0.496  0      1.06   75  
#> 8 Koyukuk 32    0.0178 0.00577 0.0309 0.0000184 0.0764 0.173  1.06   75  
#> 9 Koyukuk 33    0.0156 0.00372 0.0265 0.0000321 0.0647 0.12   1.24   76.5
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
