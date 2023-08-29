
<!-- README.md is generated from README.Rmd. Please edit that file -->

# magma.gsi <img src="man/figures/logo.png" align="center" height="140" alt="" />

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
#> Time difference of 8.702406 secs
```

Run the model:

``` r
magma_out <- magmatize_mdl(magma_data, nreps = 50, nburn = 25, thin = 1, nchains = 3)
#> Running model (and the category is... Working Girl!)
#> Time difference of 1.995209 secs
#> 2023-08-28 17:10:44.186345
```

Summarize the results:

``` r
magma_summ <- magmatize_summ(ma_out = magma_out, ma_dat = magma_data, summ_level = "district", type = "age")
#> Preparing output (patience grasshopper...)
#> Time difference of 0.5967751 secs
#> 2023-08-28 17:10:44.827379

magma_summ$age_summ[1]
#> $`(1)_D1_Koyukuk`
#> # A tibble: 9 × 10
#>   group   age     mean  median     sd       ci.05  ci.95     p0    GR n_eff
#>   <chr>   <chr>  <dbl>   <dbl>  <dbl>       <dbl>  <dbl>  <dbl> <dbl> <dbl>
#> 1 Koyukuk 11    0.0151 0.00509 0.0269 0.000107    0.0717 0.12    1.02  77.5
#> 2 Koyukuk 12    0.0176 0.00522 0.0339 0.0000499   0.0871 0.08    1.00  75  
#> 3 Koyukuk 13    0.0134 0.00280 0.0242 0.000000428 0.0663 0.16    1.01  52.7
#> 4 Koyukuk 21    0.0152 0.00323 0.0315 0.0000138   0.0684 0.187   1.04  62.6
#> 5 Koyukuk 22    0.0599 0.0483  0.0487 0.0100      0.167  0       1.00  75  
#> 6 Koyukuk 23    0.536  0.541   0.107  0.355       0.706  0       1.02  81.4
#> 7 Koyukuk 31    0.303  0.296   0.103  0.152       0.492  0       1.06  75  
#> 8 Koyukuk 32    0.0224 0.00530 0.0383 0.0000941   0.0983 0.0933  1.04  75  
#> 9 Koyukuk 33    0.0174 0.00341 0.0259 0.0000143   0.0770 0.187   1.02  68.7
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
