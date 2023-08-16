
<!-- README.md is generated from README.Rmd. Please edit that file -->

# magma.gsi

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
magma_data <- magmatize_data(wd = wd, save_data = FALSE)
#> Compiling input data, may take a minute or two...
#> No missing hatcheries
#> Time difference of 8.260864 secs
```

Run the model:

``` r
magma_out <- magmatize_mdl(magma_data, nreps = 50, nburn = 25, thin = 1, nchains = 3)
#> Running model (and the category is... Face!)
#> Time difference of 2.005737 secs
#> 2023-08-16 13:39:34.492473
```

Summarize the results:

``` r
magma_summ <- magmatize_summ(ma_out = magma_out, ma_dat = magma_data, summ_level = "district", type = "age")
#> Preparing output (patience grasshopper...)
#> Time difference of 0.4979138 secs
#> 2023-08-16 13:39:35.031269

magma_summ$age_summ[1]
#> $`(1)_D1_Koyukuk`
#> # A tibble: 9 × 10
#>   group   age     mean  median     sd       ci.05  ci.95     p0    GR n_eff
#>   <chr>   <chr>  <dbl>   <dbl>  <dbl>       <dbl>  <dbl>  <dbl> <dbl> <dbl>
#> 1 Koyukuk 11    0.0207 0.00400 0.0359 0.000000819 0.109  0.24   1.04   75  
#> 2 Koyukuk 12    0.0145 0.00434 0.0255 0.00000272  0.0740 0.227  1.04   75  
#> 3 Koyukuk 13    0.0151 0.00390 0.0249 0.0000207   0.0608 0.12   1.05   61.7
#> 4 Koyukuk 21    0.0152 0.00393 0.0287 0.0000100   0.0901 0.173  1.06   62.3
#> 5 Koyukuk 22    0.0324 0.0137  0.0428 0.00000740  0.121  0.0933 1.15   63.0
#> 6 Koyukuk 23    0.534  0.533   0.115  0.321       0.699  0      1.36  152. 
#> 7 Koyukuk 31    0.333  0.321   0.125  0.143       0.518  0      1.43   75  
#> 8 Koyukuk 32    0.0154 0.00473 0.0274 0.00000540  0.0789 0.16   0.997  62.6
#> 9 Koyukuk 33    0.0199 0.00496 0.0354 0.0000259   0.100  0.173  1.01   75
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
