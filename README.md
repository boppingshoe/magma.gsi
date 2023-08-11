
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

You can install the development version of magma.gsi from
[GitHub](https://github.com/boppingshoe/magma.gsi) with:

``` r
# install.packages("devtools")
devtools::install_github("boppingshoe/magma.gsi")
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
#> Time difference of 13.17325 secs
```

Run the model:

``` r
magma_out <- magmatize_mdl(magma_data, nreps = 50, nburn = 25, thin = 1, nchains = 3)
#> Running model (and the category is... Best Dressed!)
#> Time difference of 3.288395 secs
#> 2023-08-11 13:41:45.452511
```

Summarize the results:

``` r
magma_summ <- magmatize_summ(ma_out = magma_out, ma_dat = magma_data, summ_level = "district", type = "age")
#> Preparing output (patience grasshopper...)
#> Time difference of 0.8373759 secs
#> 2023-08-11 13:41:46.35936

magma_summ$age_summ[1]
#> $`(1)_D1_Koyukuk`
#> # A tibble: 9 × 10
#>   group   age     mean  median     sd      ci.05  ci.95     p0    GR n_eff
#>   <chr>   <chr>  <dbl>   <dbl>  <dbl>      <dbl>  <dbl>  <dbl> <dbl> <dbl>
#> 1 Koyukuk 11    0.0208 0.0120  0.0257 0.000140   0.0707 0.0667  1.09  75  
#> 2 Koyukuk 12    0.0174 0.00570 0.0251 0.0000638  0.0641 0.147   1.13  75  
#> 3 Koyukuk 13    0.0176 0.00449 0.0297 0.0000322  0.0923 0.12    1.02  59.3
#> 4 Koyukuk 21    0.0162 0.00423 0.0301 0.0000100  0.0733 0.187   1.01  75  
#> 5 Koyukuk 22    0.0531 0.0463  0.0448 0.00125    0.139  0.0267  1.35  71.0
#> 6 Koyukuk 23    0.554  0.547   0.105  0.394      0.699  0       1.01  56.6
#> 7 Koyukuk 31    0.286  0.273   0.104  0.145      0.490  0       1.01 170. 
#> 8 Koyukuk 32    0.0169 0.00571 0.0248 0.00000178 0.0726 0.12    1.01 153. 
#> 9 Koyukuk 33    0.0182 0.00520 0.0281 0.00000483 0.0832 0.173   1.10  75
```

There’s a function in the package to make trace plots and inspect mixing
of chains.

``` r
tr_plot(magma_summ$age_prop[[1]])
```

<img src="man/figures/README-example_trace_plot-1.png" width="100%" />

To see more examples on using MAGMA, you can call the manual using
`vignette("magma_vignette")` after you installed *magma.gsi*.
