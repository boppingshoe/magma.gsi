
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
#> Time difference of 8.219985 secs
```

Run the model:

``` r
magma_out <- magmatize_mdl(magma_data, nreps = 50, nburn = 25, thin = 1, nchains = 3)
#> Running model (and the category is... Best Mother!)
#> Time difference of 1.908229 secs
#> 2023-08-15 12:25:00.457544
```

Summarize the results:

``` r
magma_summ <- magmatize_summ(ma_out = magma_out, ma_dat = magma_data, summ_level = "district", type = "age")
#> Preparing output (patience grasshopper...)
#> Time difference of 0.5254169 secs
#> 2023-08-15 12:25:01.033248

magma_summ$age_summ[1]
#> $`(1)_D1_Koyukuk`
#> # A tibble: 10 × 10
#>    group   age     mean  median     sd       ci.05  ci.95     p0    GR n_eff
#>    <chr>   <chr>  <dbl>   <dbl>  <dbl>       <dbl>  <dbl>  <dbl> <dbl> <dbl>
#>  1 Koyukuk 11    0.0152 0.00277 0.0257 0.0000181   0.0548 0.16    1.04  75  
#>  2 Koyukuk 12    0.0168 0.00299 0.0323 0.000000834 0.0835 0.187   1.02  49.8
#>  3 Koyukuk 13    0.0130 0.00444 0.0177 0.00000635  0.0451 0.2     1.01  75  
#>  4 Koyukuk 21    0.0139 0.00333 0.0241 0.0000255   0.0643 0.12    1.05  60.8
#>  5 Koyukuk 22    0.0545 0.0398  0.0512 0.00781     0.131  0       1.08  75  
#>  6 Koyukuk 23    0.506  0.502   0.115  0.328       0.693  0       1.06  75  
#>  7 Koyukuk 31    0.318  0.303   0.102  0.174       0.494  0       1.13  75  
#>  8 Koyukuk 32    0.0262 0.00941 0.0403 0.000180    0.0976 0.0533  1.02  75  
#>  9 Koyukuk 33    0.0196 0.00494 0.0350 0.0000352   0.0610 0.2     1.03  75  
#> 10 Koyukuk all   0.0174 0.00856 0.0252 0.0000105   0.0774 0.147   1.01  75
```

There’s a function in the package to make trace plots and inspect mixing
of chains.

``` r
tr_plot(magma_summ$age_prop[[1]])
```

<img src="man/figures/README-example_trace_plot-1.png" width="100%" />

To see more examples on using *magma.gsi*, you can call the manual using
`vignette("magma_vignette")` after you installed *magma.gsi*.
