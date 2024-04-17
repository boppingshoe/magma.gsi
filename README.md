
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
#> FAKE is the fishery identified in the mixture.RData
#> No missing hatcheries
#> Time difference of 5.321199 secs
```

Run the model:

``` r
magma_out <- magmatize_mdl(magma_data, nreps = 50, nburn = 25, thin = 1, nchains = 3)
#> Running model... and work for a cause, not for Butch Queen Body!
#> Time difference of 1.025131 secs
#> 2024-04-17 11:50:17.629791
```

Summarize the results:

``` r
magma_summ <- magmatize_summ(ma_out = magma_out, ma_dat = magma_data, summ_level = "district")
#> Preparing output (patience grasshopper...)
#> Time difference of 0.3293231 secs
#> 2024-04-17 11:50:17.975011

magma_summ$age_summ[1]
#> $D1_Koyukuk
#> # A tibble: 9 × 10
#>   group   age     mean  median     sd     ci.05  ci.95     p0    GR n_eff
#>   <chr>   <fct>  <dbl>   <dbl>  <dbl>     <dbl>  <dbl>  <dbl> <dbl> <dbl>
#> 1 Koyukuk 11    0.0177 0.00671 0.0222 0.0000357 0.0618 0.12   1.05  166. 
#> 2 Koyukuk 12    0.0238 0.0106  0.0306 0.0000403 0.0893 0.107  1.07   62.8
#> 3 Koyukuk 13    0.0213 0.00402 0.0353 0.0000438 0.0893 0.0933 1.05   75  
#> 4 Koyukuk 21    0.0168 0.00597 0.0267 0.0000979 0.0686 0.0667 1.03   75  
#> 5 Koyukuk 22    0.0656 0.0534  0.0554 0.00614   0.166  0      1.06   61.1
#> 6 Koyukuk 23    0.536  0.546   0.115  0.301     0.698  0      1.00   75  
#> 7 Koyukuk 31    0.282  0.275   0.0994 0.128     0.422  0      1.00   75  
#> 8 Koyukuk 32    0.0125 0.00485 0.0215 0.0000138 0.0412 0.0933 1.13   75  
#> 9 Koyukuk 33    0.0243 0.00701 0.0366 0.0000318 0.0922 0.0933 0.996 145.
```

There’s a function in the package to make trace plots and inspect mixing
of chains.

``` r
magmatize_tr_plot(magma_summ$age_prop[[1]])
```

<img src="man/figures/README-example_trace_plot-1.png" width="100%" />

To see more examples on using *magma.gsi*, you can call the manual using
`vignette("magma-vignette", package = "magma.gsi")` after you installed
*magma.gsi*.
