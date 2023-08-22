
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
magma_data <- magmatize_data(wd = paste0(wd, "/vignettes"), save_data = FALSE)
#> Compiling input data, may take a minute or two...
#> No missing hatcheries
#> Time difference of 8.559501 secs
```

Run the model:

``` r
magma_out <- magmatize_mdl(magma_data, nreps = 50, nburn = 25, thin = 1, nchains = 3)
#> Running model (and the category is... Butch Queen Body!)
#> Time difference of 1.945559 secs
#> 2023-08-22 10:00:31.310867
```

Summarize the results:

``` r
magma_summ <- magmatize_summ(ma_out = magma_out, ma_dat = magma_data, summ_level = "district", type = "age")
#> Preparing output (patience grasshopper...)
#> Time difference of 0.522465 secs
#> 2023-08-22 10:00:31.87553

magma_summ$age_summ[1]
#> $`(1)_D1_Koyukuk`
#> # A tibble: 9 × 10
#>   group   age     mean  median     sd      ci.05  ci.95     p0    GR n_eff
#>   <chr>   <chr>  <dbl>   <dbl>  <dbl>      <dbl>  <dbl>  <dbl> <dbl> <dbl>
#> 1 Koyukuk 11    0.0193 0.00689 0.0394 0.0000269  0.0631 0.12   1.18   75  
#> 2 Koyukuk 12    0.0163 0.00614 0.0259 0.0000305  0.0583 0.147  1.08   75  
#> 3 Koyukuk 13    0.0144 0.00759 0.0179 0.00000981 0.0503 0.12   0.988  75  
#> 4 Koyukuk 21    0.0245 0.00900 0.0404 0.000145   0.0948 0.0533 1.07   75  
#> 5 Koyukuk 22    0.0717 0.0603  0.0595 0.00381    0.177  0      1.11   75  
#> 6 Koyukuk 23    0.498  0.500   0.131  0.307      0.717  0      0.995  75  
#> 7 Koyukuk 31    0.316  0.303   0.132  0.139      0.548  0      1.01   85.0
#> 8 Koyukuk 32    0.0214 0.00427 0.0333 0.0000144  0.101  0.16   1.01   59.9
#> 9 Koyukuk 33    0.0181 0.00610 0.0307 0.0000385  0.0626 0.0667 1.09   75
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
