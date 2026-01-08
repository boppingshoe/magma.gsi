# Individual group membership assignment summary

Individual group membership assignment summary

## Usage

``` r
magmatize_indiv(ma_out, ma_dat, out_repunit = FALSE)
```

## Arguments

- ma_out:

  MAGMA model output object name.

- ma_dat:

  MAGMA input data object name.

- out_repunit:

  Logical (default = `FALSE`). Option to output the summaries in
  populations (`FALSE`) or combined in reporting groups (`TRUE`).

  - NOTE: combining populations into reporting group recommend only for
    single district or multi-district with the same reporting groups.

## Value

A tibble contains probabilities for reporting group memberships of each
individual of MAGMA metadata.

## Examples

``` r
if (FALSE) { # \dontrun{
wd <- getwd() # path to data folder
magma_data <- magmatize_data(wd = wd, save_data = FALSE)
magma_out <- magmatize_mdl(magma_data, nreps = 50, nburn = 25, thin = 1, nchains = 2)
magma_assn <- magmatize_indiv(magma_out, magma_datt, out_repunit = TRUE)
} # }
```
