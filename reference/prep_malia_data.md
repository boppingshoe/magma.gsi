# Prepare malia data

Take a magma data object and save the items as text files in a
designated directory.

## Usage

``` r
prep_malia_data(magma_data, path)
```

## Arguments

- magma_data:

  Magma data object made using
  [`magmatize_data()`](https://boppingshoe.github.io/magma.gsi/reference/magmatize_data.md)
  function.

- path:

  A designated directory where you set up for malia data set.

## Value

Items in magma data object saved as text files

## Examples

``` r
if (FALSE) { # \dontrun{
prep_malia_data(magma_data, "D:/bobby_adfg/projects/magma/malia/data/bb2022")
} # }
```
