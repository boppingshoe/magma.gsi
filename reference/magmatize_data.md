# Preparing MAGMA input data

Preparing MAGMA input data

## Usage

``` r
magmatize_data(
  wd,
  age_classes = "all",
  fishery = NULL,
  loci_names = NULL,
  save_data = TRUE
)
```

## Arguments

- wd:

  Directory where you set up the *data* folder.

- age_classes:

  Hard code class categories for group ages.

- fishery:

  Name of the fishery. It is optional to declare the fishery name here.
  `fishery` can be included with `mixture.RData` in the data folder.

- loci_names:

  Optional. String containing loci names.

- save_data:

  Logical (with default = `TRUE`). Option to save the data in the *data*
  folder.

## Value

A list objects as the input data for `msgsi_mdl()`

## Examples

``` r
if (FALSE) { # \dontrun{
wd <- getwd() # path to data folder
magma_data <- magmatize_data(wd = wd, save_data = FALSE)
} # }
```
