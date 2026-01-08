# magma.gsi 1.0.0

* First working version as an R package.
* Functions to run MAGMA model, summarize results (age-stock composition and individual group membership assignments), and provide convergence diagnostics and trace plots.
* Additions to the non-R-package versions are: hybrid conditional GSI, "zero out" and "weak flat" configurations for age priors, documents on instructions and model descriptions, and saving posterior output in Fst format.

# magma.gsi 1.0.1
## bug fix

* An error in the code that caused mismatching of hatcheries and their proportions. The error was only applicable to analyses with hatchery groups. The error was resolved.
