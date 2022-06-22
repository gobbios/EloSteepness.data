
<!-- badges: start -->
[![R-CMD-check](https://github.com/gobbios/EloSteepness.data/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gobbios/EloSteepness.data/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This package contains the data, code and results to replicate the analyses for [*Extending Bayesian Elo-rating to quantify dominance hierarchy steepness*](https://doi.org/10.1101/2022.01.28.478016).

To install `EloSteepness.data` you require the `devtools` package (which you can get via `install.packages("devtools")`). 
You also require `EloRating` (`install.packages("EloRating")`)

At this point you need to decide whether you want to *install* the vignettes (which include the relevant code to replicate the analyses).
You could proceed without installing the vignettes and download the relevant documents [here](https://github.com/gobbios/EloSteepness.data/blob/main/vignettes) instead (which will save time during installation of the package).

Either strategy installs the package and allows access to the data.

## without vignettes

Use:

```
devtools::install_github("gobbios/EloSteepness.data", build_vignettes = FALSE)
```

That should do the trick. 
But note that in order to actually *run* the code in the vignettes, you need to install the packages listed below anyway.

## with vignettes

For this to work, it's probably a good idea to install `cmdstanr` first ([more info here](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)).

```
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan(cores = 2)
```

And also, `EloSteepness` itself should probably be installed next ([see here for instructions](https://github.com/gobbios/EloSteepness)).

```
devtools::install_github("gobbios/EloSteepness", build_vignettes = FALSE, dependencies = TRUE)
```

With this done, use:

```
devtools::install_github("gobbios/EloSteepness.data", build_vignettes = TRUE, dependencies = TRUE)
```


# Replication of analyses

In the vignette `method_evaluation` (if you installed the package this document should be available directly from R: `vignette("method_evaluation", package = "EloSteepness.data")`) you can find the code to replicate the method evaluation ([or download the file here](https://github.com/gobbios/EloSteepness.data/blob/main/documents/method_evaluation.pdf)).

In the vignette `vignette("phylogenetic_examples", package = "EloSteepness.data")` you can find the code to replicate the two examples using steepness in a comparative analysis ([or download the file here](https://github.com/gobbios/EloSteepness.data/blob/main/documents/phylogenetic_examples.pdf)).

Both vignettes use small subsets of the data for illustrative purposes, but also allow replication of the full analyses.

# References

References for the empirical data sets used are in the document `empirical_references` (`vignette("empirical_references", package = "EloSteepness.data")`, or [download the file here](https://github.com/gobbios/EloSteepness.data/blob/main/documents/empirical_references.pdf)).
