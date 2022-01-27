This package contains the data, code and results to replicate the analyses for *Extending Bayesian Elo-rating to quantify dominance hierarchy steepness*.

To install `EloSteepness.data` you require the `devtools` package (which you can get via `install.packages("devtools")`). 
You also require `EloRating` (`install.packages("EloRating")`)
Then use:

```
devtools::install_github("gobbios/EloSteepness.data")
```

More packages are required to replicate all the analyses (but *installation* of `EloSteepness.data`, and hence *access* to the data, should work without them being actually installed).
That means that you do not necessarily have to install all of the packages listed below, but only those that are listed in the vignette that you are interested in.
First you probably need/want `EloSteepness` ([see here for instructions](https://github.com/gobbios/EloSteepness)).
Furthermore, you may need:
  
  - `ape` (via `install.packages("ape")`)
  
  - `rethinking` ([see here](https://github.com/rmcelreath/rethinking))
  
  - `rstan` (via `install.packages("rstan")`)

  - `brms` (via `install.packages("brms")`)
  
  - `cmdstanr` ([see here](https://mc-stan.org/cmdstanr/articles/cmdstanr.html))


# Replication of analyses

In the vignette `vignette("method_evaluation", package = "EloSteepness.data")` you can find the code to replicate the method evaluation ([or download the file here](https://github.com/gobbios/EloSteepness.data/blob/main/vignettes/pdf_files/method_evaluation.pdf)).

In the vignette `vignette("phylogenetic_examples", package = "EloSteepness.data")` you can find the code to replicate the two examples using steepness in a comparative analysis ([or download the file here](https://github.com/gobbios/EloSteepness.data/blob/main/vignettes/pdf_files/phylogenetic_examples.pdf)).

Both vignettes use small subsets of the data for illustrative purposes, but also allow replication of the full analyses.

# References

References for the empirical data sets used can be accessed via `vignette("empirical_references", package = "EloSteepness.data")`.
