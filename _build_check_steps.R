# this is a set of instructions to check whether the package conforms to CRAN rules
# these steps here are taken from https://github.com/paul-buerkner/brms/issues/158
# stricter than necessary because this package will never go to CRAN (because of its size and its cmdstanr dependency)

devtools::clean_vignettes()
devtools::document()
devtools::clean_dll()


withr::with_envvar(c("NOT_CRAN" = "true"),
                   devtools::build_vignettes()
)

# size reductions of PDF vignettes
tools::compactPDF("doc/empirical_references.pdf", gs_quality = "ebook")
tools::compactPDF("doc/method_evaluation.pdf", gs_quality = "ebook")
tools::compactPDF("doc/phylogenetic_examples.pdf", gs_quality = "ebook")
tools::compactPDF("doc/additional_figures.pdf", gs_quality = "ebook")


withr::with_envvar(c("NOT_CRAN" = "true"),
                   devtools::build(args = c('--compact-vignettes=both'))
)


withr::with_envvar(c("NOT_CRAN" = "true"),
                   devtools::check_built("../EloSteepness.data_0.9.6.tar.gz", args = "--as-cran")
)

