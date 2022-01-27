#' result of the removal experiment
#'
#' @aliases removal_slopes
#' @format \code{removal_slopes} is a data frame with 10,020 observations, each
#'         line corresponding to one matrix subjected to one steepness algorithm
#'         multiple times (after removals)
#' @format \code{removal_global} is a data frame with 121,557 observations, each
#'         line corresponding to one matrix at a specific removal step and
#'         steepness algorithm
#' @examples
#' data(removal_global)
#' xdata <- removal_global[removal_global$data_set == "alados1992_dom_6" &
#'                           removal_global$method == "ds_pij", ]
#' xdata
#' coef(lm(steep_val ~ prunk, data = xdata))[2]
#' removal_slopes$slope_prunk[removal_slopes$data_set == "alados1992_dom_6" &
#'                              removal_slopes$method == "ds_pij"]
#'
#' removal_global[removal_global$data_set == "alados1992_dom_5" &
#'                  removal_global$step == "step04", ]
"removal_global"
