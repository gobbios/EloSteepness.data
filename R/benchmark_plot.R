#' benchmark plot
#'
#' @param pdata \code{data("removal_global")}
#' @param empirical logical, should empirical data be plotted (if \code{TRUE}),
#'        or artificial data (if \code{FALSE}). Default is \code{NULL}, which
#'        combines both sets.
#' @param prunkcut numeric, the cutoff for unknown relationships. Default
#'        is \code{0.05}.
#' @param benchmark character, either \code{"pij"} or \code{"input"}. With
#'        \code{"pij"} steepness based on Pij David's score is used as ground
#'        truth. With \code{"input"}, we use the steepness value that was used
#'        during interaction data generation as ground truth, and this is only
#'        applicable when \code{empirical = FALSE}.
#'
#' @importFrom stats lm coef cor
#' @return a plot and invisibly the results of correlation and
#'         simple regression.
#' @export
#'
#' @examples
#' data("removal_global")
#' benchmark_plot(removal_global, empirical = TRUE, benchmark = "pij")
#' x <- benchmark_plot(removal_global, empirical = FALSE, benchmark = "input")
#' x

benchmark_plot <- function(pdata,
                           empirical = NULL,
                           prunkcut = 0.05,
                           benchmark = c("pij", "input")) {
  # labels
  xlabs <- list(expression(DS[italic(P)[ij]]),
                expression(DS[italic(D)[ij]]),
                expression(DS[Bayes]),
                expression(Elo[rpt]),
                expression(Elo[paste("Bayes")]),
                expression(paste("upward steepness"))
  )

  names(xlabs) <- c("ds_pij", "ds_dij", "ds_bayes",
                    "elo_rpt", "elo_bayes_fixed", "upward")

  # select subset
  if (!is.null(empirical)) {
    if (isTRUE(empirical)) {
      pdata <- pdata[pdata$empirical, ]
    }
    if (isFALSE(empirical)) {
      pdata <- pdata[!pdata$empirical, ]
    }
  }

  if (benchmark == "pij") {
    pdata$bench_steep <- pdata$ref_steep
  }
  if (benchmark == "input") {
    pdata$bench_steep <- pdata$input_steep
  }

  # subset
  pdata <- pdata[pdata$step == "step00", ]
  pdata <- pdata[pdata$prunk <= prunkcut, ]
  pdata <- pdata[pdata$method != "elo_bayes_ori", ]

  pdata$method <- factor(as.character(pdata$method),
                         levels = c("ds_pij", "ds_dij",
                                    "ds_bayes", "elo_rpt",
                                    "elo_bayes_fixed", "upward"))

  xgrps <- levels(pdata$method)

  par(mfrow = c(2, 3), mar = c(3, 4, 3, 1))

  outres <- data.frame(method = xgrps, n = NA, rho = NA, beta = NA)

  i="ds_dij"
  for (i in xgrps) {
    if (benchmark == "pij" & i == "ds_pij") {
      plot(0, 0, axes = FALSE, ann = FALSE, type = "n")
      next()
    }
    aux <- pdata[pdata$method == i, ]
    plot(aux$bench_steep, aux$steep_val,
         xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
         las = 1, xlab = "", ylab = "observed steepness",
         lwd = 0.5, xpd = TRUE, pch = 21, bg = grey(0.2, 0.2))
    abline(0, 1)

    r <- lm(steep_val ~ bench_steep, aux)
    abline(r, lty = 2)
    mt <- sprintf("%.2f", cor(aux$steep_val, aux$bench_steep))
    outres$rho[outres$method == i] <- cor(aux$steep_val, aux$bench_steep)
    text(0.88, 0.05, bquote(rho==.(mt)))

    mt <- sprintf("%.2f", coef(r)[2])
    text(0.88, 0.12, bquote(beta==.(mt)))
    outres$beta[outres$method == i] <- coef(r)[2]

    outres$n[outres$method == i] <- nrow(aux)

    title(main = xlabs[[i]], line = 1)
    title(xlab = "ground truth", line = 2)
  }

  invisible(outres)
}
