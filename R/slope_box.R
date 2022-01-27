#' plot steepness as function of unknown relationships
#'
#' @param pdata \code{data("removal_slopes")}
#' @param empirical logical, should empirical data be plotted (if \code{TRUE}),
#'        or artificial data (if \code{FALSE}). Default is \code{NULL}, which
#'        combines both sets.
#' @param hinges numeric of length 2, determines the quantiles of the boxplot.
#'        Default is \code{c(0.05, 0.95)}.
#' @param ylims numeric of length 2, determines the y-limits.
#'
#' @importFrom stats quantile
#' @importFrom graphics segments rect points text title par abline title axis
#' @importFrom grDevices grey
#'
#' @return a plot, and invisibly the numeric data.
#' @export
#'
#' @examples
#' data("removal_slopes")
#' x <- slope_box(removal_slopes, empirical = TRUE)

slope_box <- function(pdata,
                      empirical = NULL,
                      hinges = c(0.05, 0.95),
                      ylims = c(-1.5, 1)) {
  # axis labels
  xlabs <- c(expression(DS[italic(P)[ij]]),
             expression(DS[italic(D)[ij]]),
             expression(DS[Bayes]),
             expression(Elo[rpt]),
             expression(Elo[paste("Bayes")]),
             expression(paste("upward")))

  # methods
  xgrps <- levels(pdata$method)

  # colors
  xcols <- rep("grey90", length(xgrps))


  # select subset
  if (!is.null(empirical)) {
    if (isTRUE(empirical)) {
      pdata <- pdata[pdata$empirical, ]
    }
    if (isFALSE(empirical)) {
      pdata <- pdata[!pdata$empirical, ]
    }
  }

  pd <- tapply(X = pdata$slope_prunk,
               INDEX = pdata$method,
               quantile,
               probs = sort(c(hinges, 0.25, 0.5, 0.75)),
               na.rm = TRUE)

  xlims <- c(0.5, length(xgrps) + 0.5)
  plot(0, 0, type = "n", xlim = c(0.5, length(xgrps) + 0.5), ylim = ylims,
       axes = FALSE, xlab = "", ylab = "dependence on sparseness")
  abline(h = seq(ylims[1], ylims[2], by = 0.5), col = "grey", lwd = 0.5)
  abline(h = 0, col = "grey", lty = 2, lwd = 1.5)

  for (i in 1:length(xgrps)) {
    l <- i - 0.28
    r <- i + 0.28
    p <- pd[[xgrps[i]]]

    segments(i, p[1], i, p[5], lwd = 1.5)
    rect(l, p[2], r, p[4], col = "white", border = NA)
    rect(l, p[2], r, p[4], col = xcols[i], border = "black", lwd = 0.5)
    segments(l, p[3], r, p[3], lwd = 3.5, lend = 3)
  }

  # text(x = 1:length(xgrps), y = ylims[1] * 1.02, labels = xlabs, xpd = TRUE, srt = 60, adj = c(1))
  text(x = 1:length(xgrps), y = ylims[1] * 1.1, labels = xlabs, xpd = TRUE, srt = 0, adj = c(0.5, NA))
  axis(2, las = 1)

  invisible(lapply(pd, round, digits = 5))
}
