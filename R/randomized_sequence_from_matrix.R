#' convert interaction matrix to sequence
#'
#' @param mat square matrix
#' @param n_rand number of randomized sequences to produce (default is \code{1})
#'
#' @return a list with a matrix of winner ids and a matrix with loser ids
#' @importFrom utils getFromNamespace
#' @importFrom EloRating prunk
#' @export
#'
#' @examples
#' data(bonobos, package = "EloRating")
#' bonobos[1:4, 1:4]
#' randomized_sequence_from_matrix(bonobos[1:4, 1:4])

randomized_sequence_from_matrix <- function(mat, n_rand = 1) {
  # adapted from EloSteepness:::prep_data_for_rstan
  m_foo <- getFromNamespace("mat2seqint", "EloRating")

  n <- sum(mat)
  x <- m_foo(mat)

  locmat <- matrix(rep(seq_len(n), n_rand), ncol = n_rand)
  locmat <- apply(locmat, 2, sample)

  # winner and loser index matrices
  winnermat <- matrix(x[[1]][locmat], ncol = n_rand)
  losermat <- matrix(x[[2]][locmat], ncol = n_rand)

  list(winnermat = winnermat, losermat = losermat)
}
