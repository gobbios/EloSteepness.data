n <- sample(6:15, 1)
x <- EloRating::randomsequence(nID = n, ties = 0)
ids <- colnames(x$pres)[-1] # exclude date column
y <- randomized_sequence_from_matrix(creatematrix(winners = x$seqdat$winner, losers = x$seqdat$loser), n_rand = 5)
oriints <- cbind(as.character(x$seqdat$winner), as.character(x$seqdat$loser))
x <- paste(as.character(x$seqdat$winner), as.character(x$seqdat$loser), sep = "_")
xtab <- table(as.character(x))

res <- logical(5)
for (i in seq_along(res)) {
  testints <- cbind(ids[y$winnermat[, i]], ids[y$losermat[, i]])
  testtab <- table(apply(testints, 1, paste0, collapse = "_"))
  res[i] <- identical(testtab, xtab) & !identical(testints, oriints)
}


test_that("randomized sequences are okay", {
  expect_true(all(res))
})
