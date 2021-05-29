findNearestNeighbour <- function(X, Y, Extent) {
  win <- as.owin(Extent[1:4])
  pppX <- as.ppp(X, W = win)
  pppY <- as.ppp(Y, W = win)

  nn <- data.frame(
    X.x = X[, 1L],
    X.y = X[, 2L],
    X[, 3L, drop = FALSE],
    nncross(pppX, pppY)
  )
  dfY <- data.frame(
    Y.x = Y[, 1L],
    Y.y = Y[, 2L],
    Y[, 3L, drop = FALSE],
    index = seq_len(nrow(Y))
  )
  nn <- merge(nn, dfY, by.x = "which", by.y = "index", sort = FALSE)

  nn[2:8]
}
