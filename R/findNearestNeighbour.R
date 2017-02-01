#' @export
findNearestNeighbour <- function(X, Y, Extent) {
  win <- spatstat::as.owin(Extent[1:4])
  pppX <- spatstat::as.ppp(X, W = win)
  pppY <- spatstat::as.ppp(Y, W = win)

  nn <- spatstat::nncross(pppX, pppY, what = c("dist", "which"))
  nn <- data.frame(
    X.x = X[, 1],
    X.y = X[, 2],
    X[, 3, drop = FALSE],
    nn
  )
  dfY <- data.frame(
    Y.x = Y[, 1],
    Y.y = Y[, 2],
    Y[, 3, drop = FALSE],
    index = seq_len(nrow(Y))
  )
  nn <- merge(nn, dfY, by.x = "which", by.y = "index", sort = FALSE)

  return(nn[2:8])
}
