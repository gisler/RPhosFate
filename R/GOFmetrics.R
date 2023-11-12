gmrae <- function(rae) {
  GMRAE <- exp(mean(log(rae)))

  if (qtest(GMRAE, "N1(,)")) {
    GMRAE
  } else {
    NA_real_
  }
}

kge <- function(mld, old) {
  KGE <- 1 - sqrt((cor(mld, old) - 1)^2 + (mean(mld) / mean(old) - 1)^2 +
                    (rcv(mld, old) - 1)^2)

  if (qtest(KGE, "N1(,)")) {
    KGE
  } else {
    NA_real_
  }
}

mdrae <- function(rae) {
  MdRAE <- stats::median(rae)

  if (qtest(MdRAE, "N1(,)")) {
    MdRAE
  } else {
    NA_real_
  }
}

nse <- function(mld, old, j = 2) {
  NSE <- 1 - (sum(abs(old - mld)^j) / sum(abs(old - mean(old))^j))

  if (qtest(NSE, "N1(,)")) {
    NSE
  } else {
    NA_real_
  }
}

pbias <- function(mld, old) {
  PBIAS <- sum(mld - old) / sum(old) * 100

  if (qtest(PBIAS, "N1(,)")) {
    round(PBIAS, 1)
  } else {
    NA_real_
  }
}

rcv <- function(mld, old) {
  RCV <- (sd(mld) / mean(mld)) / (sd(old) / mean(old))

  if (qtest(RCV, "N1(,)")) {
    RCV
  } else {
    NA_real_
  }
}

rmse <- function(mld, old) {
  RMSE <- sqrt(mean((mld - old)^2))

  if (qtest(RMSE, "N1(,)")) {
    RMSE
  } else {
    NA_real_
  }
}

rsr <- function(mld, old) {
  RSR <- rmse(mld, old) / sd(old)

  if (qtest(RSR, "N1(,)")) {
    RSR
  } else {
    NA_real_
  }
}
