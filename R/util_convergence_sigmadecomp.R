
Aitken <- function (ll) {
  if (any(is.infinite(ll))) {
    linf <- Inf
    a <- NA
  }
  else {
    if (ll[2] > ll[1]) {
      a <- (ll[3] - ll[2]) / (ll[2] - ll[1])
    }
    else {
      a <- 0
    }
    if (a < 1) {
      linf <- ll[3] + (ll[3] - ll[2]) / (1 - a)
    }
    else {
      linf <- Inf
    }
  }
  list(ll = ll[3], linf = linf, a = a)
}

SIGMA_COMP <-
  function(variance) {
    #variance is a list obtained by mstep in MCLUST in which information regarding the estimated variance component is contained
    #this function returns scale shape and orientation for those models for which scale shape orientation are not returned by default
    RES <- variance
    if (variance$modelName == "EII" | variance$modelName == "VII") {
      RES$shape <- rep(1, variance$d)
    }
    if (variance$modelName %in% c("EII", "VII", "EEI", "VEI", "EVI", "VVI")) {
      RES$orientation <- diag(variance$d)
    }
    if (variance$modelName == "EEE") {
      eig <- eigen(variance$Sigma)
      RES$scale <- exp(1 / variance$d * sum(log(eig$values)))
      RES$shape <- exp(log(eig$values) - log(RES$scale))
      RES$orientation <- eig$vectors
    }
    if (variance$modelName == "VVV") {
      eigenval <-
        abs(apply(variance$sigma, 3, function(x)
          eigen(x, only.values = TRUE)$val)) # abs is used for avoiding numerical problems (eigenv positive semidefinite matrix needs to always be positive)
      RES$scale <-
        tryCatch(
          expr = exp(1 / variance$d * colSums(log(eigenval))),
          warning = function(w)
            NA
        )
      RES$shape <- exp(sweep(log(eigenval), 2, log(RES$scale), "-"))
      RES$orientation <- array(apply(variance$sigma, 3, function(x)
        eigen(x)$vec),
        dim = dim(variance$sigma))
    }
    RES
  }


