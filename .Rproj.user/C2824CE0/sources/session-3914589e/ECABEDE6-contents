dirda.cv <- function(x, ina, folds = NULL, nfolds = 10, k = 2:10, stratified = FALSE,
                  type = c("vmf", "iag", "esag", "kent", "sc", "pkbd", "knn"),
                  seed = NULL, B = 1000, parallel = FALSE) {

  if ( is.null(folds) )  folds <- Directional::makefolds(ina, nfolds = nfolds, stratified = stratified, seed = seed)
  nfolds <- length(folds)

  est1 <- est2 <- est3 <- est4 <- est5 <- est6 <- est7 <- list()
  for (i in 1:nfolds) {
    est1[[ i ]] <- NA
    est2[[ i ]] <- NA
    est3[[ i ]] <- NA
    est4[[ i ]] <- NA
    est5[[ i ]] <- NA
    est6[[ i ]] <- NA
    est7[[ i ]] <- NA
  }
  per1 <- per2 <- per3 <- per4 <- per5 <- per6 <- per7 <- NA
  p <- dim(x)[2]

  if ( sum( type == "vmf") == 1 ) {
    per1 <- matrix(0, nfolds, 1)
    colnames(per1) <- "vmf"
    for (i in 1:nfolds) {
      nu <- folds[[ i ]]
      mod <- Rfast::multivmf.mle(x[-nu, ], ina[-nu], ell = FALSE)
      ki <- mod$ki
      mat <- (p/2 - 1) * log(ki) + ki * tcrossprod(mod$mi, x[nu, , drop = FALSE]) - log( besselI(ki, p/2 - 1, expon.scaled = TRUE) ) - ki
      est1[[ i ]] <- Rfast::colMaxs(mat)
      per1[i] <- mean(est1[[ i ]] == ina[nu])
    }
  }

  if ( sum( type == "iag") == 1 ) {
    per2 <- matrix(0, nfolds, 1)
    colnames(per2) <- "iag"
    g <- max(ina)
    if (p == 3) {
      for (i in 1:nfolds) {
        nu <- folds[[ i ]]
        mat <- matrix(0, length(nu), g)
        xtrain <- x[-nu, ]
        xtest <- x[nu, ]
        id <- ina[-nu]
        for (j in 1:g) {
          mod <- Rfast::iag.mle( xtrain[id == j, ] )
          a <- as.vector(xtest %*% mod$mesi[1, ])
          a2 <- a^2
          pa <- pnorm(a)
          da <- dnorm(a)
          gm <- pa + a2 * pa + a * da
          rl <- mod$param[1]
          mat[, j] <- 0.5 * a2 - 0.5 * rl + log(gm)
        }
        est2[[ i ]] <- Rfast::rowMaxs(mat)
        per2[i] <- mean(est2[[ i ]] == ina[nu])
      }
    } else {
      for (i in 1:nfolds) {
        est2[[ i ]] <- rep(NA, length(folds[[ i ]]) )
        per2[i] <- NA
      }
    }
  }

  if ( sum( type == "esag") == 1 ) {
    per3 <- matrix(0, nfolds, 1)
    colnames(per3) <- "esag"
    g <- max(ina)
    if (p == 3) {
      for (i in 1:nfolds) {
        nu <- folds[[ i ]]
        mat <- matrix(0, length(nu), g)
        xtrain <- x[-nu, ]
        xtest <- x[nu, ]
        id <- ina[-nu]
        for (j in 1:g) {
          mod <- Directional::esag.mle( xtrain[id == j, ] )
          mat[, j] <- Directional::desag(xtest, mod$mu, mod$gam, logden = TRUE )
        }
        est3[[ i ]] <- Rfast::rowMaxs(mat)
        per3[i] <- mean(est3[[ i ]] == ina[nu])
      }
    } else {
      for (i in 1:nfolds) {
        est3[[ i ]] <- rep(NA, length(folds[[ i ]]) )
        per3[i] <- NA
      }
    }
  }

  if ( sum( type == "kent") == 1 ) {

    per4 <- matrix(0, nfolds, 1)
    colnames(per4) <- "kent"
    g <- max(ina)
    if (p == 3) {
      for (i in 1:nfolds) {
        nu <- folds[[ i ]]
        mat <- matrix(0, length(nu), g)
        xtrain <- x[-nu, ]
        xtest <- x[nu, ]
        id <- ina[-nu]
        for (j in 1:g) {
          mod <- Directional::kent.mle( xtrain[id == j, ])
          mat[, j] <- Directional::dkent(xtest, G = mod$G, param = mod$param[1:2], logden = TRUE )
        }
        est4[[ i ]] <- Rfast::rowMaxs(mat)
        per4[i] <- mean(est4[[ i ]] == ina[nu])
      }
    } else {
      for (i in 1:nfolds) {
        est4[[ i ]] <- rep(NA, length(folds[[ i ]]) )
        per4[i] <- NA
      }
    }
  }

  if ( sum( type == "sc") == 1 ) {

    per5 <- matrix(0, nfolds, 1)
    colnames(per5) <- "sc"
    g <- max(ina)
    for (i in 1:nfolds) {
      nu <- folds[[ i ]]
      mat <- matrix(0, length(nu), g)
      xtrain <- x[-nu, ]
      xtest <- x[nu, ]
      id <- ina[-nu]
      for (j in 1:g) {
        mod <- Directional::spcauchy.mle( xtrain[id == j, ])
        mat[, j] <- Directional::dspcauchy(xtest, mod$mu, mod$rho, logden = TRUE)
      }
      est5[[ i ]] <- Rfast::rowMaxs(mat)
      per5[i] <- mean(est5[[ i ]] == ina[nu])
    }
  }

  if ( sum( type == "pkbd") == 1 ) {

    per6 <- matrix(0, nfolds, 1)
    colnames(per6) <- "pkbd"
    g <- max(ina)
    for (i in 1:nfolds) {
      nu <- folds[[ i ]]
      mat <- matrix(0, length(nu), g)
      xtrain <- x[-nu, ]
      xtest <- x[nu, ]
      id <- ina[-nu]
      for (j in 1:g) {
        mod <- Directional::pkbd.mle( xtrain[id == j, ])
        mat[, j] <- Directional::dpkbd(xtest, mod$mu, mod$rho, logden = TRUE)
      }
      est6[[ i ]] <- Rfast::rowMaxs(mat)
      per6[i] <- mean(est6[[ i ]] == ina[nu])
    }
  }

  if ( sum( type == "knn" ) == 1 ) {
    per7 <- matrix(0, nfolds, length(k) )
    colnames(per7) <- paste("knn", k, sep = " ")
    g <- max(ina)
    for (i in 1:nfolds) {
      nu <- folds[[ i ]]
      est7[[ i ]] <- Rfast::dirknn(x[nu, , drop = FALSE], x[-nu, ], ina[-nu], k = k, type = "C", parallel = parallel)
      per7[i, ] <- Rfast::colmeans(est7[[ i ]] == ina[nu])
    }
  }

  per <- cbind(per1, per2, per3, per4, per5, per6, per7)
  perf <- Rfast::colmeans(per)
  names(perf) <- colnames(per)
  names(perf)[1:6] <- c("vmf", "iag", "esag", "kent", "sc", "pkbd")

  if ( any( !is.na(est7) ) )  {

    est <- est7[[ 1 ]]
    for (i in 2:nfolds)   est <- rbind( est, est7[[ i ]] )
    diaf <- est - ina[unlist(folds)]
    n <- dim(est)[1]
    bf <- numeric(B)

    for (i in 1:B) {
      ind <- sample.int(n, n, TRUE)
      m <- Rfast::colmeans(diaf[ind, ] == 0)
      poio <- which.max(m)
      bf[i] <- mean(diaf[-ind, poio] == 0)
    }
    perf[7] <- mean(bf)
    names(perf)[7] <- "knn"
  }

  perf

}
