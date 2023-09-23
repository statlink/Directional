knnreg.tune2 <- function(y, x, M = 10, A = 10, ncores = 1, res = "eucl",
                        type = "euclidean", estim = "arithmetic", mat = NULL, graph = FALSE) {
  y <- as.matrix(y)
  x <- as.matrix(x)
  n <- dim(y)[1]
  d <- dim(y)[2]

  if ( is.null(mat) ) {
    nu <- sample(1:n, min( n, round(n / M) * M ) )
    ## It may be the case this new nu is not exactly the same
    ## as the one specified by the user
    ## to a matrix a warning message should appear
    options(warn = -1)
    mat <- matrix( nu, ncol = M )
  } else  mat <- mat

  M <- dim(mat)[2]
  nu <- dim(mat)[1]
  per <- matrix(nrow = M, ncol = A - 1)
  k <- 2:A
  klen <- length(k)

  if (d == 1  &  type == "euclidean") {
    folds <- list()
    for (i in 1:M)   folds[[ i ]] <- mat[, i]
    if (estim == "arithmetic") {
      method <- "average"
    } else  method <- "harmonic"
    runtime <- proc.time()
    mspe <- Rfast::knn.cv(folds = folds, nfolds = M, y = y, x = x, k = k, type = "R", method = method)$crit
    mspe <- as.vector(mspe)
    performance <- min(mspe)
    runtime <- proc.time() - runtime

  } else {
    runtime <- proc.time()

    for (vim in 1:M) {
      ytest <- y[mat[, vim], , drop = FALSE]  ## test set dependent vars
      ytrain <- y[-mat[, vim], , drop = FALSE]  ## train set dependent vars
      xtest <- x[mat[, vim], , drop = FALSE]  ## test set independent vars
      xtrain <- x[-mat[, vim], , drop = FALSE]  ## train set independent vars

      if ( type == "spher" ) {
        dis <- tcrossprod(xtest, xtrain)
        dis[ dis >= 1 ] <- 1
        disa <- acos(dis)
      } else   disa <- Rfast::dista(xtrain, xtest, trans = FALSE)

      if ( estim == "arithmetic" ) {
        for (j in 1:klen) {
          est <- matrix(nrow = nu, ncol = d)
          knn <- k[j]
          poia <- Rfast::colnth(disa, rep(knn, nu))
          for (i in 1:nu) {
            ind <- which( disa[, i] <= poia[i] )
            if ( length(ind) > knn ) {
              b <- sort( disa[ind, i], index.return = TRUE)$ix[1:knn]
              ind <- ind[b]
            }
            est[i, ] <- Rfast::colmeans( ytrain[ind, , drop = FALSE] )
          }
          if (res == "spher")  {
            est <- est / sqrt( Rfast::rowsums(est^2) )
            per[vim, j] <- 1 - sum( est * ytest )  / nu
          } else  per[vim, j] <- sum( (est - ytest)^2 ) / nu
        }  ## end for (j in k)
      } else {      ## estim = "harmonic"
        for (j in 1:klen) {
          est <- matrix(nrow = nu, ncol = d)
          knn <- k[j]
          poia <- Rfast::colnth(disa, rep(knn, nu))
          for (i in 1:nu) {
            ind <- which( disa[, i] <= poia[i] )
            if ( length(ind) > knn ) {
              b <- sort( disa[ind, i], index.return = TRUE)$ix[1:knn]
              ind <- ind[b]
            }
            est[i, ] <- Rfast::colhameans( ytrain[ind, , drop = FALSE] )
          }  ## end for (i in 1:nu)
          if (res == "spher")  {
            est <- est / sqrt( Rfast::rowsums(est^2) )
            per[vim, j] <- 1 - sum( est * ytest )  / nu
          } else  per[vim, j] <- sum( (est - ytest)^2 ) / nu
        }  ## end for (j in k)
      }   ## end if ( estim == "arithmetic" )
    }  ## end for (vim in 1:M)
    mspe <- Rfast::colmeans(per)
    bias <- per[ , which.min(mspe)] - Rfast::rowMins(per, value = TRUE) ## apply(per, 1, min)
    estb <- mean( bias )  ## TT estimate of bias
    performance <- c( min(mspe) + estb, estb)
    mspe <- Rfast::colmeans(per)
    runtime <- proc.time() - runtime

  }  ## end if (d == 1  &  res == "euclidean")

  names(mspe) <- paste("k=", 2:A, sep = "")
  if ( graph )  plot(2:c(length(mspe) + 1), mspe, xlab = "Nearest neighbours", ylab = "MSPE", type = "b")
  list(crit = mspe, best_k = k[ which.min(mspe) ], performance = performance, runtime = runtime)
}




