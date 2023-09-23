kern.reg <- function(xnew = NULL, y, x, h, res = "eucl", type = "gauss") {
  if ( !is.matrix(y)  & res == "spher" )  y <- cbind( cos(y), sin(y) )
  y <- as.matrix(y)
  d <- dim(y)[2] 
  x <- as.matrix(x)

  if ( type == "gauss" || type == "laplace" ) {
    m <- Rfast::colmeans(x)
    s <- Rfast::colVars(x, std = TRUE)
    tx <- ( t(x) - m ) / s  ## standardize the independent variables
    x <- t(tx)
    if ( !is.null(xnew) ) {
      xnew <- as.matrix(xnew)
      txnew <- ( t(xnew) - m ) / s  ## standardize the x values
      xnew <- t(txnew)
      if (type == "gauss") {
        a1 <-  - 0.5 * Rfast::dista(xnew, x )^2
      } else  a1 <-  - Rfast::dista(xnew, x, type = "manhattan" )
      z <- exp(a1/h)
      ta <- Rfast::rowsums(z)
      ta[ta == 0] <- Inf
      mhx <- ( z %*% y) / ta
    } else {
      if (type == "gauss") {
        a1 <-  - 0.5 * Rfast::Dist(x, square = TRUE )
      } else  a1 <-  - Rfast::Dist(x, method = "manhattan" )
      z <- exp(a1/h)
      ta <- Rfast::rowsums(z)
      ta[ta == 0] <- Inf
      mhx <- ( z %*% y) / ta
    }
  } else { 
    xnew <- as.matrix(xnew)
    disa <- tcrossprod(xnew, x)
    disa[ disa >= 1 ] <- 1
    f <- exp( disa / h^2 ) 
    ta <- Rfast::rowsums(f)
    ta[ta == 0] <- Inf
    est <- ( f %*% y ) / ta
  }
  if (d == 1)  est <- as.vector(est)
  if ( res == "spher" )  est <- est / sqrt( Rfast::rowsums(est) )
  est
}






