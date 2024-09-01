IAGd <- function(y, mu, logden = FALSE) {  
  y <- as.matrix(y)
  if  ( dim(y)[2] == 1 )   y <- t(y)
  d <- dim(y)[2]
  p <- d - 1
  Cd <- (2 * pi)^(-0.5 * p)
  a <- as.vector( y %*% mu )
  Mp <- matrix(1, nrow = n, ncol = p)
  Mp[, 1] <- a * pnorm(a) + dnorm(a)
  Mp[, 2] <- (1 + a^2) * pnorm(a) + a * dnorm(a)
  if ( p >= 3 )   for ( j in 3:p )  Mp[, j] <- a * Mp[, j - 1] + (j - 1) * Mp[, j - 2]
  l <- log(Cd) + 0.5*( a^2 - sum(mu^2) ) + log( Mp[, p] )
  if ( logden )  {
    return( l )
  } else   return( exp(l) )
}
