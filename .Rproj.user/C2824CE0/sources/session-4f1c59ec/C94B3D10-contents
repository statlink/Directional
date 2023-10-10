pkbd.mle <- function(x, tol = 1e-6) {

  dm <- dim(x)
  n <- dm[1]  ;  d <- dm[2]
  tx <- t(x)

  pkbd <- function(para, tx, n, d) {
    rho <- 1 / ( 1 + exp(-para[1]) )
    mu <- para[-1]
    mu <- mu / sqrt( sum(mu^2) )
    down <- sqrt( Rfast::colsums( (tx - rho * mu)^2 ) )
    - n * log( 1 - rho^2 ) + d * sum( log(down) )
  }

  m1 <- optim( c( runif(1), rnorm(d) ), pkbd, tx = tx, n = n, d = d, control = list(maxit = 5000) )
  m2 <- optim( c( runif(1), rnorm(d) ), pkbd, tx = tx, n = n, d = d, control = list(maxit = 5000) )
  while (m1$value - m2$value > tol) {
    m1 <- m2
    m2 <- optim( m1$par, pkbd, tx = tx, n = n, d = d, control = list(maxit = 5000) )
  }
  rho <- 1 / ( 1 + exp( -m2$par[1] ) )
  mu <- m2$par[-1]
  mu <- mu / sqrt( sum(mu^2) )
  loglik <- -m2$value
  list(mu = mu, rho = rho, loglik = loglik)
}
