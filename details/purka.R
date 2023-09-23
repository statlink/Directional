circ.laplace <- function(x) {
  fun <- function(pa, x, n) {
    a <- pa[1]   ;   m <- pa[2] 
    - n * log(a) + n * log( 2 - exp( -a * m ) - exp(a * m - 2 * pi * a) ) + a * sum( abs(x - m) ) 
  }
  n <- length(x)
  mod <- optim( runif(2, 0, 2 * pi), fun, x = x, n = n, method = "L-BFGS-B", lower = c(0.00001, min(x)), 
         upper = c( 500, max(x) ), control = list(maxit = 2000) )   
  param <- c(mod$pa) 
  names(param) <- c("alpha", "mu")
  list(loglik = - mod$value, param = param )
}


rcirc <- function(n, a, m) {
  u <- runif(n)
  con <- ( 2 - exp(-a * m) - exp(a * m - 2 * pi * a) ) / a
  u1 <- u[u < 0.5]
  x1 <- log(1 + a / con * u1 * exp(a * m) ) / a
  u1 <- u[u > 0.5]
  x2 <- m - log(1 - a / con * (u1 - 0.5) ) / a
  c(x1, x2)
}
 
x <- rcirc(10000, 0.6, 3)
purka.mle(x)
m2 <- purka.mle(x)$theta
( atan(m2[2] / m2[1]) + pi * I(m2[1] < 0) ) %%(2 * pi)
circ.laplace(x)
  
  

