dpurka <- function(x, a, m, logged = FALSE) {
  con <- a / ( 2 - exp(-a * m) - exp(a * m - a * 2 * pi) ) 
  den <- con * exp(-a * abs(x - m) ) 
  if (logged) den <- log(den)
  den
}


ppurka <- function(x, m, a, logged = FALSE) {
  if (x == m)  {
    prob <- 0.5
  } else if (x < m) {
    con <- a / ( 2 - exp(-a * m) - exp(a * m - a * 2 * pi) ) 
    prob <- con / a * exp(-a * m) * ( exp(a * x) - 1 )
  } else {
    con <- a / ( 2 - exp(-a * m) - exp(a * m - a * 2 * pi) ) 
    prob <- 0.5 + con / a * ( 1 - exp(a * m - a * x) ) 
  }     
  prob
}


con <- function(m) {
  pou <- function(rho) {
    f <- log(a) - log(2) - log( 1 - exp(-a * pi) ) - a * abs(x - m)
    g <-  -log(2 * pi) + log( 1 - rho^2 ) - log( 1 + rho^2 - 2 * rho * cos(x - m) )
    g - f
  }
  mod <- optimise( pou, c(0.001, 0.999) )

  M <- exp( mod$objective)
}


x <- seq(0, 2 * pi, length = 100)
m <- 3
a <- 2
rho <- 0.8

f <- log(a) - log(2) - log( 1 - exp(-a * pi) ) - a * abs(x - m)
g <-  -log(2 * pi) + log( 1 - rho^2 ) - log( 1 + rho^2 - 2 * rho * cos(x - m) )
plot(x, exp(f))
lines(x, 1.5 * exp(g),col=2)

