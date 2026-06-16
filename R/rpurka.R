rpurka <- function(n, mu, a) {
  u <- runif(n)
  theta <-  -log( 1 - u * ( 1 - exp(-a * pi) ) ) / a
  x <- ifelse(runif(n) < 0.5, 1, -1)
  (mu + x * theta) %% (2 * pi)
}