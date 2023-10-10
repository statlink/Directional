dpkbd <- function(y, mu, rho, logden = FALSE) {

  y <- as.matrix(y)
  if ( dim(y)[2] == 1 )  y <- t(y)
  d <- dim(y)[2]

  down <- sqrt( Rfast::colsums( (t(y) - rho * mu)^2 ) )
  den <- log( 1 - rho^2 ) - d * log(down)
  if ( !logden )  den <- exp(den)
  den

}
