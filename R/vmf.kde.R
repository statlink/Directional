################################
#### Kernel density estimation of circular data with a von Mises kernel
#### Tsagris Michail 2/2015
#### mtsagris@yahoo.gr
#### Garcia-Portugues E. (2013)
#### Exact risk improvement of bandwidth selectors for kernel
#### density estimation with directional data
#### Electronic Journal of Statistics
################################
vmf.kde <- function(x, h = NULL, thumb = "none") {
  ## x is the data
  ## h is the bandwidth you want
  p <- dim(x)[2]  ## dimensionality of the data
  n <- dim(x)[1]  ## sample size of the data
  ## thumb is either 'none' (defualt), or 'rot' (Garcia-Portugues, 2013)
  if ( is.null(h) ) {

    if (thumb == "rot") {
      k <- Rfast::vmf.mle(x)$kappa  ## concentration parameter
      q <- p - 1
      if (q == 2)  h <- ( 8 * sinh(k)^2 / ( k * n * ( (1 + 4 * k^2) * sinh(2 * k) -  2 * k * cosh(2 * k) ) ) )^(1/6)
      if (q >= 3) {
        up <- 4 * pi^0.5 * besselI(k, (q - 1)/2)^2
        down <- k^( (q + 1)/2) * n * (2 * q * besselI(2 * k, (q + 1)/2) +  (2 + q) * k * besselI(2 * k, (q + 3)/2) )
        h <- ( up / down ) ^ ( 1/(4 + q) )
      }

    } else if (thumb == "none")  h <- as.numeric( Directional::vmfkde.tune(x, low = 0.1, up = 1)[1] )

  } else h <- h

  d <- tcrossprod( x / h )
  diag(d) <- 1/h
  cpk <- (1/h^2)^( p/2 - 1) / ( (2 * pi) ^ (p/2) * besselI(1/h^2, p/2 - 1) )
  f <- Rfast::rowmeans( exp(d) ) * cpk
  list( h = h, f = f )
}
