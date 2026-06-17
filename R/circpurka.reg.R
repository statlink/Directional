circpurka.reg <- function(y, x, rads = TRUE, xnew = NULL) {

  tic <- proc.time()
  if ( !is.matrix(y) ) {
    if ( !rads )   y <- y * pi/180
    z <- cbind( cos(y), sin(y) )
  } else z <- y
  x <- model.matrix( ~., data.frame(x) )

  suppressWarnings({
    ini <- .reg.nr(z, x)
    mod <- optim(ini, .reg, z = z, x = x, method = "BFGS" )
    lik1 <- mod$value
    mod <- optim(mod$par, .reg, z = z, x = x, hessian = TRUE )
    lik2 <- mod$value
    while ( lik1 - lik2 > 1e-6 ) {
      lik1 <- lik2
      mod <- optim(mod$par, .reg, z = z, x = x, hessian = TRUE )
      lik2 <- mod$value
    }
  })
  be <- matrix(mod$par, ncol = 2)
  seb <- solve( mod$hessian )
  seb <- matrix( sqrt( diag(seb) ), ncol = 2)
  runtime <- proc.time() - tic

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew))
    est <- xnew %*% be
    est <- ( atan(est[, 2]/est[, 1]) + pi * I(est[, 1] < 0) ) %% (2 * pi)
    if ( !rads )  est <- est * 180 / pi
  }
  colnames(be) <- colnames(seb) <- c("Cosinus of y", "Sinus of y")
  rownames(be) <- rownames(seb) <- colnames(x)

  list( runtime = runtime, be = be, seb = seb, loglik = - mod$value - dim(x)[1] * log(2), est = est )
}




.reg <- function(be, z, x) {
  be <- matrix(be, ncol = 2)
  est <- x %*% be
  a <- sqrt( Rfast::rowsums(est^2) )
  est <- est / a
  - sum( log(a) - log(1 - exp(-pi * a) ) - a * acos( Rfast::rowsums(z * est) ) )
}

.reg.score <- function(be, z, x) {
  be  <- matrix(be, ncol = 2)
  eta <- x %*% be
  a <- sqrt( Rfast::rowsums(eta^2) )
  m <- eta / a
  c <- Rfast::rowsums(z * m)
  c <- pmin( pmax(c, -1 + 1e-10), 1 - 1e-10 )
  s <- sqrt(1 - c^2)
  d <- acos(c)
  phi <- pi / expm1(pi * a)
  w <- 1 / a - phi - d
  tang <- (z - c * m) / s
  G <- w * m + tang                     # n x 2, grad of loglik_i wrt eta_i
  cbind(-G[, 1] * x, -G[, 2] * x)          # n x 2p, grad of -loglik_i wrt be
}

.reg.nr <- function(z, x, tol = 1e-6, maxiters = 300) {

  be <- as.vector( solve( crossprod(x), crossprod(x, z) ) )

  f <- function(b)  .reg(b, z, x)
  l <- f(be)

  for ( i in 1:maxiters ) {
    com <- .reg.score(be, z, x)
    g <-  Rfast::colsums(com)
    I <- crossprod(com)
    step <- tryCatch( solve(I, g), error = function(e) g )  # fallback if singular
    a.step <- 1
    be.new <- be - a.step * step
    l.new  <- f(be.new)
    while ( !is.finite(l.new) || l.new > l ) {
      a.step <- a.step / 2
      if (a.step < 1e-14) { be.new <- be; l.new <- l; break }
      be.new <- be - a.step * step
      l.new  <- f(be.new)
    }
    conv <- abs(l.new - l) < tol
    be <- be.new ;  l <- l.new
    if (conv)  break
  }
  be
}

