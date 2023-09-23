#### Kent regression
kent.reg <- function(y, x, con = TRUE, xnew = NULL, tol = 1e-06) {

  n <- dim(y)[1]
  xbar <- Rfast::colmeans(y)
  S <- crossprod(y)/n
  xbar <- xbar/sqrt(sum(xbar^2))
  u <- c(acos(xbar[1]), (atan(xbar[3]/xbar[2]) + pi * I(xbar[2] < 0))%%(2 * pi))
  theta <- u[1]
  phi <- u[2]
  costheta <- cos(theta)
  sintheta <- sin(theta)
  cosphi <- cos(phi)
  sinphi <- sin(phi)
  H <- matrix(c(costheta, sintheta * cosphi, sintheta * sinphi, 
       -sintheta, costheta * cosphi, costheta * sinphi, 0, -sinphi, 
        cosphi), ncol = 3)
  B <- crossprod(H, S) %*% H
  psi <- 0.5 * atan(2 * B[2, 3]/(B[2, 2] - B[3, 3]))
  K <- matrix(c(1, 0, 0, 0, cos(psi), sin(psi), 0, -sin(psi), cos(psi)), ncol = 3)
  Q <- H %*% K
  y <- y %*% Q

  x <- model.matrix( ~., data.frame(x) )
  if ( !con ) x <- x[, -1]
 

  reg2 <- function(param, y, x) {
    gam1 <- param[1]
    gam2 <- param[2]
    b <- sqrt( gam1^2 + gam2^2 ) 
    be <- matrix(param[ - c( 1:2 ) ], ncol = 3)
    m <- x %*% be
    k <- sqrt( Rfast::rowsums(m) )
    m0 <- sqrt( m[, 2]^2 + m[, 3]^2 )
    rl <- rowSums(m^2)
    x1b <- cbind( -m0^2, m[, 1] * m[, 2], m[, 1] * m[, 3] ) / ( m0 * sqrt(rl) )
    x2b <- cbind( 0, -m[, 3], m[, 2] ) / m0
    Rfast::rowsums(m * y) + crossprod(x1b)


  }
  options(warn = - 1)
  on.exit( options(oop) )
  ini <- rnorm( 3 * dim(x)[2] + 2 )
  val1 <- nlm(reg2, ini, y = y, x = x, za = za, iterlim = 5000)
  val2 <- nlm(reg2, val1$estimate, y = y, x = x, za = za, iterlim = 5000)
  while (val1$minimum - val2$minimum > 1e-06) {
    val1 <- val2
    val2 <-  nlm(reg2, val1$estimate, y = y, x = x, za = za, iterlim = 5000)
  }
  da <- optim(val2$estimate, reg2, y = y, x = x, za = za, control = list(maxit = 10000), hessian = TRUE)
  gam1 <- da$par[1]
  gam2 <- da$par[2]
  rho <- sqrt( gam1^2 + gam2^2 + 1 ) - sqrt( gam1^2 + gam2^2 )
  psi <- 0.5 * acos( 2 * gam1 / ( 1/rho - rho ) )
  be <- matrix(da$par[ -c( 1:2 ) ], ncol = 3)
  m <- x %*% be
  rl <- sqrt( rowsums(m^2) )
  est <- m/rl
  di <- sum( y * est )
  se <- sqrt( diag( solve(da$hessian) ) )
  seb <- matrix(se[-c(1:2)], ncol = 3)

  if ( is.null(xnew) ) {
    mu <- x %*% be
    ki <- sqrt( Rfast::rowsums(mu^2) )
    est <- mu / ki
    fit <- sum( y * est )
  } else {
    xnew <- model.matrix( ~., data.frame(xnew) )
    if ( !con )  xnew <- xnew[, -1]
    mu <- xnew %*% be
    est <- mu / sqrt( Rfast::rowsums(mu^2) )
    fit <- NULL
  }

  if ( is.null( colnames(y) ) ) {
    colnames(est) <- colnames(be) <- colnames(seb) <- c("X", "Y", "Z")
  } else  colnames(est) <- colnames(be) <- colnames(seb) <- colnames(y)
  rownames(be) <- rownames(sb) <- c( colnames(x) )

  param <- c(di, rho, psi)
  names(param) = c("Fit measure", "Rho", "Psi")
  gamma <- c(gam1, gam2)
  names(gamma) <- c("Gamma_1", "Gamma_2")
  list(loglik = -da$value - n * log(2 * pi), param = param, rl = rl,
       gamma = gamma, beta = be, seb = seb, est = est)
}

