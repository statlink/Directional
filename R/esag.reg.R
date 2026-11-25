################################################################################
## ESAG regression (Structure 2 of Paine, Preston, Tsagris & Wood, 2020,
## Statistics and Computing 30: 153-165)
##
##   Q' y_i ~ ESAG( mu(x_i), gamma(x_i) ),   mu(x_i) = B1' x_i,
##   gamma(x_i) = B2' x_i  (gamx = TRUE)  or  gamma constant (gamx = FALSE).
##
## Q is a tuning parameter. The likelihood depends on Q only through the axis
## +-Q e1, because the frame (xi1~, xi2~) of eqs. (11)-(12) is built from
## e1 x mu: it turns with any rotation about e1, and flipping e1 -> -e1 changes
## only the signs of xi1~ and xi2~, which cancel in V^{-1}. Hence the search is
## over axes on a hemisphere (Fibonacci grid), not over the whole of O(3).
##
## Every output except gamma (and its standard errors) is given in the
## ORIGINAL coordinates of y. gamma is defined relative to the frame induced by
## the selected Q, which is returned so that submodels can be fitted with the
## same Q (argument Q).
################################################################################

esag.reg <- function(y, x, con = TRUE, gamx = FALSE, xnew = NULL,
                     n.axes = 50, Q = NULL) {
  ## y:      n x 3 matrix of unit vectors
  ## x:      covariates (vector, matrix or data.frame)
  ## con:    include an intercept (TRUE) or not
  ## gamx:   FALSE: constant gamma;  TRUE: gamma(x) = B2' x
  ## xnew:   new covariate values for prediction (optional)
  ## n.axes: approximate number of axes searched over for Q
  ## Q:      if supplied, use this Q and skip the search (e.g. for submodels)

  X <- model.matrix( ~., data.frame(x) )
  if ( !con )  X <- X[, -1, drop = FALSE]

  if ( is.null(Q) ) {
    wa <- .hemi.axes(n.axes)
    lik <- numeric( nrow(wa) )
    best <- NULL
    for ( j in 1:nrow(wa) ) {
      Qj <- .rot.e1(wa[j, ])                  ## Qj %*% e1 = wa[j, ]
      mod <- .esag.reg2(y %*% Qj, X, gamx)    ## rows are t(Qj) %*% y_i
      lik[j] <- mod$loglik
      if ( is.null(best) || mod$loglik > best$loglik ) {
        best <- mod
        Q <- Qj
      }
    }
    search <- list(axes = wa, loglik = lik)
  } else {
    best <- .esag.reg2(y %*% Q, X, gamx)
    search <- NULL
  }

  .esag.backtransform(best, Q, y, X, gamx, con, xnew, search)
}


## ---------------------------------------------------------------------------
## Fit for a given (already rotated) response z_i = Q' y_i
## ---------------------------------------------------------------------------
.esag.reg2 <- function(z, X, gamx = FALSE) {
  n <- dim(z)[1]
  p <- dim(X)[2]
  kg <- if ( gamx )  2 * p  else  2      ## number of gamma parameters
  za <- cbind(z[, 1]^2, z[, 1] * z[, 2], z[, 1] * z[, 3], z[, 1] * z[, 2],
              z[, 2]^2, z[, 2] * z[, 3], z[, 1] * z[, 3], z[, 2] * z[, 3], z[, 3]^2)

  nll <- function(param) {
    pa <- .esag.unpack(param, X, n, p, gamx)
    v <- .esag.frame(pa$m, pa$G, za)
    if ( any( !is.finite(v$g1) ) || any(v$g1 <= 0) )  return(1e+10)
    g2 <- rowSums(z * pa$m)
    a <- g2 / sqrt(v$g1)
    a2 <- a^2
    M2 <- (1 + a2) * pnorm(a) + a * dnorm(a)
    if ( any( !(M2 > 0) ) )  return(1e+10)
    - 0.5 * sum(a2) + 0.5 * sum(pa$m^2) + 1.5 * sum( log(v$g1) ) - sum( log(M2) )
  }

  ## deterministic starting values: least squares for B1, isotropy for gamma
  ini <- c( numeric(kg), as.vector( .lm.fit(X, z)$coefficients ) )
  val <- optim(ini, nll, method = "BFGS", control = list(maxit = 10000))
  repeat {
    val2 <- optim(val$par, nll, method = "Nelder-Mead", control = list(maxit = 5000))
    val2 <- optim(val2$par, nll, method = "BFGS", control = list(maxit = 10000))
    if ( val$value - val2$value < 1e-06 )  break
    val <- val2
  }
  if ( val2$value < val$value )  val <- val2

  H <- optimHess(val$par, nll)
  covb <- tryCatch( solve(H), error = function(e) matrix(NA, length(val$par), length(val$par)) )
  pa <- .esag.unpack(val$par, X, n, p, gamx)
  list( loglik = -val$value - n * log(2 * pi), par = val$par, covb = covb,
        B1 = pa$B1, B2 = pa$B2, kg = kg, za = za )
}


## ---------------------------------------------------------------------------
## Rotate everything (except gamma) back to the original coordinates
## ---------------------------------------------------------------------------
.esag.backtransform <- function(fit, Q, y, X, gamx, con, xnew, search) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  kg <- fit$kg

  ## mean regression coefficients: mu_orig = Q mu_z  =>  B1_orig = B1_z Q'
  B1 <- fit$B1 %*% t(Q)
  ## covariance: vec(B1_orig) = (Q kron I_p) vec(B1_z)
  Tm <- diag( kg + 3 * p )
  Tm[ -(1:kg), -(1:kg) ] <- kronecker( Q, diag(p) )
  covb <- Tm %*% fit$covb %*% t(Tm)
  se <- suppressWarnings( sqrt( diag(covb) ) )   ## NaN flags a non-PD Hessian
  seB1 <- matrix( se[ -(1:kg) ], ncol = 3 )

  ## gamma (relative to the frame induced by Q)
  if ( gamx ) {
    B2 <- fit$B2
    seB2 <- matrix( se[1:kg], ncol = 2 )
    G <- X %*% B2
  } else {
    B2 <- as.vector(fit$B2)
    seB2 <- se[1:2]
    G <- matrix(B2, n, 2, byrow = TRUE)
  }

  ## fitted mean directions and principal axes, in original coordinates
  m <- X %*% B1
  rl <- sqrt( Rfast::rowsums(m^2) )
  est <- m / rl
  ax <- .esag.axes(X %*% fit$B1, G)      ## computed in the rotated frame
  axis1 <- ax$xi1 %*% t(Q)               ## direction of LEAST dispersion
  axis2 <- ax$xi2 %*% t(Q)               ## direction of GREATEST dispersion
  fitm <- sum(y * est)

  nam <- if ( is.null( colnames(y) ) )  c("X", "Y", "Z")  else  colnames(y)
  colnames(B1) <- colnames(seB1) <- colnames(est) <- colnames(axis1) <- colnames(axis2) <- nam
  rownames(B1) <- rownames(seB1) <- colnames(X)
  if ( gamx ) {
    colnames(B2) <- colnames(seB2) <- c("Gamma_1", "Gamma_2")
    rownames(B2) <- rownames(seB2) <- colnames(X)
  } else  names(B2) <- names(seB2) <- c("Gamma_1", "Gamma_2")

  pred <- NULL
  if ( !is.null(xnew) ) {
    Xn <- model.matrix( ~., data.frame(xnew) )
    if ( !con )  Xn <- Xn[, -1, drop = FALSE]
    mu <- Xn %*% B1
    pred <- mu / sqrt( rowSums(mu^2) )
    colnames(pred) <- nam
  }

  list( loglik = fit$loglik, fit.measure = fitm,
        beta = B1, seb = seB1, gamma = B2, seg = seB2,
        rho = ax$rho, psi = ax$psi, rl = rl, est = est,
        axis1 = axis1, axis2 = axis2, pred = pred,
        Q = Q, covb = covb, search = search )
}


## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------
.esag.unpack <- function(param, X, n, p, gamx) {
  if ( gamx ) {
    B2 <- matrix( param[1:(2 * p)], ncol = 2 )
    G <- X %*% B2
    B1 <- matrix( param[ -(1:(2 * p)) ], ncol = 3 )
  } else {
    B2 <- param[1:2]
    G <- matrix(B2, n, 2, byrow = TRUE)
    B1 <- matrix( param[-(1:2)], ncol = 3 )
  }
  list( B1 = B1, B2 = B2, G = G, m = X %*% B1 )
}

## y' V^{-1} y for each observation, eq. (15)
.esag.frame <- function(m, G, za) {
  m0 <- sqrt( m[, 2]^2 + m[, 3]^2 )
  rl <- sqrt( rowSums(m^2) )
  x1b <- cbind( -m0^2, m[, 1] * m[, 2], m[, 1] * m[, 3] ) / (m0 * rl)
  x2b <- cbind( 0, -m[, 3], m[, 2] ) / m0
  z12 <- x1b[, 1] * x1b[, 2]  ;  z13 <- x1b[, 1] * x1b[, 3]  ;  z23 <- x1b[, 2] * x1b[, 3]
  tx1 <- cbind(x1b[, 1]^2, z12, z13, z12, x1b[, 2]^2, z23, z13, z23, x1b[, 3]^2)
  w23 <- x2b[, 2] * x2b[, 3]
  tx2 <- cbind(0, 0, 0, 0, x2b[, 2]^2, w23, 0, w23, x2b[, 3]^2)
  vinv <- cbind( x1b[, 1] * x2b, x1b[, 1] * x2b[, 2], 2 * x1b[, 2] * x2b[, 2],
                 x1b[, 3] * x2b[, 2] + x1b[, 2] * x2b[, 3], x1b[, 1] * x2b[, 3],
                 x1b[, 3] * x2b[, 2] + x1b[, 2] * x2b[, 3], 2 * x1b[, 3] * x2b[, 3] )
  gam1 <- G[, 1]  ;  gam2 <- G[, 2]
  heta <- sqrt(gam1^2 + gam2^2 + 1) - 1
  g1 <- 1 + rowSums( (gam2 * vinv + gam1 * (tx1 - tx2) + heta * (tx1 + tx2)) * za )
  list(g1 = g1, x1b = x1b, x2b = x2b)
}

## rho, psi and the principal axes xi1, xi2 of eq. (13), per observation
.esag.axes <- function(m, G) {
  m0 <- sqrt( m[, 2]^2 + m[, 3]^2 )
  rl <- sqrt( rowSums(m^2) )
  x1b <- cbind( -m0^2, m[, 1] * m[, 2], m[, 1] * m[, 3] ) / (m0 * rl)
  x2b <- cbind( 0, -m[, 3], m[, 2] ) / m0
  b <- sqrt( G[, 1]^2 + G[, 2]^2 )
  rho <- sqrt(b^2 + 1) - b
  psi <- 0.5 * atan2( G[, 2], G[, 1] )
  xi1 <- cos(psi) * x1b + sin(psi) * x2b
  xi2 <- -sin(psi) * x1b + cos(psi) * x2b
  list(rho = rho, psi = psi, xi1 = xi1, xi2 = xi2)
}

## approximately k axes, evenly spread on the hemisphere with first coordinate > 0
.hemi.axes <- function(k) {
  N <- 2 * k
  i <- seq(0.5, N - 0.5)
  z <- 1 - 2 * i / N
  ph <- pi * (1 + sqrt(5)) * i
  w <- cbind( sqrt(1 - z^2) * cos(ph), sqrt(1 - z^2) * sin(ph), z )
  w[ w[, 1] > 0, , drop = FALSE ]
}

## rotation matrix R with R %*% e1 = q (Rodrigues); stable because q[1] > 0
.rot.e1 <- function(q) {
  q <- q / sqrt( sum(q^2) )
  v <- c(0, -q[3], q[2])                 ## e1 x q
  cc <- q[1]
  K <- matrix( c(0, v[3], -v[2], -v[3], 0, v[1], v[2], -v[1], 0), 3, 3 )
  diag(3) + K + K %*% K / (1 + cc)
}