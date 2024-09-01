\name{MLE of the ESAG distribution in arbitrary dimensions}
\alias{ESAG.mle}
\title{
MLE of the ESAG distribution in arbitrary dimensions
}

\description{
MLE of the ESAG distribution in arbitrary dimensions.
}

\usage{
ESAG.mle(y)
}

\arguments{
\item{y}{
A matrix with the data expressed in Euclidean coordinates, i.e. unit vectors.
}
}

\details{
MLE of the ESAG distribution on arbitrary dimensions is implemented. 
ESAG stands for Elliptically Symmetric Angular Gaussian and it was suugested by Paine et al. (2018). 
It was then generalized to arbitrary dimensions by Yu and Huang (2024).
}

\value{
A list including:
\item{mu}{
The mean vector in \eqn{R^3}.
}
\item{gam}{
The \eqn{\gamma} parameters.
}
\item{loglik}{
The log-likelihood value.
}
\item{vinv}{
The inverse of the covariance matrix. It is returned if the argument "full" is TRUE.
}
\item{rho}{
The \eqn{rho} parameter (smallest eigenvalue of the covariance matrix). It is returned if the argument "full" is TRUE.
}
\item{psi}{
The angle of rotation \eqn{\psi} set this equal to TRUE. It is returned if the argument "full" is TRUE.
}
}

\references{
Zehao Yu and Xianzheng Huang (2024). A new parameterization for elliptically symmetric angular Gaussian distributions of arbitrary dimension.
Electronic Journal of Statististics, 18(1): 301--334.
}

\author{
Zehao Yu.

R implementation and documentation: Zehao Yu \email{zehaoy@email.sc.edu}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{dESAG}, \link{rESAG}, \link{IAG.mle}, \link{esag.mle}
}
}

\examples{
m <- colMeans( as.matrix( iris[,1:3] ) )
y <- rESAG(1000, m, c(1,0.5) )
ESAG.mle(y)
}


