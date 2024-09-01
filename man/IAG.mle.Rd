\name{MLE of the IAG distribution in arbitrary dimensions}
\alias{IAG.mle}
\title{
MLE of the IAG distribution in arbitrary dimensions
}

\description{
MLE of the IAG distribution in arbitrary dimensions.
}

\usage{
IAG.mle(y)
}

\arguments{
\item{y}{
A matrix with the data expressed in Euclidean coordinates, i.e. unit vectors.
}
}

\details{
MLE of the independent angular Gaussian (IAG) distribution, in arbitrary dimensions, is implemented. 
ESAG stands for Elliptically Symmetric Angular Gaussian and it was suugested by Paine et al. (2018). 
It was then generalized to arbitrary dimensions by Yu and Huang (2024). The IAG is a special case of the ESAG,
when the covariance matrix equals the identity matrix.
}

\value{
A list including:
\item{mu}{
The mean vector in \eqn{R^d}.
}
\item{loglik}{
The log-likelihood value.
}
}

\references{
Zehao Yu and Xianzheng Huang (2024). A new parameterization for elliptically symmetric angular Gaussian distributions of arbitrary dimension.
Electronic Journal of Statististics, 18(1): 301--334.

Paine P.J., Preston S.P., Tsagris M. and Wood A.T.A. (2018). An Elliptically Symmetric
Angular Gaussian Distribution. Statistics and Computing, 28: 689--697.
}

\author{
Zehao Yu.

R implementation and documentation: Zehao Yu \email{zehaoy@email.sc.edu}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rIAG}, \link{IAGd}, \link{ESAG.mle}
}

\examples{
m <- colMeans( as.matrix( iris[, 1:4] ) )
mu <- m / sqrt( rowsums(m^2) )
y <- rESAG(1000, mu, gam = NULL)
IAG.mle(y)
}


