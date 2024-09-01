\name{Density of the IAG distribution in arbitrary dimensions}
\alias{dESAG}
\title{
Density of the IAG distribution in arbitrary dimensions
}

\description{
Density of the IAG distribution in arbitrary dimensions.
}

\usage{
IAGd(y, mu, logden = FALSE)
}

\arguments{
\item{y}{
A matrix or a vector with the data expressed in Euclidean coordinates, i.e. unit vectors.
}
\item{mu}{
The mean vector the IAG distribution, a vector in \eqn{R^d}.
}
\item{logden}{
If you the logarithm of the density values set this to TRUE.
}
}

\details{
The density of the IAG distribution in arbitrary dimensions is computed.
}

\value{
A vector with the (log) density values of y.
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
\code{\link{rIAG}, \link{IAG.mle}, \link{ESAG.mle}
}
}

\examples{
m <- colMeans( as.matrix( iris[,1:3] ) )
mu <- m / sqrt( rowsums(m^2) )
y <- rIAG(1000, mu)
a <- IAGd(y, mu)
}

