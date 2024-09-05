\name{Simulation of random values from the ESAG distribution in arbitrary dimensions}
\alias{rESAG}
\title{
Simulation of random values from the ESAG distribution in arbitrary dimensions
}

\description{
Simulation of random values from the ESAG distribution in arbitrary dimensions.
}

\usage{
rESAG(n, mu, gam)
}

\arguments{
\item{n}{
A number, how many vectors you want to generate.
}
\item{mu}{
The mean vector the ESAG distribution, a vector in \eqn{R^d}.
}
\item{gam}{
The \eqn{\gamma} parameters of the ESAG distribution. If this is NULL, the sample 
is generated from the IAG distribution.
}
}

\details{
A random sample from the ESAG distribution is generated. In case the \eqn{\gamma_s} are zero, the sample is drawn
from the Independent Angular Gaussian (IAG) (or projected normal) distribution.
}

\value{
An \eqn{n \times d} matrix with the simulated unit vectors.
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
\code{ \link{ESAG.mle}, \link{dESAG}, \link{IAG.mle}, \link{esag.mle}, \link{resag}
}
}

\examples{
m <- colMeans( as.matrix( iris[,1:4] ) )
y <- rESAG(1000, m, c(1, 0.5) )
}

