\name{Simulation of random values from a mixture of spherical Cauchy distributions}
\alias{rmixspcauchy}
\title{
Simulation of random values from a mixture of spherical Cauchy distributions
}

\description{
The function simulates random values simulation from a given mixture of spherical Cauchy distributions.
}

\usage{
rmixspcauchy(n, probs, mu, k)
}

\arguments{
\item{n}{
The sample size.
}
\item{probs}{
This is avector with the mixing probability of each group.
}
\item{mu}{
A matrix with the mean direction of each group.
}
\item{k}{
A vector with the concentration parameter of each group.
}
}

\details{
The function simulates random values simulation from a given mixture of spherical Cauchy distributions.
}

\value{
A list including:
\item{id}{
An indicator of the group of each simulated vector.
}
\item{x}{
A matrix with the simulated data.
}
}

\references{
Tsagris M. (2024). Directional data analysis using the spherical Cauchy and the Poisson-kernel based distribution. https://arxiv.org/pdf/2409.03292
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mixspcauchy.mle}, \link{rmixpkbd}, \link{rmixvmf}
}
}

\examples{
k <- runif(3, 4, 20)
probs <- c(0.2, 0.5, 0.3)
mu <- matrix(rnorm(9), ncol = 3)
mu <- mu / sqrt( rowSums(mu^2) )
x <- rmixspcauchy(200, probs, mu, k)$x
bic.mixspcauchy(x, 5)
}

\keyword{ Random values simulation }
\keyword{ Mixtures of von Mises-Fisher distributions }
