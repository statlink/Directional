\name{Column-wise MLE of the angular Gaussian and the von Mises Fisher distributions}
\alias{colspml.mle}
\alias{colvm.mle}
\title{
Column-wise MLE of the angular Gaussian and the von Mises Fisher distributions
}

\description{
Column-wise MLE of the angular Gaussian and the von Mises Fisher distributions.
}

\usage{
colspml.mle(x ,tol = 1e-07, maxiters = 100, parallel = FALSE)
colvm.mle(x, tol = 1e-07)
}

\arguments{
\item{x}{
A numerical matrix with data. Each column refers to a different vector of observations of the same distribution.
The values of for Lognormal must be greater than zero, for the logitnormal they must by percentages, exluding 0 and 1,
whereas for the Borel distribution the x must contain integer values greater than 1.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{maxiters}{
The maximum number of iterations that can take place in each regression.
}
\item{parallel}{
Do you want this to be executed in parallel or not. The parallel takes place in C++, and the number of threads
is defined by each system's availiable cores.
}
}

\details{
For each column, spml.mle function is applied that fits the angular Gaussian distribution estimates
its parameters and computes the maximum log-likelihood.
}

\value{
A matrix with four columns. The first two are the mean vector, then the \eqn{\gamma} parameter, and the fourth
column contains maximum log-likelihood.
}

\references{
Presnell Brett, Morrison Scott P. and Littell Ramon C. (1998). Projected multivariate linear models for directional data.
Journal of the American Statistical Association, 93(443): 1068--1077.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{spml.mle}, \link{spml.reg}, \link{vmf.mle}
}
}

\examples{
x <- matrix( runif(100 * 10), ncol = 10)
a <- colspml.mle(x)
b <- colvm.mle(x)
x <- NULL
}
