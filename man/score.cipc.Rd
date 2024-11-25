\name{Score test for many simple CIPC regressions}
\alias{score.cipc}
\title{
Score test for many simple CIPC regressions
}

\description{
Score test for many simple CIPC regressions.
}

\usage{
score.cipc(y, X, rads = TRUE, tol = 1e-06)
}

\arguments{
\item{y}{
The dependent variable,  a numerical vector, it can be in radians or degrees.
}
\item{X}{
A matrix with many numerical independent variables.
}
\item{rads}{
If the dependent variable is expressed in rads, this should be TRUE and FALSE otherwise.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm in the null model (CIPC without covariates).
}
}

\details{
The score test uses the first derivative (score function) of the regression log-likelihood and it is asymptotically correct. So, this function requires sample sizes or at least 1,000 observations.
}

\value{
A matrix with two columns, the test statistic and its associated p-value.
}

\references{
Tsagris M. and Alzeley O. (2024). Circular and spherical projected Cauchy distributions: A Novel Framework for Circular and Directional Data Modeling. Australian & New Zealand Journal of Statistics (accepted for publication). https://arxiv.org/pdf/2302.02468.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{cipc.reg}, \link{cpc.mle}, \link{spher.reg}, \link{spml.regs}
}
}

\examples{
y <- rcipc(500, omega = 2, g = 5)
x <- matrix( rnorm(500 * 10), ncol = 10 )
a <- score.cipc(y, x)
}


