\name{Hypothesis test for SIPC distribution over the SESPC distribution}
\alias{pc.test}
\title{
Hypothesis test for SIPC distribution over the SESPC distribution
}

\description{
The null hypothesis is whether an SIPC distribution fits the data well, where the altenrative is that SESPC distribution is more suitable.
}

\usage{
pc.test(x, B = 1, tol = 1e-06)
}

\arguments{
\item{x}{
A numeric matrix with three columns containing the data as unit vectors in Euclidean coordinates.
}
\item{B}{
The number of bootstrap re-samples. By default is set to 999. If it is equal to 1, no bootstrap is performed and the
p-value is obtained throught the asymptotic distribution.
}
\item{tol}{
The tolerance to accept that the Newton-Raphson algorithm used in the IAG distribution has converged.
}
}

\details{
Essentially it is a test of rotational symmetry, whether the two \eqn{\theta} parameters are equal to zero. This works for spherical data only.
}

\value{
This is an "htest"class object. Thus it returns a list including:
\item{statistic}{
The test statistic value.
}
\item{parameter}{
The degrees of freedom of the test. If bootstrap was employed this is "NA".
}
\item{p.value}{
The p-value of the test.
}
\item{alternative}{
A character with the alternative hypothesis.
}
\item{method}{
A character with the test used.
}
\item{data.name}{
A character vector with two elements.
}
}

\references{
Tsagris M. and Alzeley O. (2025). Circular and spherical projected Cauchy distributions: 
A Novel Framework for Circular and Directional Data Modeling. 
Australian & New Zealand Journal of Statistics, 67(1): 77--103. https://arxiv.org/pdf/2302.02468.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{iagesag}, \link{fishkent}, \link{sespc.mle}
}
}

\examples{
x <- rvmf(100, rnorm(3), 15)
iagesag(x)
pc.test(x)
}



