\name{BIC for the model based clustering using mixtures of spherical Cauchy distributions}
\alias{bic.mixspcauchy}
\title{
BIC to choose the number of components in a model based clustering using mixtures of spherical Cauchy distributions
}

\description{
BIC to choose the number of components in a model based clustering using mixtures of spherical Cauchy distributions
}

\usage{
bic.mixspcauchy(x, G = 5, n.start = 10)
}

\arguments{
\item{x}{
A matrix containing directional data.
}
\item{G}{
The maximum number of clusters to be tested. Default value is 5.
}
\item{n.start}{
The number of random starts to try. See also R's built-in function \code{\link{kmeans}} for more information about this.
}
}

\details{
If the data are not unit vectors, they are transformed into unit vectors.
}

\value{
A plot of the BIC values and a list including:
\item{BIC}{
The BIC values for all the models tested.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
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
\code{\link{mixpkbd.mle}, \link{rmixpkbd}, \link{bic.mixspcauchy}, \link{bic.mixvmf}
}
}

\examples{
x <- as.matrix( iris[, 1:4] )
x <- x / sqrt( rowSums(x^2) )
bic.mixpkbd(x)
}
