kent.logcon3 = function (k, b) 
{
    j <- 0:40
    ka <- (2 * pi)^1.5 * k^(-0.5) * factorial(2 * j) / factorial(j)^2 * (b / k) ^ (2 * j) * 
          besselI(k, 2 * j + 0.5) / gamma(2 * j + 1)  ## the final gamma is necessary
    log( sum(ka) )
}
