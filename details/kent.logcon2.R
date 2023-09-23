kent.logcon2 = function (k, b, j = 85) 
{
    j <- 0:j
    ka <- 2 * pi * gamma(j + 0.5) / gamma(j + 1) * b^(2 * j) * 
        (k/2) ^ (-2 * j - 0.5) * besselI(k, 2 * j + 0.5) / gamma(2 * j + 1)
    log( sum(ka) )
}
