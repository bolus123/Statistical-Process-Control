normalOrderStat <- function(t, n, r) {  ## moment of normal order statistics
                                        ## t is moment, n is sample size and r is order

    integrand <- function(x, t, n, r) x ^ t * factorial(n) / factorial(r - 1) / factorial(n - r) * pnorm(x) ^ (r - 1) * dnorm(x) * (1 - pnorm(x)) ^ (n - r)
    
    integrate(integrand, -Inf, Inf, t = t, n = n, r = r)$value

}

normalPITOrderStat <- function(t, n, r) {   ## moment of normal order statistics using PIT
                                            ## t is moment, n is sample size and r is order

    integrand <- function(u, t, n, r) qnorm(u) ^ t * factorial(n) / factorial(r - 1) / factorial(n - r) * u ^ (r - 1) * (1 - u) ^ (n - r)
    
    integrate(integrand, 0, 1, t = t, n = n, r = r)$value

}

normalOrderStat(4, 5, 3)

normalPITOrderStat(4, 5, 3)