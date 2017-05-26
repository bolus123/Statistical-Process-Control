pc.f <- function(z, k, m, n, eta = 1, alternative = '2-sided'){
    
    nu.num <- (1 + z * eta * m / n) ^ 2
    nu.den <- 1 / (m - 1) + (z * eta) ^ 2 * m ^ 2 / n ^ 2 / (n - 1)
    
    nu <- nu.num / nu.den
    
    t <- k / sqrt(1 + z * eta * m / n)
    
    pt(t, nu) - pt(-t, nu)


}

ARL.f <- function(k, m, n, eta = 1, alternative = '2-sided', subdivisions = 100){

    ARLc.f <- function(z, k, m, n, eta = 1, alternative = alternative) 1 / (1 - pc.f(z, k, m, n, eta = eta, alternative = alternative)) * df(z, n - 1, m - 1)

    integrate(ARLc.f, lower = 0, upper = Inf, k = k, m = m, n = n, eta = eta, alternative = alternative, subdivisions = subdivisions)$value


}


ARL0.f <- function(k, m, n, eta = 1, subdivisions = 100, alternative = '2-sided', ARL = 500){

    ARL.f(k, m, n, eta = eta, alternative = alternative, subdivisions = subdivisions) - ARL

}

uniroot(ARL0.f, interval = c(0.0001, 9), m = 100, n = 20, ARL = 500, eta = 1)$root
