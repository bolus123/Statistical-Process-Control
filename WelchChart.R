pc.f <- function(z, k, m, n, eta = 1, alternative = '2-sided'){
    
    nu.num <- (1 + z * eta * m / n) ^ 2
    nu.den <- 1 / (m - 1) + (z * eta) ^ 2 * m ^ 2 / n ^ 2 / (n - 1)
    
    nu <- nu.num / nu.den
    
    x <- k / sqrt(1 + z * eta * m / n)
    
    pt(x, nu) - pt(-x, nu)


}

ARL.f <- function(k, m, n, eta = 1, alternative = '2-sided', subdivisions = 100, rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol){

    ARLc.f <- function(z, k, m, n, eta = 1, alternative = alternative) (1 - pc.f(z, k, m, n, eta = eta, alternative = alternative)) ^ (-1) * df(z, n - 1, m - 1)

    integrate(ARLc.f, lower = 0, upper = Inf, k = k, m = m, n = n, eta = eta, alternative = alternative, subdivisions = subdivisions, rel.tol = rel.tol, abs.tol = abs.tol)$value


}


ARL0.f <- function(k, m, n, eta = 1, subdivisions = 100, ARL = 500, alternative = '2-sided', rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol){

    ARL.f(k, m, n, eta = eta, alternative = alternative, subdivisions = subdivisions, rel.tol = rel.tol, abs.tol = abs.tol) - ARL

}

Welch.K <- function(m, n, ARL, eta, search.lower = 0.0001, search.upper = 9, subdivisions = 100, rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol){

	uniroot(ARL0.f, interval = c(search.lower, search.upper), m = m, n = n, ARL = ARL, eta = eta, subdivisions = subdivisions, rel.tol = rel.tol, abs.tol = abs.tol)$root

}

Welch.K(100, 10, 500, 1, search.lower = 0.0001, search.upper = 9, subdivisions = 1000)
##uniroot(ARL0.f, interval = c(0.0001, 9), m = 100, n = 20, ARL = 500, eta = 1)$root
