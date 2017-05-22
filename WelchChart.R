pc.f <- function(z, k, m, n, alternative = '2-sided'){

    nu.num <- (1 + z * m * (n - 1) / n / (m - 1)) ^ 2
    nu.den <- 1 / (m - 1) + z ^ 2 * m ^ 2 * (n - 1) / n ^ 2 / (n - 1) ^ 2
    
    nu <- nu.num / nu.den
    
    t <- sqrt(m) * k / sqrt(1 + z * m * (n - 1) / n / (m - 1))
    
    pt(t, nu) - pt(-t, nu)


}

ARL.f <- function(k, m, n, alternative = '2-sided'){

    ARLc.f <- function(z, k, m, n, alternative = '2-sided') 1 / (1 - pc.f(z, k, m, n, alternative = '2-sided')) * df(z, n - 1, m - 1)

    integrate(ARLc.f, lower = 0, upper = Inf, k = k, m = m, n = n, alternative = '2-sided')$value


}


ARL0.f <- function(k, m, n, alternative = '2-sided', ARL = 500){

    ARL.f(k, m, n, alternative = '2-sided') - ARL

}

uniroot(ARL0.f, interval = c(0.0001, 0.8), m = 100, n = 20, ARL = 500)$root


ARL.f(0.8, 100, , alternative = '2-sided')