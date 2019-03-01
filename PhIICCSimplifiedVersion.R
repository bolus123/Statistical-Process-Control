c4.f <- function(nu) sqrt(2 / (nu)) / beta((nu) / 2, 1 / 2) * sqrt(pi)   

tau.f <- function(m, ub.cons, lambda = 1) {

    
    m / (m + 1) / ub.cons ^ 2 * lambda ^ 2

}

qPhIICARL <- function(p, CC, tau, nu){

    inner <- qchisq(p, nu) / nu * tau * CC ^ 2
    
    (1 - pchisq(inner, 1)) ^ (-1)

}

Moment1PhIICARL <- function(CC, tau, nu) {

    integrate(qPhIICARL, 0, 1, CC = CC, tau = tau, nu = nu)$value

}

Moment1PhIICARL <- Vectorize(Moment1PhIICARL, c('CC', 'tau', 'nu'))

Moment2PhIICARL <- function(CC, tau, nu) {

    qPhIICARL2 <- function(p, CC, tau, nu) {
    
        qPhIICARL(p = p, CC = CC, tau = tau, nu = nu)^2
    
    }

    integrate(qPhIICARL2, 0, 1, CC = CC, tau = tau, nu = nu)$value

}

VarPhIICARL <- function(CC, tau, nu) {

    Moment2PhIICARL(CC = CC, tau = tau, nu = nu) - Moment1PhIICARL(CC = CC, tau = tau, nu = nu)^2

}

VarPhIICARL <- Vectorize(VarPhIICARL, c('CC', 'tau', 'nu'))


tau.f(aa[, 1], c4.f(nu.vec), lambda = 1)

Moment1PhIICARL(aa[, 2], tau.f(aa[, 1], c4.f(nu.vec), lambda = 1), nu.vec)
sqrt(VarPhIICARL(aa[, 2], tau.f(aa[, 1], c4.f(nu.vec), lambda = 1), nu.vec))


getCC.EPC <- function(p0, nu, tau, ARL0 = 370, eps = 0, lambda = 1, ub.cons = 1) {

    ARLb <- (1 - eps) * ARL0

    numer <- nu * qchisq(1 - 1/ARLb, 1)
    denom <- tau * qchisq(p0, nu)
    
    sqrt(numer / denom)

}
#
#m <- 10000
#nu <- m - 1
#
#c4 <- c4.f(nu)
#
#tau1 <- tau.f(m, ub.cons = c4, lambda = 1)
#
#getCC.EPC(p0 = 0.1, nu = nu, tau = tau1, ARL0 = 370, eps = 0, lambda = 1, ub.cons = 1) 

getCC <- function(ARL0, m, nu, ub.cons = 1, lambda = 1, interval = c(1, 3)) {

    root.finding <- function(CC, ARL0, tau, nu) {
    
        ARLin <- integrate(qPhIICARL, 0, 1, CC = CC, tau = tau, nu = nu)$value
    
        cat('CC:', CC, 'ARLin:', ARLin, '\n')
    
        ARL0 - ARLin
    
    }
    
    tau1 <- tau.f(m = m, ub.cons = ub.cons, lambda = lambda)
    
    uniroot(root.finding, interval, ARL0 = ARL0, tau = tau1, nu = nu)$root


}





m <- 10000
nu <- m - 1


c4 <- c4.f(nu)


getCC(ARL0 = 370, m = m, nu = nu, ub.cons = c4, lambda = 1, interval = c(1, 3))

