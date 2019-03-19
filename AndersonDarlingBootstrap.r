
AD.gamma.bootstrap.test <- function(X, pars, sim = 10000) {

    integrand.f <- function(u, Fn, pars) {
    
        (Fn(qgamma(u, shape = pars[1], scale = pars[2])) - u) ^ 2 / u / (1 - u)
    
    }
    
    A2.f <- function(n, Fn, pars) {
        n * integrate(integrand.f, 0, 1, Fn = Fn, pars = pars, subdivisions = 700)$value
    }

    n <- length(X)
    
    Fn0 <- ecdf(X)
    
    A20 <- A2.f(n, Fn0, pars)
    
    A2.vec <- unlist(lapply(
                    1:sim,
                    function(x) {
                        Y <- rgamma(n, shape = pars[1], scale = pars[2])
                        Fn <- ecdf(Y)
                        
                        A2.f(n, Fn, pars)
                        
                    }
                ))
                
    pval <- ifelse(A20 > median(A2.vec), 2 * mean(A2.vec > A20), 2 * mean(A2.vec < A20))
                
    list(
        pval = pval
        , A2 = A20
        , A2.vec = A2.vec
    )    
    
}

X <- as.vector(dat) + 0.5


pars <- c(mean(X) / (var(X) / mean(X)), var(X) / mean(X))

result <- AD.gamma.bootstrap.test(X, pars, sim = 10000)
