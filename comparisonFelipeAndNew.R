
c4.f <- function(nu) sqrt(2 / (nu)) / beta((nu) / 2, 1 / 2) * sqrt(pi)   

L.f <- function(m, n, alpha, eps, p) {

    nu <- m * (n - 1)

    c4 <- c4.f(nu)

    ratio <- qchisq(1 - alpha * (1 + eps), 1) / qchisq(p, nu)

    sqrt(ratio * (m + 1) * (n - 1) * c4 ^ 2)

}

m.vec <- c(15, 20, 25)
n.vec <- c(3, 5, 9)
alpha.vec <- 0.0027
eps.vec <- 0
p.vec <- 0.05

pars.mat <- expand.grid(m.vec, n.vec, alpha.vec, eps.vec, p.vec)

L.f(m = pars.mat[, 1], n = pars.mat[, 2], alpha = pars.mat[, 3], eps = pars.mat[, 4], p = pars.mat[, 5])