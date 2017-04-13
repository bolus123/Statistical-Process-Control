source('https://raw.githubusercontent.com/bolus123/Statistical-Process-Control/master/RLSim.R')


median.sim <- function(n, p, maxsim = 10000){

    unlist(lapply(1:maxsim, function(x) quantile(rnorm(n), p)))

}

prec.limits.sim <- function(m, n, percentile, ARL = 370, x.maxsim = 100000, y.maxsim = 100000){

    p <- 1 / ARL

    L <- quantile(median.sim(n, percentile, maxsim = y.maxsim), c(p / 2, 1 - p / 2))
    LCL <- L[1]
    UCL <- L[2]
    
    L.x <- lapply(
                1:x.maxsim,
                function(x){
                
                    X <- rnorm(m)
                    
                    r <- rank(X)
                    
                    if (any(X > UCL)) {
                    
                        UCL.X <- min(X[X > UCL])
                        b <- min(r[X > UCL])
                    
                    } else {
                    
                        UCL.X <- Inf
                        b <- Inf
                    
                    }
                    
                    if (any(X < LCL)) {
                    
                        LCL.X <- max(X[X < LCL])
                        a <- max(r[X < LCL])
                    
                    } else {
                    
                        LCL.X <- -Inf
                        a <- -Inf
                    
                    }

                    return(c(LCL.X, UCL.X, a, b))
                
                }
            )
    
    u.a <- rep(NA, x.maxsim)
    u.b <- rep(NA, y.maxsim)
    
    for (i in 1:x.maxsim){
    
        u.a[i] <- ifelse(is.infinite(L.x[[i]][3]), 0, L.x[[i]][3] / m)
    
        u.b[i] <- ifelse(is.infinite(L.x[[i]][4]), 1, L.x[[i]][4] / m)
    
    }
    
    u.u.a <- sum(u.a) / x.maxsim
    u.u.b <- sum(u.b) / x.maxsim

    j <- ceiling((n + 1) * percentile)
    
    h <- 0:(n - j)
    
    BB <- 1 / beta(j, n - j + 1)
    
    pc <- lapply(
            1:x.maxsim,
            function(i) {
            
                adds <- 0
            
                for (h in 0:(n - j)){
    
                    adds <- adds + (-1) ^ h / (j + h) * choose(n - j, h) * (u.b[i] ^ (j + h) - u.a[i] ^ (j + h))
                
                }
                
                BB * adds
            
            }
    )
    
    pu <- sum(unlist(pc)) / x.maxsim
    
    
    return(c(u.u.a, u.u.b, 1 - pu))
    
}

#debug(limits.sim)



prec.chart <- function(X, Y, percentile = 0.5, ARL = 370, graph = TRUE, x.maxsim = 100000, y.maxsim = 100000){

    if (is.vector(Y)) Y <- t(as.matrix(Y))

    m <- length(X)
    n <- dim(Y)[2]
    t <- dim(Y)[1]
    
    Lu <- prec.limits.sim(m, n, percentile, ARL, x.maxsim = x.maxsim, y.maxsim = y.maxsim)

    a <- round(Lu[1] * m)
    b <- m - a + 1
    
    r <- rank(X)
    
    LCL <- X[r == a]
    UCL <- X[r == b]
    
    #LCL <- quantile(X, Lu[1])
    #UCL <- quantile(X, Lu[2])
    
    CL <- (UCL - LCL) / 2

    CS <- apply(Y, 1, function(x) quantile(x, percentile))

    if (graph == TRUE){
        plot(c(1, t), c(5 / 4 * min(c(CS, LCL)), 5 / 4 * max(c(CS, UCL))), type = 'n', xlab = 't', ylab = 'CS')
        points(c(1, t), c(UCL, UCL), type = 'l', lty = 2)
        points(c(1, t), c(CL, CL), type = 'l')
        points(c(1, t), c(LCL, LCL), type = 'l', lty = 2)
        points(1:t, CS, type = 'l')
    }
    
 
    names(CL) <- NULL
    list(CL = CL, LCL = LCL, UCL = UCL, CS = CS)
    

}


RL.prec <- function(
                m, 
                n, 
                percentile, 
                ARL = 370, 
                xtype = 'norm', 
                ytype = 'norm', 
                subgroup.amt = 1000,
                shift = 0,
                x.maxsim = 100000, 
                y.maxsim = 100000, 
                maxsim = 10000
){

    Lu <- prec.limits.sim(m, n, percentile, ARL, x.maxsim = x.maxsim, y.maxsim = y.maxsim)

    rl.sim <- rep(NA, maxsim)
    
    for (sim in 1:maxsim){
    
        X <- get.data(m, 1, type = xtype)
        a <- round(Lu[1] * m)
        b <- m - a + 1
        
        r <- rank(X)
        
        LCL <- X[r == a]
        UCL <- X[r == b]
        
        rl <- 0
        
        repeat{
           
            Y <- get.data(n, subgroup.amt = subgroup.amt, type = ytype, shift = shift)
            
            CS <- apply(Y, 1, function(x) quantile(x, percentile))
            
            ll <- check.RL(CS, LCL, UCL)                                                             #check whether this iteration is out of control
                                                                                                     #If so, OOC = TRUE.
            rl <- rl + ll$RL
        
            if (ll$OOC == TRUE) {                                                                    #If OOC = TRUE, stop this iteration
                                                                                                     #
                break                                                                                #
                                                                                                     #
            }                                                                                        #
            
        }
        
        rl.sim[sim] <- rl
    
    }

    rl.sim

}

#debug(RL.prec)
#debug(prec.limits.sim)
a <- RL.prec(
                100, 
                20, 
                0.5, 
                ARL = 370, 
                xtype = 'norm', 
                ytype = 'norm', 
                subgroup.amt = 1000,
                shift = 0,
                x.maxsim = 1000, 
                y.maxsim = 1000, 
                maxsim = 1000
)

########################################
    #example
########################################
#x <- rnorm(100)
#y <- matrix(rnorm(100), ncol = 5)
#
#prec.chart(x, y)
