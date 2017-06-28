####################################################################################################################################################
    #This code need package parallel for speeding up the process of obtaining pmf of U.modified by simulations
####################################################################################################################################################

source('https://raw.githubusercontent.com/bolus123/Statistical-Process-Control/master/MKLswitch.R')

####################################################################################################################################################


get.U.Modi.Stat <- function(X, Y, type = 'exact'){               
                                                               #The purpose of this function is 
                                                               #to calculate U.modified by the simplified version of U modified
    m <- length(X)                                      
    n <- length(Y)

    X <- sort(X)    
    Y <- sort(Y)

    x.freq <- as.data.frame(table(X))                          #Obtain the empirical cdf of X
    x.prob <- x.freq$Freq / m                                  #
    Fm <- cumsum(x.prob)                                       #Fm is the empirical cdf of X
                                                               #
    y.freq <- as.data.frame(table(Y))                          #Obtain the empirical cdf of Y
    y.prob <- y.freq$Freq / n                                  #
    Gn <- cumsum(y.prob)                                       #Gn is the empirical cdf of Y
   
                
    P <- lapply(X, function(x) ifelse(sum(Y <= x) == 0, 0, Gn[max(which(Y <= x))]))  #calculate Gn(X(i))
    P <- n * unlist(P)                                                               #P = m Gn(X(i))
    
    S <- lapply(Y, function(y) ifelse(sum(X <= y) == 0, 0, Fm[max(which(X <= y))]))  #calculate Fm(Y(j))
    S <- m * unlist(S)                                                               #S = n Fm(Y(j))
    
    P.bar <- sum(P) / m
    S.bar <- sum(S) / n
    
    if (type == 'exact') {
        U.modified <- (sum(P) - sum(S)) / 2 / sqrt((m - 1) / m * sum((S - S.bar) ^ 2) + (n - 1) / n * sum((P - P.bar) ^ 2) + P.bar * S.bar)
    }
    
    else {
        U.modified <- (sum(P) - sum(S)) / 2 / sqrt(sum((S - S.bar) ^ 2) + sum((P - P.bar) ^ 2) + P.bar * S.bar)
    }                                                                                 #According to 3.2
    
    U.modified
  
}


get.U.Modi.Dist <- function(m, n, cores = 1, type = 'exact', maxiter = 1000){
                                                           #The purpose of this function is
                                                           #to obtain the pmf of U.modified under the null hypothesis 
                                                           #by keeping generating 2 samples 
                                                           #from the standard normal distribution
    require(parallel)

    cl <- makeCluster(cores)
    
    clusterExport(cl, c('m', 'n', 'get.U.Modi.Stat', 'type'), envir = environment())
    
    U.modified.dist <- parLapply(                           #This is a loop to generate samples
        cl,                                                 #and then calculate its U modified
        1:maxiter,                                          #
        function(X){                                        #
            x <- rnorm(m)                                   #
            y <- rnorm(n)                                   #
            get.U.Modi.Stat(x, y, type)                     #
        }                                                   #
                                                            #
                                                            #
    )                                                       #
    
    stopCluster(cl)
    
    U.modified.dist <- unlist(U.modified.dist)
    
    U.modified.dist <- as.data.frame(table(U.modified.dist))
    
    U.modified.dist <- cbind(as.numeric(as.character(U.modified.dist$U.modified.dist)), U.modified.dist$Freq / maxiter)
                                                            #The first column will be U modified
                                                            #The second column will be its probability
    colnames(U.modified.dist) <- c('U.modified', 'Prob')
    
    U.modified.dist

}


FP.test <- function(x, y, alternative = '2-sided', alpha = 0.05, type = 'exact', cores = 1, maxiter = 1000){
                                                            #The purpose of this function is 
                                                            #to test 2 samples by WM test with FP modification.
                                                            #The distribution of statistics will be obtained by simulations
    m <- length(x)
    n <- length(y)
    
    U.modified.pmf <- get.U.Modi.Dist(m, n, cores = cores, maxiter = maxiter)
    
    U.modified.dist <- cbind(U.modified.pmf[, 1], cumsum(U.modified.pmf[, 2]))

    U.modified.stat <- get.U.Modi.Stat(x, y, type)
    
    if (alternative == 'less') {                                                              #For testing X < Y
                                                                                              #
        lower.bound <- -Inf                                                                   #
        upper.bound <- U.modified.dist[head(which(U.modified.dist[, 2] >= 1 - alpha), 1), ]   #
        p.value <- sum(U.modified.pmf[which(U.modified.pmf[, 1] > U.modified.stat), 2])       #Obtain the lower tail probability
                                                                                              #
                                                                                              #
    } else if (alternative == 'greater'){                                                     #For testing X > Y
                                                                                              #
        lower.bound <- U.modified.dist[tail(which(U.modified.dist[, 2] <= alpha), 1), ]       #
        upper.bound <- Inf                                                                    #
        p.value <- sum(U.modified.pmf[which(U.modified.pmf[, 1] < U.modified.stat), 2])       #Obtain the upper tail probability
                                                                                              #
    } else if (alternative == '2-sided'){                                                          #For testing X != Y
                                                                                                   #
        upper.bound <- U.modified.dist[head(which(U.modified.dist[, 2] >= 1 - alpha / 2), 1), ]    #
        lower.bound <- U.modified.dist[tail(which(U.modified.dist[, 2] <= alpha / 2), 1), ]        #
        p.value <- sum(U.modified.pmf[which(abs(U.modified.pmf[, 1]) > abs(U.modified.stat)), 2])  #Obtain the 2-tail probability
                                                                                                   #
    } else {
    
        cat('warning:', '\n')
    
    } 

    list(
        lower = lower.bound, 
        upper = upper.bound, 
        stat = U.modified.stat,
        p.value = p.value
    )

}
