####################################################################################################################################################   
    #Based on Champ and Jones(2004), Rose and Does(1995)
####################################################################################################################################################
library(mvtnorm)
library(adehabitatLT)
####################################################################################################################################################

source('https://raw.githubusercontent.com/bolus123/Statistical-Process-Control/master/MKLswitch.R')

####################################################################################################################################################
    #parts of getting L
####################################################################################################################################################

c4.f <- function(nu) sqrt(2 / nu) * 1 / beta(nu / 2, 1 / 2) * sqrt(pi)             #c4.function

corr.f <- function(m, off.diag = - 1 / (m - 1)){
                                                                                   #correlation matrix
    crr <- diag(m)
    crr[which(crr == 0)] <- off.diag

    crr
    
}

####################################################################################################################################################
    #get L by multivariate T
####################################################################################################################################################

get.L.mvt <- function(
                 m 
                 ,nu
                 ,FAP = 0.1
                 #,ARL = 370
                 ,Phase1 = TRUE
                 ,off.diag = NULL
                 ,alternative = '2-sided'
                 ,maxiter = 10000
                 #,MCMC = FALSE
                 #,MCMC.search.interval = c(1, 5)
                 #,MCMC.maxsim = 10000
                 #,MCMC.burn = 1000

){
                                                                              #The purpose of this function is 
    MCMC <- FALSE                                                             #to obtain L and K based on
                                                                              #multivariate T.
                                                                              #MCMC part is not available now.

    #if (off.diag == NULL) off.diag <- ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))
    
    corr.P <- corr.f(m = m, off.diag = off.diag)                              #get correlation matrix with equal correlations

    pu <- 1 - FAP
    
    if (MCMC == TRUE) {
        #MVN.Q.Gibbs.Sampling(
        #    pu, 
        #    MCMC.maxsim, 
        #    corr.P, 
        #    tails = alternative, 
        #    burn = MCMC.burn, 
        #    search.interval = MCMC.search.interval
        #)
    
        
    
    } else {
 
        L <- ifelse(
                alternative == '2-sided',
                qmvt(pu, df = nu, sigma = corr.P, tail = 'both.tails', maxiter = maxiter)$quantile,
                qmvt(pu, df = nu, sigma = corr.P, maxiter = maxiter)$quantile
            )
                                                      #get L by multivariate T
        

    }
    
    
    K <- L * c4.f(nu) * sqrt((m - 1) / m)             #get K
   
    list(L = L, K = K)

}

####################################################################################################################################################
    #get L by multivariate Normal
####################################################################################################################################################

joint.pdf.mvn.chisq <- function(Y, K, m, nu, sigma, alternative = '2-sided') {

    s <- length(Y)

    #L <- K / sqrt((m - 1) / m * nu) * Y * c4.f(nu)
    L <- K / sqrt((m - 1) / m * nu) * Y / c4.f(nu)
    
    dpp <- lapply(
            1:s,
            function(i){

                LL <- rep(L[i], m)
                
                ifelse(
                    alternative == '2-sided',
                    pmvnorm(lower = -LL, upper = LL, sigma = sigma),
                    pmvnorm(lower = rep(-Inf, m), upper = LL, sigma = sigma)
                )
            
            } 
    
    )    

    dpp <- unlist(dpp)
    
    dpp * dchi(Y, df = nu)


}



root.mvn.F <- function(K, m, nu, sigma, pu, alternative = '2-sided', subdivisions = 2000, rel.tol = 1e-2) {
                                                            #The purpose of this function is
    #s <- length(Y)                                          #to obtain appropriate L
    #                                                        #by multivariate normal
    #L <- K / sqrt((m - 1) / m * nu) * Y * c4.f(nu)          #
    #
    #pp <- lapply(
    #        1:s,
    #        function(i){
    #        
    #            LL <- rep(L[i], m)
    #            ifelse(
    #                alternative == '2-sided',
    #                pmvnorm(lower = -LL, upper = LL, sigma = sigma),
    #                pmvnorm(lower = rep(-Inf, m), upper = LL, sigma = sigma)
    #            )
    #        
    #        } 
    #
    #)
    #
    #pp <- mean(unlist(pp))
    
    pp <- integrate(
            joint.pdf.mvn.chisq, 
            lower = 0, 
            upper = Inf, 
            K = K, 
            m = m, 
            nu = nu, 
            sigma = sigma, 
            alternative = alternative, 
            subdivisions = subdivisions, 
            rel.tol = rel.tol
        )$value
    
    pu - pp


}

#sigma = corr.f(20)
#joint.pdf.mvn.chisq(1, 3, 20, 80, sigma, alternative = '2-sided')
#root.mvn.F(3, 20, 80, sigma, 0.7)

get.L.mvn <- function(
                 m
                 ,nu
                 ,FAP = 0.1
                 #,ARL = 370
                 ,Phase1 = TRUE
                 ,off.diag = NULL
                 ,alternative = '2-sided'
                 #,maxsim = 10000
                 ,maxiter = 10000
                 ,subdivisions = 2000
                 ,tol = 1e-2   
                 #,MCMC = FALSE
                 #,MCMC.search.interval = c(1, 5)
                 #,MCMC.maxsim = 10000
                 #,MCMC.burn = 1000

){                                                          #The purpose of this function is 
                                                            #to obtain L and K based on
    MCMC <- FALSE                                           #multivariate normal.
                                                            #MCMC part is not available now.
    if (is.null(off.diag)) off.diag <- ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))
    
    
    corr.P <- corr.f(m = m, off.diag = off.diag)  

    pu <- 1 - FAP
    
    #Y <- sqrt(rchisq(maxsim, df = nu))
    
    if (MCMC == TRUE) {
    
        #MVN.Q.Gibbs.Sampling(
        #    pu, 
        #    MCMC.maxsim, 
        #    corr.P, 
        #    tails = alternative, 
        #    burn = MCMC.burn, 
        #    search.interval = MCMC.search.interval
        #)
    
        
    
    } else {

        #L <- ifelse(
        #        alternative == '2-sided',
        #        qmvnorm(pu, sigma = corr.P, tail = 'both.tails', maxiter = maxiter)$quantile,
        #        qmvnorm(pu, sigma = corr.P, maxiter = maxiter)$quantile
        #    )

        K <- uniroot(
                root.mvn.F,
                interval = c(1, 7),
                m = m,
                nu = nu,
                #Y = Y,
                sigma = corr.P,
                pu = pu,
                alternative = alternative,
                subdivisions = 2000,
                tol = tol,
                rel.tol = tol,
                maxiter = maxiter
        )$root

    }
    
    
    #K <- mean(L / Y * sqrt((m - 1) / m * nu)) #/ c4.f(nu))
    
    L <- K / c4.f(nu) * sqrt(m / (m - 1))
    
    list(L = L, K = K)
    

}

#get.L.mvn(20, 80)

####################################################################################################################################################

get.L <- function(
            m
            ,nu
            ,FAP = 0.1
            #,ARL = 370
            #,Phase1 = TRUE
            ,off.diag = NULL
            ,alternative = '2-sided'
            ,maxiter = 10000
            ,method = 'direct'
            ,indirect_subdivisions = 100L
            ,indirect_tol = .Machine$double.eps^0.25
            #,indirect.maxsim = 10000
            #,MCMC = FALSE
            #,MCMC.search.interval = c(1, 5)
            #,MCMC.maxsim = 10000
            #,MCMC.burn = 1000

){                                                  #The purpose of this function is to obtain L and K
                                                    #by multivariate T or multivariate normal.
                                                    #Multivariate normal is so time-consuming
                                                    #that I do not recommend.
    Phase1 <- TRUE
    
    
    if (is.null(off.diag)) off.diag <- ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))

    is.int <- ifelse(nu == round(nu), 1, 0)

    if (method == 'direct' & is.int == 1) {                       #using multivariate T to obtain L and K
    
        get.L.mvt(
            m = m
            ,nu = nu
            ,FAP = FAP
            ,Phase1 = Phase1
            ,off.diag = off.diag
            ,alternative = alternative
            ,maxiter = maxiter
        )
    
    } else {                                       #using multivariate normal to obtain L and K
    
        if (is.int == 0) cat('Nu is not an integer. The indirect method will be conducted.', '\n')
    
        get.L.mvn(
            m = m
            ,nu = nu
            ,FAP = FAP
            ,Phase1 = Phase1
            ,off.diag = off.diag
            ,alternative = alternative
            #,maxsim = indirect.maxsim
            ,subdivisions = indirect_subdivisions
            ,maxiter = maxiter
            ,tol = indirect_tol
        )
    
    }


}

####################################################################################################################################################
    #reverse the process
####################################################################################################################################################
get.FAP0 <- function(
            K
            ,m
            ,nu
            #,ARL = 370
            #,Phase1 = TRUE
            ,off.diag = NULL
            ,alternative = '2-sided'
            ,maxiter = 10000
            #,method = 'direct'
            ,indirect_subdivisions = 100L
            ,indirect_tol = .Machine$double.eps^0.25
            #,indirect.maxsim = 10000
            #,MCMC = FALSE
            #,MCMC.search.interval = c(1, 5)
            #,MCMC.maxsim = 10000
            #,MCMC.burn = 1000
){

    method <- 'direct'

    Phase1 <- TRUE
    
    is.int <- ifelse(nu == round(nu), 1, 0)
    
    if (is.null(off.diag)) off.diag <- ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))
    
    corr.P <- corr.f(m = m, off.diag = off.diag)
    
    L <- K / c4.f(nu) * sqrt(m / (m - 1))
    

    if (method == 'direct' & is.int == 1) {                       #using multivariate T to obtain L and K
    
        ll <- rep(L, m)
        
        ifelse(
            alternative == '2-sided',
            1 - pmvt(lower = -ll, upper = ll, df = nu, corr = corr.P, algorithm = TVPACK, abseps = 1e-12),
            1 - pmvt(lower = -Inf, upper = ll, df = nu, corr = corr.P, algorithm = TVPACK, abseps = 1e-12)
        )
    
    } else{
    
    
    
    }




}


####################################################################################################################################################
    #Example about how to get K and L by FAP0
####################################################################################################################################################
#get.L(
#            20
#            ,80
#            ,FAP = 0.1
#            ,off.diag = NULL
#            ,alternative = '2-sided'
#            ,maxiter = 10000
#            ,method = 'direct'
#)
#get.L(
#            20
#            ,80
#            ,FAP = 0.1
#            ,off.diag = NULL
#            ,alternative = '2-sided'
#            ,maxiter = 10000
#            ,method = 'indirect'
#)

####################################################################################################################################################
    #Example about how to get FAP by K
####################################################################################################################################################
#
#k <- c(
#22.84794
#,10.52163
#,7.386245
#,6.057107
#,5.340985
#,4.898076
#,4.598652
#,4.383316
#,3.842428
#,3.619871
#,3.498807
#,3.422752
#,3.280788
#,3.182391
#,3.151007
#,3.135567
#,3.126384
#,3.120294
#
#
#
#
#
#)
#
#m.seq <- c(
#    3
#    ,4
#    ,5
#    ,6
#    ,7
#    ,8
#    ,9
#    ,10
#    ,15
#    ,20
#    ,25
#    ,30
#    ,50
#    ,100
#    ,150
#    ,200
#    ,250
#    ,300
#
#)
#
#
#i <- 0
#
#record <- rep(NA, 18)
#
#for (m in m.seq){
#
#    i <- i + 1
#
#    nu <- m - 1
#
#    record[i] <- get.FAP0(k[i], m, nu)
#
#}
#
#as.matrix(record)




