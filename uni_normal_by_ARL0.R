#library(adehabitatLT)

#indirect_method <- function(k, m, nu, subdivisions = 100){
#
#
#    pcc.f <- function(Y, k, m, nu) {
#    
#        int <- 1 / (1 - pnorm(k / sqrt((m + 1) / m) * Y)) * dchi(Y, df = nu)
#        int[is.infinite(int)] <- NaN
#        int[is.nan(int)] <- 0
#        
#        int
#    
#    }
#    
#    #debug(pcc.f)
#    
#    integrate(pcc.f, lower = 0, upper = Inf, k = k, m = m, nu = nu, subdivisions = subdivisions)$value
#
#    
#
#}

CC_norm <- function(
                m, 
                nu, 
                ARL = 370, 
                method = 'direct', 
                alternative = '2-sided'
){

    if (method == 'direct'){
    
        L <- ifelse(alternative == '2-sided', qt(1 - 1 / ARL / 2, df = nu),qt(1 - 1 / ARL, df = nu))
        k <- L * sqrt((m + 1) / m) * c4.f(nu)
    
    }   else    {
    
        
    
    }
    
    k
    

}

ARL_norm <- function(
                k,
                m, 
                nu, 
                delta = 0,
                method = 'direct', 
                alternative = '2-sided'
){

    if (method == 'direct'){
    
        L <- sqrt(m / (m + 1)) * k / c4.f(nu) - delta
        p <- pt(L, df = nu)
        ARL <- ifelse(alternative == '2-sided', 1 / ((1 - p) * 2), 1 / (1 - p))
    
    }

    ARL

}


####################################################################################################################################################
    #Example about how to get K by ARL0
####################################################################################################################################################
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
#ARL <- 500
#
#record <- rep(NA, 18)
#
#i <- 0
#
#for (m in m.seq){
#
#    i <- i + 1
#
#    nu <- m - 1
#
#    record[i] <- CC_norm(m, nu, ARL)
#
#}
#
#as.matrix(record)





####################################################################################################################################################
    #Example about how to get ARL0
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
#delta <- 2
#
#for (m in m.seq){
#
#    i <- i + 1
#
#    nu <- m - 1
#
#    record[i] <- ARL_norm(k[i], m, nu, delta)
#
#}
#
#as.matrix(record)




