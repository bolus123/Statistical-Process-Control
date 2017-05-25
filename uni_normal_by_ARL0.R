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
        k <- L / sqrt((m + 1) / m)
    
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
    
        L <- sqrt((m + 1) / m) * k - delta
        p <- pt(L, df = nu)
        ARL <- ifelse(alternative == '2-sided', 1 / ((1 - p) * 2), 1 / (1 - p))
    
    }

    ARL

}