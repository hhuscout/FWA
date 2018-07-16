library(MASS)  ; library(nloptr)

# y  : univariate response matrix input
# X  : covariates matrix input
# Kn : oga path number, less than sample size
# varinat: different type of choosing gamma_k
###############################################################
FWA = function(y, X, Kn, variant)
{
  n=nrow(X) ; p=ncol(X) ; y=y-mean(y) ; for(i in 1:p) X[,i]=X[,i] - mean(X[,i])     # 標準化
  beta0.rowsum = rowSums(beta0)   
  #-OGA ------------------------------------------------------------
  jhat = ehat = ehat1 = integer(Kn) ; u=y 
  for(k in 1:Kn)   
  {
    #derivative of Lost function
    LF = colSums( abs(t(u)%*%X) ) 
    #get gamma_k
    if(variant == 0) { gamma = 2 / (k+2) } 
    if(k > 1) { LF[jhat[1:(k-1)]] = 0 }  ; 
    #find
    jhat[k] = which( abs(LF) == max(abs(LF)) )
    
    if(k == 1) { 
      #get gamma_k by minimization
     if(variant == 1) 
     {
      fn = function(gamma) { t(gamma * X[,jhat[1]]) %*% (gamma * X[,jhat[1]]) }
      gamma = optim(0.5, fn, method = "L-BFGS-B", lower = 0, upper = 1) $ par
     }
      #update
     y_hat = gamma * X[,jhat[1]] }
    else    { 
      #get gamma_k by minimization
     if(variant == 1) 
     {
       fn = function(gamma) { t((1-gamma) * y_hat + gamma * X[,jhat[k]]) %*% 
                              ((1-gamma) * y_hat + gamma * X[,jhat[k]]) }
       gamma = optim(0.5, fn, method = "L-BFGS-B", lower = 0, upper = 1) $ par
     }
      #update
     y_hat = (1-gamma) * y_hat + gamma * X[,jhat[k]] }
    u = u - y_hat
    ehat.temp = u ; ehat1[k]=abs(det( t(ehat.temp)%*%ehat.temp ))/ncol(y)
  }
  return(jhat)
}

#--test------------
n=100 ; p=2000 ; eta = 0 ; rho = 0 ; y.col = 1 ; sigma.set=array( rho,c(y.col, y.col)) ; diag(sigma.set)=1  ;  xx=array(0,c(n,p)) 
#beta0=array( c( 3,-3.5,4,-2.8,3.2, rep(0,p-5) ), c(p,1) ) 
beta0=array( c( 3,-3.5,4,-2.8,3.2, 1.3, rep(0,p-6) ), c(p,1) ) 
Kn=min( floor(5*sqrt(n/log(p))), p) 
xx = array( rnorm(n*p, mean=1, sd=sqrt(2.25)), c( n, p) )  ; for(i in 1:p) xx[,i]=xx[,i] + eta*array( rnorm(n), c(n,1) )
yy = xx%*%beta0+array( rnorm(n), c(n,1) )
variant = 0

y = yy ; X = xx 

FWA( yy, xx, Kn, variant)
