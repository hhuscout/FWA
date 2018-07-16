library(MASS)  ; library(mvtnorm) ; library(nloptr)

# y  : univariate response matrix input
# X  : covariates matrix input
# Kn : oga path number, less than sample size
###############################################################
FRA = function(y, X, Kn, variant)
{
  n=nrow(X) ; p=ncol(X) ; y=y-mean(y) ; for(i in 1:p) X[,i]=X[,i] - mean(X[,i])     # 標準化
  beta0.rowsum = rowSums(beta0)   
  #-OGA ------------------------------------------------------------
  jhat = ehat = ehat1 = integer(Kn) ; u=y 
  for(k in 1:Kn)   
  {
    LF = colSums( abs(t(u)%*%X) ) 
    if(k >  1) { LF[jhat[1:(k-1)]] = 0 }  ; jhat[k] = which( abs(LF) == max(abs(LF)) )
    if(k == 1) { fn = function(x) { t(x * X[,jhat[1]]) %*% (x * X[,jhat[1]]) }
      y_hat = optim( par = 0.5, fn, method = "L-BFGS-B", lower = 0, upper = 1 ) $ value }
    else    { 
        fn = function( x ) t((1-x[1]) * y_hat + x[2] * X[,jhat[k]]) %*% 
          ((1-x[1]) * y_hat + x[2] * X[,jhat[k]])
        gamma = optim(par = c(0.5,0.5), fn, method = "L-BFGS-B", lower = 0, upper = 1) $ par [2]
        w     = optim(par = c(0.5,0.5), fn, method = "L-BFGS-B", lower = 0, upper = 1) $ par [1]
      ; y_hat = (1-w) * y_hat + gamma * X[,jhat[k]] }
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

y = yy ; X = xx 

FRA( yy, xx, Kn, variant)
