# some functions that might be used. 
fun1stat <- function(idx,predDat,names) round(mean(( (predDat[idx,names]) -  (predDat[idx,"Concentration"]))^2,na.rm=T),2)
fun2stat <- function(idx,predDat,names) round(mean(( (predDat[idx,names]) -  (predDat[idx,"Concentration"]))),2)
fun3stat <- function(idx,predDat,names) {if(!is.na(mean((predDat[idx,names])))){(round(mean((predDat[idx,names])),2))}else{return(NA)}}
fun4stat <- function(idx,predDat,names) round(cor( (predDat[idx,names]), (predDat[idx,"Concentration"]), use="pairwise.complete.obs"),2)
fun5stat <- function(idx,predDat,names) round(mean(abs( (predDat[idx,names]) -  (predDat[idx,"Concentration"]))),2)

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p)))){stop("Dimension problem!")}
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

K2A<-function(k,q){
  K_fill <- matrix(NA,q,q)
  K_fill[,1] <- k[1:q];K_fill[2:3,2] <- k[(q+1):(q+2)];K_fill[3,3] <- k[(q+3)]
  K_fill[1,2] <- K_fill[2,1];K_fill[1:2,3]<- K_fill[3,1:2]
  A_fill <- t(chol(K_fill))
  A <-A_fill[A_fill!=0]
  return(A)
}


A2K<-function(k,q){
  A <- matrix(0,q,q)
  A[lower.tri(A,TRUE)] <- k
  K <- A%*%t(A)
  return(K[lower.tri(A,TRUE)])
}


A2cor<-function(k,q){
  A <- matrix(0,q,q)
  A[lower.tri(A,TRUE)] <- k
  K <- A%*%t(A)
  cor <- cov2cor(K)
  return(cor)
}

K2cor<-function(k,q){
  K_fill <- matrix(NA,q,q)
  K_fill[,1] <- k[1:q];K_fill[2:3,2] <- k[(q+1):(q+2)];K_fill[3,3] <- k[(q+3)]
  K_fill[1,2] <- K_fill[2,1];K_fill[1:2,3]<- K_fill[3,1:2]
  cor <- cov2cor(K_fill)
  return(cor)
}

forecastMultiSpMisalign <- function(n,pred.covars,pred.coords,beta,A.each,psi,rho,mkCov,dist=rdist){
  N = sum(n);q = length(pred.covars);p=ncol(pred.covars[[1]]);nltr <- q*(q+1)/2
  pred.coord=do.call("rbind",pred.coords)
  
  n.p = unlist(lapply(pred.covars,ncol))
  X = matrix(0,nrow=N,ncol=sum(n.p))
  for(i in 1:q) {if(n[[i]]!=0) X[(cumsum(n)[i]-n[i]+1):cumsum(n)[i],(cumsum(n.p)[i]-n.p[i]+1):cumsum(n.p)[i]] = as.matrix(pred.covars[[i]])}
  niter = nrow(A.each)
  ysave = matrix(NA,niter,N)
  
  for(iter in 1:niter){
  A = matrix(0,q,q)
  A[lower.tri(A,TRUE)] <- A.each[iter,]
  V_sp  = mkCov(A,rho[iter,],pred.coords,psi[iter,],n,dist)
  mu = X%*%beta[iter,]
  ysave[iter,]<- mu + t(chol(V_sp))%*%rnorm(length(mu))
  }
  return(list(forecast = ysave))
}

forecastUnivSpMisalign  <- function(sp,Xp,betas,params){
  
  library(mvtnorm)
  library(geoR)
  
  corfx <- function(d,params){
    r   <- params[2] # sigma-to-noise ratio
    rho <- params[3]        # range
    nu  <- params[4]        # smoothness
    COR <- r*matern(d,rho,nu)
    COR[d==0]<-1
    return(COR)}
  
  np        <- nrow(sp)
  d11       <- rdist(sp,sp)
  diag(d11) <- 0
  pred.y    <- matrix(0,nrow(betas),np)
  
  for(iters in 1:nrow(betas)){
    beta = betas[iters,]
    sigma2 = params[iters,1] 
    par  = params[iters,] 
    
    # 
    S11      <- corfx(d11,par)*sigma2
    MMM      <- Xp%*%beta
    
    pred.y[iters,] <- rmvnorm(1,MMM,S11)
    #MMM + t(chol(VVV))%*%rnorm(np)
  }
  return(pred.y)
}

predMultiSpMisalign <- function(y,x,coords,n,pred.covars,pred.coords,beta,A.each,psi,rho,mkCov,mkVcros,mkVar,dist=rdist){
  q = length(x);p=ncol(x[[1]]);nltr <- q*(q+1)/2
  Y = NULL;pred.y=NULL
  for(i in 1:q){
    if(n[i]>0){
      Y = c(Y,y[[i]]); 
    }
  }
  N = length(Y)
  n.p = unlist(lapply(x,ncol))
  X = matrix(0,nrow=N,ncol=sum(n.p))
  for(i in 1:q) {if(n[[i]]!=0) X[(cumsum(n)[i]-n[i]+1):cumsum(n)[i],(cumsum(n.p)[i]-n.p[i]+1):cumsum(n.p)[i]] = as.matrix(x[[i]])}
  niter = nrow(A.each)
  
  for(i in 1:length(pred.coords)){
    if(!any(is.na(pred.coords[[i]]))&nrow(pred.coords[[i]])>0){
    ysave = matrix(NA,niter,nrow(pred.coords[[i]]))
    for(t in 1:niter){
      A = matrix(0,q,q)
      A[lower.tri(A,TRUE)] <- A.each[t,]
      V_sp  = mkCov(A,rho[t,],coords,psi[t,],n,dist)
      r<-Y-X%*%beta[t,]
      
      Var   = mkVar(A,rho[t,i],pred.coords[[i]],psi[t,i],i,dist)
      Vcros = mkVcros(A,rho[t,],pred.coords[[i]],i,coords,psi[t,],n,dist)
      mu = as.matrix(pred.covars[[i]])%*%beta[t,(cumsum(n.p)[i]-n.p[i]+1):cumsum(n.p)[i]] + Vcros%*%solve(V_sp,r)
      var = Var-Vcros%*%solve(V_sp,t(Vcros))
      ysave[t,]<- mu + t(chol(var))%*%rnorm(length(mu))
      # cat("mcmciter",t,"\n")
    }
    pred.y[[i]] <-ysave
    }else{
      pred.y[[i]] <- NA
    }
  }
  return(pred.y)
}

predmapMultiSpMisalign <- function(y,x,coords,n,pred.covars,pred.coords,beta.hat,A.each.hat,psi.hat,rho.hat,mkCov,mkVcros,mkVar,dist=rdist){
  cat(head(beta.hat))
  q = length(x);p=ncol(x[[1]]);nltr <- q*(q+1)/2
  Y = NULL;pred.y=NULL
  for(i in 1:q){
    if(n[i]>0){
      Y = c(Y,y[[i]]); 
    }
  }
  N = length(Y)
  n.p = unlist(lapply(x,ncol))
  X = matrix(0,nrow=N,ncol=sum(n.p))
  for(i in 1:q) {if(n[[i]]!=0) X[(cumsum(n)[i]-n[i]+1):cumsum(n)[i],(cumsum(n.p)[i]-n.p[i]+1):cumsum(n.p)[i]] = as.matrix(x[[i]])}
  
  
  for(i in 1:length(pred.coords)){
    if(!any(is.na(pred.coords[[i]]))&nrow(pred.coords[[i]])>0){
    
    
	  # A.each.hat = colMeans(A.each)
	  # rho.hat = colMeans(rho)
	  # psi.hat = colMeans(psi)
	  # beta.hat = colMeans(beta)
	  
      A = matrix(0,q,q)
      A[lower.tri(A,TRUE)] <- A.each.hat
      V_sp  = mkCov(A,rho.hat,coords,psi.hat,n,dist)
      r<-Y-X%*%beta.hat
      
      Vcros = mkVcros(A,rho.hat,pred.coords[[i]],i,coords,psi.hat,n,dist)
      mu = as.matrix(pred.covars[[i]])%*%beta.hat[(cumsum(n.p)[i]-n.p[i]+1):cumsum(n.p)[i]] + Vcros%*%solve(V_sp,r)
	  
    pred.y[[i]] <-mu
    }else{
      pred.y[[i]] <- NA
    }
  }
  return(pred.y)
}

predmapMultiSpMisalign_mean <- function(y,x,coords,n,pred.covars,pred.coords,beta,A.each,psi,rho,mkCov,mkVcros,mkVar,dist=rdist){
  q = length(x);p=ncol(x[[1]]);nltr <- q*(q+1)/2
  Y = NULL;pred.y=NULL
  for(i in 1:q){
    if(n[i]>0){
      Y = c(Y,y[[i]]); 
    }
  }
  N = length(Y)
  n.p = unlist(lapply(x,ncol))
  X = matrix(0,nrow=N,ncol=sum(n.p))
  for(i in 1:q) {if(n[[i]]!=0) X[(cumsum(n)[i]-n[i]+1):cumsum(n)[i],(cumsum(n.p)[i]-n.p[i]+1):cumsum(n.p)[i]] = as.matrix(x[[i]])}
  niter = nrow(A.each)
  
  for(i in 1:length(pred.coords)){
    if(!any(is.na(pred.coords[[i]]))&nrow(pred.coords[[i]])>0){
      
      mu = rep(0,nrow(pred.coords[[i]]))      
      for(t in 1:niter){
      
      A = matrix(0,q,q)
      A[lower.tri(A,TRUE)] <- A.each[t,]
      V_sp  = mkCov(A,rho[t,],coords,psi[t,],n,dist)
      r<-Y-X%*%beta[t,]
      
      Vcros = mkVcros(A,rho[t,],pred.coords[[i]],i,coords,psi[t,],n,dist)
      mu = mu+as.matrix(pred.covars[[i]])%*%beta[t,(cumsum(n.p)[i]-n.p[i]+1):cumsum(n.p)[i]] + Vcros%*%solve(V_sp,r)

      }
      pred.y[[i]] <-mu/niter
    }else{
      pred.y[[i]] <- NA
    }
  }
  return(pred.y)
}

forecastmapMultiSpMisalign <- function(n,pred.covars,pred.coords,beta,A.each,psi,rho,mkCov,dist=rdist){
  N = sum(n);q = length(pred.covars);p=ncol(pred.covars[[1]]);nltr <- q*(q+1)/2
  pred.y=NULL
  pred.coord=do.call("rbind",pred.coords)
  
  n.p = unlist(lapply(pred.covars,ncol))
  X = matrix(0,nrow=N,ncol=sum(n.p))
  for(i in 1:q) {if(n[[i]]!=0) X[(cumsum(n)[i]-n[i]+1):cumsum(n)[i],(cumsum(n.p)[i]-n.p[i]+1):cumsum(n.p)[i]] = as.matrix(pred.covars[[i]])}
  
  mu = X%*%colMeans(beta)
  return(list(forecast = mu))
}


predUnivSpMisalign  <- function(y,X,s,sp,Xp,betas,params){
  
  library(mvtnorm)
  library(geoR)
  
  corfx <- function(d,params){
    r   <- params[2] # sigma-to-noise ratio
    rho <- params[3]        # range
    nu  <- params[4]        # smoothness
    COR <- r*matern(d,rho,nu)
    COR[d==0]<-1
    return(COR)}
  
  np        <- nrow(sp)
  d         <- rdist(s,s)
  d12       <- rdist(sp,s)
  d11       <- rdist(sp,sp)
  diag(d11) <- 0
  pred.y    <- matrix(0,nrow(betas),np)
  
  for(iters in 1:nrow(betas)){
    beta = betas[iters,]
    sigma2 = params[iters,1] 
    par  = params[iters,] 
    
    # 
    S11      <- corfx(d11,par)*sigma2
    S12      <- corfx(d12,par)*sigma2
    Sigma    <- corfx(d,par)
    Siginv   <- solve(Sigma)
    S22inv   <- Siginv/sigma2
    VVV      <- S11-S12%*%S22inv%*%t(S12)
    MMM      <- Xp%*%beta + S12%*%S22inv%*%(y-X%*%beta)
    
    pred.y[iters,] <- rmvnorm(1,MMM,VVV)
      #MMM + t(chol(VVV))%*%rnorm(np)
  }
  return(pred.y)
}

logit <- function(theta,a,b){
  return(log((theta-a)/(b-theta)))
}

logitinv <- function(z,a,b){
  return(b-(b-a)/(1+exp(z)))
} 

rinvgamma<-function (n, shape, scale = 1) 
{
  return(1/rgamma(n = n, shape = shape, rate = scale))
}

quants <- function(x){if(length(x)>0){c(quantile(x, prob=c(0.025,0.5,0.975)),mean(x))}else{NA}}

getmode <- function(x){foo = density(x);return(foo$x[which.max(foo$y)])}

library(scales)
library(fields)
cols = tim.colors(64)

plotpoints = function(loc,obs,zrange=NULL,xlim=range(loc[,1]),ylim=range(loc[,2]),alpha = 0.5,cex=0.1,main="",add=F,legend = F,...){
  if(length(obs)!=nrow(loc)) print("warningloc")
  if(is.null(zrange)) zrange = quantile(obs,prob=c(0.025,0.975),na.rm=T)
  obs[obs>zrange[2]|obs<zrange[1]] <- NA
  col=cols[findInterval(obs,seq(zrange[1],zrange[2],len=64))]
  col[is.na(col)] = "gray"
  if(add){
    points(x=loc[,1],y=loc[,2],col=scales::alpha(col, alpha),pch=19,cex=cex)
  }else{
    plot(x=loc[,1],y=loc[,2],col=scales::alpha(col, alpha),xlab="",ylab="",axes=F,
         xlim=xlim,ylim=ylim,pch=19,cex=cex,main=main)
    if(legend)image.plot(legend.only = TRUE,zlim=zrange,xlab=NULL,ylab=NULL,
               midpoint=F,col=cols,...)
  }
}

