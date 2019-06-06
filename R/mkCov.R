# # make covariance matrix; only works for equal range; otherwise need to fix the diagonal block
# use Euclidean distance, does not allow for missing outputs
mkCov <- function(A,rho,coords,psi,n=NULL,dist = rdist){
  # assume rho are the same for all outputs. can change this later
  K = A%*%t(A)
  q = nrow(A)
  if(is.null(n)){
  n = rep(0,q)
  for (i in 1:q) n[i] <- nrow(coords[[i]])
  }
  N = sum(n)
  V = matrix(NA,N,N)
  
  for (i in 1:q){   
    if(n[i]!=0) V[(1+cumsum(c(0,n))[i]):(cumsum(n)[i]),(1+cumsum(c(0,n))[i]):(cumsum(n)[i])] = 
      K[i,i]*exp(-dist(coords[[i]],coords[[i]])*rho[i])
  }
  
  for(i in 1:(q-1)){
    for(j in (i+1):q){
      if(n[i]!=0 & n[j]!=0){
      Dcross <- dist(coords[[i]],coords[[j]])
      
      if(!all(rho==rho[1])){
        temp = matrix(0,nrow = n[i],ncol= n[j])
        for(k in 1:q) {
          temp <-temp+A[i,k]*t(A)[k,j]*exp(-Dcross*rho[k])
        }
      }else{
        temp <- K[i,j]*exp(-Dcross*rho[1])
      }
      
      V[(1+cumsum(c(0,n))[i]):(cumsum(n)[i]),(1+cumsum(c(0,n))[j]):(cumsum(n)[j])] <- temp
      V[(1+cumsum(c(0,n))[j]):(cumsum(n)[j]),(1+cumsum(c(0,n))[i]):(cumsum(n)[i])] <- t(temp)
      }
    }
  }
  
  for(i in 1:q){
    if(n[i]!=0) diag(V)[(1+cumsum(c(0,n))[i]):(cumsum(n)[i])] <- K[i,i]+psi[i]
  }
  return(V)
}

mkVcros <- function(A,rho,predcoords,qpol,coords,psi,n=NULL,dist=rdist){ # for one pollutant. can change this for multiple pollutant later
  # assume rho are the same for all outputs. can change this later
  K = A%*%t(A)
  q = nrow(A)
  N = sum(n)
  Npred = sum(nrow(predcoords))
  V = matrix(NA,Npred,N)
  
  for(j in 1:q){
    if(n[j]!=0){
      
      Dcross  <- dist(predcoords,coords[[j]])
      
      if(!all(rho==rho[1])){
      temp = matrix(0,nrow = Npred,ncol= n[j])
      for(k in 1:q) {
        temp <-temp+A[qpol,k]*t(A)[k,j]*exp(-rho[k]*Dcross)
      }
      }else{
      temp <- K[qpol,j]*exp(-rho[1] * Dcross)
      }
      
      V[,(1+cumsum(c(0,n))[j]):(cumsum(n)[j])] <- temp
    }
  }
  
  return(V)
}
  
mkVar <- function(A,rho,coords,psi,qpol,dist=rdist){
  # assume rho are the same for all outputs. can change this later
  K = A%*%t(A)
  V = K[qpol,qpol]*exp(-dist(coords,coords)*rho)
  diag(V) = diag(V) + psi
  return(V)
}

# mkCov1 <- function(A,rho,coords,psi){
#   foo = NULL
#   K = A%*%t(A)
#   q = nrow(A)
#   for (i in 1:q) foo[[i]] = K[i,i]*exp(-rdist(coords[[i]])*rho[i])
# 
#     D12  = rdist(coords[[1]],coords[[2]]);  D13  = rdist(coords[[1]],coords[[3]]);D23  = rdist(coords[[2]],coords[[3]])
#   foo23 = foo13 = foo12 = NULL
#   for (i in 1:q) foo12[[i]] = exp(-D12*rho[i])
#   for (i in 1:q) foo13[[i]] = exp(-D13*rho[i])
#   for (i in 1:q) foo23[[i]] = exp(-D23*rho[i])
# 
#   V12 = matrix(0,nrow=nrow(D12),ncol=ncol(D12))
#   V13 = matrix(0,nrow=nrow(D13),ncol=ncol(D13))
#   V23 = matrix(0,nrow=nrow(D23),ncol=ncol(D23))
# 
#   for (i in 1:q){
#   V12 = V12 + A[1,i]*t(A)[i,2]*foo12[[i]]
#   V13 = V13 + A[1,i]*t(A)[i,3]*foo13[[i]]
#   V23 = V23 + A[2,i]*t(A)[i,3]*foo23[[i]]
#   }
# 
#   V = rbind(cbind(foo[[1]],V12,V13),cbind(t(V12),foo[[2]],V23),cbind(t(V13),t(V23),foo[[3]]))
#   diag(V)[1:nrow(coords[[1]])] <- K[1,1]+psi[1]
#   diag(V)[nrow(coords[[1]])+1:nrow(coords[[2]])] <- K[2,2]+psi[2]
#   diag(V)[nrow(coords[[2]])+nrow(coords[[1]])+1:nrow(coords[[3]])] <- K[3,3]+psi[3]
#   return(V)
# }


# coords = cbind(runif(4),runif(4))
# coords=list(coords[1:2,],coords[3:4,],coords[c(1,4),])
# C1 = mkCov(A,rho,coords=coords,psi=Psi)
# # # compare with spBayes
# sl = sk = 1;
# coordsD = rdist(rbind(coords[[1]],coords[[2]],coords[[3]])); n= c(nrow(coords[[1]]),nrow(coords[[2]]),nrow(coords[[3]]))
# N = sum(n)
# C = matrix(NA,N,N)
# for(k in 1:q){
#   sl = 1;
#   for(l in 1:q){
#     for(kk in 1: n[k]){
#       for(jj in 1: n[l]){
#         C[(sl+jj-2)*N+(sk+kk-2)+1] = 0.0;
#         for(ii in 1:q){
#           C[(sl+jj-2)*N+(sk+kk-2)+1] = C[(sl+jj-2)*N+(sk+kk-2)+1] + A[k+q*(ii-1)]*A[l+q*(ii-1)]*exp(-coordsD[(sl+jj-2)*N+(sk+kk-2)+1]*rho[ii])
#         }
#       }
#     }
#     sl = sl + n[l]
#   }
#   sk = sk + n[k]
# }
# 
# 
#   sl = 1;
#   for(l in 1:q){
#     for(k in 1: n[l]){
#       C[(sl+k-2)*N+(sl+k-2)+1] =C[(sl+k-2)*N+(sl+k-2)+1]+ Psi[l];
#     }
#     sl =sl + n[l];
#   }
# 
# C==C1
