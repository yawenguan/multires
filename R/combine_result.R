# combine mcmc samples from subposterior
filename=paste0(resultdir,"realdat_Rcpp_batch",1,"_sea",sea,"_nbasis",nbasis,"_P",p,"_fold",fold,".Rda")
load(filename)

# for storage
size = ncol(out$beta)
keep = (size*0.2+1):size
size = length(keep)
beta_array = array(NA,dim=c(p*6,size,120))
A_array = array(NA,dim=c(21,size,120))
psi_array = array(NA,dim=c(6,size,120))
rho_array = array(NA,dim=c(1,size,120))
rm(out)

# load each mcmc sample
for(rep in 1:batchnum){
  filename=paste0(resultdir,"realdat_Rcpp_batch",rep,"_sea",sea,"_nbasis",nbasis,"_P",p,"_fold",fold,".Rda")
  if(file.exists(filename)){
    load(file= filename)
    beta_array[,,rep] = out$beta[,keep]
    A_array[,,rep] = out$A[,keep]
    psi_array[,,rep] = out$psi[,keep]
    rho_array[,,rep] = out$rho[keep]
    rm(out)
  }else{cat(sea,rep)}
}

# index of each parameter
betaidx = 1:dim(beta_array)[1]
Aidx = 1:dim(A_array)[1] + max(betaidx)
psiidx = 1:dim(psi_array)[1] + max(Aidx)
rhoidx = max(psiidx) + 1

# create matrix to store the comined MCMC
comb = array(NA,dim=c(rhoidx,size,1))
subchain = array(NA,dim=c(rhoidx,size,120))
subchain[betaidx,,] <- beta_array
subchain[Aidx,,] <- A_array
subchain[psiidx,,] <- psi_array
subchain[rhoidx,,] <- rho_array

comb[,,1] <- consensusMCcov(subchain[,,which(!is.na(beta_array[1,1,]))])
comb[betaidx,,] -> beta_comb
comb[Aidx,,] -> A_comb
comb[psiidx,,] -> psi_comb
comb[rhoidx,,] -> rho_comb
save(beta_comb,A_comb,psi_comb,rho_comb,
     file=paste0(resultdir,"combchain_Rcpp_sea",sea,"_nbasis",nbasis,"_P",p,"_fold",fold,".Rda"))