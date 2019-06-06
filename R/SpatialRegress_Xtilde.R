# fit multivariate spatial model to 3 days as a batch
# load library and source functions
library(Rcpp)
library(fields)
source(paste0(fndir,'mkCov.R'))

# setup
sea = "JAS"; fold = 1;nbasis = 5; batch = 3; report = 0.1*iters

# load prediction site for five fold CV
load(paste0(datadir,"predsite.Rda"))
predsite = PREDSITE[[fold]]
load(paste0(datadir,"proc",nbasis,".RData")); 
df$one = 1
# scale both spectral covariates and station data
foonames = unlist(lapply(poll_names, function(pol) c(paste0("Xtilde.",pol,1:nbasis),paste0("X.",pol))))
df[,colnames(df)%in%foonames] = scale(df[,colnames(df)%in%foonames])
df_bypol = split(df,df$Pol)
df_bypol = lapply(df_bypol, function(x){
  x[,"Concentration"] = scale(x[,"Concentration"])
  return(x)}
  )
df = do.call(rbind, df_bypol)
predDat = subset(df, Site %in% predsite)
df = subset(df, !Site %in% predsite)

# make three days as a batch for model fitting ----------------------------
days =  range(df$Sample.Date)
set = seq(min(days),max(days),by = batch)
for (rep in 1:batchnum){
dateset = set[rep]:(set[rep+1]-1)
cat("Fitting model to date",dateset,"\n")

# # Build training data set
traindf = subset(df, Sample.Date%in%dateset)
traindf$Pol = factor(traindf$Pol, levels=poll_names)
traindf = split(traindf,traindf$Sample.Date)
y.train = lapply(traindf,function(x) split(x[,"Concentration"],f=x$Pol))
x.train = lapply(traindf,function(x) lapply(1:q, function(i) subset(x,Pol==poll_names[i])[,covar.names[[i]]]))
coords.train = lapply(traindf,function(x) split(x[,c("Longitude","Latitude")],f=x$Pol))
n.train = lapply(y.train,function(x) lapply(1:q, function(i) length(x[[i]])))
do.call(rbind,n.train)

maxd = 3000 #max(rdist.earth(coords[[1]][[1]])) # distance in miles
##call spMisalignLM
n.samples <- iters
psi.init=rep(0.2,q)
nt = length(y.train)

#set up for Rcpp
nltr <- q*(q+1)/2
Aidx <- cumsum(c(1,q:2))
mAidx <-(1:nltr)[-Aidx]
AMindx <- matrix(1:q^2,q)
AMindx <- AMindx[lower.tri(AMindx,T)]

Xlist <- Ylist <- Nlist <- nlist <- coordslist<-Dlist<- coordslistUTM<-NULL
for(t in 1:nt){
  Y = NULL;n = NULL
  ytemp = y.train[[t]]; xtemp = x.train[[t]]
  for(i in 1:q) n[[i]]<- sum(!is.na(ytemp[[i]]))
  N = sum(n)
  X = matrix(0,nrow=N,ncol=length(unlist(covar.names)))
  n.p = unlist(lapply(covar.names,length))
  for(i in 1:q){
    if(n[[i]]!=0) {
      X[(cumsum(n)[i]-n[i]+1):cumsum(n)[i],(cumsum(n.p)[i]-n.p[i]+1):cumsum(n.p)[i]] = as.matrix(xtemp[[i]])
      Y = c(Y,ytemp[[i]])
    }
  }
  Ylist[[t]] = Y
  Xlist[[t]] = X
  Nlist[[t]] = N
  coordslist[[t]] = as.matrix(do.call(rbind, coords.train[[t]])); 
  coordslist[[t]] = coordslist[[t]][!is.na(coordslist[[t]][,1]),]
  # coordslistUTM[[t]] = as.matrix(do.call(rbind, coordsUTM.train[[t]]))
  # coordslistUTM[[t]] = coordslistUTM[[t]][!is.na(coordslistUTM[[t]][,1]),]
  Dlist[[t]]  = rdist.earth(coordslist[[t]])
  nlist = rbind(nlist,n)
  rm(ytemp);rm(xtemp)
}

rho_step = psi_step=0.2
rho_init = rep(3/100,q); psi_init = rep(0.1,q); beta_init = rep(1,ncol(X))
sd_beta = 100; a_psi = 2; b_psi = 0.1; SK = diag(0.1,q);
IW_df = q+1; a_rho = 3/(700); b_rho = 3/(50) # changed from 3/(0.75*maxd), 3/(0.1*maxd)
n_report=10; fix_psi=F
A_step = 0.08*c(0.028, 0.036, 0.572, 0.156, 0.174, 0.169, 0.001, 0.885, 0.199, 0.194, 0.327,
0.382, 0.062, 0.300, 0.225, 0.010, 0.354, 0.070, 0.276, 0.111, 0.008)

fooidx = c(2,7,9,12,13,16,19,21)
A_init = c( 0.51,  0.44,  1.92,  2.61,  0.83,  1.56,  0.07,  1.57,  2.50,  0.46,
            1.90,  1.02,  1.33,  0.07,  0.75,  0.31,  0.16,  0.22,  0.48, -0.08,
            0.27)

Rcpp::sourceCpp(paste0(rcppdir,'MultiSpMisalign_all_RCPP.cpp'))
set.seed(1)
tick=proc.time()
set.seed(1)
out = MultiSpMisaligncpp_all(Ylist, Xlist, Dlist, Nlist,nlist,
                             nt, q, p,(Aidx-1),(mAidx-1),(AMindx-1),
                             A_step, rho_step, psi_step,
                             rho_init,  psi_init, beta_init, A_init,
                             sd_beta,  a_psi,  b_psi,SK,
                             IW_df, a_rho, b_rho,
                             iters,  n_report, fix_psi,report)
tock=proc.time()-tick
tock
out$time <- tock

# thin it
keep =  1:iters

out$beta = out$beta[,keep]
out$A    = out$A[,keep]
out$psi  = out$psi[,keep]
out$rho  = out$rho[keep]
out$iters = iters
out$Nlist = Nlist
out$nlist = nlist
out$fold = fold
print(out$time)
out$dateset = dateset
save(out,file=paste0(resultdir,"realdat_Rcpp_batch",rep,"_sea",sea,"_nbasis",nbasis,"_P",p,"_fold",fold,".Rda"))
}