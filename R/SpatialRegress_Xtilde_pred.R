# setupd
source(paste0(fndir,"utils.R"))
source(paste0(fndir,'mkCov.R'))

prediction = T; forecast = T; mappredict = T
sea = "JAS"; fold = 1;nbasis = 5; batch = 3; report = 0.1*iters

# load prediction site for five fold CV
load(paste0(datadir,"predsite.Rda"))
load(paste0(datadir,"proc",nbasis,".RData")); 
predsite = PREDSITE[[fold]]
df$one = 1
df$idx = 1:nrow(df)

# scale both spectral covariates and station data
foonames = unlist(lapply(poll_names, function(pol) c(paste0("Xtilde.",pol,1:nbasis),paste0("X.",pol))))
sc = scale(df[,colnames(df)%in%foonames])
sc.center= attr(sc,"scaled:center")
sc.scale = attr(sc,"scaled:scale")

df[,colnames(df)%in%foonames] = scale(df[,colnames(df)%in%foonames])
df_bypol = split(df,df$Pol)
df_bypol = lapply(df_bypol, function(x){
  x[,"Concentration"] = scale(x[,"Concentration"])
  return(x)}
)
df = do.call(rbind, df_bypol)
predDat  = subset(df, Site %in% predsite)
trainDat = subset(df, !Site %in% predsite)

days =  range(df$Sample.Date)
set = seq(min(days),max(days),by = batch)
dateset = set[1]:(set[10]-1)
cat("Training dates",range(dateset),"\n")

# # Build training data set
traindf = subset(trainDat, Sample.Date%in%dateset)
traindf$Pol = factor(traindf$Pol, levels=poll_names)
traindf = split(traindf,traindf$Sample.Date)
y.train = lapply(traindf,function(x) split(x[,"Concentration"],f=x$Pol))
x.train = lapply(traindf,function(x) lapply(1:q, function(i) subset(x,Pol==poll_names[i])[,covar.names[[i]]]))
coords.train = lapply(traindf,function(x) split(x[,c("Longitude","Latitude")],f=x$Pol))
n.train = lapply(y.train,function(x) lapply(1:q, function(i) length(x[[i]])))

# # Build prediction data set
preddf = subset(predDat, Sample.Date%in%dateset)
preddf$Pol = factor(preddf$Pol, levels=poll_names)
preddf = split(preddf,preddf$Sample.Date)
y.pred = lapply(preddf,function(x) split(x[,"Concentration"],f=x$Pol))
x.pred = lapply(preddf,function(x) lapply(1:q, function(i) subset(x,Pol==poll_names[i])[,covar.names[[i]]]))
coords.pred = lapply(preddf,function(x) split(x[,c("Longitude","Latitude")],f=x$Pol))
n.pred = lapply(y.pred,function(x) lapply(1:q, function(i) length(x[[i]])))

# # Build forecast data set
foredays <- (max(dateset)+1):max(days)
cat("Forecast dates",range(foredays),"\n")
foredf = subset(df, Sample.Date%in%foredays)
foredf$Pol = factor(foredf$Pol, levels=poll_names)
foredf = split(foredf,foredf$Sample.Date)
y.fore = lapply(foredf,function(x) split(x[,"Concentration"],f=x$Pol))
x.fore = lapply(foredf,function(x) lapply(1:q, function(i) subset(x,Pol==poll_names[i])[,covar.names[[i]]]))
coords.fore = lapply(foredf,function(x) split(x[,c("Longitude","Latitude")],f=x$Pol))
n.fore = lapply(y.fore,function(x) lapply(1:q, function(i) length(x[[i]])))
idx.fore = lapply(foredf,function(x) split(x[,"idx"],f=x$Pol))

# # prediction ####
if(prediction){
  cat("Making spatial prediction","\n")
  load(file=paste0(resultdir,"combchain_Rcpp_sea",sea,"_nbasis",nbasis,"_P",p,"_fold",fold,".Rda"))
  predcombsave = NULL
  for(t in 1:length(dateset)){
    day = dateset[t]; predtemp = NULL
    tick= proc.time()
    predout <- predMultiSpMisalign(y = y.train[[t]],x = x.train[[t]], coords = coords.train[[t]],
                                   n=unlist(n.train[[t]]),pred.covars=x.pred[[t]],
                                   pred.coords=coords.pred[[t]],
                                   beta=t(beta_comb),A.each=t(A_comb),psi = t(psi_comb),
                                   rho = cbind(rho_comb,rho_comb,rho_comb,
                                               rho_comb,rho_comb,rho_comb),
                                   mkCov,mkVcros,mkVar,dist=rdist.earth)
    tock = proc.time()-tick
    #loc,day,outcome,stat
    for(j in 1:length(poll_names)) if(!any(is.na(coords.pred[[t]][[j]]))&nrow(coords.pred[[t]][[j]])>0) predtemp[[j]] <- 
      t(rbind(apply(predout[[j]], 2, quants),colMeans(predout[[j]])))
    
    predcombsave[[t]] <- predtemp
    print(day)
  }
  cat("Complete spatial prediction","\n")
  save(predcombsave,file=paste0(resultdir,"predcombchain_Rcpp_sea",sea,"_nbasis",nbasis,"_P",p,"_fold",fold,".Rda"))
}

# # Forecast ####
cat("Making forecast day ahead","\n")
if(forecast){
  load(file=paste0(resultdir,"combchain_Rcpp_sea",sea,"_nbasis",nbasis,"_P",p,"_fold",fold,".Rda"))
  predoutput = NULL
  for(t in 1:length(foredays)){
    predout <- forecastMultiSpMisalign(n=unlist(n.fore[[t]]),pred.covars=x.fore[[t]],pred.coords=coords.fore[[t]],
                                       beta=t(beta_comb),A.each=t(A_comb),psi = t(psi_comb),
                                       rho = cbind(rho_comb,rho_comb,rho_comb,
                                                   rho_comb,rho_comb,rho_comb),
                                       mkCov,dist=rdist.earth)
    CI = apply(predout[[1]], 2, quantile, prob = c(0.025,0.975))
    Mean = colMeans(predout[[1]])
    predoutput[[t]] = data.frame(Mean, CI[1,], CI[2,],unlist(idx.fore[[t]]), foredays[t])
    colnames(predoutput[[t]]) = c("Mean","2.5%","97.5%","idx","Sample.Date")
  }
  cat("Complete forecast","\n")
  save(predoutput,file=paste0(resultdir,"forecastcombchain_Rcpp_sea",sea,"_nbasis",nbasis,"_P",p,"_fold",fold,".Rda"))
}

# # mapprediction ####
if(mappredict){
  load(file=paste0(resultdir,"combchain_Rcpp_sea",sea,"_nbasis",nbasis,"_P",p,"_fold",fold,".Rda"))
  load(paste0(datadir,"CMAQ.RData"))
  
  mapcoords = layout[,2:3]
  mapcoords.pred = lapply(1:q,function(x) mapcoords)
  
  for(mapT in 1:length(dateset)){
    day = dateset[mapT]; predtemp = NULL
    mappredcombsave = NULL
    cat("Making map prediction for day", day,"\n")
    mappreddf <- data.frame(one=rep(1, 137241))
    
    for(pol in poll_names){
      foovar = (read.csv(file=paste0(cmaq_spectdir, pol,"/basis",nbasis,"_day",day,".csv"))[,-1])
      colnames(foovar) <- c(paste0("Xtilde.",pol,1:nbasis),paste0("X.",pol))
      mappreddf<- cbind(mappreddf, foovar)
    }
    # scale each covariates
    mappreddf[,colnames(mappreddf)%in%foonames]=sapply(foonames, function(x) (mappreddf[,x]-sc.center[x])/sc.scale[x])
    x.mappred = lapply(1:q, function(i) mappreddf[,covar.names[[i]]])
    idx = seq(1,ncol(beta_comb),l=10)
    tick= proc.time()
    mappredcombsave <- predmapMultiSpMisalign_mean(y = y.train[[mapT]],x = x.train[[mapT]], coords = coords.train[[mapT]],
                                              n=unlist(n.train[[mapT]]),pred.covars=x.mappred,
                                              pred.coords=mapcoords.pred,
                                              beta=(t(beta_comb[,idx])),A.each=(t(A_comb[,idx])),psi= (t(psi_comb[,idx])),
                                              rho = (cbind(rho_comb[idx],rho_comb[idx],rho_comb[idx],
                                                           rho_comb[idx],rho_comb[idx],rho_comb[idx])),
                                              mkCov,mkVcros,mkVar,dist=rdist.earth)
    tock = proc.time()-tick
    cat("Complete map prediction for day", day,"\n")
    
    save(mappredcombsave,file=paste0(resultdir,"mappredcombchain_Rcpp_sea",sea,"Day",day,"_nbasis",nbasis,"_P",p,"_fold",fold,".Rda"))
  }
}