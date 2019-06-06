# plot coefficients by km -------------------------------------------------
load(paste0(datadir,"proc",nbasis,".RData"))
source(paste0(fndir,"utils.R"))

label    = c(60,120,240,600)
delta    = seq(12*2*pi/600,12*2*pi/48,l=50)
period   = 12*2*pi/delta
xtick    = log10(period)
cols     = tim.colors(6)
cols[4]  ="goldenrod"
ylimsave = NULL

foonames = unlist(lapply(poll_names,function(x) paste0("Xtilde.",x,1:nbasis)))
foo      = scale(df[,colnames(df)%in%foonames])
center   = attr(foo,"scaled:center")
scale    = attr(foo,"scaled:scale")
bmat     = bs(delta,nbasis,intercept=T)
yscale   = unlist(lapply(poll_names, function(x) sd(subset(df, Pol==x)$Concentration)))
xscale   = apply(df[,paste0("X.",poll_names)],2,sd)

{pdf(paste0(figdir,"Figure 6.pdf"))
  load(paste0(resultdir,"combchain_Rcpp_seaJAS_nbasis",nbasis,"_P",p,"_fold1.Rda"))
  for(pol in 1:6){
    betaoutput = lapply(1:30, function(x) beta_comb[1:p+(pol-1)*p,][x+1,]*yscale[pol])#*scale[x])
    betaoutput = do.call(rbind,betaoutput)
    betaUmatrix = betaLmatrix = betamatrix = matrix(NA,nrow =6, nc=length(delta))
    polj = pol
    betaoutputj = betaoutput[1:nbasis+(polj-1)*nbasis,]*xscale[polj]
    btmean = rowMeans(betaoutputj)
    betamatrix[polj,] =  bmat%*%btmean
    foo1 = apply(betaoutputj,2,function(x) bmat%*%x)
    CI = apply(foo1,1,quants)
    betaLmatrix[polj,] <- CI[1,]
    betaUmatrix[polj,] <- CI[3,]
    ylim = range(rbind(betaUmatrix,betaLmatrix),na.rm=T)
    par(mar=c(5.1, 5.1, 5.1, 2.1), xpd=TRUE)
    plot(log10(period),betamatrix[polj,],xlab="Period (km)", xlim = log10(c(48,600)),
         ylim = ylim,ylab="",type="l",xaxt="n",main=poll_names[pol])
    axis(1,at=log10(label),labels=label)
    legend("bottom",legend = poll_names[polj], horiz=T,col=cols[polj],lty=1,lwd=2,cex=0.8)
    lines(log10(period),rep(0,length(period)),lty=2,lwd=2)
    polygon(c(xtick, rev(xtick)),
            c(betaLmatrix[polj,], rev(betaUmatrix[polj,])),
            col=alpha(cols[polj],0.2),border = NA)
  }
  dev.off()}
