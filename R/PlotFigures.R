# plot Figure 1
source(paste0(fndir,"utils.R"))
load(paste0(datadir,"proc",nbasis,".RData"))

pdf(paste0(figdir,"Figure 1.pdf"))
par(mar=c(0,1,0,0))
PMloc = unique(df[df$Pol=="PM25",c("Longitude","Latitude")])
Speciesloc = unique(df[df$Pol%in%poll_names[-1],c("Longitude","Latitude")])
id1 = match(Speciesloc[,1],PMloc[,1])
id2 = match(Speciesloc[,2],PMloc[,2])
overlapid <- id1 == id2
PM_SP = Speciesloc[which(overlapid==TRUE),]
Speciesloconly = Speciesloc[-which(overlapid==TRUE),]
PMonly = PMloc[-na.omit(id1),]
xlim = range(df$Longitude); ylim = range(df$Latitude)
cols  <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","magenta")
plot(PMonly,xlim = xlim, ylim = ylim,col = cols[3], pch = 3,axes=F,xlab="",ylab="")
points(Speciesloconly,xlim = xlim, ylim = ylim,col = cols[2], pch = 0)
points(PM_SP,xlim = xlim, ylim = ylim,col = cols[9], pch = 2)
US(add=T,col="black")
legend("bottomleft",legend=c("Total PM2.5","Speciated PM2.5", "Both"),pch = c(3,0,2),col=cols[c(3,2,9)],bty="n")
dev.off()


# Figure 2 ----------------------------------------------------------------
pdf(paste0(figdir,"Figure 2.pdf"),width = 6.5, height = 4.5)
par(mar=c(0,1,2,1))
load(paste0(datadir,"CMAQ.RData"))
m1 <- dim(cmaq)[1];m2 <- dim(cmaq)[2]

VAR = read.csv(paste0(cmaq_spectdir,"EC","/basis",nbasis,"_day",182,".csv"))
cmaq = matrix(VAR[,7],ncol=m2,byrow=T)
image.plot(x.proj,y.proj,t(apply(cmaq,1,rev)),main="CMAQ EC (log) on Jul 1, 2011",axes=FALSE,xlab = "",ylab="",
           zlim= c(-7,2.5),col = viridis(64),xlim = range(us.map[,1],na.rm=T),ylim = range(us.map[,2],na.rm=T))
lines(us.map)

VAR = read.csv(paste0(cmaq_spectdir,"EC","/basis",nbasis,"_day",183,".csv"))
cmaq = matrix(VAR[,7],ncol=m2,byrow=T)
image.plot(x.proj,y.proj,t(apply(cmaq,1,rev)),main="CMAQ EC (log) on Jul 2, 2011",axes=FALSE,xlab = "",ylab="",
           zlim= c(-7,2.5),col = viridis(64),xlim = range(us.map[,1],na.rm=T),ylim = range(us.map[,2],na.rm=T))
lines(us.map)

VAR = read.csv(paste0(cmaq_spectdir,"PM25","/basis",nbasis,"_day",182,".csv"))
cmaq = matrix(VAR[,7],ncol=m2,byrow=T)
image.plot(x.proj,y.proj,t(apply(cmaq,1,rev)),main="CMAQ PM25 (log) on Jul 1, 2011",axes=FALSE,xlab = "",ylab="",
           zlim=  c(-1.5,4.5),col = viridis(64),xlim = range(us.map[,1],na.rm=T),ylim = range(us.map[,2],na.rm=T))
lines(us.map)

VAR = read.csv(paste0(cmaq_spectdir,"PM25","/basis",nbasis,"_day",183,".csv"))
cmaq = matrix(VAR[,7],ncol=m2,byrow=T)
image.plot(x.proj,y.proj,t(apply(cmaq,1,rev)),main="CMAQ PM25 (log) on Jul 2, 2011",axes=FALSE,xlab = "",ylab="",
           zlim=  c(-1.5,4.5),col = viridis(64),xlim = range(us.map[,1],na.rm=T),ylim = range(us.map[,2],na.rm=T))
lines(us.map)


xlim = range(df$Longitude); ylim = range(df$Latitude)
foo = subset(df, Sample.Date == 182 & Pol == "PM25")
cols = viridis(64)
plotpoints(foo[,c("Longitude","Latitude")],xlim = xlim, ylim = ylim,
           foo$Concentration,cex = 1, alpha = 1,legend=T,zrange = c(-1.5,4.5),
           main="Subsampled station PM25 (log) on Jul 1, 2011")
US(add=T)
par(mar=c(0,1,2,1))
foo = subset(df, Sample.Date == 183 & Pol == "PM25")
plotpoints(foo[,c("Longitude","Latitude")],xlim = xlim, ylim = ylim,
           foo$Concentration,cex = 1, alpha = 1,legend=T,zrange = c(-1.5,4.5),
           main="Subsampled station PM25 (log) on Jul 2, 2011")
US(add=T)
par(mar=c(0,1,2,1))
foo = subset(df, Sample.Date == 182 & Pol == "EC")
plotpoints(foo[,c("Longitude","Latitude")],xlim = xlim, ylim = ylim,
           foo$Concentration,cex = 1, alpha = 1,legend=T,zrange = c(-7,2.5),
           main="Station EC (log) on Jul 1, 2011")
US(add=T)
par(mar=c(0,1,2,1))
foo = subset(df, Sample.Date == 183 & Pol == "EC")
plotpoints(foo[,c("Longitude","Latitude")],xlim = xlim, ylim = ylim,
           foo$Concentration,cex = 1, alpha = 1,legend=T,zrange = c(-7,2.5),
           main="Station EC (log) on Jul 2, 2011")
US(add=T)
dev.off()


# Figure 3 ----------------------------------------------------------------
pdf(paste0(figdir,"Figure 3.pdf"),width = 12, height = 4)
par(mar=c(0,0,3,1))
VAR = read.csv(paste0(cmaq_spectdir,"PM25","/basis",nbasis,"_day",182,".csv"))[,7]
cmaq = matrix(NA,nrow=max(layout[,4]),ncol=max(layout[,5]),byrow=T)
cmaq[cbind(layout[,4],layout[,5])] = VAR

#Fourier frequencies
nvec = c(m1,m2)
w1 <- 2*pi*as.matrix(c(seq(0, ceiling((nvec[1] - 1)/2)), 
                       seq(floor((nvec[1] - 1)/2), 1, by = -1)))/nvec[1]
w2 <- 2*pi*as.matrix(c(seq(0, ceiling((nvec[2] - 1)/2)),
                       seq(floor((nvec[2] - 1)/2), 1, by = -1)))/nvec[2]
w1 <- w1[, rep(1, nvec[2])]
w2 <- t(w2[, rep(1, nvec[1])])
delta  <- sqrt(w1^2+w2^2)
scaling = m1*m2

#fft
z <- fft(cmaq,inverse=TRUE)

#Create the constructed covariates
Xtilde <- matrix(NA,nrow=m1*m2,ncol=10)
freqseg = seq(0,2*pi,by=pi/5)
# pdf(paste0("PM25filter10.pdf"),width = 7.5, height = 4.5)
# pdf(paste0("PM25filter10_zoom_fixed.pdf"),width = 12, height = 4)
for(k in 1:(length(freqseg)-1)){
  W           <- delta
  W[which(delta < freqseg[k] | delta >= freqseg[k+1])]=0
  W[which(delta >= freqseg[k] & delta < freqseg[k+1])]=1
  x           <- fft(z*W)
  xtilde <- Re(x)/scaling
  Xtilde[,k] =  xtilde[cbind(layout[,4],layout[,5])]
}

# # zoom to nc region
foox = matrix(layout$Lon,nrow=max(layout[,4]),ncol=max(layout[,5]),byrow=T)
idx = foox <= -75.5 & foox >= -84.5
IDX.X = apply(idx,1 ,any)
fooy = matrix(layout$Lat,nrow=max(layout[,4]),ncol=max(layout[,5]),byrow=T)
idx = fooy <= 36.6 & fooy >= 34
IDX.Y = apply(idx,2 ,any)

image.plot(matrix(layout$Lon,nrow=max(layout[,4]),ncol=max(layout[,5]),byrow=T)[IDX.X,IDX.Y],
           matrix(layout$Lat,nrow=max(layout[,4]),ncol=max(layout[,5]),byrow=T)[IDX.X,IDX.Y],
           cmaq[IDX.X,IDX.Y],col = viridis(64),ylim = c(34,36.6),xlim = c(-84.5,-75.5),
           zlim = range(cmaq[IDX.X,IDX.Y]),main="CMAQ PM2.5 (log) on Jul 1, 2011",
           xlab="",ylab="",axes=T)
US(add=T)

for(k in c(1,2,5,8)){
  foo = matrix(Xtilde[,k],nrow=max(layout[,4]),ncol=max(layout[,5]),byrow=T)
  image.plot(matrix(layout$Lon,nrow=max(layout[,4]),ncol=max(layout[,5]),byrow=T)[IDX.X,IDX.Y],
             matrix(layout$Lat,nrow=max(layout[,4]),ncol=max(layout[,5]),byrow=T)[IDX.X,IDX.Y],
             foo[IDX.X,IDX.Y],col = viridis(64),ylim = c(34,36.6),xlim = c(-84.5,-75.5),
             zlim = range(foo[IDX.X,IDX.Y]),main=paste0("CMAQ PM2.5 Xtilde",k),
             xlab="",ylab="",axes=T)
  US(add=T)
}
dev.off()


# Figure 4 ----------------------------------------------------------------
pdf(paste0(figdir,"Figure 4.pdf"),width = 6.5, height = 4.5)
par(mar=c(0,0,3,1))
yscale = lapply(poll_names, function(x) return(c(m=mean(subset(df, Pol==x)$Concentration),
                                                 sd=sd(subset(df, Pol==x)$Concentration))))
yscale = do.call(rbind,yscale)

subdf = subset(df,Sample.Date==183&Pol=="PM25")
plotpoints(subdf[,c("Longitude","Latitude")],cex=2,alpha=1,exp(subdf$Concentration),add=F,legend=T,zrange=c(0,35))
points(subdf[,"Longitude"],subdf[,"Latitude"],cex=2,col="black")
US(add=T,lwd=2)

load(file=paste0(resultdir,"mappredcombchain_Rcpp_seaJASDay183_nbasis5_P31_fold1.Rda"))
cmaq = matrix(NA,nrow=max(layout[,4]),ncol=max(layout[,5]),byrow=T)
cmaq[cbind(layout[,4],layout[,5])] = exp(mappredcombsave[[1]]*yscale[1,"sd"]+yscale[1,"m"])
image.plot(x.proj,y.proj,t(apply(cmaq,1,rev)),axes=FALSE,xlab = "",ylab="",
           main="Predictied PM2.5 on Jul 2, 2011",col = viridis(64),
           xlim = range(us.map[,1],na.rm=T),ylim = range(us.map[,2],na.rm=T),
           zlim= c(0,35))
lines(us.map,lwd=2)

VAR = read.csv(paste0(cmaq_spectdir,"PM25","/basis",nbasis,"_day",183,".csv"))[,7]
cmaq[cbind(layout[,4],layout[,5])] = exp(VAR)
image.plot(x.proj,y.proj,t(apply(cmaq,1,rev)),axes=FALSE,xlab = "",ylab="",
           main="Predictied PM2.5 on Jul 2, 2011",col = viridis(64),
           xlim = range(us.map[,1],na.rm=T),ylim = range(us.map[,2],na.rm=T))
lines(us.map,lwd=2)
dev.off()


# Figure 5 ----------------------------------------------------------------
# plot results from model fit
# plot prediction results by region
library(rgeos)
library(sp)
library(dplyr)
library(rgdal)
library(ggplot2)
library(raster)          
load(paste0(datadir,"proc5.RData"))
load(paste0(datadir,"CMAQ.RData"))
yscale = lapply(poll_names, function(x) return(c(m=mean(subset(df, Pol==x)$Concentration),
                                                 sd=sd(subset(df, Pol==x)$Concentration))))
yscale = do.call(rbind,yscale)

# The palette with grey:
mycol <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
df$idx = 1:nrow(df)
sitedf = uniquesite = unique(df[,c("Site","Longitude","Latitude")])
uniquesite$siteid = 1:nrow(uniquesite)
uniquesite$region <- NA
coordinates(uniquesite) <- ~ Longitude + Latitude

us<-getData('GADM', country='USA', level=1)  #Get the County Shapefile for the US
# Set the projection of the SpatialPointsDataFrame using the projection of the shapefile
proj4string(uniquesite) <- proj4string(us)

# The five regions in cmaq maps
northeast <- c("Maine", "New Hampshire", "Vermont", "Massachusetts",
               "New York", "New Jersey", "Maryland", "Delaware", "Connecticut",
               "Rhode Island", "Pennsylvania", "District of Columbia", "Virginia", "West Virginia")
great_lakes <- c("Ohio", "Michigan", "Indiana","Illinois","Wisconsin")
atlantic <- c("North Carolina", "South Carolina","Georgia", "Florida")
south <- c("Kentucky", "Tennessee", "Mississippi", "Alabama", "Louisiana", "Missouri", "Oklahoma","Arkansas")
west <- c("California", "Oregon", "Washington", "Arizona", "Nevada", "New Mexico")
regionid <- list(northeast,great_lakes,atlantic,south,west)
regionnames <- c("Northeast","Great Lakes","Atlantic","South","West")

for(idx in 1:length(regionnames)){
  subdat <-subset(us,NAME_1%in%regionid[[idx]])
  foo = (over(uniquesite,subdat))
  uniquesite$region[which(!is.na(foo[,1]))] <- regionnames[idx]
}

sitedf$siteid = uniquesite$siteid
sitedf$region = uniquesite$region
df <- merge(df,uniquesite)

cmaqdf = layout
layout$region <- NA
coordinates(layout) <- ~ Longitude + Latitude
# Set the projection of the SpatialPointsDataFrame using the projection of the shapefile
proj4string(layout) <- proj4string(us)

for(idx in 1:length(regionnames)){
  subdat <-subset(us,NAME_1%in%regionid[[idx]])
  foo = (over(layout,subdat))
  layout$region[which(!is.na(foo[,1]))] <- regionnames[idx]
}
cmaqdf$region = layout$region
pred = do.call(rbind,lapply(poll_names,function(x) cbind(cmaqdf,x)))
colnames(pred)[7] <- "Pol"
sea = "JAS"
sum = 0
numday = NULL
for(day in 182:208){
  file=paste0(resultdir,"mappredcombchain_Rcpp_sea",sea,"Day",day,"_nbasis",nbasis,"_P",p,"_fold1.Rda")
  if(file.exists(file)){
    load(file)
    x=unlist(lapply(1:6, function(x) exp(mappredcombsave[[x]]*yscale[x,"sd"]+yscale[x,"m"])))
    sum = sum + x
    cat(day)
    numday = c(numday,day)
  }
}
pred$Concentration = sum/length(numday)

sum = 0
for(day in numday){
  foovar = lapply(poll_names,function(pol) as.matrix(fread(file=paste0(cmaq_spectdir, pol,"/basis5_day",day,".csv")))[,7])
  x=unlist(lapply(foovar,exp))
  sum = sum + x
}
pred$cmaq = sum/length(numday)

# average mapprediciton and cmaq
pred$Pol <- factor(pred$Pol,poll_names)
foo <- pred %>% 
  group_by(region,Pol) %>%
  summarise(mean = mean(Concentration),cmaqmean = mean(cmaq))

# average dataset
df$Pol <- factor(df$Pol,poll_names)
foo1 <- df %>%
  filter(Sample.Date%in%numday) %>%
  group_by(region,Pol) %>%
  summarise(mean = mean(exp(Concentration)))
foo$mean2 = foo1$mean
assign(paste0("foo",sea),foo)


pdf(paste0(figdir,"Figure 5.pdf"),height=3, width=13)
par(mfrow = c(1,6),mar=c(2,4,3,1))
for(name in regionnames){
  foo = get(paste0("foo",sea))
  tab = subset(foo,region==name)
  tab = as.matrix(tab[,3:5])
  tab[1,] <- tab[1,] - colSums(tab[-1,])
  if(name==regionnames[1]){
    barplot(tab,main=name,
            names.arg = c("SD + Cross","CMAQ","Station"),ylab="Concentration",col=(mycol[1:6]),cex.axis=1.5)
  }else{
    barplot(tab,main=name,
            names.arg = c("SD + Cross","CMAQ","Station"),col=(mycol[1:6]),cex.axis=1.5)
  }
}
plot(NA, xlim = c(1,5),ylim=c(1,3),xlab="",ylab="",axes=F)
legend(1,3,legend=rev(c("Other",poll_names[-1])), fill=rev(mycol[1:6]),cex = 1.5)
dev.off()



# Figure 6  ---------------------------------------------------------------
# plot coefficients by km -------------------------------------------------
load(paste0(datadir,"proc",nbasis,".RData"))
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

