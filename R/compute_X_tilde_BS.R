#  ------------------------------------------------------------------------
# Create spectral covariates using K number of basis function
# 1. fft to get coef
# 2. filter freq (using B-spline basis)
# 3. inverse fft to get spectral cov
# 4. match grid cell to station location
# saved as file=paste0("proclog",K,".Rda"
#  ------------------------------------------------------------------------
library(splines)
library(fields)

poll_names<- c("PM25","EC","NO3","NH4","OC","SO4")
K <- nbasis #Number of basis function

# create directory to save all spectral covariates
for(pol in poll_names){
  ifelse(!dir.exists(paste0(cmaq_spectdir, pol)), dir.create(paste0(cmaq_spectdir, pol),recursive = T),FALSE)
}

#load station data
load(paste0(datadir,"AQdata.RData"))
df = AQdata
#log transform station data
df$Concentration <- log(df$Concentration)

#Load the CMAQ layout
load(paste0(datadir,"CMAQ.RData")) #"cmaq"   "us.map" "x.proj" "y.proj" "layout"
#dimensions CMAQ
m1 <- dim(cmaq)[1]
m2 <- dim(cmaq)[2]

# create spectral covariates ----------------------------------------------
load(paste0(datadir,"CMAQdata_SO4.Rdata"))
load(paste0(datadir,"CMAQdata_EC.Rdata"))
load(paste0(datadir,"CMAQdata_NH4.Rdata"))
load(paste0(datadir,"CMAQdata_NO3.Rdata"))
load(paste0(datadir,"CMAQdata_OC.Rdata"))
load(paste0(datadir,"CMAQdata_PM25.Rdata"))
days  = unique(df$Sample.Date)
dates = unique(df$Orig.Date)

for(pol in poll_names){
  VAR = get(paste0("VAR_",pol))
  VAR = log(VAR) # changed to log scale
  cat(range(VAR),"\n")
  
  Xtildekeep <- matrix(NA,dim(df)[1],K)
  Xkeep <- matrix(NA,dim(df)[1],1)
  
  for(dayidx in 1:length(days)){
    cmaq = matrix(NA,nrow=m1,ncol=m2,byrow=T)
    cmaq[cbind(layout[,4],layout[,5])] = VAR[,colnames(VAR)==dates[dayidx]]
    scaling <-m1*m2
    
    #Fourier frequencies
    nvec = c(m1,m2)
    w1 <- 2*pi*as.matrix(c(seq(0, ceiling((nvec[1] - 1)/2)), 
                           seq(floor((nvec[1] - 1)/2), 1, by = -1)))/nvec[1]
    w2 <- 2*pi*as.matrix(c(seq(0, ceiling((nvec[2] - 1)/2)),
                           seq(floor((nvec[2] - 1)/2), 1, by = -1)))/nvec[2]
    w1 <- w1[, rep(1, nvec[2])]
    w2 <- t(w2[, rep(1, nvec[1])])
    delta  <- sqrt(w1^2+w2^2)
    
    #fft
    z <- fft(cmaq,inverse=TRUE)

    #Create the constructed covariates
    Xtilde <- matrix(NA,nrow=m1*m2,ncol=K)
    foo = bs(delta,K,intercept = T)
    for(k in 1:K){
      x      <- fft(z*foo[,k])
      xtilde <- Re(x)/scaling
      Xtilde[,k] = xtilde[cbind(layout[,4],layout[,5])]
    }

    write.csv(cbind(Xtilde,VAR[,colnames(VAR)==dates[dayidx]]),
              file=paste0(cmaq_spectdir, pol,"/basis",K,"_day",days[dayidx],".csv"))

    # extract grid cells corresponding to stations
    extidx = df$Sample.Date == days[dayidx] # get day
    Xtildekeep[extidx,] = Xtilde[df[extidx,"cell"],] # get cells for the day
    Xkeep[extidx,] = VAR[df[extidx,"cell"],colnames(VAR)==dates[dayidx]]
  }
  colnames(Xtildekeep) <- paste0("Xtilde.",pol,1:K)
  colnames(Xkeep) <- paste0("X.",pol)
  df = cbind(df,Xtildekeep,Xkeep)
}
save(df,file=paste0(datadir,"proc",K,".RData"))
