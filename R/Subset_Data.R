#  ------------------------------------------------------------------------
# Create subset of the original data for code submission
# Sample 100 out of total 830 PM stations
# Subset date July 1 - Aug 
#  ------------------------------------------------------------------------
rm(list=ls())
library(data.table)

obsdir  <- "../Observation 2011/"
cmaqdir <- "../cmaq/"
poll_names <- c("PM25","EC","NO3","NH4","OC","SO4")

#load station data
df <- data.frame()
for(pol in poll_names){
  Y = read.csv(paste0(obsdir,"2011 ",pol,".csv"))
  Y$Pol = pol
  df = rbind(df,Y)
}
rm(Y)
dim(df)

#sample stations
df_pm  = subset(df,Pol =="PM25")
pmsite = sample(unique(df_pm$Site),100)
df_pm  = subset(df_pm, Site%in%pmsite)

df_species  = subset(df,Pol !="PM25")
df     = rbind(df_pm,df_species)
#subset by date
dayrange = as.numeric(as.Date(c("7/1/2011","8/1/2011"),"%m/%d/%Y")-as.Date("2011-01-01"))+1

Date <- df$Sample.Date
date     <- as.numeric(as.Date(df$Sample.Date,"%m/%d/%Y")-as.Date("2011-01-01"))+1
df$Sample.Date <- date
df$Orig.Date <- Date
df = subset(df, Sample.Date %in% dayrange[1]:dayrange[2])

# match station location to cmaq grid
layout = read.csv(paste0(cmaqdir,"CMAQ.layout.csv"))
locid = unique(df[,c("Longitude","Latitude")])
for(i in 1:nrow(locid)){
  cellidx = which.min(rdist.earth(locid[i,],layout[,c("Longitude","Latitude")]))
  df$cell[df$Longitude == locid$Longitude[i] & df$Latitude == locid$Latitude[i]] <- cellidx
}
AQdata = df
save(AQdata, file="AQdata.Rdata")

# subset cmaq data
for(pol in poll_names){
  VAR = as.matrix(fread(paste0(cmaqdir,"CMAQ 2011 ",pol,".csv")))
  VAR = VAR[,dayrange[1]:dayrange[2]]
  assign(paste0("VAR_",pol),VAR)
}
save(VAR_SO4, file="CMAQdata_SO4.Rdata")
save(VAR_EC, file="CMAQdata_EC.Rdata")
save(VAR_NH4, file="CMAQdata_NH4.Rdata")
save(VAR_NO3, file="CMAQdata_NO3.Rdata")
save(VAR_OC, file="CMAQdata_OC.Rdata")
save(VAR_PM25, file="CMAQdata_PM25.Rdata")
