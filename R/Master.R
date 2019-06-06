# Perform multivariate multirest downscaling for PM2.5
# A subset of the original data is taken
#     subsample PM station data to 100 locations
#     subset date from 7/1/2011 - 7/31/2011

# load libraries
library(parallelMCMCcombine)
library(splines)
library(fields)
library(Rcpp)
library(viridis)
library(data.table)

# set some directory ------------------------------------------------------
datadir       = "../Data/"
fndir         = ""
resultdir     = "../Result/"
rcppdir       = "../Rcpp/"
cmaq_spectdir = "../Data/spectral_cov/"
figdir        = "../Figure/"
# Create specral covariates -----------------------------------------------
nbasis = 5 # number of basis
poll_names = c("PM25","EC","NO3","NH4","OC","SO4")
q = length(poll_names)

if(!file.exists(paste0(datadir,"proc",nbasis,".RData"))) 
  source(paste0(fndir,"compute_X_tilde_BS.R"))


# Fit models --------------------------------------------------------------
# Each model is fit to 10 batches with batch size 3
# This can be parallelized easily on a cluster

# Fit spectral downscaler with cross-species
covar.names <- lapply(1:q,function(x) unlist(lapply(poll_names,function(x) paste0("Xtilde.",x,1:nbasis))))
covar.names = lapply(covar.names,function(x) c("one",x))
p = nbasis*6+1; iters = 1e3; batchnum = 9# mcmc iterations
source(paste0(fndir,"SpatialRegress_Xtilde.R"))
source(paste0(fndir,"combine_result.R"))

fit_all_model = F # Change to T to fit other models
if(fit_all_model){
# Fit spectral downscaler without cross-species
covar.names <- lapply(poll_names,function(x) paste0("Xtilde.",x,1:nbasis))
covar.names = lapply(covar.names,function(x) c("one",x))
p = nbasis+1; iters = 5e2; batchnum = 9 # mcmc iterations
source(paste0(fndir,"SpatialRegress_Xtilde.R"))
source(paste0(fndir,"combine_result.R"))

# Fit spatial downscaler with cross-species
covar.names <- lapply(1:q,function(x) paste0("X.",poll_names))
covar.names = lapply(covar.names,function(x) c("one",x))
p = 7; iters = 5e2; batchnum = 9 # mcmc iterations
source(paste0(fndir,"SpatialRegress_Xtilde.R"))
source(paste0(fndir,"combine_result.R"))

# Fit spatial downscaler without cross-species
covar.names <- lapply(poll_names,function(x) paste0("X.",x))
covar.names = lapply(covar.names,function(x) c("one",x))
p = 2; iters = 5e2; batchnum = 9 # mcmc iterations
source(paste0(fndir,"SpatialRegress_Xtilde.R"))
source(paste0(fndir,"combine_result.R"))
}

# Make prediction and forecast --------------------------------------------
# spectral downscaler with cross-species
covar.names <- lapply(1:q,function(x) unlist(lapply(poll_names,function(x) paste0("Xtilde.",x,1:nbasis))))
covar.names = lapply(covar.names,function(x) c("one",x))
p = nbasis*6+1; 
source(paste0(fndir,"SpatialRegress_Xtilde_pred.R"))

fit_all_model = F # Change to T to fit other models
if(fit_all_model){
  # Fit spectral downscaler without cross-species
  covar.names <- lapply(poll_names,function(x) paste0("Xtilde.",x,1:nbasis))
  covar.names = lapply(covar.names,function(x) c("one",x))
  p = nbasis+1; 
  source(paste0(fndir,"SpatialRegress_Xtilde_pred.R"))
  
  # Fit spatial downscaler with cross-species
  covar.names <- lapply(1:q,function(x) paste0("X.",poll_names))
  covar.names = lapply(covar.names,function(x) c("one",x))
  p = 7; 
  source(paste0(fndir,"SpatialRegress_Xtilde_pred.R"))
  
  # Fit spatial downscaler without cross-species
  covar.names <- lapply(poll_names,function(x) paste0("X.",x))
  covar.names = lapply(covar.names,function(x) c("one",x))
  p = 2; 
  source(paste0(fndir,"SpatialRegress_Xtilde_pred.R"))
}

# Make figures
source(paste0(fndir, "PlotFigures.R"))
