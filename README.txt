Example data and code for spectral downscaler.

Note that the example data is subset to dates from July 1-July 30, 2011
and 100 PM stations locations are sampled from the total 840 stations. 

Therefore the analysis based on the example data can be different from 
the manuscript, which is based on the data for the entire 2011. 

To perform the analysis, source("Master.R") in R. 
This will perform the following
  1. create spectral covariates
  2. fit spectral downscaler with cross-species for 9 batches with 3 days in a batch
  3. combine MCMC samples from subposteriors using parallelMCMCcombine
  4. make spatial prediction and forecast at hold out locations. And make map predicton for the contiguous US. 
  5. Make figures 
