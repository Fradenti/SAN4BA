library(SANvimvt)
library(parallel)
library(crayon)
source("Multivariate/A2_MVT_Hyperpar_VI.R")

sample_sizes <- c(50, 500, 1000) # group size: we have 6 groups of equal size, hence the total sample size ranges from 300 to 6000
dimensions <- c(2,5,10)

sims <- 50
root <- 23442

for(dd in 1:3){
  for(ss in 1:3){
    ################################################################################
    SAMPLE_TO_RUN <- ss
    DIM <- dd
    ################################################################################
    cat(paste("*************************************************\n"))
    cat(paste("*********** SCENARIO",SAMPLE_TO_RUN," ***********\n"))
    cat(paste("*************************************************\n"))
    
    y      <- readRDS(paste0("Multivariate/Datasets/samplesize",sample_sizes[SAMPLE_TO_RUN],
                             "_dim",dimensions[DIM],"/data_list",sims,".RDS"))
    group  <- readRDS(paste0("Multivariate/Datasets/samplesize",sample_sizes[SAMPLE_TO_RUN],
                             "_dim",dimensions[DIM],"/group_list",sims,".RDS"))
    name   <- paste0("Multivariate/Runs/VI_RDS/samplesize",sample_sizes[SAMPLE_TO_RUN],
                     "_dim",dimensions[DIM])
    
    ################################################################################
    
  
  if(is.null(name)){stop("missing variable: name")}
  if(is.null(root)){stop("missing variable: root - sure you do not want a seed?")}
  
  source("Multivariate/B2_MVT_Aux_RUN_VI.R")
  
}
}
