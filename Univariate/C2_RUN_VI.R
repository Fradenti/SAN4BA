library(tidyverse)
library(crayon)
library(SANvi)
source("Univariate/A2_Hyperpar_VI.R")

################################################################################
step <-  50 # every step iterations, I save the results
runs <-  50
sample_sizes <- c(10,50,500,2500)
sims <- 50

for(SAMPLE_TO_RUN in 1:4){
  gc()
  cat(paste("*************************************************\n"))
  cat(paste("*********** SCENARIO",SAMPLE_TO_RUN," ***********\n"))
  cat(paste("*************************************************\n"))
  
  
  name <- paste0("Univariate/Runs/VI_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN])
  
  y      <- readRDS(paste0("Univariate/Datasets/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"/data_list",sims,".RDS"))
  group  <- readRDS(paste0("Univariate/Datasets/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"/group_list",sims,".RDS"))
  
  
  if(is.null(name)){stop("missing variable: name")}
  if(is.null(root)){stop("missing variable: root - sure you do not want a seed?")}
  
  source("Univariate/B2_Aux_RUN_VI.R")
  
}

