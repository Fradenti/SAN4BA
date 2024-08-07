library(tidyverse)
library(parallel)
library(mclust)
Rcpp::sourceCpp("Univariate/Helpers/estimate_postdens.cpp")
source("Univariate/Helpers/postprocessing_output.R")
sample_sizes = c(10,50,500,2500)#0)
##################################
##################################
ncores = 10
sims   = 50

for(SAMPLE_TO_RUN in 1:4){
  
  cat(paste("*************************************************\n"))
  cat(paste("*********** SCENARIO",SAMPLE_TO_RUN," ***********\n"))
  cat(paste("*************************************************\n"))
  
  name <- paste0("Univariate/Runs/GIBBS_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_GIBBS")
  
  simu_CAM   <- readRDS(paste0(name,"_CAM.RDS"))
  dens_cam <- mclapply(1:sims,function(z) 
    post_mean_dens_i(i = z, interval = c(-10,10),
                     object = simu_CAM,
                     type = "CAM"),
    mc.cores = ncores)
  rm(simu_CAM)
  cat("\n Done with CAM!")
  
  simu_fCAM   <- readRDS(paste0(name,"_fCAM.RDS"))
  dens_fcam <- mclapply(1:sims,function(z) 
    post_mean_dens_i(i = z, interval = c(-10,10),
                     object = simu_fCAM,
                     type = "fCAM"),
    mc.cores = ncores)
  rm(simu_fCAM)
  cat("\n Done with fCAM!")
  
  simu_fSAN_OF   <- readRDS(paste0(name,"_fSAN_OF.RDS"))
  dens_fsan <- mclapply(1:sims,function(z) 
    post_mean_dens_i(i = z, interval = c(-10,10),
                     object = simu_fSAN_OF,
                     type = "fSAN_OF"),
    mc.cores = ncores)
  rm(simu_fSAN_OF)
  
  cat("\n Done with fSAN_OF!")
  
  simu_fiSAN_OF <- readRDS(paste0(name,"_fiSAN_OF.RDS"))
  dens_fisan <- mclapply(1:sims,function(z) 
    post_mean_dens_i(i = z, interval = c(-10,10),
                     object = simu_fiSAN_OF,type = "fiSAN_OF"),
    mc.cores = ncores)
  rm(simu_fiSAN_OF)
  
  cat("\n Done with fiSAN_OF!\n")
  
  dd1 = do.call(rbind,dens_cam)
  dd2 = do.call(rbind,dens_fcam)
  dd3 = do.call(rbind,dens_fisan)
  dd4 = do.call(rbind,dens_fsan)
  
  all_densities = rbind(dd1,dd2,dd3,dd4)
  name_res <- paste0("Univariate/Results/GIBBS_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_GIBBS")
  saveRDS(all_densities, paste0(name_res,"_posterior_mean_dens_over_grid_GIBBS.RDS"))
  cat(SAMPLE_TO_RUN)
  
}

