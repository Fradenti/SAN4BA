library(SANplemvt)
library(parallel)

source("Multivariate/A1_MVT_Hyperpar_Gibbs.R")
source("Multivariate/B1_MVT_Aux_Gibbs.R")

nrep <- 10000
burn <- 5000
Ncores <- 10

sample_sizes <- c(50, 500, 1000) # group size: we have 6 groups of equal size, hence the total sample size ranges from 300 to 6000
dimensions <- c(2,5,10)

sims <- 50
root <- 23432

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
  name   <- paste0("Multivariate/Runs/GIBBS_RDS/samplesize",sample_sizes[SAMPLE_TO_RUN],
                   "_dim",dimensions[DIM])
  
  ################################################################################
  
  
  if(is.null(name)){stop("missing variable: name")}
  if(is.null(root)){stop("missing variable: root - sure you do not want a seed?")}
  
  print(Sys.time())
  cat("\n start the running\n")
  
  #fiSAN
  simu_fiSAN <- parallel::mclapply(1:sims,function(x){
    one_sweep_mvtfiSAN_of(i = x,
                    nrep = nrep,
                    burn = burn,
                    root = root,
                    y = y, 
                    group = group,Hyperpar_list = Hyperpar_list)
  },
  mc.cores = Ncores)
  print(Sys.time())
  cat("\n Finished running - now saving\n")
  saveRDS(simu_fiSAN,paste0(name,"_GIBBS_MVT_fiSAN.RDS"))
  rm(simu_fiSAN)
  print(Sys.time())
  cat(paste("*************************************************\n"))
  cat(paste0("Done with sample size ",ss,"\n"))
    
  }
  
  cat(paste0("Done with dimension ",dd,"\n"))
  
}
