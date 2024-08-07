library(SANple)
library(parallel)

source("Univariate/A1_Hyperpar_GIBBS.R")
source("Univariate/B1_Aux_GIBBS.R")

nrep <- 10000
burn <- 5000
Ncores <- 10

sample_sizes <- c(10,50,500,2500)
sims <- 50
root <- 23442


for(ss in 1:4){
  
  simu_CAM <- simu_fSAN_of <- simu_fiSAN_of <- list()
  ################################################################################
  SAMPLE_TO_RUN <- ss
  ################################################################################
  cat(paste("*************************************************\n"))
  cat(paste("*********** SCENARIO",SAMPLE_TO_RUN," ***********\n"))
  cat(paste("*************************************************\n"))
  
  y      <- readRDS(paste0("Univariate/Datasets/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"/data_list",sims,".RDS"))
  group  <- readRDS(paste0("Univariate/Datasets/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"/group_list",sims,".RDS"))
  name   <- paste0("Univariate/Runs/GIBBS_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN])
  
  ################################################################################
  
  
  if(is.null(name)){stop("missing variable: name")}
  if(is.null(root)){stop("missing variable: root - sure you do not want a seed?")}
  root <- 9241 * ss
  
  #CAM
  cat(paste("Starting with CAM\n"))
  simu_CAM <- parallel::mclapply(1:sims,function(x){
      one_sweep_CAM(i = x,
                    Hyperpar_list = Hyperpar_list,
                    nrep = nrep,
                    burn = burn,
                    root = root,
                    y = y,
                    group = group,
                    burnin = BI )
  },mc.cores = Ncores)
  cat(paste("---- Finished running, now saving...\n"))
  saveRDS(simu_CAM, paste0(name,"_GIBBS_CAM2.RDS"))
  rm(simu_CAM)
  cat(paste(" --- Done with CAM \n"))
  cat(paste("*************************************************\n"))
  
  #fCAM
  cat(paste("Starting with fCAM\n"))
  simu_fCAM <- parallel::mclapply(1:sims,function(x){
    one_sweep_fCAM(i = x,
                  Hyperpar_list = Hyperpar_list,
                  nrep = nrep,
                  burn = burn,
                  root = root,
                  y = y,
                  group = group,
                  burnin = BI )
  },mc.cores = Ncores)
  cat(paste("---- Finished running, now saving...\n"))
  saveRDS(simu_fCAM, paste0(name,"_GIBBS_fCAM2.RDS"))
  rm(simu_fCAM)
  cat(paste(" --- Done with fCAM \n"))
  cat(paste("*************************************************\n"))
  
  # fisan of
  cat(paste("Starting with fiSAN\n"))
  simu_fiSAN_of <- parallel::mclapply(1:sims,function(x){
    one_sweep_fiSAN_of(i = x,
                  Hyperpar_list = Hyperpar_list,
                  nrep = nrep,
                  burn = burn,
                  root = root,
                  y = y,
                  group = group,
                  burnin = BI )
  },mc.cores = Ncores)
  cat(paste("---- Finished running, now saving...\n"))
  saveRDS(simu_fiSAN_of, paste0(name,"_GIBBS_fiSAN_OF2.RDS"))
  rm(simu_fiSAN_of)
  cat(paste("---- Done with fiSAN\n"))
  cat(paste("*************************************************\n"))
  
  # fsan of
  cat(paste("Starting with fSAN\n"))
  simu_fSAN_of <- parallel::mclapply(1:sims,function(x){
    one_sweep_fSAN_of(i = x,
                  Hyperpar_list = Hyperpar_list,
                  nrep = nrep,
                  burn = burn,
                  root = root,
                  y = y,
                  group = group,
                  burnin = BI )
  },mc.cores = Ncores)
  cat(paste("---- Finished running, now saving...\n"))
  saveRDS(simu_fSAN_of,paste0(name,"_GIBBS_fSAN_OF2.RDS"))
  rm(simu_fSAN_of)
  cat(paste("---- Done with fSAN\n"))
  cat(paste("*************************************************\n"))
  
}  
