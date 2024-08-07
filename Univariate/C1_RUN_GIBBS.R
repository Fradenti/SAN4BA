library(SANple)
library(parallel)

source("Univariate/A1_Hyperpar_GIBBS.R")
source("Univariate/B1_Aux_GIBBS.R")

nrep <- 10000
burn <- 5000
Ncores <- 1

sample_sizes <- c(10,50,500,2500)
sims <- 50
root <- 23442


for(ss in 1:4){
  
  simu_CAM <- simu_fCAM <- simu_fSAN_of <- simu_fiSAN_of <- list()
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
  for(i0 in 1:sims){
  simu_CAM[[i0]] <-
    one_sweep_CAM(i = i0,
                  Hyperpar_list = Hyperpar_list,
                  nrep = nrep,
                  burn = burn,
                  root = root,
                  y = y,
                  group = group,
                  burnin = BI )
  cat(i0)
  }
  cat(paste("--- Finished running CAM - Now Saving\n"))
  saveRDS(simu_CAM, paste0(name,"_GIBBS_CAM.RDS"))
  rm(simu_CAM)
  cat(paste(" --- Done with CAM \n"))
  cat(paste("*************************************************\n"))

  #fCAM
  cat(paste("Starting with fCAM\n"))
  for(i1 in 1:sims){
    simu_fCAM[[i1]] <-
      one_sweep_fCAM(i = i1,
                    Hyperpar_list = Hyperpar_list,
                    nrep = nrep,
                    burn = burn,
                    root = root,
                    y = y,
                    group = group,
                    burnin = BI )
    cat(i1)
  }
  cat(paste("--- Finished running fCAM - Now Saving\n"))
  saveRDS(simu_fCAM, paste0(name,"_GIBBS_fCAM.RDS"))
  rm(simu_fCAM)
  cat(paste(" --- Done with fCAM \n"))
  cat(paste("*************************************************\n"))
  
  
  # fisan of
  cat(paste("Starting with fiSAN\n"))
  for(i2 in 1:sims){
  simu_fiSAN_of[[i2]] <- 
    one_sweep_fiSAN_of(i = i2,
                       Hyperpar_list = Hyperpar_list,
                       nrep = nrep, burn = burn,
                       root = root,
                       y = y,
                       group = group,
                       burnin = BI )
  cat(i2)
  }
  cat(paste("--- Finished running fiSAN - Now Saving\n"))
  saveRDS(simu_fiSAN_of, paste0(name,"_GIBBS_fiSAN_OF.RDS"))
  rm(simu_fiSAN_of)
  cat(paste(" --- Done with fiSAN \n"))
  cat(paste("*************************************************\n"))
  
  cat(paste("Starting with fSAN\n"))
  # fsan of
  for(i3 in 1:sims){
    simu_fSAN_of[[i3]] <-   one_sweep_fSAN_of(i = i3,
                      Hyperpar_list = Hyperpar_list,
                      nrep = nrep,burn = burn,
                      root = root,
                      y = y, group = group,
                      burnin = BI )
    cat(i3)
  }
  cat(paste("--- Finished running fSAN - Now Saving\n"))
  saveRDS(simu_fSAN_of,paste0(name,"_GIBBS_fSAN_OF.RDS"))
  rm(simu_fSAN_of)
  cat(paste("---- Done with fSAN\n"))
  cat(paste("*************************************************\n"))
  
}  
