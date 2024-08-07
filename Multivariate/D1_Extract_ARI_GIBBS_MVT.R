library(mclust)
source("Multivariate/Helpers/postprocessing_output.R")
sample_sizes <- c(50, 500, 1000) # group size: we have 6 groups of equal size, hence the total sample size ranges from 300 to 6000
dimensions <- c(2,5,10)

sims <- 50

for(dd in 1:3){
  for(ss in 1:3){
  ###############################################################################
  SAMPLE_TO_RUN=ss
  DIM = dd
  ###############################################################################
  cat(paste("*************************************************\n"))
  cat(paste("*********** SCENARIO",SAMPLE_TO_RUN," ***********\n"))
  cat(paste("*************************************************\n"))
  
  
  oc_groundtruth <- readRDS(paste0("Multivariate/Datasets/samplesize",sample_sizes[SAMPLE_TO_RUN],
                     "_dim",dimensions[DIM],"/trueOC_list50.RDS"))
  dc_groundtruth <- readRDS(paste0("Multivariate/Datasets/samplesize",sample_sizes[SAMPLE_TO_RUN],
                     "_dim",dimensions[DIM],"/trueDC_list50.RDS"))
  
  Clustering_results <- list()
  name   <- paste0("Multivariate/Runs/GIBBS_RDS/samplesize",sample_sizes[SAMPLE_TO_RUN],
                   "_dim",dimensions[DIM])
  
  simu_fiSAN <- readRDS(paste0(name,"_GIBBS_MVT_fiSAN.RDS"))
  Clustering_results$fiSAN_OCL <- parallel_OCL_extract(res = simu_fiSAN)
  Clustering_results$fiSAN_DCL <- parallel_DCL_extract(res = simu_fiSAN)
  time_fiSAN <- as.numeric(lapply(simu_fiSAN,function(x)    as.numeric(x$time,units="secs")))
  memo_fiSAN <- as.numeric(lapply(simu_fiSAN,function(x)    log(as.numeric(x$size/10^6))))
  rm(simu_fiSAN)
  
  ###############################################################################
  name_res <- paste0("Multivariate/Results/GIBBS_RDS/samplesize",sample_sizes[SAMPLE_TO_RUN],
                     "_dim",dimensions[DIM],"_GIBBS")
  ###############################################################################
  saveRDS(Clustering_results,paste0(name_res,"_estimated_clustering_results.RDS"))
  ###############################################################################
  ###############################################################################
  
  # Rand_indexes ------------------------------------------------------------
  
  Obs_rand <- matrix(NA,sims,1)
  Obs_rand[,1] <- apply(Clustering_results$fiSAN_OCL, 1, function(x) mclust::adjustedRandIndex(x, oc_groundtruth))
  colnames(Obs_rand) = c("fiSAN") 
  
  
  Dis_rand <- matrix(NA,sims,1)
  Dis_rand[,1] <- apply(Clustering_results$fiSAN_DCL, 1, function(x) mclust::adjustedRandIndex(x, dc_groundtruth))
  colnames(Dis_rand) = c("fiSAN")
  
  par(mfrow=c(1,2))
  boxplot(Obs_rand,main = "OC - Rand Indexes")
  boxplot(Dis_rand,main = "DC - Rand Indexes")
  par(mfrow=c(1,1))
  
  All_rands <- list(Dis_rand,Obs_rand)
  saveRDS(All_rands,paste0(name_res,"_Rand_indexes.RDS"))
  
  ###############################################################################
  ###############################################################################
  ###############################################################################
  
  # time --------------------------------------------------------------------
  Dis_time <- cbind(time_fiSAN)
  colnames(Dis_time) = c("fiSAN")
  boxplot(Dis_time, main = "Elapsed Seconds")
  saveRDS(Dis_time,paste0(name_res,"_elapsed_seconds.RDS"))
  
  # memo --------------------------------------------------------------------
  Dis_memo <- cbind(memo_fiSAN)
  colnames(Dis_memo) = c("fiSAN")
  boxplot(Dis_memo, main = "log_Megabytes")
  saveRDS(Dis_memo,paste0(name_res,"_log_Megabytes.RDS"))
  
  
  ###############################################################################
  cat(paste("*********************************"))
  cat(paste("done with ss",ss,'and dd',dd,"\n"))
  cat(paste("*********************************"))
  ###############################################################################
  
}
}
