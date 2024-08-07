library(mclust)
Rcpp::sourceCpp("Univariate/Helpers/estimate_postdens.cpp")
source("Univariate/Helpers/postprocessing_output.R")
sample_sizes = c(10,50,500,2500)#5000)
sims <- 50
ss   <- length(sample_sizes)
Ncores <- 10

for(SAMPLE_TO_RUN in 1:ss){
  ###############################################################################
  name <- paste0("Univariate/Runs/GIBBS_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_GIBBS")
  ###############################################################################
  cat(paste("*************************************************\n"))
  cat(paste("*********** SCENARIO",SAMPLE_TO_RUN," ***********\n"))
  cat(paste("*************************************************\n"))
  
  
  oc_groundtruth <- readRDS(paste0("Univariate/Datasets/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"/trueOC_list50.RDS"))
  dc_groundtruth <- readRDS(paste0("Univariate/Datasets/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"/trueDC_list50.RDS"))
  
  Clustering_results <- list()
  simu_CAM   <- readRDS(paste0(name,"_CAM.RDS"))
  Clustering_results$CAM_OCL   <- parallel_OCL_extract(res = simu_CAM)
  Clustering_results$CAM_DCL   <- parallel_DCL_extract(res = simu_CAM)
  time_CAM <- as.numeric(lapply(simu_CAM,function(x)      as.numeric(x$time,units="secs")))
  memo_CAM <- as.numeric(lapply(simu_CAM,function(x)      log(as.numeric(x$size/10^6))))
  rm(simu_CAM)
  
  simu_fCAM  <- readRDS(paste0(name,"_fCAM.RDS"))
  Clustering_results$fCAM_OCL  <- parallel_OCL_extract(res = simu_fCAM)
  Clustering_results$fCAM_DCL  <- parallel_DCL_extract(res = simu_fCAM)
  time_fCAM <- as.numeric(lapply(simu_fCAM,function(x)     as.numeric(x$time,units="secs")))
  memo_fCAM <- as.numeric(lapply(simu_fCAM,function(x)     log(as.numeric(x$size/10^6))))
  rm(simu_fCAM)
  
  simu_fiSAN <- readRDS(paste0(name,"_fiSAN_OF.RDS"))
  Clustering_results$fiSAN_OCL <- parallel_OCL_extract(res = simu_fiSAN)
  Clustering_results$fiSAN_DCL <- parallel_DCL_extract(res = simu_fiSAN)
  time_fiSAN <- as.numeric(lapply(simu_fiSAN,function(x)    as.numeric(x$time,units="secs")))
  memo_fiSAN <- as.numeric(lapply(simu_fiSAN,function(x)    log(as.numeric(x$size/10^6))))
  rm(simu_fiSAN)

  simu_fSAN <- readRDS(paste0(name,"_fSAN_OF.RDS"))
  Clustering_results$fSAN_OCL <- parallel_OCL_extract(res = simu_fSAN)
  Clustering_results$fSAN_DCL <- parallel_DCL_extract(res = simu_fSAN)
  time_fSAN <- as.numeric(lapply(simu_fSAN,function(x)    as.numeric(x$time,units="secs")))
  memo_fSAN <- as.numeric(lapply(simu_fSAN,function(x)    log(as.numeric(x$size/10^6))))
  rm(simu_fSAN)
  
  ###############################################################################
  name_res <- paste0("Univariate/Results/GIBBS_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_GIBBS")
  ###############################################################################
  saveRDS(Clustering_results,paste0(name_res,"_estimated_clustering_results.RDS"))
  ###############################################################################
  ###############################################################################

  # Rand_indexes ------------------------------------------------------------
  
  Obs_rand <- matrix(NA,sims,4)
  Obs_rand[,1] <- apply(Clustering_results$CAM_OCL,   1, function(x) mclust::adjustedRandIndex(x, oc_groundtruth))
  Obs_rand[,2] <- apply(Clustering_results$fCAM_OCL,  1, function(x) mclust::adjustedRandIndex(x, oc_groundtruth))
  Obs_rand[,3] <- apply(Clustering_results$fiSAN_OCL, 1, function(x) mclust::adjustedRandIndex(x, oc_groundtruth))
  Obs_rand[,4] <- apply(Clustering_results$fSAN_OCL,  1, function(x) mclust::adjustedRandIndex(x, oc_groundtruth))
  colnames(Obs_rand) = c("CAM","fCAM","fiSAN","fSAN") 
  

  Dis_rand <- matrix(NA,sims,4)
  Dis_rand[,1] <- apply(Clustering_results$CAM_DCL,   1, function(x) mclust::adjustedRandIndex(x, dc_groundtruth))
  Dis_rand[,2] <- apply(Clustering_results$fCAM_DCL,  1, function(x) mclust::adjustedRandIndex(x, dc_groundtruth))
  Dis_rand[,3] <- apply(Clustering_results$fiSAN_DCL, 1, function(x) mclust::adjustedRandIndex(x, dc_groundtruth))
  Dis_rand[,4] <- apply(Clustering_results$fSAN_DCL,  1, function(x) mclust::adjustedRandIndex(x, dc_groundtruth))
  colnames(Dis_rand) = c("CAM","fCAM","fiSAN","fSAN") 
  
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
  Dis_time <- cbind(time_CAM,
                    time_fCAM,
                    time_fiSAN ,
                    time_fSAN )
  colnames(Dis_time) = c("CAM","fCAM","fiSAN","fSAN")
  boxplot(Dis_time, main = "Elapsed Seconds")
  saveRDS(Dis_time,paste0(name_res,"_elapsed_seconds.RDS"))
  
  # memo --------------------------------------------------------------------
  Dis_memo <- cbind(memo_CAM,
                    memo_fCAM,
                    memo_fSAN,
                    memo_fiSAN)
  colnames(Dis_memo) = c("CAM","fCAM","fiSAN","fSAN")
  boxplot(Dis_memo, main = "log_Megabytes")
  saveRDS(Dis_memo,paste0(name_res,"_log_Megabytes.RDS"))
  
  
  cat(SAMPLE_TO_RUN)
}

