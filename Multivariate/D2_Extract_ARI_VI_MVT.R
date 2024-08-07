library(mclust)
# Rcpp::sourceCpp("Postprocessing_functions/estimate_postdens.cpp")
# source("Postprocessing_functions/postprocessing_output.R")

sample_sizes <- c(50, 500, 1000) # group size: we have 6 groups of equal size, hence the total sample size ranges from 300 to 6000
dimensions <- c(2,5,10)
sims <- 50

for(dd in 1:3){
  for(ss in 1:3){
  ###############################################################################
  SAMPLE_TO_RUN=ss
  DIM = dd
    
  ################################################################################
  name <- paste0("Multivariate/Runs/VI_RDS/samplesize",sample_sizes[SAMPLE_TO_RUN],
                 "_dim",dimensions[DIM],"_VI")
  ###############################################################################
  
  res_fiSAN <- readRDS(paste0(name,"_MVT_fiSAN.RDS"))
  
  oc_groundtruth <- readRDS(paste0("Multivariate/Datasets/samplesize",sample_sizes[SAMPLE_TO_RUN],
                                     "_dim",dimensions[DIM],"/trueOC_list50.RDS"))
  dc_groundtruth <- readRDS(paste0("Multivariate/Datasets/samplesize",sample_sizes[SAMPLE_TO_RUN],
                                     "_dim",dimensions[DIM],"/trueDC_list50.RDS"))
    
  
  L           <- lapply(res_fiSAN, function(x) x$results)
  OC <- DC <- list()
  for(i in 1:sims){
    OC[[i]] <- as.numeric(factor(unlist(lapply(L[[i]]$sim$XI,
                                                 function(x) apply(x, 1, which.max)))))
    DC[[i]] <- as.numeric(factor(apply(L[[i]]$sim$RHO,1,which.max)))
  }  
    
  clock_fisan   <- unlist(lapply(res_fiSAN, function(z) median(as.numeric(z$time_50,unit="secs"))))
  clock_fisan2   <- unlist(lapply(res_fiSAN, function(z) sum(as.numeric(z$time_50,unit="secs"))))
  clock_fisan_max   <- unlist(lapply(res_fiSAN, function(z) max(as.numeric(z$time_50,unit="secs"))))
  
  logmemo_fisan <- log(unlist(lapply(res_fiSAN, function(z) as.numeric(z$size/10^6))))
  est_dis     <- do.call(rbind,DC)
  est_obs     <- do.call(rbind,OC)
  fisan_oc      <- apply((est_obs),1,function(r) mcclust::arandi(r,oc_groundtruth)) 
  fisan_dc      <- apply(est_dis,1,function(r)   mcclust::arandi(r,dc_groundtruth)) 
  
  par(mfrow=c(2,2))
  ari_oc = (cbind(fisan_oc))
  colnames(ari_oc) = c("fiSAN")
  boxplot(ari_oc)
  ari_dc = (cbind(fisan_dc))
  colnames(ari_dc) = c("fiSAN")
  boxplot(ari_dc)
  clock = (cbind(clock_fisan))
  colnames(clock) = c("fiSAN")
  boxplot(clock)
  clock_sum = (cbind(clock_fisan2))
  colnames(clock_sum) = c("fiSAN")
  boxplot(clock_sum)
  clock_max = (cbind(clock_fisan_max))
  colnames(clock_max) = c("fiSAN")
  boxplot(clock_max)
  logmega = (cbind(logmemo_fisan))
  colnames(logmega) = c("fiSAN")
  boxplot(logmega)
  par(mfrow=c(1,1))
  name <- paste0("Multivariate/Results/VI_RDS/samplesize",sample_sizes[SAMPLE_TO_RUN],
         "_dim",dimensions[DIM],"_VI")
  cat(SAMPLE_TO_RUN) 
  saveRDS(list(ari_dc=ari_dc, 
               ari_oc = ari_oc,
               clock=clock,
               clock_sum=clock_sum,
               clock_max=clock_max,
               logmega = logmega),
          file = paste0(name,"_Results_parallel_times_sum_max.RDS"))
  
}
}
