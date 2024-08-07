source("Univariate/Helpers/postprocessing_output.R")
sample_sizes = c(10,50,500,2500)#

for(SAMPLE_TO_RUN in 1:4){

  ################################################################################
  name <- paste0("Univariate/Runs/VI_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_VI")
  ###############################################################################

  res_fSAN  <- readRDS(paste0(name,"_fSAN.RDS"))
  res_fiSAN <- readRDS(paste0(name,"_fiSAN.RDS"))
  
  oc_groundtruth <- readRDS(paste0("Univariate/Datasets/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"/trueOC_list50.RDS"))
  dc_groundtruth <- readRDS(paste0("Univariate/Datasets/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"/trueDC_list50.RDS"))

  clock_fsan   <- unlist(lapply(res_fSAN, function(z) median(as.numeric(z$time_50,unit="secs"))))
  clock_fsan2   <- unlist(lapply(res_fSAN, function(z) sum(as.numeric(z$time_50,unit="secs"))))
  clock_fsan_max   <- unlist(lapply(res_fSAN, function(z) max(as.numeric(z$time_50,unit="secs"))))
  logmemo_fsan <- log(unlist(lapply(res_fSAN, function(z) as.numeric(z$size/10^6))))
  est_dis     <- do.call(rbind,lapply(res_fSAN, function(z) SANvi::estimate_clustering_vi(z$results)$dis_level))
  est_obs     <- do.call(rbind,lapply(res_fSAN, function(z) SANvi::estimate_clustering_vi(z$results)$obs_level$OC))
  fsan_oc      <- apply((est_obs),1,function(r) mcclust::arandi(r,oc_groundtruth)) 
  fsan_dc      <- apply(est_dis,1,function(r)   mcclust::arandi(r,dc_groundtruth)) 
  
  clock_fisan   <- unlist(lapply(res_fiSAN, function(z) median(as.numeric(z$time_50,unit="secs"))))
  clock_fisan2   <- unlist(lapply(res_fiSAN, function(z) sum(as.numeric(z$time_50,unit="secs"))))
  clock_fisan_max   <- unlist(lapply(res_fiSAN, function(z) max(as.numeric(z$time_50,unit="secs"))))
  logmemo_fisan <- log(unlist(lapply(res_fiSAN, function(z) as.numeric(z$size/10^6))))
  est_dis     <- do.call(rbind,lapply(res_fiSAN, function(z) SANvi::estimate_clustering_vi(z$results)$dis_level))
  est_obs     <- do.call(rbind,lapply(res_fiSAN, function(z) SANvi::estimate_clustering_vi(z$results)$obs_level$OC))
  fisan_oc      <- apply((est_obs),1,function(r) mcclust::arandi(r,oc_groundtruth)) 
  fisan_dc      <- apply(est_dis,1,function(r)   mcclust::arandi(r,dc_groundtruth)) 
  
  par(mfrow=c(1,2))
  ari_oc = (cbind(fsan_oc,fisan_oc))
  colnames(ari_oc) = c("fSAN","fiSAN")
  boxplot(ari_oc)
  ari_dc = (cbind(fsan_dc,fisan_dc))
  colnames(ari_dc) = c("fSAN","fiSAN")
  boxplot(ari_dc)
  par(mfrow=c(1,3))
  clock = (cbind(clock_fsan,clock_fisan))
  colnames(clock) = c("fSAN","fiSAN")
  boxplot(clock)
  clock_sum = (cbind(clock_fsan2,clock_fisan2))
  colnames(clock_sum) = c("fSAN","fiSAN")
  boxplot(clock_sum, main = "time sum")
  clock_max = (cbind(clock_fsan_max,clock_fisan_max))
  colnames(clock_max) = c("fSAN","fiSAN")
  boxplot(clock_max,main = "MAX")
  par(mfrow=c(1,1))
  logmega = (cbind(logmemo_fsan,logmemo_fisan))
  colnames(logmega) = c("fSAN","fiSAN")
  boxplot(logmega, main = "memo")
  par(mfrow=c(1,1))
  
  name <- paste0("Univariate/Results/VI_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_VI")
 cat(SAMPLE_TO_RUN) 
saveRDS(list(ari_dc=ari_dc, 
             ari_oc = ari_oc,
             clock=clock,
             clock_sum=clock_sum,
             clock_max=clock_max,
             logmega = logmega),
        file = paste0(name,"_Results_median_sum_max_times.RDS"))
}

