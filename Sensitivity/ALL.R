library(tidyverse)
library(crayon)
library(SANvi)
runs  = 50
cores = 18
maxL  = 25
maxK  = 20
warmstart = F
hyp_alpha1 = 1
hyp_alpha2 = 1
hyp_beta1 = 1
hyp_beta2 = 1
alpha_bar = .05
beta_bar = .05
root = 250892

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------

sample_sizes <- c(10,50,500,2500)
LL <-  c(15,25,35)
KK <-  c(10,20,30)
BB <- c(0.001,0.05,1)
ss <- 2

sims <- 50
root <- 23442
runs <-  50

SAMPLE_TO_RUN <- ss

y      <- readRDS(paste0("Univariate/Datasets/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"/data_list",sims,".RDS"))
group  <- readRDS(paste0("Univariate/Datasets/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"/group_list",sims,".RDS"))
name   <- paste0("Sensitivity/Runs/VI_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN])

if(is.null(name)){stop("missing variable: name")}
if(is.null(root)){stop("missing variable: root - sure you do not want a seed?")}

Ofisan <- list()
res_fisan <- list()
ALL_FISAN <- list()
for(ll in LL){
  for(kk in KK){
    for(bb in BB){
      
      ################################################################################
      root <- 9241 * (ll+kk+bb)
      
      for(i in 1:sims){

        for(rr in 1:runs){
        
        ALL_FISAN[[rr]] <-  SANvi::variational_fiSAN(y = y[[i]],
                                                    group = group[[i]],
                                                    maxL = ll,
                                                    maxK = kk,
                                                    maxSIM = 500,
                                                    epsilon = 1e-4,
                                                    warmstart = warmstart,
                                                    conc_hyperpar = c(hyp_alpha1,
                                                                      hyp_alpha2),
                                                    beta_bar = bb,
                                                    seed = rr*root)
        if(rr %% 10 == 0 ) cat(i);   
        
      }
      
      Ofisan$size <- object.size(ALL_FISAN)
      LT = unlist(lapply(ALL_FISAN, function(x) as.numeric(x$time,units="secs" )))
      Ofisan$time_50 <- LT
      ind_best_fisan <- which.max(unlist(lapply(ALL_FISAN, function(u) max(u$sim$Elbo_val))))
      Ofisan$results <- (ALL_FISAN[[ind_best_fisan]])
      cat(green(" \u2713 \n"))
      res_fisan[[i]] <-  Ofisan
      
      if(i == sims){
        cat("Saving progress...\n")
        saveRDS(res_fisan,paste0(name,"_VI_fiSAN_OF_L_",
                                 ll,
                                 "_K_",kk,
                                 "_b_",bb,".RDS"))
      }
      }
    }
  }
}


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 
library(mclust)
source("Univariate/Helpers/postprocessing_output.R")
oc_groundtruth <- readRDS(paste0("Univariate/Datasets/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"/trueOC_list50.RDS"))
dc_groundtruth <- readRDS(paste0("Univariate/Datasets/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"/trueDC_list50.RDS"))
# 




i=0
for(ll in LL){
  for(kk in KK){
    for(bb in BB){
      name   <- paste0("Sensitivity/Runs/VI_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN])
      
      resnam <- paste0(name,"_VI_fiSAN_OF_L_",
                       ll,
                       "_K_",kk,
                       "_b_",bb,".RDS")
      
      Clustering_results <- list()
      res_fiSAN <- readRDS(resnam)
      cat("----\n")
      cat(unlist(lapply(res_fiSAN,function(x) c(x$results$params$L)))[1])
      cat("...\n")
      cat(unlist(lapply(res_fiSAN,function(x) c(x$results$params$K)))[1])
      cat("...\n")
      cat(unlist(lapply(res_fiSAN,function(x) c(x$results$params$beta_bar)))[1])
      cat("----\n")
      Los           <- lapply(res_fiSAN, function(x) x$results)
      # now max time!
      clock_fisan   <- unlist(lapply(res_fiSAN, function(z) max(as.numeric(z$time_50,unit="secs"))))
      logmemo_fisan <- log(unlist(lapply(res_fiSAN, function(z) as.numeric(z$size/10^6))))
      est_dis     <- do.call(rbind,lapply(Los, function(z) SANvi::estimate_clustering_vi(z)$dis_level))
      est_obs     <- do.call(rbind,lapply(Los, function(z) SANvi::estimate_clustering_vi(z)$obs_level$OC))
      fisan_oc    <- apply((est_obs),1,function(r) mcclust::arandi(r,oc_groundtruth)) 
      fisan_dc    <- apply(est_dis,1,function(r)   mcclust::arandi(r,dc_groundtruth)) 
      
      name_res <- paste0("Sensitivity/Results/VI_RDS/VI_fiSAN_OF_L_",
                         ll,
                         "_K_",kk,
                         "_b_",bb)
      
      ###############################################################################
      saveRDS(list(est_dis,est_obs),paste0(name_res,"_estimated_clustering_results.RDS"))
      
      All_rands <- list(fisan_oc,fisan_dc)
      saveRDS(All_rands,paste0(name_res, "_Rand_indexes.RDS"))
      i = i +1
      cat(paste(i,"\n"))
    }
  }
}

RESU <- c()
for(ll in LL){
  for(kk in KK){
    for(bb in BB){
      
      name_res <- paste0("Sensitivity/Results/VI_RDS/VI_fiSAN_OF_L_",
                         ll,
                         "_K_",kk,
                         "_b_",bb)
      
      ari = readRDS(paste0(name_res,"_Rand_indexes.RDS"))
      RESU1 <- data.frame(ari = ari[[2]], type = "DC", L = ll, b = bb, K = kk)
      RESU2 <- data.frame(ari = ari[[1]], type = "OC", L = ll, b = bb, K = kk)
      RESU <-  rbind(RESU, RESU1,RESU2)    
      
      
    }
  }
}

table(RESU$L,RESU$K)
table(RESU$L,RESU$b)
# 
library(tidyverse)

RESU <- RESU %>% mutate(L = paste("L =",L),K = paste("T =",K),b = paste(b))

a1 = ggplot(RESU %>% filter(type=="DC"))+
  theme_bw()+
  geom_boxplot(aes(x=factor(b),y=ari))+
  facet_grid(L~K) + theme(legend.position = "none", text = element_text(size=16))+
  ylab("ARI") + xlab("b") + ggtitle("VI - Distributional clustering")
a1
#ggview::ggview(h=8,w=8)
ggsave("Sensitivity/Output/VI_DC_ari.pdf",h=8,w=8)
ggsave("Sensitivity/Output/VI_DC_ari.eps",h=8,w=8)

a2 = ggplot(RESU %>% filter(type=="OC"))+
  theme_bw()+
  geom_boxplot(aes(x=factor(b),y=ari))+
  facet_grid(L~K) + theme(legend.position = "none", text = element_text(size=16))+
  ylab("ARI") + xlab("b") + ggtitle("VI - Observational clustering")
a2
#ggview::ggview(h=8,w=8)
ggsave("Sensitivity/Output/VI_OC_ari.pdf",h=8,w=8)
ggsave("Sensitivity/Output/VI_OC_ari.eps",h=8,w=8)
# 
# 
# 
# 
RESUcl <- c()
for(ll in LL){
  for(kk in KK){
    for(bb in BB){
      
      name_res <- paste0("Sensitivity/Results/VI_RDS/VI_fiSAN_OF_L_",
                         ll,
                         "_K_",kk,
                         "_b_",bb)
      
      clr = readRDS(paste0(name_res,"_estimated_clustering_results.RDS"))
      RESU1cl <- data.frame(cl = apply((clr[[1]]),1,function(x) length(unique(x))), 
                            type = "DC", L = ll, b = bb, K = kk)
      
      RESU2cl <- data.frame(cl =  apply((clr[[2]]),1,function(x) length(unique(x))), type = "OC", L = ll, b = bb, K = kk)
      RESUcl <-  rbind(RESUcl, RESU1cl,RESU2cl)    
      
      
    }
  }
}
# 
# 
# RESUcl
# 
RESUcl <- RESUcl %>% mutate(L = paste("L =",L),K = paste("T =",K),b = paste(b))

b1 = ggplot(RESUcl %>% filter(type=="DC"))+
  theme_bw()+
  geom_hline(yintercept = 3,col=2)+
  geom_boxplot(aes(x=factor(b),y=cl))+
  facet_grid(L~K) + theme(legend.position = "none", text = element_text(size=16))+
  ylab("Estimated number of clusters") + xlab("b") + ggtitle("VI - Distributional clustering")
b1
#ggview::ggview(h=8,w=8)
ggsave("Sensitivity/Output/VI_CLnum_DC.pdf",h=8,w=8)
ggsave("Sensitivity/Output/VI_CLnum_DC.eps",h=8,w=8)
#
b2 = ggplot(RESUcl %>% filter(type=="OC"))+
  theme_bw()+
  geom_hline(yintercept = 5,col=2)+
  geom_boxplot(aes(x=factor(b),y=cl))+
  facet_grid(L~K) + theme(legend.position = "none", text = element_text(size=16))+
  ylab("Estimated number of clusters") + xlab("b") + ggtitle("VI - Observational clustering")
b2
#ggview::ggview(h=8,w=8)
ggsave("Sensitivity/Output/VI_CLnum_OC.pdf",h=8,w=8)
ggsave("Sensitivity/Output/VI_CLnum_OC.eps",h=8,w=8)
#
#
library(patchwork)
a1+b1
#ggview::ggview(h=8,w=16)
ggsave("Sensitivity/Output/VI_ARICLnum_DC.pdf",h=8,w=16)
ggsave("Sensitivity/Output/VI_ARICLnum_DC.eps",h=8,w=16)
#
#
a2+b2
#ggview::ggview(h=8,w=16)
ggsave("Sensitivity/Output/VI_ARICLnum_OC.pdf",h=8,w=16)
ggsave("Sensitivity/Output/VI_ARICLnum_OC.eps",h=8,w=16)






# time and memo -----------------------------------------------------------

RESU <- c()
for(ll in LL){
  for(kk in KK){
    for(bb in BB){
      
      name   <- paste0("Sensitivity/Runs/VI_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN])
      
      resnam <- paste0(name,"_VI_fiSAN_OF_L_",
                       ll,
                       "_K_",kk,
                       "_b_",bb,".RDS")
      
      res_fiSAN <- readRDS(resnam)
      
      RESU1 <- data.frame(value = unlist(lapply(res_fiSAN, function(z) sum(as.numeric(z$time_50,unit="secs")))), 
                          type = "seconds", L = ll, b = bb, K = kk)

      RESU2 <- data.frame(value = unlist(lapply(res_fiSAN,function(x) x$size/10^6)), 
                          type = "mega", L = ll, b = bb, K = kk)
      RESU <-  rbind(RESU, RESU1,RESU2)    
      
      
    }
  }
}

table(RESU$L,RESU$K)
table(RESU$L,RESU$b)
# 
RESU <- RESU %>% mutate(L = paste("L =",L),K = paste("T =",K),b = paste(b))


b1 = ggplot(RESU %>% filter(type=="seconds"))+
  theme_bw()+
  #geom_hline(yintercept = 3,col=2)+
  geom_boxplot(aes(x=factor(b),y=value))+
  #  scale_y_log10()+
  facet_grid(L~K) + theme(legend.position = "none", text = element_text(size=16))+
  ylab("Elapsed seconds") + xlab("b") + ggtitle("VI - Computational time")
b1

ggview::ggview(h=9, w=9)
ggsave("Sensitivity/Output/VI_LKtime_max.pdf",h=9,w=9)
ggsave("Sensitivity/Output/VI_LKtime_max.eps",h=9,w=9)
#
b2 = ggplot(RESU %>% filter(type=="mega"))+
  theme_bw()+
  #geom_hline(yintercept = 5,col=2)+
  geom_boxplot(aes(x=factor(b),y=value))+
  facet_grid(L~K) + theme(legend.position = "none", text = element_text(size=16))+
  ylab("Megabytes") + xlab("b") + ggtitle("VI - Memory usage")
b2
#ggview::ggview(h=8,w=8)
ggsave("Sensitivity/Output/VI_LKmega.eps",h=9,w=9)
ggsave("Sensitivity/Output/VI_LKmega.pdf",h=9,w=9)




