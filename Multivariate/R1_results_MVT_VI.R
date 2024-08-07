
library(parallel)
library(tidyverse)

source("Multivariate/A1_MVT_Hyperpar_Gibbs.R")
source("Multivariate/B1_MVT_Aux_Gibbs.RR")

nrep <- 10000
burn <- 5000
Ncores <- 10

sample_sizes <- c(50, 500, 1000) # group size: we have 6 groups of equal size, hence the total sample size ranges from 300 to 6000
dimensions <- c(2,5,10)
oc = dc = c()
for(dd in 1:3){
  for(ss in 1:3){
    ###############################################################################
    SAMPLE_TO_RUN=ss
    DIM = dd
    name <- paste0("Multivariate/Results/VI_RDS/samplesize",sample_sizes[SAMPLE_TO_RUN],
                                              "_dim",dimensions[DIM],"_VI_Results_parallel_times_sum_max.RDS")
  d = readRDS(name)  
  dc = cbind(dc,d$ari_dc)  
  oc = cbind(oc,d$ari_oc)
  }
}    
boxplot(dc)
boxplot(oc)



RES_vi = c()
for(dd in 1:3){
  for(ss in 1:3){
    ###############################################################################
    SAMPLE_TO_RUN=ss
    DIM = dd
    
  name_vi <- paste0("Multivariate/Results/VI_RDS/samplesize",sample_sizes[SAMPLE_TO_RUN],
               "_dim",dimensions[DIM],"_VI_Results_parallel_times_sum_max.RDS")

  h <- readRDS(name_vi)

  q <- cbind(DCari = h$ari_dc, 
        OCari = h$ari_oc,
        #time = h$clock, #  median times
        time = h$clock_max, #  MAX times
        #time = h$clock_sum,
        mega = exp(h$logmega)/50, # verify
        dim = dimensions[DIM],
        samsize = sample_sizes[SAMPLE_TO_RUN])
    colnames(q)= c("DCari","OCari","time","mega","dim","samplesize")
    RES_vi = rbind(RES_vi, q)

  }
  }
RES_VI <- as_tibble(RES_vi) %>% mutate(alg = "VI")


RES_gib = c()
for(dd in 1:3){
  for(ss in 1:3){
    ###############################################################################
    SAMPLE_TO_RUN=ss
    DIM = dd
    
    Ari<- readRDS(paste0("Multivariate/Results/GIBBS_RDS/samplesize",sample_sizes[SAMPLE_TO_RUN],
                      "_dim",dimensions[DIM],"_GIBBS_Rand_indexes.RDS"))    
    
    Clo<- readRDS(paste0("Multivariate/Results/GIBBS_RDS/samplesize",sample_sizes[SAMPLE_TO_RUN],
                         "_dim",dimensions[DIM],"_GIBBS_elapsed_seconds.RDS"))    
    
    mega<- readRDS(paste0("Multivariate/Results/GIBBS_RDS/samplesize",sample_sizes[SAMPLE_TO_RUN],
                         "_dim",dimensions[DIM],"_GIBBS_log_Megabytes.RDS"))    
    
    h <- readRDS(name_vi)
    
    q <- cbind(DCari =Ari[[1]], 
               OCari =Ari[[2]],
               time = Clo,
               mega = exp(mega),
               dim = dimensions[DIM],
               samsize = sample_sizes[SAMPLE_TO_RUN])
    colnames(q)= c("DCari","OCari","time","mega","dim","samplesize")
    RES_gib = rbind(RES_gib, q)
    
  }
}
RES_GI <- as_tibble(RES_gib) %>% mutate(alg = "MCMC")
RES_ALL <- rbind(RES_GI,RES_VI)

################################ 1

cols <- c("sienna3","cyan4")
melted <- reshape2::melt(RES_ALL,c("dim","samplesize","alg"))
melted2 <- melted

melted2 %>% filter(variable == "DCari" | variable == "OCari") %>% 
  mutate(variable = case_when(variable == "DCari"~"DC",
                              variable == "OCari"~"OC"
  )) %>%
  mutate(samplesizes = case_when(samplesize == 50 ~ "N = 300",
                                samplesize == 500 ~ "N = 3000",
                                samplesize == 1000 ~ "N = 6000",
  )) %>% 
  
ggplot()+theme_bw()+
  geom_boxplot(aes(x=factor(dim),y=value,col=alg))+
  facet_grid(variable~samplesizes)+
  theme(legend.position = "bottom", text = element_text(size=16),
        axis.text.x = element_text(angle = 0))+
  scale_color_manual("Algorithm",values=cols)+
  xlab("Dimension")+ylab("Adjusted Rand Index")#+coord_flip()

ggview::ggview(h=6,w=9)
ggsave("Multivariate/Output/SIM_MVT_ARI.pdf",h=6,w=9)
ggsave("Multivariate/Output/SIM_MVT_ARI.png",h=6,w=9)
ggsave("Multivariate/Output/SIM_MVT_ARI.eps",h=6,w=9)



melted2 %>% filter(variable == "mega" | variable == "time" ) %>%
  mutate(samplesizes = case_when(samplesize == 50 ~ "N = 300",
                                 samplesize == 500 ~ "N = 3000",
                                 samplesize == 1000 ~ "N = 6000",
  )) %>% mutate(varib2 = case_when(variable == "time" ~ "Computational time (sec)",
                                   variable == "mega" ~ "Memory usage (MB)")) %>% 
  ggplot()+theme_bw()+
  geom_boxplot(aes(x=factor(dim),y=value,col=alg))+
  facet_grid(varib2~samplesizes,scales = "free_y")+
  theme(legend.position = "bottom", text = element_text(size=16),
        axis.text.x = element_text(angle = 0))+
  scale_color_manual("Algorithm",values=cols)+scale_y_log10()+
  xlab("Dimension")+ylab(NULL)#+coord_flip()

ggview::ggview(h=6,w=9)
ggsave("Multivariate/Output/SIM_MVT_timemax_Memo_flipped_logYAxes.pdf",h=6,w=9)
ggsave("Multivariate/Output/SIM_MVT_timemax_Memo_flipped_logYAxes.png",h=6,w=9)
ggsave("Multivariate/Output/SIM_MVT_timemax_Memo_flipped_logYAxes.eps",h=6,w=9)
