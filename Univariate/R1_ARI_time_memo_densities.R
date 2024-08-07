library(tidyverse)
library(reshape2)
library(patchwork)
library(LaplacesDemon)


# Useful values -----------------------------------------------------------
sample_sizes <- c(10,50,500,2500) #5000
sims         <- 50

# ARI ---------------------------------------------------------------------
ARI_data = c()
for(k in 1:4){
  
  res = readRDS(file = paste0("Univariate/Results/VI_RDS/scen1_samplesize",sample_sizes[k],"_VI_Results_median_and_sum_times.RDS"))
  VB_ari_DC <- reshape2::melt(res$ari_dc) %>% mutate(alg = "VI",type = "DC" )
  VB_ari_OC <- reshape2::melt(res$ari_oc) %>% mutate(alg = "VI",type = "OC" )
  
  gibbs_ari <- readRDS(paste0("Univariate/Results/GIBBS_RDS/scen1_samplesize",sample_sizes[k],"_GIBBS_Rand_indexes.RDS"))
  
  GB_ari_DC <- reshape2::melt(gibbs_ari[[1]]) %>% mutate(alg = "Gibbs",type = "DC")
  GB_ari_OC <- reshape2::melt(gibbs_ari[[2]]) %>% mutate(alg = "Gibbs",type = "OC")
  
  ARI_data <- rbind(ARI_data,rbind(VB_ari_DC,VB_ari_OC,
                                   GB_ari_DC,GB_ari_OC
                                   ) %>% mutate(sim = k))
}

ARI_data$Var2 <- as.character(ARI_data$Var2)
table(ARI_data$Var2)
ARI_data <- ARI_data %>% 
  mutate(sim = paste("Configuration",sim))
table(ARI_data$Var2)
table(ARI_data$Var2,ARI_data$alg)

## Plotting ARI data -----------------------------------------------------------


###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
ARI2 <- ARI_data  %>% mutate(sim2 = case_when(sim == "Configuration 1" ~"N = 60",
                                              sim == "Configuration 2" ~"N = 300",
                                              sim == "Configuration 3" ~"N = 3000",
                                              sim == "Configuration 4" ~"N = 15000"),
                             Var3=Var2) %>% # old code used Var3 renaming - I put this line to not disrupt everything that follows
  mutate(alg2 = case_when(alg == "Gibbs-OF"~"MCMC",
                          alg == "Gibbs" ~ "MCMC",
                          alg == "VI" ~ "VI",
                          TRUE ~ alg) )%>% 
  mutate(sim2 = factor(sim2,levels=c("N = 60","N = 300","N = 3000","N = 15000"))) %>% filter(!(alg == "VI" & Var2 =="CAM"))

table(ARI2$Var3)

# MCMC - VB
cols <- c("sienna3","cyan4")
ggplot(ARI2)+theme_bw()+
  geom_boxplot(aes(x=Var3,y=value,col=alg2))+
  facet_grid(type~sim2)+
  theme(legend.position = "bottom", text = element_text(size=16),
        axis.text.x = element_text(angle = 0))+
  scale_color_manual("Algorithm",values=cols)+
  xlab(NULL)+ylab("Adjusted Rand Index")#+coord_flip()

ggview::ggview(h=6,w=12)
ggsave("Univariate/Output/SIM1_ARI.pdf",h=6,w=12)
ggsave("Univariate/Output/SIM1_ARI.png",h=6,w=12)
ggsave("Univariate/Output/SIM1_ARI.eps",h=6,w=12)

cols <- c("sienna3","cyan4")
ggplot(ARI2)+theme_bw()+
  geom_boxplot(aes(x=Var3,y=value,col=alg2))+
  facet_grid(type~sim2,scales = "free_x")+
  theme(legend.position = "bottom", text = element_text(size=14),
        axis.text.x = element_text(angle = 0))+
  scale_color_manual("Algorithm",values=cols)+
  xlab(NULL)+ylab("Adjusted Rand Index")+coord_flip()


#ggview::ggview(h=7,w=10)
# ggsave("Univariate/Output/SIM1_ARI_freex.pdf",h=7,w=15)
# ggsave("Univariate/Output/SIM1_ARI_freex.png",h=7,w=15)
# ggsave("Univariate/Output/SIM1_ARI_freex.eps",h=7,w=15)

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------

# TIME ------------------------------------------------------------------------

TIME_data = c()
for(k in 1:4){
  
  res = readRDS(file = paste0("Univariate/Results/VI_RDS/scen1_samplesize",sample_sizes[k],"_VI_Results_median_sum_max_times.RDS"))
  
  VB_sec_DC <- reshape2::melt(res$clock_max) %>% mutate(alg = "VI") %>% mutate(value = value)
  
  gibbs_sec <- readRDS(paste0("Univariate/Results/GIBBS_RDS/scen1_samplesize",sample_sizes[k],"_GIBBS_elapsed_seconds.RDS"))
  GB_sec <- reshape2::melt(gibbs_sec) %>% mutate(alg = "Gibbs")
  
  TIME_data <- rbind(TIME_data,rbind(VB_sec_DC,GB_sec) %>% mutate(sim = k))
  
}

plot(TIME_data$value)
table(TIME_data$alg)
table(TIME_data$Var2,TIME_data$alg)
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
Time2 <- TIME_data %>% mutate(Var3 = Var2,
                                                                              sim2 = case_when(sim == "1" ~"N = 60",
                                                                                               sim == "2" ~"N = 300",
                                                                                               sim == "3" ~"N = 3000",
                                                                                               sim == "4" ~"N = 15000")) %>%
  mutate(alg2 = case_when(alg == "Gibbs-OF"~"MCMC",
                          alg == "Gibbs" ~ "MCMC",
                          alg == "VI" ~ "VI",
                          TRUE ~ alg) )%>% 
  mutate(sim2 = factor(sim2,levels=c("N = 60","N = 300","N = 3000","N = 15000"))) %>% filter(!(alg == "VI" & Var2 =="CAM"))

ggplot(Time2)+theme_bw()+
  geom_boxplot(aes(x=Var3, y=value, col=alg2))+
  facet_grid(~sim2)+
  theme(legend.position = "bottom", text = element_text(size=16),
        axis.text.x = element_text(angle = 90))+
  scale_color_manual("Algorithm",values=cols)+
  xlab(NULL)+ylab("Computational time in seconds\n(log10 scale)")+
  scale_y_log10()

# ggsave("Univariate/Output/SIM1_Timelog_max.pdf",h=5,w=12)
# ggsave("Univariate/Output/SIM1_Timelog_max.png",h=5,w=12)
# ggsave("Univariate/Output/SIM1_Timelog_max.eps",h=5,w=12)


###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------


# MEMORY ------------------------------------------------------------------------

MEMO_data = c()
for(k in 1:4){
  res = readRDS(file = paste0("Univariate/Results/VI_RDS/scen1_samplesize",sample_sizes[k],"_VI_Results_median_sum_max_times.RDS"))
  VB_mega_DC <- reshape2::melt(res$logmega) %>% mutate(alg = "VI")
  
  gibbs_mega <- readRDS(paste0("Univariate/Results/GIBBS_RDS/scen1_samplesize",sample_sizes[k],"_GIBBS_log_Megabytes.RDS"))
  
  GB_mega <- reshape2::melt(gibbs_mega) %>% mutate(alg = "Gibbs")
  
  MEMO_data <- rbind(MEMO_data,rbind(VB_mega_DC,GB_mega) %>% mutate(sim = k))
  
}

MEMO_data$Var2 <- as.character(MEMO_data$Var2)
table(MEMO_data$Var2)

MEMO_data = MEMO_data %>% mutate(memo_per_run = case_when(alg=="VI" ~ exp(value)/50, # xk posso sempre scartare e salvarne una!
                                                          alg=="Gibbs" ~ exp(value),
                                                          alg=="Gibbs-OF" ~ exp(value)))



Memo2 <- MEMO_data  %>% mutate(Var3 = Var2,
                                                                              sim2 = case_when(sim == "1" ~"N = 60",
                                                                                               sim == "2" ~"N = 300",
                                                                                               sim == "3" ~"N = 3000",
                                                                                               sim == "4" ~"N = 15000")) %>%
  mutate(alg2 = case_when(alg == "Gibbs-OF"~"MCMC",
                          alg == "Gibbs" ~ "MCMC",
                          alg == "VI" ~ "VI",
                          TRUE ~ alg) )%>% 
  mutate(sim2 = factor(sim2,levels=c("N = 60","N = 300","N = 3000","N = 15000"))) %>% filter(!(alg == "VI" & Var2 =="CAM"))

ggplot(Memo2)+theme_bw()+
  geom_boxplot(aes(x=Var3, y=memo_per_run, col=alg2))+
  facet_grid(~sim2)+
  theme(legend.position = "bottom", text = element_text(size=16),
        axis.text.x = element_text(angle = 90))+
  scale_color_manual("Algorithm",values=cols)+
  xlab(NULL)+ylab("Memory usage in MB\n(log10 scale)")+
  scale_y_log10()

#ggview::ggview(h=5,w=12)
# ggsave("Univariate/Output/SIM1_Memolog_log.pdf",h=5,w=12)
# ggsave("Univariate/Output/SIM1_Memolog_log.png",h=5,w=12)
# ggsave("Univariate/Output/SIM1_Memolog_log.eps",h=5,w=12)


# -------------------------------------------------------------------------


Time3 = Time2 %>% mutate(VVar = "Computational time (sec)")
Memo3 = Memo2 %>% mutate(VVar = "Memory usage (MB)", value= memo_per_run) %>% select(-memo_per_run)
TiMe = rbind(Time3,Memo3)

TiMe <- TiMe %>% mutate(Var3 = factor(Var3,levels = c("CAM","fCAM","fiSAN","fSAN")))

ggplot(TiMe)+theme_bw()+
  geom_boxplot(aes(x=Var3, y=value, col=alg2))+
  facet_grid(VVar~sim2)+#, scales = "free_y")+
  theme(legend.position = "bottom", text = element_text(size=16))+
  scale_color_manual("Algorithm",values=cols)+
  xlab(NULL)+
  ylab(NULL)+
  scale_y_log10() #+coord_flip()

ggview::ggview(h=6,w=12)
ggsave("Univariate/Output/SIM1_maxtimeMega_log.pdf", h=6, w=12)
ggsave("Univariate/Output/SIM1_maxtimeMega_log.png", h=6, w=12)
ggsave("Univariate/Output/SIM1_maxtimeMega_log.eps", h=6, w=12)

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------

