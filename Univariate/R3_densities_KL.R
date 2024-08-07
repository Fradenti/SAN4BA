library(tidyverse)
library(reshape2)
library(patchwork)
library(LaplacesDemon)

# Useful values -----------------------------------------------------------
sample_sizes <- c(10,50,500,2500) #5000
sims         <- 50

# DENSITIES --------------------------------------------------------------------

## TRUTH -----------------------------------------------------------
sample_sizes <- c(10,50,500,2500)
replica_group = 2
s  = .6
mu = c(-5,-2,2,5,0)
sims = 50

f1 = function(x) .5*dnorm(x,mu[1],s)+.5*dnorm(x,mu[2],s)
f2 = function(x) .5*dnorm(x,mu[3],s)+.5*dnorm(x,mu[4],s)
f3 = function(x) .5*dnorm(x,mu[5],s)+.5*dnorm(x,mu[5],s)

Q = seq(-10,10,by=.05)
QQ <- rbind(
  cbind(x=Q,y=f1(Q),X2=1),
  cbind(x=Q,y=f2(Q),X2=2),
  cbind(x=Q,y=f3(Q),X2=3)
)
cols <- c("sienna3","cyan4")

## ESTIMATED ----------------------------------------------------------------
### Configuration 1 ---------------------------------------------------------

ADof = AD = list()
ADvb = list()

for(SAMPLE_TO_RUN in 1:4){
  

  name_gibbs <- paste0("Univariate/Results/GIBBS_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_GIBBS_posterior_mean_dens_over_grid_GIBBS.RDS")
  all_densities_gibbs <- readRDS(name_gibbs)
  all_densities_gibbs <- all_densities_gibbs %>% mutate(X4 = paste("Group",X2))
  
  
  name_vi <- paste0("Univariate/Results/VI_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_VI_posterior_mean_dens_VI.RDS")
  all_d_vb <- readRDS(name_vi)
  all_dens_vb <- all_d_vb %>% rename(X1 = Var1, X2 = Var2, sim = i) %>% 
    mutate(Y2 = case_when(X2==2  ~ 3,
                          X2==3  ~ 2,
                          TRUE ~ X2)) %>%  
    mutate(X4 = paste("Group",Y2))
  
  AD[[SAMPLE_TO_RUN]] <- all_densities_gibbs %>% filter(X2<4) %>% mutate(truth = case_when(X4 == "Group 1" ~ f1(x),
                                                                                     X4 == "Group 2" ~ f2(x),
                                                                                     X4 == "Group 3" ~ f3(x))) %>% 
    group_by(type, sim, X4) %>% summarise(KLdist = KLD(value,truth)$mean.sum.KLD)  %>% mutate(Sample = SAMPLE_TO_RUN, algo = "GIBBS")
  
  ADvb[[SAMPLE_TO_RUN]] <- all_dens_vb %>% filter(X2<4) %>% mutate(truth = case_when(X4 == "Group 1" ~ f1(x),
                                                                                     X4 == "Group 2" ~ f2(x),
                                                                                     X4 == "Group 3" ~ f3(x))) %>% 
    group_by(type, sim, X4) %>% summarise(KLdist = KLD(value,truth)$mean.sum.KLD) %>% mutate(Sample = SAMPLE_TO_RUN, algo = "VI")  %>% mutate(type = case_when(type == "fsan" ~ "fSAN", 
                                                                                                                                                               type == "fisan"~ "fiSAN",
                                                                                                                                                               TRUE ~ type))
}

KL_gibbs <- do.call("rbind",AD)
table(KL_gibbs$type)
KL_vinf <- do.call("rbind",ADvb)
table(KL_vinf$type)
#saveRDS(KL_vinf,"Univariate/Results/VI_RDS/scen1_all_KL.RDS")

KL_all <- rbind(KL_gibbs,KL_vinf)

colnames(KL_all) = c("Var2","run","Group","value","sim","alg")
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
KL2 <- KL_all %>% filter(!(Var2 == "fiSAN" & alg == "Gibbs")) %>% mutate(Var3 = case_when(Var2 ==  "fSAN_OF" ~ "fSAN",
                                                                                          Var2 == "fSAN" & alg == "GIBBS" ~ "fCAM",
                                                                                          Var2 == "fiSAN_OF" ~ "fiSAN",
                                                                                          TRUE ~ Var2),
                                                                         sim2 = case_when(sim == "1" ~"N = 60",
                                                                                          sim == "2" ~"N = 300",
                                                                                          sim == "3" ~"N = 3000",
                                                                                          sim == "4" ~"N = 15000")) %>%
  mutate(alg2 = case_when(alg == "GIBBS-OF"~"MCMC",
                          alg == "GIBBS" ~ "MCMC",
                          alg == "VI" ~ "VI",
                          TRUE ~ alg) )%>% 
  mutate(sim2 = factor(sim2,levels=c("N = 60","N = 300","N = 3000","N = 15000"))) %>% filter(!(alg == "VI" & Var2 =="CAM"))


a = factor(KL2$Var3,
           levels= rev(levels(as.factor(KL2$Var3))))
table(a)
table(a,KL2$alg)

ggplot(KL2)+theme_bw()+
  geom_boxplot(aes(y=value, x=Var3, col=alg2))+
  facet_grid(Group~sim2)+
  theme(legend.position = "bottom", text = element_text(size=16),
        axis.text.x = element_text(angle = 0))+
  scale_color_manual("Algorithm",values=cols)+
  ylab(NULL)+xlab("KL divergence")

ggsave("Univariate/Output/SIM1_KL.pdf",h=9,w=12)
ggsave("Univariate/Output/SIM1_KL.png",h=9,w=12)
ggsave("Univariate/Output/SIM1_KL.eps",h=9,w=12)


###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------



# ---------------------
# SAMPLE TO RUN
SAMPLE_TO_RUN=1

name <- paste0("Univariate/Results/GIBBS_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_GIBBS_posterior_mean_dens_over_grid_GIBBS.RDS")
all_densities = readRDS(name)
all_densities = all_densities %>% mutate(X4 = paste("Group",X2))

table(all_densities$type)

name_vi <- paste0("Univariate/Results/VI_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_VI_posterior_mean_dens_VI.RDS")
all_d_vb <- readRDS(name_vi)
all_dens_vb <- all_d_vb %>% rename(X1 = Var1, X2 = Var2, sim = i) %>% 
  mutate(Y2 = case_when(X2==2  ~ 3,
                        X2==3  ~ 2,
                        TRUE ~ X2)) %>%  
  mutate(X4 = paste("Group",Y2))

AD[[SAMPLE_TO_RUN]] <- all_densities %>% filter(X2<4) %>% mutate(truth = case_when(X4 == "Group 1" ~ f1(x),
                                                                                   X4 == "Group 2" ~ f2(x),
                                                                                   X4 == "Group 3" ~ f3(x))) %>% 
  group_by(type, sim, X4) %>% summarise(KLdist = KLD(value,truth)$mean.sum.KLD)  %>% mutate(Sample = SAMPLE_TO_RUN, algo = "GIBBS") 

head(all_densities)
head(all_d_vb)

## CAM sim 4 and 11 are rightfully labeled as G2, no need to relabel
all_dens_vb <- all_d_vb %>% rename(X1 = Var1, X2 = Var2, sim = i) %>% 
  mutate(Y2 = case_when(X2==2  ~ 3,
                        X2==3  ~ 2,
                        TRUE ~ X2)) %>%  
  mutate(X4 = paste("Group",Y2))
# all_dens_vb = all_dens_vb %>% mutate(Y2 = case_when((type=="CAM" & sim ==4  & X2 ==2) ~ 2,
#                                                     (type=="CAM" & sim ==11 & X2 ==3) ~ 3))

table(all_dens_vb$type)
all_dens_vb = all_dens_vb %>% mutate(type = case_when(type == "fisan"~"fiSAN",
                                                      type == "fsan"~"fSAN",
                                                      TRUE ~ "CAM"))

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------

alld2 <- all_densities %>% 
  filter(!(type == "fiSAN")) %>% mutate(Var3 = case_when(type ==  "fSAN_OF" ~ "fSAN",
                                                         type == "fSAN"   ~ "fCAM",
                                                         type == "fiSAN_OF" ~ "fiSAN",
                                                         TRUE ~ type))
table(alld2$Var3)


R = ggplot()+ theme_bw()+
  geom_path(data = alld2 %>% filter( sim == 1 & X2<=3),
            aes(x=x,y=value, group = rev(X2) ),col="lightgray",alpha=.5) + facet_grid(X4~Var3,scales = "free_y")
P = ggplot()+ theme_bw()+
  geom_path(data = all_dens_vb %>% filter( sim == 1 & X2<=3),
            aes(x=x,y=value, group = rev(X2) ),col="lightgray",alpha=.5) + facet_grid(X4~type,scales = "free_y")

# Q = ggplot()+ theme_bw()+
#   geom_path(data = all_densities %>% filter( sim == 1 & X2<=3),
#             aes(x=x,y=value, group = rev(X2) ),col="lightgray",alpha=.5) +
#   geom_path(data = all_dens_vb %>% filter( sim == 1 & X2<=3),
#             aes(x=x,y=value, group = rev(X2) ),col="blue",alpha=.5) + facet_grid(X4~type,scales = "free_y")

for(i in 2:50){
  R = R +
    geom_path(data = alld2 %>% filter( sim ==i& X2<=3),
              aes(x=x,y=value, group = rev(X2)),col="lightgray",alpha=.5) + 
    facet_grid(X4~Var3,scales = "free_y")
  P = P +
    geom_path(data = all_dens_vb %>% filter( sim ==i& X2<=3),
              aes(x=x,y=value, group = rev(X2)),col="lightgray",alpha=.5) +
    facet_grid(X4~type,scales = "free_y")
  cat(i)
}
R = R +
  scale_color_manual("Group", values = c(2,4,3))+
  ggtitle(paste("MCMC"))+
  theme(text = element_text(size=16) )+
  ylab("Density") +xlab("y")
R = R +
  geom_path(data=as_tibble(QQ) %>% mutate(X4 = paste("Group",X2)), aes(x=x,y=y,group = X4))
R

ggview::ggview(h=6,w=12)
ggsave("Univariate/Output/SIM1_dens_conf1.pdf",h=6,w=12)
ggsave("Univariate/Output/SIM1_dens_conf1.png",h=6,w=12)
ggsave("Univariate/Output/SIM1_dens_conf1.eps",h=6,w=12,device = cairo_ps)


###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
