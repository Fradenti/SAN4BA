library(LaplacesDemon)

# ---------------------
# SAMPLE TO RUN
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

name_vi1 <- paste0("Univariate/Results/VI_RDS/scen1_samplesize",sample_sizes[1],"_VI_posterior_mean_dens_VI.RDS")
name_vi2 <- paste0("Univariate/Results/VI_RDS/scen1_samplesize",sample_sizes[2],"_VI_posterior_mean_dens_VI.RDS")
name_vi3 <- paste0("Univariate/Results/VI_RDS/scen1_samplesize",sample_sizes[3],"_VI_posterior_mean_dens_VI.RDS")
name_vi4 <- paste0("Univariate/Results/VI_RDS/scen1_samplesize",sample_sizes[4],"_VI_posterior_mean_dens_VI.RDS")

# there is a potential problem with CAM in Scenario 1 - it is removed anyway in the paper...
all_d_vb1 <- cbind(readRDS(name_vi1),run = 1)
all_d_vb2 <- cbind(readRDS(name_vi2),run = 2)
all_d_vb3 <- cbind(readRDS(name_vi3),run = 3)
all_d_vb4 <- cbind(readRDS(name_vi4),run = 4)
dim(all_d_vb1)
dim(all_d_vb2)
dim(all_d_vb3)
dim(all_d_vb4)

# View(all_d_vb1)

all_VB <- rbind(all_d_vb1,all_d_vb2,all_d_vb3,all_d_vb4)

dim(all_d_vb1)
dim(all_d_vb2)
dim(all_d_vb3)
dim(all_d_vb4)

all_VB <- all_VB %>% rename(X1 = Var1, X2 = Var2, sim = i) %>% 
  mutate(Y2 = case_when(X2==2  ~ 3,
                        X2==3  ~ 2,
                        TRUE ~ X2)) %>%  
  mutate(X4 = paste("Group",Y2))

table(all_VB$run,all_VB$type)

table(all_VB$type)
all_VB = all_VB %>% mutate(type = case_when(type == "fisan"~"fiSAN",
                                            type == "fsan"~"fSAN",
                                            TRUE ~ NA),
                           run =  case_when(run == "1" ~"N = 60",
                                            run == "2" ~"N = 300",
                                            run == "3" ~"N = 3000",
                                            run == "4" ~"N = 15000"),
                           run = factor(run,levels=c("N = 60","N = 300","N = 3000","N = 15000")))

table(all_VB$run)

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------

GGG = all_VB %>% filter( sim == 1 & X2<=3)

table(GGG$X2)

P = ggplot()+ theme_bw()+
  geom_path(data = all_VB %>% filter( sim == 1 & X2<=3),
            aes(x=x,y=value, group = rev(X2) ),col="lightgray",alpha=.5) + 
  facet_grid(type+X4~run,scales = "free_y")

for(i in 2:50){
  P = P +
    geom_path(data = all_VB %>% filter( sim == i & X2<=3),
              aes(x=x,y=value, group = rev(X2)),col="lightgray",alpha=.5) +
    facet_grid(type+X4~run,scales = "free_y")
  cat(i)
}


P = P +
  scale_color_manual("Group", values = c(2,4,3))+
  ggtitle(paste("VI"))+
  theme(text = element_text(size=16) )+
  ylab("Density") +xlab("y")
P = P +
  geom_path(data=as_tibble(QQ) %>% mutate(X4 = paste("Group",X2)), aes(x=x,y=y,group = X4))
P

ggview::ggview(height =12,w=12)
ggsave("Univariate/output/dens_est_vb_ALL.png",height = 12,width = 12)
ggsave("Univariate/output/dens_est_vb_ALL.pdf",height = 12,width = 12)
ggsave("Univariate/output/dens_est_vb_ALL.eps",
       height = 12, width = 12, device = cairo_ps)


###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------




