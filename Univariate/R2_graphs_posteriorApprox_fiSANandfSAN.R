library(tidyverse)
library(patchwork)
library(reshape2)
library(label.switching)
library(ggExtra)

SAMPLE_TO_RUN <- 2
sample_sizes= c(10,50,500,2500)#)0)


# Read runs ---------------------------------------------------------------

name <- paste0("Univariate/Runs/VI_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_VI")
vif_fiSAN <- readRDS(paste0(name,"_fiSAN.RDS"))

name <- paste0("Univariate/Runs/GIBBS_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_GIBBS")
gib_fiSAN  <- readRDS(paste0(name,"_fiSAN_OF.RDS"))


### functions to sample VI posterior--------------------------------------------
sampling_vi_G <- function(output, nsim){
  
  check  <- output$results$params$gamma0
  ind    <- which(!output$results$sim$theta_l[,4]==check)
  active <- as_tibble(output$results$sim$theta_l[ind,]) %>% arrange(V1)
  active <- as.matrix(active)
  RES  = array(NA, c(nsim,2,nrow(active)))
  for(i in 1:nrow(active)){
    sig2 <- 1/rgamma(nsim, active[i,3], active[i,4])
    mu   <- rnorm(nsim, active[i,1], sqrt(sig2/active[i,2]))
    RES[,,i] <- cbind(mu = mu, sig2 = sig2)
  }
  return(RES)
}


RUN <- 2
# sample from variational posterior
vi_sampled = sampling_vi_G(output =vif_fiSAN[[RUN]], nsim = 5000)

# remove label switching
one_model <- gib_fiSAN[[RUN]]
l1 <- ecr.iterative.1(one_model$sim$obs_cluster,K=one_model$params$maxL)

par(mfrow=c(1,1))
disentangled_mu <- t(sapply(1:nrow(one_model$sim$mu),function(x) one_model$sim$mu[x,l1$permutations[x,]]))
head(disentangled_mu)
disentangled_sig2 <- t(sapply(1:nrow(one_model$sim$sigma2),function(x) one_model$sim$sigma2[x,l1$permutations[x,]]))

# five chains has the lowest sd: they are the one that we are looking for!
inds = sort(apply(disentangled_mu,2,sd),ind = T)


par(mfrow=c(1,2))
# not really representative
matplot(one_model$sim$mu[,inds$ix[1:5]],type="l")
matplot(disentangled_mu[,inds$ix[1:5]],type="l")
matplot(one_model$sim$sigma2[,inds$ix[1:5]],type="l")
matplot(disentangled_sig2[,inds$ix[1:5]],type="l")
par(mfrow=c(1,1))

# ---------------------------------------------------------------------------------------
plot(apply(one_model$sim$mu,2,sd))
inds_orig = sort(apply(one_model$sim$mu,2,sd),ind = T)

matplot(one_model$sim$mu[,inds_orig$ix[1:5]],type="l")
matplot(disentangled_mu[,inds$ix[1:5]],type="l")

Omu = reshape2::melt(one_model$sim$mu[,inds_orig$ix[1:5]])
Dmu = reshape2::melt(disentangled_mu[,inds$ix[1:5]])

plot(ts(one_model$sim$mu[,inds_orig$ix[1:5]]))

Omu = as_tibble(Omu) %>% mutate(type = "Original chains")
Dmu = as_tibble(Dmu) %>% mutate(type = "Post-processed chains")
ALLmu = rbind(Omu,Dmu)

ggplot(ALLmu)+theme_bw()+
  geom_line(aes( x = Var1, y= value, col=factor(Var2), group = factor(Var2)))+
  facet_wrap(~type,nrow = 2,scales = "free_y")+
  scale_color_brewer(palette = "Set1")  +
  ylab("Posterior means") + xlab("MCMC iterations")+
  theme(text = element_text(size=16), legend.position = "none")
################################################################################
#ggview::ggview(h=8,w=12)
ggsave("Univariate/Output/labelsw.pdf",h=8,w=12)
ggsave("Univariate/Output/labelsw.eps",h=8,w=12)
################################################################################


library(coda)
geweke.plot(as.mcmc(one_model$sim$alpha))
MU <- disentangled_mu[,inds$ix[1:5]]
D <- cbind(MU,one_model$sim$alpha)
colnames(D) = c("mu1","mu2","mu3","mu4","mu5","alpha")
geweke.plot(as.mcmc(D),pvalue = .05,nbins = 10,frac1 = .25,frac2 = .25)
geweke.diag(as.mcmc(disentangled_mu[,inds$ix[1:5]]),frac1 = .25,frac2 = .25)


D = as_tibble(D)
mD = D %>% mutate(mu1 = dplyr::cummean(mu1),
                  mu2 = dplyr::cummean(mu2),
                  mu3 = dplyr::cummean(mu3),
                  mu4 = dplyr::cummean(mu4),
                  mu5 = dplyr::cummean(mu5),
                  alpha = dplyr::cummean(alpha))

rD =reshape2::melt(D) %>% mutate(x = rep(1:5000,6))
rmD =reshape2::melt(mD) %>% mutate(x = rep(1:5000,6))

ggplot()+
  geom_line(data=rD,aes(x=x,y=value))+
  geom_line(data=rmD,aes(x=x,y=value),col=4,lwd=1)+
  facet_wrap(~variable,scales = "free_y")+
  ylab("") + xlab("MCMC iterations")+theme_bw()+
  theme(text = element_text(size=16), legend.position = "none")

###############################################################################
#ggview::ggview(h=8,w=12)
ggsave("Univariate/Output/chains_stable.pdf",h=8,w=12)
ggsave("Univariate/Output/chains_stable.eps",h=8,w=12)
###############################################################################
pdf(file = "Univariate/Output/geweke.pdf", w=10,h=5)
geweke.plot(as.mcmc(D),pvalue = .05,nbins = 10,frac1 = .25,frac2 = .25)
dev.off()
geweke.diag(as.mcmc(D),frac1 = .25,frac2 = .25)
################################################################################

sub_disent_mu <- disentangled_mu[,inds$ix[1:5]]
sub_disent_s2 <- disentangled_sig2[,inds$ix[1:5]]

# ok, I keep the first five columns, and I reorder them oin terms of magnitude
idx = sort(apply(sub_disent_mu,2,median),index=T)
sub_disent_mu = sub_disent_mu[,idx$ix]
sub_disents2 = sub_disent_s2[,idx$ix]

par(mfrow=c(1,2))
matplot(sub_disent_mu,type="l")
matplot(sub_disent_s2,type="l")
par(mfrow=c(1,1))


m <- melt(sub_disent_mu)
v <- melt(sub_disent_s2)
D <- cbind(m,var = v[,3])
D <- tibble(D) %>% mutate(algo = "gibbs")
colnames(D)[3] <- "mean"

D1 = c()
for(j in 1:5){
  D1 = rbind(D1,cbind(1:5000,j,vi_sampled[,,j]))
}
str(D1)
colnames(D1) = c("Var1","Var2","mean","var")
D1 = as_tibble(D1) %>% mutate(algo="VB")
D2 = rbind(D,D1) %>% mutate(Var2=paste("Mixture component",Var2))

ggplot(D2)+
  geom_point(aes(x=mean,y=var,col=algo),alpha=.2,pch=".")+
  geom_density2d(aes(x=mean,y=var,col=algo),lineend = "butt",linemitre = 25)+
  scale_color_manual(values = c("turquoise4","coral1")) + theme_bw()+
  facet_wrap(~Var2,scales = "free")


ggplot(D2)+
  geom_point(aes(x=mean,y=var,col=algo),alpha=.1,pch=".")+
  geom_density2d(aes(x=mean,y=var,col=algo),lineend = "butt",linemitre = 25)+
  scale_color_manual(values = c("turquoise4","coral1")) + theme_bw()+
  facet_wrap(~Var2,scales = "free")+
  theme(legend.position = c(.85,.25))

sumD2 <- D2 %>% group_by(Var2,algo) %>% mutate(mmed=median(mean),
                                          vmed=median(var))

ggplot(sumD2)+
  geom_density(aes(x=mean,col=algo))+
  geom_vline(aes(xintercept = mmed,col=algo,lty=algo))+
  scale_color_manual(values = c(1,2)) + theme_bw()+
  facet_wrap(~Var2,scales = "free")+
  theme(legend.position = c(.85,.25))

ggplot(sumD2)+
  geom_density(aes(x=var,col=algo))+
  geom_vline(aes(xintercept = vmed,col=algo,lty=algo))+
  scale_color_manual(values = c(1,2)) + theme_bw()+
  facet_wrap(~Var2,scales = "free")+
  theme(legend.position = c(.85,.25))

# MCMC - VB
cols <- c("orangered","turquoise4")

##------------------------------------------------------------------------------
for(i in 1:5){
  mix1=D2 %>% filter(Var2 == paste0("Mixture component ",i))
  D = ggplot(mix1)+
  geom_point(aes(x=mean,y=var,col=algo),alpha=.1,size=.25)+
  stat_density2d(aes(x=mean,y=var,col=algo),
                 h = .12,n = 200,lwd=.8,bins=5)+
  scale_color_manual(values = cols) + theme_bw()+
  theme(legend.position = "none",text = element_text(size=15))

D <- D + 
  xlab("Posterior mean") + 
  ylab("Posterior variance")

R <- ggExtra::ggMarginal(
  p = D,
  type = 'density',
  margins = 'both',
  size = 5,lwd=.8,
  groupColour = TRUE,
  fill = 'transparent'
)
 R
 ggsave(plot = R,paste0("Univariate/Output/fiSAN_RUN",RUN,"_scatter",i,"_conf",SAMPLE_TO_RUN,".pdf"),height =6,width =7)
 ggsave(plot = R,filename = paste0("Univariate/Output/fiSAN_RUN",RUN,"_scatter",i,"_conf",SAMPLE_TO_RUN,".eps"),h=5,w=6,device = cairo_ps)
}

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------



SAMPLE_TO_RUN <- 3
sample_sizes= c(10,50,500,2500)#)0)


# Read runs ---------------------------------------------------------------

name <- paste0("Univariate/Runs/VI_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_VI")
vif_fiSAN <- readRDS(paste0(name,"_fiSAN.RDS"))

name <- paste0("Univariate/Runs/GIBBS_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_GIBBS")
gib_fiSAN  <- readRDS(paste0(name,"_fiSAN_OF.RDS"))


RUN <- 44 #
# sample from variational posterior
vi_sampled = sampling_vi_G(output =vif_fiSAN[[RUN]], nsim = 5000)

# remove label switching
one_model <- gib_fiSAN[[RUN]]
l1 <- ecr.iterative.1(one_model$sim$obs_cluster,K=one_model$params$maxL)

par(mfrow=c(1,1))
disentangled_mu <- t(sapply(1:nrow(one_model$sim$mu),function(x) one_model$sim$mu[x,l1$permutations[x,]]))
head(disentangled_mu)
disentangled_sig2 <- t(sapply(1:nrow(one_model$sim$sigma2),function(x) one_model$sim$sigma2[x,l1$permutations[x,]]))

# five chains has the lowest sd: they are the one that we are looking for!
inds = sort(apply(disentangled_mu,2,sd),ind = T)
plot(inds$x)


par(mfrow=c(1,2))
matplot(one_model$sim$mu[,inds$ix[1:5]],type="l")
matplot(disentangled_mu[,inds$ix[1:5]],type="l")
matplot(one_model$sim$sigma2[,inds$ix[1:5]],type="l")
matplot(disentangled_sig2[,inds$ix[1:5]],type="l")
par(mfrow=c(1,1))

# ---------------------------------------------------------------------------------------
plot(apply(one_model$sim$mu,2,sd))
inds_orig = sort(apply(one_model$sim$mu,2,sd),ind = T)
plot(inds_orig$x)

matplot(one_model$sim$mu[,inds_orig$ix[1:5]],type="l")
matplot(disentangled_mu[,inds$ix[1:5]],type="l")

Omu = reshape2::melt(one_model$sim$mu[,inds_orig$ix[1:5]])
Dmu = reshape2::melt(disentangled_mu[,inds$ix[1:5]])


Omu = as_tibble(Omu) %>% mutate(type = "Original chains")
Dmu = as_tibble(Dmu) %>% mutate(type = "Post-processed chains")
ALLmu = rbind(Omu,Dmu)

ggplot(ALLmu)+
  theme_bw()+
  geom_line(aes( x = Var1, y= value, col=factor(Var2), group = factor(Var2)))+
  facet_wrap(~type,nrow = 2,scales = "free_y")+
  scale_color_brewer(palette = "Set1")  +
  ylab("Posterior means") + xlab("MCMC iterations")+
  theme(text = element_text(size=16), legend.position = "none")
################################################################################
#ggview::ggview(h=8,w=12)
ggsave("Univariate/Output/labelsw_fiSAN_conf3.pdf",h=8,w=12)
ggsave("Univariate/Output/labelsw_fiSAN_conf3.eps",h=8,w=12)
################################################################################

# ----------------------------------------------------------------------------------------


sub_disent_mu <- disentangled_mu[,inds$ix[1:5]]
sub_disent_s2 <- disentangled_sig2[,inds$ix[1:5]]

# ok, I keep the first five columns, in order of magnitude
idx = sort(apply(sub_disent_mu,2,median),index=T)
sub_disent_mu = sub_disent_mu[,idx$ix]
sub_disents2 = sub_disent_s2[,idx$ix]

par(mfrow=c(1,2))
matplot(sub_disent_mu,type="l")
matplot(sub_disent_s2,type="l")
par(mfrow=c(1,1))

m = melt(sub_disent_mu)
v = melt(sub_disent_s2)
D <- cbind(m,var = v[,3])
D <- tibble(D) %>% mutate(algo = "gibbs")
colnames(D)[3] = "mean"

D1 = c()
for(j in 1:5){
  D1 = rbind(D1,cbind(1:5000,j,vi_sampled[,,j]))
}
str(D1)
colnames(D1) = c("Var1","Var2","mean","var")
D1 = as_tibble(D1) %>% mutate(algo="VB")
D2 = rbind(D,D1) %>% mutate(Var2=paste("Mixture component",Var2))

ggplot(D2)+
  geom_point(aes(x=mean,y=var,col=algo),alpha=.2,pch=".")+
  geom_density2d(aes(x=mean,y=var,col=algo),lineend = "butt",linemitre = 25)+
  scale_color_manual(values = c("turquoise4","coral1")) + theme_bw()+
  facet_wrap(~Var2,scales = "free")



# MCMC - VB
cols <- c("orangered","turquoise4")

##------------------------------------------------------------------------------
i=1
##------------------------------------------------------------------------------
mix1=D2 %>% filter(Var2 == paste0("Mixture component ",i))
D = ggplot(mix1)+
  geom_point(aes(x=mean,y=var,col=algo),alpha=.1,size=.25)+
  stat_density2d(aes(x=mean,y=var,col=algo),
                 h = .05,n = 200,lwd=.8,bins=5)+
  scale_color_manual(values = cols) + theme_bw()+
  theme(legend.position = "bottom",text = element_text(size=15))
D



##------------------------------------------------------------------------------
for(i in 1:5){
  mix1=D2 %>% filter(Var2 == paste0("Mixture component ",i))
  D = ggplot(mix1)+
    geom_point(aes(x=mean,y=var,col=algo),alpha=.1,size=.25)+
    stat_density2d(aes(x=mean,y=var,col=algo),
                   h = .02,n = 200,lwd=.8,bins=5)+
    scale_color_manual(values = cols) + theme_bw()+
    theme(legend.position = "none",text = element_text(size=15))
  
  #  facet_wrap(~Var2,scales = "free")
  #ggMarginalGadget(D)
  D <- D + 
    xlab("Posterior mean") + 
    ylab("Posterior variance")
  
  R <- ggExtra::ggMarginal(
    p = D,
    type = 'density',
    margins = 'both',
    size = 5,lwd=.8,
    groupColour = TRUE,
    fill = 'transparent'
  )
  R
  ggsave(plot = R,paste0("Univariate/Output/fiSAN_RUN",RUN,"_scatter",i,"_conf",SAMPLE_TO_RUN,".pdf"),height =6,width =7)
  ggsave(plot = R,filename = paste0("Univariate/Output/fiSAN_RUN",RUN,"_scatter",i,"_conf",SAMPLE_TO_RUN,".eps"),h=5,w=6,device = cairo_ps)
}
















################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

# Same analysis with fSAN -------------------------------------------------

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


library(tidyverse)
library(patchwork)
library(reshape2)
library(label.switching)
library(ggExtra)

SAMPLE_TO_RUN <- 2
sample_sizes= c(10,50,500,2500)#)0)


# Read runs ---------------------------------------------------------------

name <- paste0("Univariate/Runs/VI_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_VI")
vif_fSAN <- readRDS(paste0(name,"_fSAN.RDS"))

name <- paste0("Univariate/Runs/GIBBS_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_GIBBS")
gib_fSAN  <- readRDS(paste0(name,"_fSAN_OF.RDS"))


### functions to sample VI posterior--------------------------------------------
sampling_vi_G <- function(output, nsim){
  
  check  <- output$results$params$gamma0
  ind    <- which(!output$results$sim$theta_l[,4]==check)
  active <- as_tibble(output$results$sim$theta_l[ind,]) %>% arrange(V1)
  active <- as.matrix(active)
  RES  = array(NA, c(nsim,2,nrow(active)))
  for(i in 1:nrow(active)){
    sig2 <- 1/rgamma(nsim, active[i,3], active[i,4])
    mu   <- rnorm(nsim, active[i,1], sqrt(sig2/active[i,2]))
    RES[,,i] <- cbind(mu = mu, sig2 = sig2)
  }
  return(RES)
}


plot(unlist(lapply(vif_fSAN,function(x) max(x$results$sim$Elbo_val))))
which.max(unlist(lapply(vif_fSAN,function(x) max(x$results$sim$Elbo_val))))

RUN <- 1
# sample from variational posterior
vi_sampled = sampling_vi_G(output =vif_fSAN[[RUN]], nsim = 5000)

# remove label switching
one_model <- gib_fSAN[[RUN]]
plot(ts(one_model$sim$mu)[,1:10])
plot(ts(one_model$sim$mu)[,1:10+10])
l1 <- ecr.iterative.1(one_model$sim$obs_cluster,K=one_model$params$maxL)

par(mfrow=c(1,1))
disentangled_mu <- t(sapply(1:nrow(one_model$sim$mu),function(x) one_model$sim$mu[x,l1$permutations[x,]]))
head(disentangled_mu)
disentangled_sig2 <- t(sapply(1:nrow(one_model$sim$sigma2),function(x) one_model$sim$sigma2[x,l1$permutations[x,]]))




# five chains has the lowest sd: they are the one that we are looking for!
inds = sort(apply(disentangled_mu,2,sd),ind = T)

par(mfrow=c(1,2))
# not really representative
matplot(one_model$sim$mu[,inds$ix[1:5]],type="l")
matplot(disentangled_mu[,inds$ix[1:5]],type="l")
matplot(one_model$sim$sigma2[,inds$ix[1:5]],type="l")
matplot(disentangled_sig2[,inds$ix[1:5]],type="l")
par(mfrow=c(1,1))

# ---------------------------------------------------------------------------------------
plot(apply(one_model$sim$mu,2,sd))
inds_orig = sort(apply(one_model$sim$mu,2,sd),ind = T)

matplot(one_model$sim$mu[,inds_orig$ix[1:5]],type="l")
matplot(disentangled_mu[,inds$ix[1:5]],type="l")

Omu = reshape2::melt(one_model$sim$mu[,inds_orig$ix[1:5]])
Dmu = reshape2::melt(disentangled_mu[,inds$ix[1:5]])

plot(ts(one_model$sim$mu[,inds_orig$ix[1:5]]))

Omu = as_tibble(Omu) %>% mutate(type = "Original chains")
Dmu = as_tibble(Dmu) %>% mutate(type = "Post-processed chains")
ALLmu = rbind(Omu,Dmu)

ggplot(ALLmu)+
  theme_bw()+
  geom_line(aes( x = Var1, y= value, col=factor(Var2), group = factor(Var2)))+
  facet_wrap(~type,nrow = 2,scales = "free_y")+
  scale_color_brewer(palette = "Set1")  +
  ylab("Posterior means") + xlab("MCMC iterations")+
  theme(text = element_text(size=16), legend.position = "none")
################################################################################
#ggview::ggview(h=8,w=12)
ggsave("Univariate/Output/labelsw_fsan_conf2.pdf",h=8,w=12)
ggsave("Univariate/Output/labelsw_fsan_conf2.eps",h=8,w=12)
################################################################################

# ----------------------------------------------------------------------------------------

sub_disent_mu <- disentangled_mu[,inds$ix[1:5]]
sub_disent_s2 <- disentangled_sig2[,inds$ix[1:5]]

# ok, I keep the first five columns, in order of magnitude
idx = sort(apply(sub_disent_mu,2,median),index=T)
sub_disent_mu = sub_disent_mu[,idx$ix]
sub_disents2 = sub_disent_s2[,idx$ix]

par(mfrow=c(1,2))
matplot(sub_disent_mu,type="l")
matplot(sub_disent_s2,type="l")
par(mfrow=c(1,1))


m = melt(sub_disent_mu)
v = melt(sub_disent_s2)
D <- cbind(m,var = v[,3])
D <- tibble(D) %>% mutate(algo = "gibbs")
colnames(D)[3] = "mean"

D1 = c()
for(j in 1:5){
  D1 = rbind(D1,cbind(1:5000,j,vi_sampled[,,j]))
}
str(D1)

colnames(D1) = c("Var1","Var2","mean","var")
plot(D1[,3])
D1 = as_tibble(D1) %>% mutate(algo="VB")
D2 = rbind(D,D1) %>% mutate(Var3=paste("Mixture component",Var2))

plot(D2$Var2)


RESULTS <- list(mu = one_model$sim$mu,
                sig2 = one_model$sim$sigma2,
                dis_mu = disentangled_mu,
                dis_sig2 = disentangled_sig2,
                vi_sample = vi_sampled,
                for_plot = D2)

cols <- c("orangered","turquoise4")

ggplot(D2)+
  geom_point(aes(x=mean,y=var,col=algo),alpha=.2,pch=".")+
  geom_density2d(aes(x=mean,y=var,col=algo),lineend = "butt",linemitre = 25)+
  scale_color_manual(values = cols) + theme_bw()+
  facet_wrap(~Var2,scales = "free")


ggplot(D2)+
  geom_point(aes(x=mean,y=var,col=algo),alpha=.1,pch=".")+
  geom_density2d(aes(x=mean,y=var,col=algo),lineend = "butt",linemitre = 25)+
  scale_color_manual(values =  cols) + theme_bw()+
  facet_wrap(~Var2,scales = "free")+
  theme(legend.position = c(.85,.25))

sumD2 <- D2 %>% group_by(Var2,algo) %>% mutate(mmed=median(mean),
                                               vmed=median(var))


# MCMC - VB
cols <- c("orangered","turquoise4")
dev.off()
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
# trial
mix1=D2 %>% filter(Var3 == paste0("Mixture component ",Var2))
D = ggplot(mix1)+
  geom_point(aes(x=mean,y=var,col=algo),alpha=.1,size=.25)+
  stat_density2d(aes(x=mean,y=var,col=algo),
                 h = .1,n = 200,lwd=.8,bins=5)+
  scale_color_manual(values = cols) + theme_bw()+
  theme(legend.position = "bottom",text = element_text(size=15))
D


##------------------------------------------------------------------------------
for(i in 1:5){
  mix1=D2 %>% filter(Var3 == paste0("Mixture component ",i))
  D = ggplot(mix1)+
    geom_point(aes(x=mean,y=var,col=algo),alpha=.1,size=.25)+
    stat_density2d(aes(x=mean,y=var,col=algo),
                   h = .1,n = 200,lwd=.8,bins=5)+
    scale_color_manual(values = cols) + theme_bw()+
    theme(legend.position = "none",text = element_text(size=15))
  
  #  facet_wrap(~Var2,scales = "free")
  #ggMarginalGadget(D)
  D <- D + 
    xlab("Posterior mean") + 
    ylab("Posterior variance")
  
  R <- ggExtra::ggMarginal(
    p = D,
    type = 'density',
    margins = 'both',
    size = 5,lwd=.8,
    groupColour = TRUE,
    fill = 'transparent'
  )
  R
  ggsave(plot = R,paste0("Univariate/Output/fSAN_RUN",RUN,"_scatter",i,"_conf",SAMPLE_TO_RUN,".pdf"),height =6,width =7)
  ggsave(plot = R,filename = paste0("Univariate/Output/fSAN_RUN",RUN,"_scatter",i,"_conf",SAMPLE_TO_RUN,".eps"),h=5,w=6,
         device = cairo_ps)
}

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------




SAMPLE_TO_RUN <- 3
sample_sizes= c(10,50,500,2500)#)0)


# Read runs ---------------------------------------------------------------

name <- paste0("Univariate/Runs/VI_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_VI")
vif_fSAN <- readRDS(paste0(name,"_fSAN.RDS"))

name <- paste0("Univariate/Runs/GIBBS_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_GIBBS")
gib_fSAN  <- readRDS(paste0(name,"_fSAN_OF.RDS"))



RUN <- 10
# sample from variational posterior
# vi_sampled = sampling_vi_G(output =vif_fSAN[[RUN]], nsim = 5000)
## ATTENTION the function works only if the number of estimated OC is exact
nsim = 5000
output =vif_fSAN[[RUN]]
check  <- output$results$params$gamma0
ind    <- which(!output$results$sim$theta_l[,4]==check)
active <- as_tibble(output$results$sim$theta_l[ind,]) %>% arrange(V1) %>% filter(V2>20)
active <- as.matrix(active)
RES  = array(NA, c(nsim,2,nrow(active)))
for(i in 1:nrow(active)){
  sig2 <- 1/rgamma(nsim, active[i,3], active[i,4])
  mu   <- rnorm(nsim, active[i,1], sqrt(sig2/active[i,2]))
  RES[,,i] <- cbind(mu = mu, sig2 = sig2)
}
vi_sampled <- RES


# remove label switching
one_model <- gib_fSAN[[RUN]]
l1 <- ecr.iterative.1(one_model$sim$obs_cluster,K=one_model$params$maxL)

par(mfrow=c(1,1))
disentangled_mu <- t(sapply(1:nrow(one_model$sim$mu),function(x) one_model$sim$mu[x,l1$permutations[x,]]))
head(disentangled_mu)
disentangled_sig2 <- t(sapply(1:nrow(one_model$sim$sigma2),function(x) one_model$sim$sigma2[x,l1$permutations[x,]]))



# five chains has the lowest sd: they are the one that we are looking for!
inds = sort(apply(disentangled_mu,2,sd),ind = T)

par(mfrow=c(1,2))
# not really representative
matplot(one_model$sim$mu[,inds$ix[1:5]],type="l")
matplot(disentangled_mu[,inds$ix[1:5]],type="l")
matplot(one_model$sim$sigma2[,inds$ix[1:5]],type="l")
matplot(disentangled_sig2[,inds$ix[1:5]],type="l")
par(mfrow=c(1,1))

# ---------------------------------------------------------------------------------------
plot(apply(one_model$sim$mu,2,sd))
inds_orig = sort(apply(one_model$sim$mu,2,sd),ind = T)

matplot(one_model$sim$mu[,inds_orig$ix[1:5]],type="l")
matplot(disentangled_mu[,inds$ix[1:5]],type="l")

Omu = reshape2::melt(one_model$sim$mu[,inds_orig$ix[1:5]])
Dmu = reshape2::melt(disentangled_mu[,inds$ix[1:5]])

plot(ts(one_model$sim$mu[,inds_orig$ix[1:5]]))

Omu = as_tibble(Omu) %>% mutate(type = "Original chains")
Dmu = as_tibble(Dmu) %>% mutate(type = "Post-processed chains")
ALLmu = rbind(Omu,Dmu)

ggplot(ALLmu)+
  theme_bw()+
  geom_line(aes( x = Var1, y= value, col=factor(Var2), group = factor(Var2)))+
  facet_wrap(~type,nrow = 2,scales = "free_y")+
  scale_color_brewer(palette = "Set1")  +
  ylab("Posterior means") + xlab("MCMC iterations")+
  theme(text = element_text(size=16), legend.position = "none")
################################################################################
#ggview::ggview(h=8,w=12)
ggsave("Univariate/Output/labelsw_fsan_conf3.pdf",h=8,w=12)
ggsave("Univariate/Output/labelsw_fsan_conf3.eps",h=8,w=12)
################################################################################
# ----------------------------------------------------------------------------------------

sub_disent_mu <- disentangled_mu[,inds$ix[1:5]]
sub_disent_s2 <- disentangled_sig2[,inds$ix[1:5]]

# sub_disent_mu <- one_model$sim$mu[,inds$ix[1:5]]
# sub_disent_s2 <- one_model$sim$sigma2[,inds$ix[1:5]]

# ok, I keep the first five columns, in order of magnitude
idx = sort(apply(sub_disent_mu,2,median),index=T)
sub_disent_mu = sub_disent_mu[,idx$ix]
sub_disents2 = sub_disent_s2[,idx$ix]

par(mfrow=c(1,2))
matplot(sub_disent_mu,type="l")
matplot(sub_disent_s2,type="l")
par(mfrow=c(1,1))


m = melt(sub_disent_mu)
v = melt(sub_disent_s2)
D <- cbind(m,var = v[,3])
D <- tibble(D) %>% mutate(algo = "gibbs")
colnames(D)[3] = "mean"

D1 = c()
for(j in 1:5){
  D1 = rbind(D1,cbind(1:5000,j,vi_sampled[,,j]))
}
str(D1)
colnames(D1) = c("Var1","Var2","mean","var")
D1 = as_tibble(D1) %>% mutate(algo="VB")
D2 = rbind(D,D1) %>% mutate(Var2=paste("Mixture component",Var2))

RESULTS <- list(mu = one_model$sim$mu,
                sig2 = one_model$sim$sigma2,
                dis_mu = disentangled_mu,
                dis_sig2 = disentangled_sig2,
                vi_sample = vi_sampled,
                for_plot = D2)

ggplot(D2)+
  geom_point(aes(x=mean,y=var,col=algo),alpha=.2,pch=".")+
  geom_density2d(aes(x=mean,y=var,col=algo),lineend = "butt",linemitre = 25)+
  scale_color_manual(values = cols) + theme_bw()+
  facet_wrap(~Var2,scales = "free")


# MCMC - VB
dev.off()
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
i=2
mix1=D2 %>% filter(Var2 == paste0("Mixture component ",i))
D = ggplot(mix1)+
  geom_point(aes(x=mean,y=var,col=algo),alpha=.1,size=.25)+
  stat_density2d(aes(x=mean,y=var,col=algo),
                 h = .08,n = 200,lwd=.8,bins=5)+
  scale_color_manual(values = cols) + theme_bw()+
  theme(legend.position = "bottom",text = element_text(size=15))
D

##------------------------------------------------------------------------------
for(i in 1:5){
  mix1=D2 %>% filter(Var2 == paste0("Mixture component ",i))
  D = ggplot(mix1)+
    geom_point(aes(x=mean,y=var,col=algo),alpha=.1,size=.25)+
    stat_density2d(aes(x=mean,y=var,col=algo),
                   h = .08,n = 200,lwd=.8,bins=5)+
    scale_color_manual(values = cols) + theme_bw()+
    theme(legend.position = "none",text = element_text(size=15))
  
  #  facet_wrap(~Var2,scales = "free")
  #ggMarginalGadget(D)
  D <- D + 
    xlab("Posterior mean") + 
    ylab("Posterior variance")
  
  R <- ggExtra::ggMarginal(
    p = D,
    type = 'density',
    margins = 'both',
    size = 5,lwd=.8,
    groupColour = TRUE,
    fill = 'transparent'
  )
  R
  ggsave(plot = R,paste0("Univariate/Output/fSAN_RUN",RUN,"_scatter",i,"_conf",SAMPLE_TO_RUN,".pdf"),height =6,width =7)
  ggsave(plot = R,filename = paste0("Univariate/Output/fSAN_RUN",RUN,"_scatter",i,"_conf",SAMPLE_TO_RUN,".eps"),h=5,w=6,
         device = cairo_ps)
}

