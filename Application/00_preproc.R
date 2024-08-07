
# Libraries ----------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(SANvimvt)
theme_set(theme_bw())

# Useful functions -------------------------------------------------------------
logit = function(x)  {
  x = ifelse(x==1,.999,ifelse(x==0,1e-3,x))
  log(x/(1-x))
}

qn <- function(x){
  x = ifelse(x==1,.999,ifelse(x==0,1e-3,x))
  qnorm(x)
}

normalize <- function(x){
  (x-min(x))/(max(x) + .01 -min(x)) # avoid exactly 1
}

# Hyperpar ---------------------------------------------------------------------

L = 35
K = 30
beta_dir = 0.05

# Import and transform data ----------------------------------------------------

data_sub <- readRDS("Application/data_sub.RDS")

y_group <- as.numeric(as.factor(data_sub$FA))

# I keep only numeric variables
num <- select_if(data_sub, is.numeric) %>%                
  dplyr::select(-mode,-explicit,-key,-year,-valence,-count)
pheatmap::pheatmap(cor(num))

# scale the data in 0-1
norm_num <- apply(num,2,normalize)
colnames(norm_num)
summary(norm_num)
pheatmap::pheatmap(cor(norm_num))
inds = sample(1:20270,5000)
# pairs(norm_num[inds,],pch=".",col=y_group)
colnames(norm_num)

# keep the three interesting variables: duration, energy, speechiness
sub_norm <- norm_num[,c(3,4,9)]
### remove outliers
ind_problem1 <- which(sub_norm[,1]>.99)
data_sub[which(sub_norm[,1]>.99),]
sub_norm = sub_norm[-ind_problem1,]
y_group = y_group[-ind_problem1]
## adjust zero speech
sub_norm2 <- apply(sub_norm,2,function(x) ifelse(x==0,min(x[x>0]),x))

# map to normal
y_obser <- apply(sub_norm2,2,qn)
set.seed(1234)
inds_for_plotting = sample(1:nrow(y_obser),5000)
# pairs(y_obser[inds_for_plotting,],pch=".",col=y_group)

# removed duplicated entries
inds_problem2 <- which(duplicated(y_obser))
y_obser2 <- y_obser[-inds_problem2,]
y_group2 <- y_group[-inds_problem2]

# Marginal distributions 

hist(y_obser2 [,1],breaks = 100)
hist(y_obser2 [,2],breaks = 100)
hist(y_obser2 [,3],breaks = 100)
plot(y_obser2 [,1],pch=".")
plot(y_obser2 [,2],pch=".")
plot(y_obser2 [,3],pch=".")

pairs(y_obser2,col=y_group2,pch=".")

library(GGally)

ggpairs(as_tibble(y_obser2),upper = list("box_no_facet"),diag = list("histogram"),mapping = aes(col=factor(y_group2)))


Yg = list()
Nj = numeric(max(y_group2))
for(i in 1:max(y_group2)){
  inds = which(y_group2==i)
  Yg[[i]] = y_obser2[inds,]
  Nj[i] = length(inds)
}

# here, it is a good moment to run SANvi! Source script 01. Otherwise, skip to script 02!
