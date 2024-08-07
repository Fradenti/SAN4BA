source("Univariate/Helpers/postprocessing_output.R")
library(SANvi)
library(tidyverse)
sims <- 50
sample_sizes = c(10,50,500,2500)#0)


summarize_densities_expectations = function (output, filter = TRUE, ...) 
{
  if ((output$results$model)[1] == "CAM") {
    post_weights = SANvi:::post_sb_weight(alk = output$results$sim$a_bar_lk, blk = output$results$sim$b_bar_lk)
  }
  else if (output$results$model == "fSAN") {
    post_weights = apply(output$results$sim$beta_bar_lk, 2, function(x) x/sum(x))
  }
  else if (output$results$model == "fiSAN") {
    post_weights = apply(output$results$sim$beta_bar_lk, 2, function(x) x/sum(x))
  }
  else {
    stop("Provide valid VB output")
  }
  ######################### sovrascrivo da pacchetto per fare match di intervalli  
  seq <- seq(-10, 10, length.out = 500)
  ##########################
  L <- output$results$params$L
  if (filter) {
    inds <- which(colSums(output$results$sim$RHO) >= 0.5)
  }else {
    inds <- rep(TRUE, output$results$params$K)
  }
  exp_sigma2_star <- output$results$sim$theta_l[, 4]/(output$results$sim$theta_l[, 
                                                                 3] - 1)
  exp_mu <- output$results$sim$theta_l[, 1]
  Dens = apply(as.matrix(post_weights[, inds]), 2, function(w) {
    sapply(seq, function(g) {
      sum(w * dnorm(g, exp_mu, sqrt(exp_sigma2_star)))
    })
  })
  return(list(x = seq, Dens = Dens))
}


stickBre <- function(v){
  p1 <- v
  p2 <- c(1,1-v[-length(v)])
  p3 <- log(p1) + cumsum(log(p2))
exp(p3)
}


ord_dens <- function(temp){
  L <-  apply(temp$Dens,2, function(r) mean(temp$x * r))
  ran <- sort(L,index=T)  
  temp$Dens[,ran$ix]
}

# dens --------------------------------------------------------------------

for(SAMPLE_TO_RUN in 1:4){
  name <- paste0("Univariate/Runs/VI_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_VI")
  
  res_fsan  <- readRDS(paste0(name,"_fSAN.RDS"))
  res_fisan <- readRDS(paste0(name,"_fiSAN.RDS"))
  dens_fisan=dens_fsan=c()
  

  for(i in 1:sims){
    
    temp_fsan = summarize_densities_expectations(res_fsan[[i]])
    ord = ord_dens(temp_fsan)
    dens_fsan = rbind(dens_fsan, cbind(reshape2::melt(ord),
                                       i,
                                       x=rep(temp_fsan$x,ncol(ord))) )
    
    temp_fisan = summarize_densities_expectations(res_fisan[[i]])
    ord = ord_dens(temp_fisan)
    
    dens_fisan = rbind(dens_fisan, cbind(reshape2::melt(ord),
                                         i,
                                         x=rep(temp_fisan$x,ncol(ord))) )
    
    if(i==sims){
      dens_fsan =  as.data.frame(dens_fsan)
      dens_fsan[,"type"] <- "fsan"
      dens_fisan =  as.data.frame(dens_fisan)
      dens_fisan[,"type"] <- "fisan"
      D_all <- rbind(dens_fsan,dens_fisan)
    }
    cat(i)
  }
  name_res <- paste0("Univariate/Results/VI_RDS/scen1_samplesize",sample_sizes[SAMPLE_TO_RUN],"_VI")
  saveRDS(D_all,paste0(name_res,"_posterior_mean_dens_VI.RDS"))
cat(i)
}
