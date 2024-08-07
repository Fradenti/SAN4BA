res_fsan = res_fisan = list()
Ofisan <- Ofsan <- Ocam <- list()
# no parallel

for(i in 1:sims){
  
  cat(paste0("fsan: simulation number ",i))

  tt1 = Sys.time();
  ALL_FSAN <- list()
  for(simsim1 in 1:50){
    ALL_FSAN[[simsim1]] <-  SANvi::variational_fSAN(y = y[[i]],
                                      group = group[[i]],
                                      maxL = maxL,
                                      maxK = maxK,
                                      warmstart = warmstart,
                                      alpha_bar = alpha_bar,
                                      beta_bar = beta_bar,
                                      epsilon = 1e-4,
                                      verbose = F,
                                      maxSIM = 500,
                                      seed = root*simsim1
                                      )
    if(simsim1 %% 10 == 0 ) cat(paste("-",simsim1,"-"))   
  }
  tt2 = Sys.time() ;
  
  Ofsan$size <- object.size(ALL_FSAN)
  LT = unlist(lapply(ALL_FSAN, function(x) as.numeric(x$time, units="secs" )))
  Ofsan$time_50 <- LT
  Ofsan$time_parallel <- tt2-tt1
  ind_best_fsan <- which.max(unlist(lapply(ALL_FSAN, function(u) max(u$sim$Elbo_val))))
  
  Ofsan$results <- (ALL_FSAN[[ind_best_fsan]])
  cat(green(" \u2713 \n"))
  res_fsan[[i]]  <-  Ofsan
  
  if(i %% step == 0 | i == sims){
    cat("Saving progress...\n")
    saveRDS(res_fsan, paste0(name,"_VI_fSAN2.RDS"))
  }  
  
  cat(paste0("fisan: simulation number ",i))
  tti1 = Sys.time(); #tic()
  ALL_FISAN <- list()
  for(simsim2 in 1:50){
    ALL_FISAN[[simsim2]] <-  SANvi::variational_fiSAN(y = y[[i]],
                                                   group = group[[i]],
                                                   maxL = maxL,
                                                   maxK = maxK,
                                                   warmstart = warmstart,
                                                   conc_hyperpar = c(hyp_alpha1,
                                                                     hyp_alpha2),
                                                   beta_bar = beta_bar,
                                                   epsilon = 1e-4,
                                                   verbose = F,
                                                   maxSIM = 500,
                                                   seed = root*simsim2)
    if(simsim2 %% 10 ==0 ) cat(paste("-",simsim2,"-"))
  }
  tti2 = Sys.time() ; #time3 = toc()
  
  Ofisan$size <- object.size(ALL_FISAN)
  LT = unlist(lapply(ALL_FISAN, function(x) as.numeric(x$time,units="secs" )))
  Ofisan$time_50 <- LT
  Ofisan$time_parallel <- tti2-tti1
  ind_best_fisan <- which.max(unlist(lapply(ALL_FISAN, function(u) max(u$sim$Elbo_val))))
  
  Ofisan$results <- (ALL_FISAN[[ind_best_fisan]])
  cat(green(" \u2713 \n"))
  res_fisan[[i]] <-  Ofisan
  
  if(i %% step == 0 | i == sims){
    cat("Saving progress...\n")
    saveRDS(res_fisan,paste0(name,"_VI_fiSAN2.RDS"))
  }
  
  
  cat(paste0("------------------------------\n"))
  
  
}



