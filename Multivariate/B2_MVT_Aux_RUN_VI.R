res_fisan = list()
Ofisan <- list()

for(i in 1:sims){
  
  
  cat(paste0("fisan: dataset number ",i,"\n"))
  t1 = Sys.time()
  ALL_FISAN <-  list()
  for(simsim in 1:50){
    ALL_FISAN[[simsim]] = SANvimvt::variational_mvt_fiSAN(y = y[[i]],
                                                          group = group[[i]],
                                                          maxL = maxL,
                                                          maxK = maxK,
                                                          m0 = rep(0,ncol(y[[i]])),
                                                          beta0 = 0.01,
                                                          nu0 = ncol(y[[i]])+5,
                                                          W0 = diag(ncol(y[[i]])), 
                                                          warmstart = warmstart,
                                                          ###############################
                                                          verbose = FALSE,
                                                          maxSIM = 500,
                                                          epsilon = 1e-4,
                                                          ###############################
                                                          conc_hyperpar = c(hyp_alpha1,
                                                                            hyp_alpha2),
                                                          beta_bar = 0.05,
                                                          seed = root*i*simsim)
    if(simsim %% 10 == 0) cat(paste("Done with starting point #",simsim,"\n"))
  }
  
  t2 = Sys.time()
  Ofisan$size <- object.size(ALL_FISAN)
  LT = unlist(lapply(ALL_FISAN, function(x) as.numeric(x$time, units="secs" )))
  Ofisan$time_50 <- LT
  Ofisan$time_parallel <- t2-t1
  class(ALL_FISAN) = "multistart"
  inds <- which.max(unlist(lapply(ALL_FISAN,function(x) max(x$sim$Elbo_val))))
  Ofisan$results <- ALL_FISAN[[inds]]
  cat(green(" \u2713 \n"))
  res_fisan[[i]] <-  Ofisan
  

  
  cat(paste0("------------------------------\n"))
  
}

cat("Saving progress...\n")
saveRDS(res_fisan,paste0(name,"_VI_MVT_fiSAN.RDS"))
