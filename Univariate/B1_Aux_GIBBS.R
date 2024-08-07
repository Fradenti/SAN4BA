
one_sweep_CAM <- function(i, nrep, burn, root, y, group, Hyperpar_list,
                          burnin,
                          remove_bi = FALSE){
  
  run_CAM = SANple::sample_CAM(nrep,burn, y[[i]], group[[i]],
                       maxK = Hyperpar_list$maxK,
                       maxL = Hyperpar_list$maxL,
                       m0 = Hyperpar_list$m0,
                       tau0 = Hyperpar_list$tau0,
                       lambda0 = Hyperpar_list$lambda0,
                       gamma0  = Hyperpar_list$gamma0,
                       hyp_alpha1 = Hyperpar_list$hyp_alpha1_DP,
                       hyp_alpha2 = Hyperpar_list$hyp_alpha2_DP,
                       hyp_beta1 = Hyperpar_list$hyp_beta1_DP,
                       hyp_beta2 = Hyperpar_list$hyp_beta2_DP,
                       seed = i * root,
                       S_start = c(1:length(unique(group[[i]])))-1, 
                       nclus_start = 5, # Hyperpar_list$maxL,
                       progress = FALSE)
  run_CAM[["size"]] <- object.size(run_CAM)
  if(remove_bi){
    run_CAM <- discard_burnin(run_CAM,burnin = burnin)
  }
  
  return(run_CAM)
  
}

discard_burnin <- function(output, burnin){
  
  output_reduced <- output
  output_reduced$sim$mu <- output_reduced$sim$mu[-burnin,]
  output_reduced$sim$sigma2 <- output_reduced$sim$sigma2[-burnin,]
  output_reduced$sim$obs_cluster <- output_reduced$sim$obs_cluster[-burnin,]
  output_reduced$sim$distr_cluster <- output_reduced$sim$distr_cluster[-burnin,]
  output_reduced$sim$pi    <- output_reduced$sim$pi[-burnin,]
  output_reduced$sim$omega <- output_reduced$sim$omega[,,-burnin]
  output_reduced$sim$alpha <- output_reduced$sim$alpha[-burnin]
  output_reduced$sim$beta  <- output_reduced$sim$beta[-burnin]
  output_reduced$sim$maxK  <- output_reduced$sim$maxK[-burnin]
  output_reduced$sim$maxL  <- output_reduced$sim$maxL[-burnin]
  output_reduced$sim$K  <- output_reduced$sim$K[-burnin]
  output_reduced$sim$L  <- output_reduced$sim$L[-burnin]
  
  return(output_reduced)
}


one_sweep_fCAM <- function(i, nrep, burn, root, y, group, Hyperpar_list,
                           burnin,
                           remove_bi = FALSE){
  
  run_fCAM = SANple::sample_fSAN(nrep,burn, y[[i]], group[[i]],
                              maxK = Hyperpar_list$maxK,
                              maxL = Hyperpar_list$maxL,
                              m0 = Hyperpar_list$m0,
                              tau0 = Hyperpar_list$tau0,
                              lambda0 = Hyperpar_list$lambda0,
                              gamma0  = Hyperpar_list$gamma0,
                              hyp_alpha1 = Hyperpar_list$hyp_alpha1_Dir,
                              hyp_alpha2 = Hyperpar_list$hyp_alpha2_Dir,
                              hyp_beta1 = Hyperpar_list$hyp_beta1_Dir,
                              hyp_beta2 = Hyperpar_list$hyp_beta2_Dir,
                              eps_beta  = Hyperpar_list$eps_beta,
                              eps_alpha = Hyperpar_list$eps_alpha,
                              seed = i * root,warmstart = T,
                              S_start = c(1:length(unique(group[[i]])))-1, 
                              nclus_start =  5,
                              progress = FALSE )
  run_fCAM[["size"]] <- object.size(run_fCAM)
  
  if(remove_bi){
    run_fCAM <- discard_burnin(run_fCAM, burnin)
  }
  
  return(run_fCAM)
  
}


one_sweep_fiSAN_of <- function(i, nrep, burn,root, y, group, Hyperpar_list,
                            burnin,
                            remove_bi = FALSE){
  
  run_fiSAN <- SANple::sample_fiSAN_sparsemix(nrep = nrep,burn = burn,
                            y = y[[i]],
                            group = group[[i]],
                            maxK = Hyperpar_list$maxK,
                            maxL = Hyperpar_list$maxL,
                            m0 = Hyperpar_list$m0,
                            tau0 = Hyperpar_list$tau0,
                            lambda0 = Hyperpar_list$lambda0,
                            gamma0  = Hyperpar_list$gamma0,
                            hyp_alpha1 = Hyperpar_list$hyp_alpha1_DP,
                            hyp_alpha2 = Hyperpar_list$hyp_alpha2_DP,
                            beta = .05,
                            eps_beta  = Hyperpar_list$eps_beta,
                            seed = i * root,
                            S_start = c(1:length(unique(group[[i]])))-1, 
                            nclus_start =  5,
                            progress = T)
  run_fiSAN[["size"]] <- object.size(run_fiSAN)
  
  if(remove_bi){
    run_fiSAN <- discard_burnin(run_fiSAN, burnin = burnin)
  }
  
  return(run_fiSAN)
  
}


one_sweep_fSAN_of <- function(i, nrep, burn, root, y, group, Hyperpar_list,
                               burnin,
                               remove_bi = FALSE){
  
    run_fSAN <- SANple::sample_fSAN_sparsemix(nrep = nrep, burn =burn,
                                      y = y[[i]],
                                      group = group[[i]],
                                      maxK = Hyperpar_list$maxK,
                                      maxL = Hyperpar_list$maxL,
                                      m0 = Hyperpar_list$m0,
                                      tau0 = Hyperpar_list$tau0,
                                      lambda0 = Hyperpar_list$lambda0,
                                      gamma0  = Hyperpar_list$gamma0,
                                      eps_beta  = Hyperpar_list$eps_beta,
                                      seed = i * root,
                                      alpha = .05,beta = .05,
                                      S_start = c(1:length(unique(group[[i]])))-1, 
                                      nclus_start =  5,
                                      progress = T)
  run_fSAN[["size"]] <- object.size(run_fSAN)
  
  if(remove_bi){
    run_fSAN <- discard_burnin(run_fSAN, burnin = burnin)
  }
  
  return(run_fSAN)
  
}
