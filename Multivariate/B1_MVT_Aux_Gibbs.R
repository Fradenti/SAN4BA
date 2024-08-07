
one_sweep_mvtfiSAN_of <- function(i, nrep, burn, root, y, group, Hyperpar_list){
  
  mvt_overfiSAN <- SANplemvt::sample_mvt_overfiSAN(nrep = nrep,
                                              burn = burn, 
                                              y = y[[i]], group = group[[i]], 
                                              maxK = Hyperpar_list$maxK, 
                                              maxL = Hyperpar_list$maxL, 
                                              m0 = rep(0,ncol(y[[i]])),
                                              beta0 = 0.01,
                                              nu0 = ncol(y[[i]])+5,
                                              W0 = diag(ncol(y[[i]])), 
                                              eps_beta = 1,
                                              hyp_alpha1 = 1, 
                                              hyp_alpha2 = 1,
                                              hyp_beta = 10, 
                                              beta =  0.05, 
                                              seed = i * root,
                                              warmstart = T,
                                              S_start = c(1:length(unique(group[[i]])))-1, 
                                              nclus_start = 5)
  mvt_overfiSAN[["size"]] <- object.size(mvt_overfiSAN)

  return(mvt_overfiSAN)
  
}
