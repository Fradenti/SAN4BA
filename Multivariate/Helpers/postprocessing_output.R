parallel_NFD_oc_extract <- function(res,ocgt){
  unlist(
    lapply(res,function(x)
    {  
      psm <- salso::psm(x$sim$obs_cluster,nCores = Ncores)
      gtpsm <- salso::psm(rbind(ocgt, ocgt))
      mean((psm-gtpsm)^2)
    }
    )
  )
}

parallel_NFD_dc_extract <- function(res, dcgt){
  unlist(
    lapply(res,function(x)
    {  
      psm <- salso::psm(x$sim$distr_cluster,nCores = Ncores)
      gtpsm <- salso::psm(rbind(dcgt, dcgt))
      mean((psm-gtpsm)^2)
    }
    )
  )
}



parallel_OCL_extract <- function(res){
  do.call(rbind,
          lapply(res,function(x)
          {  
            obs <- salso::salso(x$sim$obs_cluster,nCores = 10, maxNClusters = 10)
            obs  
          }
          )
  )
}


parallel_DCL_extract <- function(res){
  do.call(rbind,
          lapply(res,function(x)
          {  
            dis <- salso::salso(x$sim$distr_cluster,nCores = Ncores,
                                maxNClusters = max(x$params$group))
            dis
          }
          )
  )
}



est_dens_gibbs <- function(out, 
                           type="CAM",
                           interval = c(-10,10),
                           Lout = 500){
  
  nsim   <- length(out$sim$alpha)
  x_grid <- seq(interval[1],interval[2],length.out = Lout)
  
  if(type == "CAM"){
    out$sim$L = out$sim$maxL
  }
  
  if(type == "fSAN_OF" | type == "fiSAN_OF"){
    out$sim$L = rep(out$params$maxL,nsim)
  }
  
  dens <- loops_post_dens(mu = out$sim$mu,
                          sig = sqrt(out$sim$sigma2),
                          x_grid = x_grid, 
                          est_L = out$sim$L,
                          omega = out$sim$omega,
                          S = out$sim$distr_cluster)  
  
  
  return(list(x_grid = x_grid,
              dens = dens))
  
} 



post_mean_dens_i <- function(i, interval, object, type){
  
  d_xxx <- est_dens_gibbs(out = object[[i]],
                          interval = interval,type = type)
  
  d1 = data.frame(reshape::melt(apply(d_xxx$dens,c(1,2),mean)),  
                  x = d_xxx$x_grid, sim = i,type = type)
  return(d1)
#  return(rbind(d1,d2,d3))
}


parallel_NFD_VI_extract <- function(VI_M_S,dcgt,ocgt){
  FR_OC <-   unlist(
    lapply(VI_M_S,function(x)
    {  
      psm <- salso::psm(x$M_clust,nCores = Ncores)
      gtpsm <- salso::psm(rbind(ocgt, ocgt))
      mean((psm-gtpsm)^2)
    }
    )
  )
  FR_DC <- unlist(
    lapply(VI_M_S,function(x)
    {  
      psm <- salso::psm(x$S_clust,nCores = Ncores)
      gtpsm <- salso::psm(rbind(dcgt, dcgt))
      mean((psm-gtpsm)^2)
    }
    )) 
  return(cbind(FR_OC=FR_OC,FR_DC=FR_DC))
}


variational_sample_M_S <- function(vb_results,nsim,seed){
  
  
  S_clust <- t(replicate(nsim,
                         apply(vb_results$results$sim$RHO,1,function(r)
                           sample(1:vb_results$results$params$K,1,FALSE,r))))
  
  # is this ok?
  XI <- do.call(rbind,vb_results$results$sim$XI)
  
  
  M_clust <-  t(replicate(nsim,
                          apply(XI,1,function(w)
                            sample(1:vb_results$results$params$L,1,FALSE,w))))
  
  return(list(M_clust=M_clust,S_clust=S_clust))
}

ALLvariational_sample_M_S <- function(vb_list,nsim){
  
  IND = length(vb_list)
  parallel::mclapply(1:IND, function(x) variational_sample_M_S(vb_list[[x]],nsim),mc.cores=Ncores)
  
}
