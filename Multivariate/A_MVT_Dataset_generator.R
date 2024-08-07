group_sizes <-  sample_sizes <- c(50, 500, 1000) # group size: we have 6 groups of equal size, hence the total sample size ranges from 300 to 6000
dimensions <- c(2,5,10)
sims = 50

y <- group <- list()

# dataset 1
for(xx in 1:length(sample_sizes)){
  for(dd in 1:length(dimensions)){

  name <- paste0("Multivariate/Datasets/samplesize",sample_sizes[xx],"_dim",dimensions[dd])
  dir.create(name)

  set.seed(1221*xx)
  Y <- c()
  Gr <- c()
  n <- group_sizes[xx]
  D <- dimensions[dd]
  variance <- 0.2 ############# qui varianza da ridurre se troppo vicini
  
  
  Cor1 <- diag(D)
  for(i in 1:D){
    for(j in 1:D){
      if(abs(i-j)<2){
        Cor1[i,j] <- .25
      }
    }
  }
  diag(Cor1) <- 1  
  
  Cor2 <- diag(D)
  Cor2[upper.tri(Cor2)] <- Cor2[lower.tri(Cor2)] <-  +.50
  
  Cor3 <- diag(D)
  Cor3[upper.tri(Cor1)] <- Cor3[lower.tri(Cor3)] <-  +.85
  
  
  for(i in 1:sims){
    
    
    y1 <- rbind( matrix( rnorm((n/2)*D, -5, sqrt(variance)), (n/2),D,byrow = T),
                 mvnfast::rmvn((n/2), rep(-2,D), sigma = Cor1 %*% (diag(D)*variance)))
    
    y2 <- rbind( matrix( rnorm((n/2)*D, -5, sqrt(variance)), (n/2),D,byrow = T),
                 mvnfast::rmvn((n/2), rep(-2,D), sigma = Cor1 %*% (diag(D)*variance)))
    
    y3 <- rbind(matrix(rnorm((n/2)*D, 2, sqrt(variance)), (n/2), D),
                mvnfast::rmvn((n/2), rep(5,D), sigma = Cor2 %*% (diag(D)*variance)))
    
    y4 <- rbind(matrix(rnorm((n/2)*D, 2, sqrt(variance)), (n/2), D),
                mvnfast::rmvn((n/2), rep(5,D), sigma = Cor2 %*% (diag(D)*variance)))
    
    
    y5 <- mvnfast::rmvn(n,rep(0,D), sigma = Cor3 %*% (diag(D)*variance))
    y6 <- mvnfast::rmvn(n,rep(0,D), sigma = Cor3 %*% (diag(D)*variance))
    
    ym     <- rbind(y1,y2,y3,y4,y5,y6)
    yl    <- list(y1,y2,y3,y4,y5,y6)
    gr <- rep(1:6,unlist(lapply(yl, nrow)))
   cat(i) 
   y[[i]] = ym
   group[[i]] = gr
   cat(i)
   
  }

  oc1 = rep(1:2,each = n/2)
  dc1 = rep(1,n)
  oc2 = rep(1:2,each = n/2)
  dc2 = rep(1,n)
  
  oc3 = rep(3:4,each = n/2)
  dc3 = rep(2,n)
  oc4 = rep(3:4,each = n/2)
  dc4 = rep(2,n)
  
  oc5 = rep(5,each = n)
  dc5 = rep(3,n)
  oc6 = rep(5,each = n)
  dc6 = rep(3,n)

  
  oc_groundtruth <- c(oc1,oc2,oc3,oc4,oc5,oc6)
  dc_groundtruth <- c(1,1,2,2,3,3)

  
  saveRDS(y, paste0(name,"/data_list",sims,".RDS"))
  saveRDS(group, paste0(name,"/group_list",sims,".RDS"))
  saveRDS(dc_groundtruth, paste0(name,"/trueDC_list",sims,".RDS"))
  saveRDS(oc_groundtruth, paste0(name,"/trueOC_list",sims,".RDS"))
}
}
