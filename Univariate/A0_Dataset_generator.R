sample_sizes = c(10,50,500,2500)#,5000)
replica_group = 2
s  = .6
mu = c(-5,-2,2,5,0)
sims = 50

# dataset 1 (simil-jasa)
for(xx in 1:length(sample_sizes)){
  
  name <- paste0("Simulation1/Datasets/scen1_samplesize",sample_sizes[xx])
  dir.create(name)
  
  y = list()
  group = list()
  group_ss  = sample_sizes[xx]
  
  
  set.seed(123454321*xx)
  
  
  for(i in 1:sims){
    Y = c()
    Gr = c()
    
    for(j in 1:replica_group){
      Y = c(Y, rnorm(group_ss/2,mu[1],s), rnorm(group_ss/2,mu[2],s),
            rnorm(group_ss/2,mu[3],s), rnorm(group_ss/2,mu[4],s),
            rnorm(group_ss/1,mu[5],s))
      
      Gr = c(Gr,
             c(rep(1. +(3*(j-1)) ,  group_ss),
               rep(2. +(3*(j-1)),   group_ss),
               rep(3. +(3*(j-1)),   group_ss) ))
    }
    
    y[[i]] = Y
    group[[i]] = Gr
  }
  
  oc_groundtruth <- c(rep(1,group_ss/2), rep(2,group_ss/2),
                      rep(3,group_ss/2), rep(4,group_ss/2),
                      rep(5,group_ss))
  
  oc_groundtruth <- rep(oc_groundtruth,replica_group)
  dc_groundtruth <- rep(c(1,2,3),replica_group)
  
  plot(y[[2]],col=group[[2]])
  plot(y[[2]],col=oc_groundtruth)
  
  saveRDS(y, paste0(name,"/data_list",sims,".RDS"))
  saveRDS(group, paste0(name,"/group_list",sims,".RDS"))
  saveRDS(dc_groundtruth, paste0(name,"/trueDC_list",sims,".RDS"))
  saveRDS(oc_groundtruth, paste0(name,"/trueOC_list",sims,".RDS"))
}

f1 = function(x) .5*dnorm(x,mu[1],s)+.5*dnorm(x,mu[2],s)
f2 = function(x) .5*dnorm(x,mu[3],s)+.5*dnorm(x,mu[4],s)
f3 = function(x) .5*dnorm(x,mu[5],s)+.5*dnorm(x,mu[5],s)
