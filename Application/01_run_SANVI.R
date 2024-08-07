
# Run the algorithm 1000 times and extract the best run ------------------------
# (takes a while! + remember to change the number of cores you want to use)
NSIM <- 1000
VIall = parallel::mclapply(1:NSIM, function(i)   variational_mvt_fiSAN(y_obser2, y_group2,
                                                                      seed = 250*i,
                                                                      beta_bar = beta_dir,
                                                                      maxK = K,
                                                                      maxL = L,
                                                                      epsilon = 1e-4,
                                                                      maxSIM = 500),
                           mc.cores = 25)
saveRDS(VIall,"Application/Runs/Spot_1000runs_v2.RDS")
elbos = lapply(VIall, function(x) x$sim$Elbo)
plot(elbos[[1]][-c(1:3)],type="l",ylim = c(-16000,-10000),xlim = c(1,500))
for(i in 2:NSIM){
  lines(elbos[[i]],type="l")
}
ind.max <- which.max(unlist(lapply(elbos,max)))
res = VIall[[ind.max]]
lines(elbos[[ind.max]],type="b",col=4)
max(elbos[[ind.max]])
ind.max
lines(elbos[[906]],type="b",col=2)
saveRDS(res,"Application/Spot_BESTof1000runs_v2.RDS")



plot(diff(elbos[[1]][-c(1:3)]),type="l",ylim = c(-1,10),xlim = c(1,500))
for(i in 2:NSIM){
  lines(diff(elbos[[i]]),type="l")
}
abline(h=0)


plot(diff(elbos[[1]][-c(1:3)]),type="l",ylim = c(-.001,.001),xlim = c(1,500))
for(i in 2:NSIM){
  lines(diff(elbos[[i]]),type="l")
}
abline(h=0)

