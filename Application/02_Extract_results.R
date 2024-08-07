library(tidyverse)
library(ggrepel)
# source script 00 and 01 (only if you still need to run SANvi)


# Load the best run and estimate the clustering --------------------------------
res <- readRDS("Application/Spot_BESTof1000runs_v2.RDS")
plot(res$sim$Elbo_val,type="l")

data_sub <- readRDS("Application/data_sub.RDS")
y_obser2 <- res$params$y


# Estimate clusters ------------------------------------------------------------
vi_dcl = unlist(apply(res$sim$RHO,1,function(c)which.max(c)))
vi_ocl = unlist(lapply(res$sim$XI,function(x)
  apply(x,1,function(c)which.max(c))))
dcobs = rep(vi_dcl,Nj)
ucl = unique(vi_dcl)



table(vi_ocl)
length(table(vi_ocl))
table(vi_dcl)
length(table(vi_dcl))
pairs(y_obser2,col=vi_ocl,pch=".")

# filter the overall data
data_sub2 <- data_sub[-ind_problem1,]
data_sub2 <- data_sub2[-inds_problem2,]

artist_and_custer = data.frame(auth = unique(data_sub2$FA),dcl = vi_dcl) %>% arrange(dcl)
artist_and_custer
# interesting DCL
# run 1000
rappers   <- 1 
hard_rock <- 14
Classical <- 10
speech    <- 18

ART_LIST <- sapply(sort(ucl),function(x) artist_and_custer[artist_and_custer$dcl==x,1])

doc = lapply(ART_LIST,function(r) paste(c(r),collapse = " - "))
doc2 <- do.call(c,doc)

# list of artist for supplementary material
DOC <- as_tibble(cbind(ind = 1:length(ucl),art = doc2 ))
knitr::kable(DOC,format = "latex",)

data <- data_sub2 %>% mutate(ocl=vi_ocl, dcl = dcobs) %>% 
  dplyr::select(FA,ocl,dcl,name) %>% 
  mutate(dur = y_obser2[,1],ene = y_obser2[,2], spe =  y_obser2[,3], group = y_group2)

# description of data
art1 = "AC/DC" 
art2 = "2Pac"
art3 = "Sergei Rachmaninoff"
art4 = "Dale Carnegie"
COLS=c(1,"steelblue3","forestgreen","orangered4")

g1 <- ggplot()+
  geom_point(data=data,aes(y = ene,
                           x = dur), col="gray",size=.1,alpha=1)+
  theme_bw()+
  geom_point(data=data %>% filter(FA == art1 | FA == art2 | FA == art3 | FA == art4),
             aes(y = ene,
                 x = dur, col = FA),size=.9)+
  theme(text = element_text(size=18)) + xlab("Duration") + ylab("Energy")+
  scale_color_manual("",values = COLS)
g2 <- ggplot()+
  geom_point(data=data,aes(y = spe,
                           x = dur), col="gray",size=.1, pch=".",alpha=1)+
  theme_bw()+
  geom_point(data=data %>% filter(FA == art1 | FA == art2 | FA == art3 | FA == art4),
             aes(y = spe,
                 x = dur, col = FA),size=.9)+
  theme(text = element_text(size=18)) + xlab("Duration") + ylab("Speechiness")+
  scale_color_manual("",values = COLS)
g3 <- ggplot()+
  geom_point(data=data,aes(y = spe,
                           x = ene), col="gray",size=.1, pch=".",alpha=1)+
  theme_bw()+
  geom_point(data=data %>% filter(FA == art1 | FA == art2 | FA == art3 | FA == art4),
             aes(y = spe,
                 x = ene, col = FA),size=.9)+
  theme(text = element_text(size=18)) + xlab("Energy") + ylab("Speechiness")+
  scale_color_manual("",values = COLS)

g1+g2+g3  +
  plot_layout(guides = 'collect') & theme(legend.position = "bottom")

#ggviewggview(h=6,w=18)
#################################################################################
ggsave("Application/Output/descriptive.pdf",h=5,w=15)
ggsave("Application/Output/descriptive.eps",h=5,w=15)
#################################################################################



# overview of DC  -------------------------------------------------------------
dev.off()

boxplot(data$dur~data$dcl)
boxplot(data$ene~data$dcl)
boxplot(data$spe~data$dcl)

dbx = data %>% filter(dcl == hard_rock | dcl == rappers | dcl == Classical | dcl == speech) %>% 
  dplyr::select(ocl,dcl,dur,ene,spe) %>% mutate(d = case_when(dcl == rappers ~ "Rap",
                                                              dcl == Classical~ "Classical",
                                                              dcl == hard_rock~"Hard\nrock",
                                                              dcl == speech~"Audio\nlectures"))
dbxm <- dbx %>% reshape2::melt(c("d","ocl","dcl")) %>% 
  mutate(tit = case_when(variable == "ene" ~"Energy",
                         variable == "dur" ~"Duration",
                         variable == "spe" ~"Speechiness"))
ggplot(dbxm)+
  geom_boxplot(aes(x=value,y=d,group=d))+
  theme(text = element_text(size=18))+ 
  facet_wrap(~tit)+ylab("")+xlab("Score value")+
  scale_y_discrete(limits=rev)

#ggviewggview(h=5,w=15)
#################################################################################
ggsave("Application/Output/Genre_boxplots_1000_v2.eps",h=5,w=15)
ggsave("Application/Output/Genre_boxplots_1000_v2.pdf",h=5,w=15)
#################################################################################

# overview of DC from OC standpoint  -------------------------------------------------------------
D <- tibble(y_du=y_obser2[,1],
            y_en=y_obser2[,2],
            y_sp=y_obser2[,3],
            g=y_group2,
            ocl = vi_ocl,
            dcl = dcobs)

mus <- t(res$sim$ml)
y_obser <- y_obser2
y_group <- y_group2
inds <- c(2,3)

post_weights <- res$sim$beta_bar_lk
WW <-  apply(post_weights,2,function(r) r/sum(r))

plot_atoms <- function(inds,distr,main = NULL,thre=0.05){
  
  plot(y_obser2[,inds],
       pch=".",col="gray",main = main)
  points(y_obser2[data$dcl==distr,inds],
         pch=".",col="blue",cex=2)
  for(i in 1:L){
    if(WW[i,distr]<thre){
      next 
    }
    Sig = solve(res$sim$Wl[,,i] * res$sim$nul[i]) #/ res$sim$betal[i]
    lines(ellipse::ellipse(Sig[inds,inds],
                           center=c(mus[i,inds[1]],mus[i,inds[2]])),
          lwd=WW[i,distr]*5)
    points(mus[i,inds[1]],
           mus[i,inds[2]],
           pch=21,
           cex=WW[i,distr]*2,
           bg="red")
  }
}


par(mfrow=c(2,3))
plot_atoms(inds = c(2,3),hard_rock,main = "Hard Rock")
plot_atoms(c(2,1), hard_rock,main = "Hard Rock")
plot_atoms(c(3,1), hard_rock,main = "Hard Rock")
plot_atoms(c(2,3), Classical,main = "Classical")
plot_atoms(c(2,1), Classical,main = "Classical")
plot_atoms(c(3,1), Classical,main = "Classical")


par(mfrow=c(2,3))
plot_atoms(c(2,3), rappers,main = "Rap")
plot_atoms(c(2,1), rappers,main = "Rap")
plot_atoms(c(3,1), rappers,main = "Rap")
plot_atoms(c(2,3), speech,main = "Audio lectures")
plot_atoms(c(2,1), speech,main = "Audio lectures")
plot_atoms(c(3,1), speech,main = "Audio lectures")

# nice idea, let's make it in ggplot!

DD = reshape2::melt(D,id = c("g", "ocl", "dcl"))

all <-  combn(3,2)
Dall <-  tibble(y_du=data$dur,
              y_en=data$ene,
              y_sp=data$spe,
              g = data$group,
              ocl = data$ocl,
              dcl = data$dcl)
INDS = sample(1:nrow(Dall))
D = Dall[INDS,]


colnames(y_obser2) <- c("Duration","Energy","Speechiness")

GG_spot_major <- function(distr, thresh = .01,
                          col="cyan4",key=1 ){
  H =list()
  for(i in 1:ncol(all)){
    h =  local({
      i <- i
      used = apply(WW,2,function(x) all(round(x-1/L,5)==0))
      act = which(WW[,distr]>thresh)
      DDD = as_tibble(D[,c(all[,i],6)])
      colnames(DDD) = c("X1","X2","dcl")
      MM1= EE1 = c()
      inds  = all[,i]
      
      for(jj in 1:L){
        
        if(res$sim$Wl[inds,inds,jj][1,1]==1){
          next
        }
        Sig = solve(res$sim$Wl[inds,inds,jj] * res$sim$nul[jj]) #/ res$sim$betal[i]
        E1 = ellipse::ellipse(Sig,
                              center=c(mus[jj,inds[1]],mus[jj,inds[2]]))
        
        E1 = cbind(E1,xx = jj,w = WW[jj,distr], active = jj %in% act, w2 = WW[jj,distr]*key)
        EE1 =rbind(EE1,E1)  
        MM1 = rbind(MM1, c(mus[jj,inds],xx = jj,w = WW[jj,distr], 
                           active = jj %in% act, w2 = WW[jj,distr]*key))
      }
      colnames(MM1)[1:2] = c("x","y")
      ggplot()+
        # all points
        geom_point(data=DDD,aes(x = X1,
                                y = X2), col="gray", pch=".",alpha=.2)+
        theme_bw()+
        # all ellipses
        geom_path(data=as_tibble(EE1),aes(x=x,y=y,group = xx),alpha = .1,lty=3)+#alpha=w
        
        geom_point(data=DDD %>% filter(dcl==distr),aes(x = X1,
                                                       y = X2),
                   fill=col,pch=21,col="transparent")+
        
        xlab(colnames(y_obser2)[all[1,i]])+
        ylab(colnames(y_obser2)[all[2,i]])+
        
        geom_path(data=as_tibble(EE1),aes(x=x,y=y,group = xx),alpha = as_tibble(EE1)$active,
                  lwd = as_tibble(EE1)$w2)+#alpha=w
        geom_point(data=as_tibble(MM1),aes(x=x,y=y,
                                           group = xx,alpha = active+.001*w2,
                                           size=w2/2+.25),pch="+")+
        theme(legend.position = "none",text = element_text(size=18))
    })
    H[[i]] = h
  }
  (H[[1]]|H[[2]]|H[[3]])
}


U1 <- wrap_elements((GG_spot_major(rappers,.01,"cyan4",key = 1.2))+plot_annotation(title = 'Rap') & theme(text = element_text(size=18)))
U1

U2 <- wrap_elements((GG_spot_major(hard_rock,.01,"orangered3",key = 1.2))+plot_annotation(title = 'Hard rock') & theme(text = element_text(size=18)))
U2

U3 <- wrap_elements((GG_spot_major(Classical,.01,"orange3",key = 1.2))+plot_annotation(title = 'Classical') & theme(text = element_text(size=18)))
U3

U4 <- wrap_elements((GG_spot_major(speech,.05,"green4",key = 1.2))+plot_annotation(title = 'Audio lectures') & theme(text = element_text(size=18)))
U4

dev.off()
wrap_plots(U2, U1, nrow = 2)
#ggviewggview(h=10,w=15)
################################################################################
ggsave("Application/Output/p1_1000_v2.eps",h=10,w=15,device = cairo_ps)
ggsave("Application/Output/p1_1000_v2.pdf",h=10,w=15,device = cairo_pdf)
ggsave("Application/Output/p1_1000_v2.png",h=10,w=15)
 ################################################################################


wrap_plots(U4, U3, nrow = 2)
#ggviewggview(h=10,w=15)
################################################################################
ggsave("Application/Output/p2_1000_v2.eps",h=10,w=15,device = cairo_ps)
ggsave("Application/Output/p2_1000_v2.pdf",h=10,w=15)
ggsave("Application/Output/p2_1000_v2.png",h=10,w=15)
################################################################################


# OC analysis -------------------------------------------------------------
#View(data)
act <- res$sim$nul>8
atom <- 32

Club32 <- data %>% filter(ocl == atom )
length(unique(Club32$FA))
nrow(Club32)
length(unique(Club32$dcl))


s1 = (data %>% filter(ocl == atom ) %>% filter(FA == "AC/DC", name=="Thunderstruck") )[1,]
s2 = (data %>% filter(ocl == atom ) %>% filter(FA == "Bruce Springsteen", 
                                               name=="Born in the U.S.A.") )[1,]
s3 = (data %>% filter(ocl == atom ) %>% filter(FA == "Michael Jackson", name=="Thriller - 2003 Edit") )[1,] %>% 
  mutate(name = "Thriller")
s4 = (data %>% filter(ocl == atom ) %>% filter(FA == "Madonna", name=="Like a Prayer") )[1,]

D1 = rbind(s1,s2,s3,s4)

temp = data %>% filter( ocl == atom )
table(temp$ocl,temp$dcl)  
table(temp$ocl)  
length(table(temp$FA)  )

MM1= EE1 = c()
inds  = 1:2
for(jj in atom){
  
  if(res$sim$Wl[inds,inds,jj][1,1]==1){
    next
  }
  Sig = solve(res$sim$Wl[inds,inds,jj] * res$sim$nul[jj]) #/ res$sim$betal[i]
  E1 = ellipse::ellipse(Sig,
                        center=c(mus[jj,inds[1]],mus[jj,inds[2]]))
  
  E1 = cbind(E1,xx = jj,#w = WW[jj,distr], 
             active = 1)#jj %in% act)
  EE1 =rbind(EE1,E1)  
  MM1 = rbind(MM1, c(mus[jj,inds],xx = jj,#w = WW[jj,distr],
                     active = jj %in% act))
}


g1 <- ggplot()+
  geom_point(data=data,aes(y = ene,
                           x = dur), col="gray",
             size=.1,alpha=1)+
  theme_bw()+
  geom_point(data=data %>% filter( ocl ==atom ),
             aes(y = ene,
                 x = dur, col = factor(ocl)),alpha=.8,
             size=.5)+
  theme(text = element_text(size=18), legend.position = "none") + xlab("Duration") + ylab("Energy")+
  #scale_color_viridis_d()+
  scale_color_manual("",values = "royalblue3")+
  geom_path(data=as_tibble(EE1),aes(x=x,y=y,group = xx),
            alpha = as_tibble(EE1)$active )+#alpha=w
  geom_point(data=as_tibble(MM1),aes(x=V1,y=V2,group = xx,
                                     alpha = active+.001,size=active+.001),
             pch="+")

g11= g1 +
  geom_label_repel(data=D1,aes(x=dur,y=ene,label=name),seed=42, min.segment.length = 0, box.padding = 5)+
  geom_point(data=D1,aes(x=dur,y=ene),col="red2")
g11
#-----
MM1= EE1 = c()
inds  = c(1,3)
for(jj in c(atom)){
  
  Sig = solve(res$sim$Wl[inds,inds,jj] * res$sim$nul[jj]) #/ res$sim$betal[i]
  E1 = ellipse::ellipse(Sig,
                        center=c(mus[jj,inds[1]],mus[jj,inds[2]]))
  
  E1 = cbind(E1,xx = jj,#w = WW[jj,distr], 
             active = 1)#jj %in% act)
  EE1 =rbind(EE1,E1)  
  MM1 = rbind(MM1, c(mus[jj,inds],xx = jj,#w = WW[jj,distr], 
                     active = jj %in% act))
}
g2 <- ggplot()+
  geom_point(data=data,aes(y = spe,
                           x = dur), col="gray",
             size=.1,alpha=1)+
  theme_bw()+
  geom_point(data=data %>% filter( ocl ==atom ),
             aes(y = spe,
                 x = dur, col = factor(ocl)),alpha=.8,
             size=.5)+
  theme(text = element_text(size=18), legend.position = "none") + 
  xlab("Duration") + ylab("Speechiness")+
  #scale_color_viridis_d()+
  scale_color_manual("",values = "royalblue3")+
  geom_path(data=as_tibble(EE1),aes(x=x,y=y,group = xx),alpha = as_tibble(EE1)$active )+#alpha=w
  geom_point(data=as_tibble(MM1),aes(x=V1,y=V2,group = xx,
                                     alpha = active+.001,size=active+.001),
             pch="+")

g22=g2 +
  geom_label_repel(data=D1,aes(x=dur,y=spe,label=name),seed=42, min.segment.length = 0, box.padding = 5)+
  geom_point(data=D1,aes(x=dur,y=spe),col="red2")

#-----
MM1= EE1 = c()
inds  = c(2,3)
for(jj in c(atom)){
  
  Sig = solve(res$sim$Wl[inds,inds,jj] * res$sim$nul[jj]) #/ res$sim$betal[i]
  E1 = ellipse::ellipse(Sig,
                        center=c(mus[jj,inds[1]],mus[jj,inds[2]]))
  
  E1 = cbind(E1,xx = jj,#w = WW[jj,distr], 
             active =1)# jj %in% act)
  EE1 =rbind(EE1,E1)  
  MM1 = rbind(MM1, c(mus[jj,inds],xx = jj,#w = WW[jj,distr], 
                     active = jj %in% act))
}
g3 <- ggplot()+
  geom_point(data=data,aes(y = spe,
                           x = ene), col="gray",
             size=.1,alpha=1)+
  theme_bw()+
  geom_point(data=data %>% filter( ocl ==atom ),
             aes(y = spe,
                 x = ene, col = factor(ocl)),alpha=.8,
             size=.5)+
  theme(text = element_text(size=18), legend.position = "none") + 
  xlab("Energy") + ylab("Speechiness")+
  #scale_color_viridis_d()+
  scale_color_manual("",values = "royalblue3")+
  geom_path(data=as_tibble(EE1),aes(x=x,y=y,group = xx),alpha = as_tibble(EE1)$active )+#alpha=w
  geom_point(data=as_tibble(MM1),aes(x=V1,y=V2,group = xx,
                                     alpha = active+.001,size=active+.001),
             pch="+")

g33 = g3 +
  geom_label_repel(data=D1,aes(x=ene,y=spe,label=name),seed=42, 
                   min.segment.length = 0, box.padding = 3)+
  geom_point(data=D1,aes(x=ene,y=spe),col="red2")


g11+g22+g33  

#ggviewggview(h=5,w=15)
################################################################################
ggsave("Application/Output/OC_analisis_1000_v2.eps",h=5,w=15, device = cairo_ps)
ggsave("Application/Output/OC_analisis_1000_v2.pdf",h=5,w=15)
ggsave("Application/Output/OC_analisis_1000_v2.png",h=5,w=15)


