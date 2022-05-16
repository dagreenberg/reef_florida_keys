#Output analysis
library(here); library(dplyr);library(rlist)

###Functions###
source(here('code','functions.R'))

####Species' model summaries####
files_kl<- list.files(here('outputs','species model summary','Key Largo'))
key_largo<- cbind(read.csv(here('outputs','species model summary','Key Largo',files_kl[1])))
for(i in 2:length(files_kl)){
  key_largo<- rbind(key_largo,read.csv(here('outputs','species model summary','Key Largo',files_kl[i])))
}

files_kw<- list.files(here('outputs','species model summary','Key West'))
key_west<- read.csv(here('outputs','species model summary','Key West',files_kw[1]))
for(i in 2:length(files_kw)){
  key_west<- rbind(key_west,read.csv(here('outputs','species model summary','Key West',files_kw[i])))
}

####Fish trait data####
fish_traits<- read.csv(here('data','REEF','Caribbean_fish_trait_matrix.csv'))

#Match together datasets
m_kl<- match(key_largo$SP,fish_traits$commonname)
key_largo$size<- fish_traits$size[m_kl] #species length (cm)
key_largo[,13:14]<- fish_traits[m_kl,8:9] #colour & behaviour
key_largo[,15:16]<- fish_traits[m_kl,11:12] #Family, schooling

m_kw<- match(key_west$SP,fish_traits$commonname)
key_west$size<- fish_traits$size[m_kw] #species length (cm)
key_west[,13:14]<- fish_traits[m_kw,8:9] #colour & behaviour
key_west[,15:16]<- fish_traits[m_kw,11:12] #Family & schooling 

####1 - Mean abundance for each survey####
files_kl_par_2<- list.files(here('outputs','species parameter estimates','Key Largo','model 2'))
####RVC vs. REEF abundance correlations####
mean_abund_list_kl<- list()
for(i in 1:nrow(key_largo)){
  sp<- read.csv(here('outputs','species parameter estimates','Key Largo','model 2',files_kl_par_2[i]))
  x1=sp[,(gsub('\\..*','',colnames(sp))=='x1')]
  x2=sp[,(gsub('\\..*','',colnames(sp))=='x2')]
  c=sp[,2:5][(gsub('\\..*','',names(sp[,2:5]))=='cut')]
  m_abund_rvc<- exp(apply(x1,1,mean))
  m_abund_reef<- exp(apply(log(mean_ord_to_n(x=x2,c=c)),1,mean))
  mean_abund_list_kl[[i]]<- data.frame(mean.rvc.kl=sort(m_abund_rvc),mean.reef.kl=sort(m_abund_reef),sp=rep(key_largo$SP[i],3000))
 key_largo[i,17]=median(m_abund_rvc)
 key_largo[i,18]=quantile(m_abund_rvc,0.025)
 key_largo[i,19]=quantile(m_abund_rvc,0.975)
 key_largo[i,20]=median(m_abund_reef)
 key_largo[i,21]=quantile(m_abund_reef,0.025)
 key_largo[i,22]=quantile(m_abund_reef,0.975)
  }
colnames(key_largo)[17:22]=c('mean.x.reef','mean.x.reef.l95','mean.x.reef.u95','mean.x.rvc','mean.x.rvc.l95',
                             'mean.x.rvc.u95')

cor_kl_list<- list()
cor_kl<- NA
for(q in 1:nrow(mean_abund_list_kl[[1]])){
  cor_kl_list[[q]]<- log10(mean_abund_list_kl[[1]][q,1:2])
  for(z in 2:length(mean_abund_list_kl)){
    cor_kl_list[[q]]<- rbind(cor_kl_list[[q]],log10(mean_abund_list_kl[[z]][q,1:2]))
  }   
  cor_kl[q]<- cor.test(cor_kl_list[[q]][,1],cor_kl_list[[q]][,2])$estimate
}
median(cor_kl)
quantile(cor_kl,0.025)
quantile(cor_kl,0.975)

files_kw_par_2<- list.files(here('outputs','species parameter estimates','Key West','model 2'))
mean_abund_list_kw<- list()
for(i in 1:nrow(key_west)){
  sp<- read.csv(here('outputs','species parameter estimates','Key West','model 2',files_kw_par_2[i]))
  x1=sp[,(gsub('\\..*','',colnames(sp))=='x1')]
  x2=sp[,(gsub('\\..*','',colnames(sp))=='x2')]
  c=sp[,2:5][(gsub('\\..*','',names(sp[,2:5]))=='cut')]
  m_abund_rvc<- exp(apply(x1,1,mean))
  m_abund_reef<- exp(apply(log(mean_ord_to_n(x=x2,c=c)),1,mean))
  mean_abund_list_kw[[i]]<- data.frame(mean.rvc.kw=sort(m_abund_rvc),mean.reef.kw=sort(m_abund_reef),sp=rep(key_west$SP[i],3000))
  key_west[i,17]=median(m_abund_rvc)
  key_west[i,18]=quantile(m_abund_rvc,0.025)
  key_west[i,19]=quantile(m_abund_rvc,0.975)
  key_west[i,20]=median(m_abund_reef)
  key_west[i,21]=quantile(m_abund_reef,0.025)
  key_west[i,22]=quantile(m_abund_reef,0.975)
}
colnames(key_west)[17:22]=c('mean.x.reef','mean.x.reef.l95','mean.x.reef.u95','mean.x.rvc','mean.x.rvc.l95',
                             'mean.x.rvc.u95')

cor_kw_list<- list()
cor_kw<- NA
for(q in 1:nrow(mean_abund_list_kw[[1]])){
  cor_kw_list[[q]]<- log10(mean_abund_list_kw[[1]][q,1:2])
  for(z in 2:length(mean_abund_list_kw)){
    cor_kw_list[[q]]<- rbind(cor_kw_list[[q]],log10(mean_abund_list_kw[[z]][q,1:2]))
  }   
  cor_kw[q]<- cor.test(cor_kw_list[[q]][,1],cor_kw_list[[q]][,2])$estimate
}
median(cor_kw)
quantile(cor_kw,0.025)
quantile(cor_kw,0.975)

###Regional abundance correlations####
matched_species<-key_west$SP[key_west$SP %in% key_largo$SP]

mean_abund_combined<- list()
for(i in 1:length(matched_species)){
  mean_abund_combined[[i]]<-  data.frame(mean_abund_list_kl[[match(matched_species,key_largo$SP)[i]]],mean_abund_list_kw[[match(matched_species,key_west$SP)[i]]]) 
}

#Correlation in species' mean abundance between Key Largo & Key West in RVC surveys
cor_rvc_list<- list()
cor_rvc<- NA
for(q in 1:nrow(mean_abund_combined[[1]])){
  cor_rvc_list[[q]]<- c(log10(mean_abund_combined[[1]][q,1]),log10(mean_abund_combined[[1]][q,4]))
  for(z in 2:length(mean_abund_combined)){
    cor_rvc_list[[q]]<- rbind(cor_rvc_list[[q]],c(log10(mean_abund_combined[[z]][q,1]),log10(mean_abund_combined[[z]][q,4])))
  }   
  cor_rvc[q]<- cor.test(cor_rvc_list[[q]][,1],cor_rvc_list[[q]][,2])$estimate
}
median(cor_rvc)
quantile(cor_rvc,0.025)
quantile(cor_rvc,0.975)

#Correlation in species' mean abundance between Key Largo & Key West in REEF surveys
cor_reef_list<- list()
cor_reef<- NA
for(q in 1:nrow(mean_abund_combined[[1]])){
  cor_reef_list[[q]]<- c(log10(mean_abund_combined[[1]][q,2]),log10(mean_abund_combined[[1]][q,5]))
  for(z in 2:length(mean_abund_combined)){
    cor_reef_list[[q]]<- rbind(cor_reef_list[[q]],c(log10(mean_abund_combined[[z]][q,2]),log10(mean_abund_combined[[z]][q,5])))
  } 
  cor_reef[q]<- cor.test(cor_reef_list[[q]][,1],cor_reef_list[[q]][,2])$estimate
}
median(cor_reef)
quantile(cor_reef,0.025)
quantile(cor_reef,0.975)

key_largo_sub<- subset(key_largo, SP %in% matched_species)
key_west_sub<- subset(key_west, SP %in% matched_species)
key_west_sub<- key_west_sub[match(key_largo_sub$SP,key_west_sub$SP),]

#Figure 1 ####
par(mfrow=c(2,2))
plot(log10(mean.x.reef)~log10(mean.x.rvc),data=key_largo,bty='l',xaxt='n',yaxt='n',xlab='Mean count (RVC)',ylab='Mean count (REEF)',type='n',ylim=c(-3,1.5),xlim=c(-3,1.5),main='Key Largo')
points(log10(mean.x.reef)~log10(mean.x.rvc),data=subset(key_largo,mod=='model1'),pch=21,bg=adjustcolor('turquoise4',alpha.f=0.7),col='transparent',cex=1.5)
for(i in 1:nrow(subset(key_largo,mod=='model1'))){
  lines(c(log10(subset(key_largo,mod=='model1')$mean.x.reef.l95[i]),log10(subset(key_largo,mod=='model1')$mean.x.reef.u95[i]))~rep(c(log10(subset(key_largo,mod=='model1')$mean.x.rvc[i])),2),lwd=1,col=adjustcolor('turquoise4',alpha.f=0.4))
}
for(i in 1:nrow(subset(key_largo,mod=='model1'))){
  lines(rep(c(log10(subset(key_largo,mod=='model1')$mean.x.reef[i])),2)~c(log10(subset(key_largo,mod=='model1')$mean.x.rvc.l95[i]),log10(subset(key_largo,mod=='model1')$mean.x.rvc.u95[i])),lwd=1,col=adjustcolor('turquoise4',alpha.f=0.4))
}
points(log10(mean.x.reef)~log10(mean.x.rvc),data=subset(key_largo,mod=='model2'),pch=21,bg=adjustcolor('goldenrod',alpha.f=0.7),col='transparent',cex=1.5)
for(i in 1:nrow(subset(key_largo,mod=='model2'))){
  lines(c(log10(subset(key_largo,mod=='model2')$mean.x.reef.l95[i]),log10(subset(key_largo,mod=='model2')$mean.x.reef.u95[i]))~rep(c(log10(subset(key_largo,mod=='model2')$mean.x.rvc[i])),2),lwd=1,col=adjustcolor('goldenrod',alpha.f=0.4))
}
for(i in 1:nrow(subset(key_largo,mod=='model2'))){
  lines(rep(c(log10(subset(key_largo,mod=='model2')$mean.x.reef[i])),2)~c(log10(subset(key_largo,mod=='model2')$mean.x.rvc.l95[i]),log10(subset(key_largo,mod=='model2')$mean.x.rvc.u95[i])),lwd=1,col=adjustcolor('goldenrod',alpha.f=0.4))
}
axis(2, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(2, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)

axis(1, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
median(cor_kl)
quantile(cor_kl,0.025)
quantile(cor_kl,0.975)
text(x=log10(5),y=log10(0.001),expression(paste(italic(r),' = 0.76 [0.75 - 0.78]')))
text(-3,par("usr")[2] + 0.05,'A',cex=1,pos=3,xpd=T,font=2)

plot(log10(key_west_sub$mean.x.rvc)~log10(key_largo_sub$mean.x.rvc),data=key_largo_sub,bty='l',xaxt='n',yaxt='n',xlab='Mean count (Key Largo)',ylab='Mean count (Key West)',type='n',ylim=c(-3,1.5),xlim=c(-3,1.5),main='RVC')
points(log10(key_west_sub$mean.x.rvc)~log10(key_largo_sub$mean.x.rvc),pch=21,bg=adjustcolor('darkblue',alpha.f=0.7),col='transparent',cex=1.5)
for(i in 1:nrow(key_west_sub)){
  lines(c(log10(key_west_sub$mean.x.rvc.l95[i]),log10(key_west_sub$mean.x.rvc.u95[i]))~rep(c(log10(key_largo_sub$mean.x.rvc[i])),2),lwd=1,col=adjustcolor('darkblue',alpha.f=0.4))
}
for(i in 1:nrow(key_largo_sub)){
  lines(rep(c(log10(key_west_sub$mean.x.rvc[i])),2)~c(log10(key_largo_sub$mean.x.rvc.l95[i]),log10(key_largo_sub$mean.x.rvc.u95[i])),lwd=1,col=adjustcolor('darkblue',alpha.f=0.4))
}
axis(2, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(2, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)

axis(1, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
median(cor_rvc)
quantile(cor_rvc,0.025)
quantile(cor_rvc,0.975)
text(x=log10(5),y=log10(0.001),expression(paste(italic(r),' = 0.94 [0.92 - 0.94]')))
text(-3,par("usr")[2] + 0.05,'B',cex=1,pos=3,xpd=T,font=2)

plot(log10(key_west_sub$mean.x.reef)~log10(key_largo_sub$mean.x.reef),data=key_largo_sub,bty='l',xaxt='n',yaxt='n',xlab='Mean count (Key Largo)',ylab='Mean count (Key West)',type='n',ylim=c(-3,1.5),xlim=c(-3,1.5),main='REEF')
points(log10(key_west_sub$mean.x.reef)~log10(key_largo_sub$mean.x.reef),pch=21,bg=adjustcolor('darkred',alpha.f=0.7),col='transparent',cex=1.5)
for(i in 1:nrow(key_west_sub)){
  lines(c(log10(key_west_sub$mean.x.reef.l95[i]),log10(key_west_sub$mean.x.reef.u95[i]))~rep(c(log10(key_largo_sub$mean.x.reef[i])),2),lwd=1,col=adjustcolor('darkred',alpha.f=0.4))
}
for(i in 1:nrow(key_largo_sub)){
  lines(rep(c(log10(key_west_sub$mean.x.reef[i])),2)~c(log10(key_largo_sub$mean.x.reef.l95[i]),log10(key_largo_sub$mean.x.reef.u95[i])),lwd=1,col=adjustcolor('darkred',alpha.f=0.4))
}
axis(2, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(2, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)

axis(1, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
median(cor_reef)
quantile(cor_reef,0.025)
quantile(cor_reef,0.975)
text(x=log10(5),y=log10(0.001),expression(paste(italic(r),' = 0.90 [0.89 - 0.91]')))
text(-3,par("usr")[2] + 0.05,'C',cex=1,pos=3,xpd=T,font=2)

plot(log10(mean.x.reef)~log10(mean.x.rvc),data=key_west,bty='l',xaxt='n',yaxt='n',xlab='Mean count (RVC)',ylab='Mean count (REEF)',type='n',ylim=c(-3,1.5),xlim=c(-3,1.5),main='Key West')
points(log10(mean.x.reef)~log10(mean.x.rvc),data=subset(key_west,mod=='model1'),pch=21,bg=adjustcolor('turquoise4',alpha.f=0.7),col='transparent',cex=1.5)
for(i in 1:nrow(subset(key_west,mod=='model1'))){
  lines(c(log10(subset(key_west,mod=='model1')$mean.x.reef.l95[i]),log10(subset(key_west,mod=='model1')$mean.x.reef.u95[i]))~rep(c(log10(subset(key_west,mod=='model1')$mean.x.rvc[i])),2),lwd=1,col=adjustcolor('turquoise4',alpha.f=0.4))
}
for(i in 1:nrow(subset(key_west,mod=='model1'))){
  lines(rep(c(log10(subset(key_west,mod=='model1')$mean.x.reef[i])),2)~c(log10(subset(key_west,mod=='model1')$mean.x.rvc.l95[i]),log10(subset(key_west,mod=='model1')$mean.x.rvc.u95[i])),lwd=1,col=adjustcolor('turquoise4',alpha.f=0.4))
}

points(log10(mean.x.reef)~log10(mean.x.rvc),data=subset(key_west,mod=='model2'),pch=21,bg=adjustcolor('goldenrod',alpha.f=0.7),col='transparent',cex=1.5)
for(i in 1:nrow(subset(key_west,mod=='model2'))){
  lines(c(log10(subset(key_west,mod=='model2')$mean.x.reef.l95[i]),log10(subset(key_west,mod=='model2')$mean.x.reef.u95[i]))~rep(c(log10(subset(key_west,mod=='model2')$mean.x.rvc[i])),2),lwd=1,col=adjustcolor('goldenrod',alpha.f=0.4))
}
for(i in 1:nrow(subset(key_west,mod=='model2'))){
  lines(rep(c(log10(subset(key_west,mod=='model2')$mean.x.reef[i])),2)~c(log10(subset(key_west,mod=='model2')$mean.x.rvc.l95[i]),log10(subset(key_west,mod=='model2')$mean.x.rvc.u95[i])),lwd=1,col=adjustcolor('goldenrod',alpha.f=0.4))
}
axis(2, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(2, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
axis(1, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
median(cor_kw)
quantile(cor_kw,0.025)
quantile(cor_kw,0.975)
text(x=log10(5),y=log10(0.001),expression(paste(italic(r),' = 0.81 [0.80 - 0.86]')))
text(-3,par("usr")[2] + 0.05,'D',cex=1,pos=3,xpd=T,font=2)

####Figure 1. w/ species labels:####
par(mfrow=c(2,2))
plot(log10(mean.x.reef)~log10(mean.x.rvc),data=key_largo,bty='l',xaxt='n',yaxt='n',xlab='Mean count (RVC)',ylab='Mean count (REEF)',type='n',ylim=c(-3,1.5),xlim=c(-3,1.5),main='Key Largo')
points(log10(mean.x.reef)~log10(mean.x.rvc),data=subset(key_largo,mod=='model1'),pch=21,bg=adjustcolor('turquoise4',alpha.f=0.7),col='transparent',cex=1.5)
for(i in 1:nrow(subset(key_largo,mod=='model1'))){
  lines(c(log10(subset(key_largo,mod=='model1')$mean.x.reef.l95[i]),log10(subset(key_largo,mod=='model1')$mean.x.reef.u95[i]))~rep(c(log10(subset(key_largo,mod=='model1')$mean.x.rvc[i])),2),lwd=1,col=adjustcolor('turquoise4',alpha.f=0.4))
}
for(i in 1:nrow(subset(key_largo,mod=='model1'))){
  lines(rep(c(log10(subset(key_largo,mod=='model1')$mean.x.reef[i])),2)~c(log10(subset(key_largo,mod=='model1')$mean.x.rvc.l95[i]),log10(subset(key_largo,mod=='model1')$mean.x.rvc.u95[i])),lwd=1,col=adjustcolor('turquoise4',alpha.f=0.4))
}
points(log10(mean.x.reef)~log10(mean.x.rvc),data=subset(key_largo,mod=='model2'),pch=21,bg=adjustcolor('goldenrod',alpha.f=0.7),col='transparent',cex=1.5)
for(i in 1:nrow(subset(key_largo,mod=='model2'))){
  lines(c(log10(subset(key_largo,mod=='model2')$mean.x.reef.l95[i]),log10(subset(key_largo,mod=='model2')$mean.x.reef.u95[i]))~rep(c(log10(subset(key_largo,mod=='model2')$mean.x.rvc[i])),2),lwd=1,col=adjustcolor('goldenrod',alpha.f=0.4))
}
for(i in 1:nrow(subset(key_largo,mod=='model2'))){
  lines(rep(c(log10(subset(key_largo,mod=='model2')$mean.x.reef[i])),2)~c(log10(subset(key_largo,mod=='model2')$mean.x.rvc.l95[i]),log10(subset(key_largo,mod=='model2')$mean.x.rvc.u95[i])),lwd=1,col=adjustcolor('goldenrod',alpha.f=0.4))
}
axis(2, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(2, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)

axis(1, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
median(cor_kl)
quantile(cor_kl,0.025)
quantile(cor_kl,0.975)
text(x=jitter(log10(key_largo$mean.x.rvc),0.1),y=jitter(log10(key_largo$mean.x.reef),0.1),key_west$SP,cex=0.5)
text(x=log10(5),y=log10(0.001),expression(paste(italic(r),' = 0.82 (95% CI: 0.78-0.84)')))
text(-3,par("usr")[2] + 0.05,'A',cex=1,pos=3,xpd=T,font=2)

plot(log10(key_west_sub$mean.x.rvc)~log10(key_largo_sub$mean.x.rvc),data=key_largo_sub,bty='l',xaxt='n',yaxt='n',xlab='Mean count (Key Largo)',ylab='Mean count (Key West)',type='n',ylim=c(-3,1.5),xlim=c(-3,1.5),main='RVC')
points(log10(key_west_sub$mean.x.rvc)~log10(key_largo_sub$mean.x.rvc),pch=21,bg=adjustcolor('darkblue',alpha.f=0.7),col='transparent',cex=1.5)
for(i in 1:nrow(key_west_sub)){
  lines(c(log10(key_west_sub$mean.x.rvc.l95[i]),log10(key_west_sub$mean.x.rvc.u95[i]))~rep(c(log10(key_largo_sub$mean.x.rvc[i])),2),lwd=1,col=adjustcolor('darkblue',alpha.f=0.4))
}
for(i in 1:nrow(key_largo_sub)){
  lines(rep(c(log10(key_west_sub$mean.x.rvc[i])),2)~c(log10(key_largo_sub$mean.x.rvc.l95[i]),log10(key_largo_sub$mean.x.rvc.u95[i])),lwd=1,col=adjustcolor('darkblue',alpha.f=0.4))
}
axis(2, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(2, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)

axis(1, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
median(cor_rvc)
quantile(cor_rvc,0.025)
quantile(cor_rvc,0.975)
text(x=jitter(log10(key_largo_sub$mean.x.rvc),0.1),y=jitter(log10(key_west_sub$mean.x.rvc),0.1),key_west$SP,cex=0.5)
text(x=log10(5),y=log10(0.001),expression(paste(italic(r),' = 0.93 (95% CI: 0.88-0.95)')))
text(-3,par("usr")[2] + 0.05,'B',cex=1,pos=3,xpd=T,font=2)

plot(log10(key_west_sub$mean.x.reef)~log10(key_largo_sub$mean.x.reef),data=key_largo_sub,bty='l',xaxt='n',yaxt='n',xlab='Mean count (Key Largo)',ylab='Mean count (Key West)',type='n',ylim=c(-3,1.5),xlim=c(-3,1.5),main='REEF')
points(log10(key_west_sub$mean.x.reef)~log10(key_largo_sub$mean.x.reef),pch=21,bg=adjustcolor('darkred',alpha.f=0.7),col='transparent',cex=1.5)
for(i in 1:nrow(key_west_sub)){
  lines(c(log10(key_west_sub$mean.x.reef.l95[i]),log10(key_west_sub$mean.x.reef.u95[i]))~rep(c(log10(key_largo_sub$mean.x.reef[i])),2),lwd=1,col=adjustcolor('darkred',alpha.f=0.4))
}
for(i in 1:nrow(key_largo_sub)){
  lines(rep(c(log10(key_west_sub$mean.x.reef[i])),2)~c(log10(key_largo_sub$mean.x.reef.l95[i]),log10(key_largo_sub$mean.x.reef.u95[i])),lwd=1,col=adjustcolor('darkred',alpha.f=0.4))
}
axis(2, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(2, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)

axis(1, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
median(cor_reef)
quantile(cor_reef,0.025)
quantile(cor_reef,0.975)
text(y=jitter(log10(key_west_sub$mean.x.reef),0.1),x=jitter(log10(key_largo_sub$mean.x.reef),0.1),key_west_sub$SP,cex=0.5)
text(x=log10(5),y=log10(0.001),expression(paste(italic(r),' = 0.86 (95% CI: 0.82-0.88)')))
text(-3,par("usr")[2] + 0.05,'C',cex=1,pos=3,xpd=T,font=2)

plot(log10(mean.x.reef)~log10(mean.x.rvc),data=key_west,bty='l',xaxt='n',yaxt='n',xlab='Mean count (RVC)',ylab='Mean count (REEF)',type='n',ylim=c(-3,1.5),xlim=c(-3,1.5),main='Key West')
points(log10(mean.x.reef)~log10(mean.x.rvc),data=subset(key_west,mod=='model1'),pch=21,bg=adjustcolor('turquoise4',alpha.f=0.7),col='transparent',cex=1.5)
for(i in 1:nrow(subset(key_west,mod=='model1'))){
  lines(c(log10(subset(key_west,mod=='model1')$mean.x.reef.l95[i]),log10(subset(key_west,mod=='model1')$mean.x.reef.u95[i]))~rep(c(log10(subset(key_west,mod=='model1')$mean.x.rvc[i])),2),lwd=1,col=adjustcolor('turquoise4',alpha.f=0.4))
}
for(i in 1:nrow(subset(key_west,mod=='model1'))){
  lines(rep(c(log10(subset(key_west,mod=='model1')$mean.x.reef[i])),2)~c(log10(subset(key_west,mod=='model1')$mean.x.rvc.l95[i]),log10(subset(key_west,mod=='model1')$mean.x.rvc.u95[i])),lwd=1,col=adjustcolor('turquoise4',alpha.f=0.4))
}

points(log10(mean.x.reef)~log10(mean.x.rvc),data=subset(key_west,mod=='model2'),pch=21,bg=adjustcolor('goldenrod',alpha.f=0.7),col='transparent',cex=1.5)
for(i in 1:nrow(subset(key_west,mod=='model2'))){
  lines(c(log10(subset(key_west,mod=='model2')$mean.x.reef.l95[i]),log10(subset(key_west,mod=='model2')$mean.x.reef.u95[i]))~rep(c(log10(subset(key_west,mod=='model2')$mean.x.rvc[i])),2),lwd=1,col=adjustcolor('goldenrod',alpha.f=0.4))
}
for(i in 1:nrow(subset(key_west,mod=='model2'))){
  lines(rep(c(log10(subset(key_west,mod=='model2')$mean.x.reef[i])),2)~c(log10(subset(key_west,mod=='model2')$mean.x.rvc.l95[i]),log10(subset(key_west,mod=='model2')$mean.x.rvc.u95[i])),lwd=1,col=adjustcolor('goldenrod',alpha.f=0.4))
}
axis(2, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(2, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
axis(1, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
median(cor_kw)
quantile(cor_kw,0.025)
quantile(cor_kw,0.975)
text(x=log10(5),y=log10(0.001),expression(paste(italic(r),' = 0.73 (95% CI: 0.65-0.79)')))
text(-3,par("usr")[2] + 0.05,'D',cex=1,pos=3,xpd=T,font=2)
text(x=jitter(log10(key_west$mean.x.rvc),0.1),y=jitter(log10(key_west$mean.x.reef),0.1),key_west$SP,cex=0.5)

####2. Divergent population trajectories ####
##Taxonomic probabilities for divergent trajectories###
summary(as.factor(key_largo$mod)) #65 model 1, 22 model 2
summary(as.factor(key_west$mod)) #73 model 1, 23 model 2

key_largo$mod.bin=ifelse(key_largo$mod=='model2',1,0)
key_west$mod.bin=ifelse(key_west$mod=='model2',1,0)

#Grouped by family
kl_fam<- key_largo %>% group_by(Family) %>% summarize(n=n(),n.d=sum(mod.bin))
kl_fam$binom.p=dbinom(kl_fam$n.d,size=kl_fam$n,prob=sum(key_largo$mod.bin)/nrow(key_largo)) #Probability of observing based on binom distribution
kw_fam<- key_west %>% group_by(Family) %>% summarize(n=n(),n.d=sum(mod.bin))
kw_fam$binom.p=dbinom(kw_fam$n.d,size=kw_fam$n,prob=sum(key_west$mod.bin)/nrow(key_west)) #Probability of observing based on binom distribution
kl_fam<- subset(kl_fam,n>2)
kw_fam<- subset(kw_fam,n>2)

print(as.data.frame(kl_fam)) #Table s3
print(as.data.frame(kw_fam)) #Table s3

#2 Predictors of divergence####
#By size
par(mfrow=c(2,1))
plot(log10(size)~rep(1,nrow(key_largo)),ylim=c(0.5,3),xlim=c(1,1.2),bty='l',xaxt='n',ylab='',type='n',xlab='',data=key_largo,yaxt='n',main='Key Largo')
points(log10(size)~jitter(rep(1.05,nrow(subset(key_largo,mod=='model1'))),amount=0.02),data=subset(key_largo,mod=='model1'),cex=2,col='white',bg=adjustcolor('turquoise4',alpha.f=0.7),pch=21)
points(log10(size)~jitter(rep(1.15,nrow(subset(key_largo,mod=='model2'))),amount=0.02),data=subset(key_largo,mod=='model2'),cex=2,col='white',bg=adjustcolor('goldenrod',alpha.f=0.7),pch=21)
mtext(side=1,at=c(1.05,1.15),c('Congruent','Divergent'),line=1)
axis(2, col="black", at=seq(0,3,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(1),expression(10),expression(100),expression(1000)))
pow <- 0:3
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(2, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
mtext(expression(paste("Total Length (cm)")),side=2,line = 3,cex=1)

plot(log10(size)~rep(1,nrow(key_west)),ylim=c(0.5,3),xlim=c(1,1.2),bty='l',xaxt='n',ylab='',type='n',xlab='',data=key_west,yaxt='n',main='Key West')
points(log10(size)~jitter(rep(1.05,nrow(subset(key_west,mod=='model1'))),amount=0.02),data=subset(key_west,mod=='model1'),cex=2,col='white',bg=adjustcolor('turquoise4',alpha.f=0.7),pch=21)
points(log10(size)~jitter(rep(1.15,nrow(subset(key_west,mod=='model2'))),amount=0.02),data=subset(key_west,mod=='model2'),cex=2,col='white',bg=adjustcolor('goldenrod',alpha.f=0.7),pch=21)
mtext(side=1,at=c(1.05,1.15),c('Congruent','Divergent'),line=1)
axis(2, col="black", at=seq(0,3,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(1),expression(10),expression(100),expression(1000)))
pow <- 0:3
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(2, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
mtext(expression(paste("Total Length (cm)")),side=2,line = 3,cex=1)

summary(lm(log10(size)~as.factor(mod)-1,data=key_largo))
confint(lm(log10(size)~as.factor(mod)-1,data=key_largo))
summary(lm(log10(size)~as.factor(mod)-1,data=key_west))
confint(lm(log10(size)~as.factor(mod)-1,data=key_west))
#By mean abundance
par(mfrow=c(2,2))
plot(log10(mean.x.rvc)~rep(1,nrow(key_largo)),xlim=c(1,1.2),ylim=c(-3,1.6),bty='l',xaxt='n',ylab='Mean Abundance (RVC)',type='n',xlab='',data=key_largo,main='Key Largo',yaxt='n')
x.1<- jitter(rep(1.05,nrow(subset(key_largo,mod=='model1'))),amount=0.02)
x.2<- jitter(rep(1.15,nrow(subset(key_largo,mod=='model2'))),amount=0.02)
points(log10(mean.x.rvc)~x.1,data=subset(key_largo,mod=='model1'),cex=2,col='white',bg=adjustcolor('turquoise4',alpha.f=0.7),pch=21)
points(log10(mean.x.rvc)~x.2,data=subset(key_largo,mod=='model2'),cex=2,col='white',bg=adjustcolor('goldenrod',alpha.f=0.7),pch=21)
for(i in 1:nrow(subset(key_largo,mod=='model1'))){
  lines(c(log10(subset(key_largo,mod=='model1')$mean.x.rvc.l95[i]),log10(subset(key_largo,mod=='model1')$mean.x.rvc.u95[i]))~rep(x.1[i],2),lwd=1,col=adjustcolor('turquoise4',alpha.f=0.4))
}
for(i in 1:nrow(subset(key_largo,mod=='model2'))){
  lines(c(log10(subset(key_largo,mod=='model2')$mean.x.rvc.l95[i]),log10(subset(key_largo,mod=='model2')$mean.x.rvc.u95[i]))~rep(x.2[i],2),lwd=1,col=adjustcolor('goldenrod',alpha.f=0.4))
}
axis(2, col="black", at=seq(-3,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(''),expression(0.01),expression(1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(2, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
mtext(side=1,at=c(1.05,1.15),c('Congruent','Divergent'),line=1)

plot(log10(mean.x.reef)~rep(1,nrow(key_largo)),xlim=c(1,1.2),ylim=c(-3,1.6),bty='l',xaxt='n',ylab='Mean Abundance (REEF)',type='n',xlab='',data=key_largo,main='Key Largo',yaxt='n')
x.1<- jitter(rep(1.05,nrow(subset(key_largo,mod=='model1'))),amount=0.02)
x.2<- jitter(rep(1.15,nrow(subset(key_largo,mod=='model2'))),amount=0.02)
points(log10(mean.x.reef)~x.1,data=subset(key_largo,mod=='model1'),cex=2,col='white',bg=adjustcolor('turquoise4',alpha.f=0.7),pch=21)
points(log10(mean.x.reef)~x.2,data=subset(key_largo,mod=='model2'),cex=2,col='white',bg=adjustcolor('goldenrod',alpha.f=0.7),pch=21)
for(i in 1:nrow(subset(key_largo,mod=='model1'))){
  lines(c(log10(subset(key_largo,mod=='model1')$mean.x.reef.l95[i]),log10(subset(key_largo,mod=='model1')$mean.x.reef.u95[i]))~rep(x.1[i],2),lwd=1,col=adjustcolor('turquoise4',alpha.f=0.4))
}
for(i in 1:nrow(subset(key_largo,mod=='model2'))){
  lines(c(log10(subset(key_largo,mod=='model2')$mean.x.reef.l95[i]),log10(subset(key_largo,mod=='model2')$mean.x.reef.u95[i]))~rep(x.2[i],2),lwd=1,col=adjustcolor('goldenrod',alpha.f=0.4))
}
axis(2, col="black", at=seq(-3,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(''),expression(0.01),expression(1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(2, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
mtext(side=1,at=c(1.05,1.15),c('Congruent','Divergent'),line=1)

plot(log10(mean.x.rvc)~rep(1,nrow(key_west)),xlim=c(1,1.2),ylim=c(-3,1.6),bty='l',xaxt='n',ylab='Mean Abundance (RVC)',type='n',xlab='',data=key_west,main='Key West',yaxt='n')
x.1<- jitter(rep(1.05,nrow(subset(key_west,mod=='model1'))),amount=0.02)
x.2<- jitter(rep(1.15,nrow(subset(key_west,mod=='model2'))),amount=0.02)
points(log10(mean.x.rvc)~x.1,data=subset(key_west,mod=='model1'),cex=2,col='white',bg=adjustcolor('turquoise4',alpha.f=0.7),pch=21)
points(log10(mean.x.rvc)~x.2,data=subset(key_west,mod=='model2'),cex=2,col='white',bg=adjustcolor('goldenrod',alpha.f=0.7),pch=21)
for(i in 1:nrow(subset(key_west,mod=='model1'))){
  lines(c(log10(subset(key_west,mod=='model1')$mean.x.rvc.l95[i]),log10(subset(key_west,mod=='model1')$mean.x.rvc.u95[i]))~rep(x.1[i],2),lwd=1,col=adjustcolor('turquoise4',alpha.f=0.4))
}
for(i in 1:nrow(subset(key_west,mod=='model2'))){
  lines(c(log10(subset(key_west,mod=='model2')$mean.x.rvc.l95[i]),log10(subset(key_west,mod=='model2')$mean.x.rvc.u95[i]))~rep(x.2[i],2),lwd=1,col=adjustcolor('goldenrod',alpha.f=0.4))
}
axis(2, col="black", at=seq(-3,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(''),expression(0.01),expression(1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(2, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
mtext(side=1,at=c(1.05,1.15),c('Congruent','Divergent'),line=1)

plot(log10(mean.x.reef)~rep(1,nrow(key_west)),xlim=c(1,1.2),ylim=c(-3,1.6),bty='l',xaxt='n',ylab='Mean Abundance (REEF)',type='n',xlab='',data=key_west,main='Key West',yaxt='n')
x.1<- jitter(rep(1.05,nrow(subset(key_west,mod=='model1'))),amount=0.02)
x.2<- jitter(rep(1.15,nrow(subset(key_west,mod=='model2'))),amount=0.02)
points(log10(mean.x.reef)~x.1,data=subset(key_west,mod=='model1'),cex=2,col='white',bg=adjustcolor('turquoise4',alpha.f=0.7),pch=21)
points(log10(mean.x.reef)~x.2,data=subset(key_west,mod=='model2'),cex=2,col='white',bg=adjustcolor('goldenrod',alpha.f=0.7),pch=21)
for(i in 1:nrow(subset(key_west,mod=='model1'))){
  lines(c(log10(subset(key_west,mod=='model1')$mean.x.reef.l95[i]),log10(subset(key_west,mod=='model1')$mean.x.reef.u95[i]))~rep(x.1[i],2),lwd=1,col=adjustcolor('turquoise4',alpha.f=0.4))
}
for(i in 1:nrow(subset(key_west,mod=='model2'))){
  lines(c(log10(subset(key_west,mod=='model2')$mean.x.reef.l95[i]),log10(subset(key_west,mod=='model2')$mean.x.reef.u95[i]))~rep(x.2[i],2),lwd=1,col=adjustcolor('goldenrod',alpha.f=0.4))
}
axis(2, col="black", at=seq(-3,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(''),expression(0.01),expression(1),expression(1),expression(10),expression(100)))
pow <- -3:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(2, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
mtext(side=1,at=c(1.05,1.15),c('Congruent','Divergent'),line=1)


#Key Largo - mean abundance comparison between congruent (model 1) and divergent (model 2) species 
m1_kl<- mean_abund_list_kl[which(key_largo$mod=='model1')]
m2_kl<- mean_abund_list_kl[which(key_largo$mod=='model2')]

m1_rvc_m<- NA
m1_reef_m<- NA
mean_x_kl1<- list()
for(q in 1:nrow(m1_kl[[1]])){
  mean_x_kl1[[q]]<- log10(m1_kl[[1]][q,1:2])
  for(z in 2:length(m1_kl)){
    mean_x_kl1[[q]]<- rbind(mean_x_kl1[[q]],log10(m1_kl[[z]][q,1:2]))
  }   
  m1_rvc_m[q]<- mean(mean_x_kl1[[q]][,1])
  m1_reef_m[q]<- mean(mean_x_kl1[[q]][,2])
}

m2_rvc_m<- NA
m2_reef_m<- NA
mean_x_kl2<- list()
for(q in 1:nrow(m2_kl[[1]])){
  mean_x_kl2[[q]]<- log10(m2_kl[[1]][q,1:2])
  for(z in 2:length(m2_kl)){
    mean_x_kl2[[q]]<- rbind(mean_x_kl2[[q]],log10(m2_kl[[z]][q,1:2]))
  }   
  m2_rvc_m[q]<- mean(mean_x_kl2[[q]][,1])
  m2_reef_m[q]<- mean(mean_x_kl2[[q]][,2])
}

median(10^m1_rvc_m)
quantile(10^m1_rvc_m,0.025)
quantile(10^m1_rvc_m,0.975)

median(10^m2_rvc_m)
quantile(10^m2_rvc_m,0.025)
quantile(10^m2_rvc_m,0.975)

median(10^m1_reef_m)
quantile(10^m1_reef_m,0.025)
quantile(10^m1_reef_m,0.975)

median(10^m2_reef_m)
quantile(10^m2_reef_m,0.025)
quantile(10^m2_reef_m,0.975)


m1_kw<- mean_abund_list_kw[which(key_west$mod=='model1')]
m2_kw<- mean_abund_list_kw[which(key_west$mod=='model2')]
m1_kw_abund<- do.call(rbind, lapply(m1_kw, data.frame, stringsAsFactors=FALSE))

m1_rvc_m<- NA
m1_reef_m<- NA
mean_x_kw1<- list()
for(q in 1:nrow(m1_kw[[1]])){
  mean_x_kw1[[q]]<- log10(m1_kw[[1]][q,1:2])
  for(z in 2:length(m1_kw)){
    mean_x_kw1[[q]]<- rbind(mean_x_kw1[[q]],log10(m1_kw[[z]][q,1:2]))
  }   
  m1_rvc_m[q]<- mean(mean_x_kw1[[q]][,1])
  m1_reef_m[q]<- mean(mean_x_kw1[[q]][,2])
}

m2_rvc_m<- NA
m2_reef_m<- NA
mean_x_kw2<- list()
for(q in 1:nrow(m2_kw[[1]])){
  mean_x_kw2[[q]]<- log10(m2_kw[[1]][q,1:2])
  for(z in 2:length(m2_kw)){
    mean_x_kw2[[q]]<- rbind(mean_x_kw2[[q]],log10(m2_kw[[z]][q,1:2]))
  }   
  m2_rvc_m[q]<- mean(mean_x_kw2[[q]][,1])
  m2_reef_m[q]<- mean(mean_x_kw2[[q]][,2])
}

median(10^m1_rvc_m)
quantile(10^m1_rvc_m,0.025)
quantile(10^m1_rvc_m,0.975)

median(10^m2_rvc_m)
quantile(10^m2_rvc_m,0.025)
quantile(10^m2_rvc_m,0.975)

median(10^m1_reef_m)
quantile(10^m1_reef_m,0.025)
quantile(10^m1_reef_m,0.975)

median(10^m2_reef_m)
quantile(10^m2_reef_m,0.025)
quantile(10^m2_reef_m,0.975)

#By schooling
kl_school<- key_largo %>% group_by(as.factor(Schooling)) %>% summarize(n=n(),n.d=sum(mod.bin))
kl_school$binom.p=dbinom(kl_school$n.d,size=kl_school$n,prob=sum(key_largo$mod.bin)/nrow(key_largo))
kw_school<- key_west %>% group_by(as.factor(Schooling)) %>% summarize(n=n(),n.d=sum(mod.bin))
kw_school$binom.p=dbinom(kw_school$n.d,size=kw_school$n,prob=sum(key_west$mod.bin)/nrow(key_west))
print(as.data.frame(kl_school)) #Table Sx
print(as.data.frame(kw_school)) #Table Sx


#By colour
kl_col<- key_largo %>% group_by(as.factor(color)) %>% summarize(n=n(),n.d=sum(mod.bin))
kl_col$binom.p=dbinom(kl_col$n.d,size=kl_col$n,prob=sum(key_largo$mod.bin)/nrow(key_largo))
kw_col<- key_west %>% group_by(as.factor(color)) %>% summarize(n=n(),n.d=sum(mod.bin))
kw_col$binom.p=dbinom(kw_col$n.d,size=kw_col$n,prob=sum(key_west$mod.bin)/nrow(key_west))
print(as.data.frame(kl_col)) #Table Sx
print(as.data.frame(kw_col)) #Table Sx

#By crypsis
kl_cryp<- key_largo %>% group_by(as.factor(behavior)) %>% summarize(n=n(),n.d=sum(mod.bin))
kl_cryp$binom.p=dbinom(kl_col$n.d,size=kl_col$n,prob=sum(key_largo$mod.bin)/nrow(key_largo))
kw_cryp<- key_west %>% group_by(as.factor(behavior)) %>% summarize(n=n(),n.d=sum(mod.bin))
kw_cryp$binom.p=dbinom(kw_col$n.d,size=kw_col$n,prob=sum(key_west$mod.bin)/nrow(key_west))
print(as.data.frame(kl_cryp)) #Table Sx
print(as.data.frame(kw_cryp)) #Table Sx

###3. Measurement error####
files_kl_par_2<- list.files(here('outputs','species parameter estimates','Key Largo','model 2'))
temp_var_r_q_list_kl<- list()
temp_var_r_q_long_kl<- list()
for(i in 1:nrow(key_largo)){
  sp<- read.csv(here('outputs','species parameter estimates','Key Largo','model 2',files_kl_par_2[i]))
  q1=sp[,(gsub('\\..*','',colnames(sp))=='sd_q1')]
  q2=sp[,(gsub('\\..*','',colnames(sp))=='sd_q2')]
  r1=sp[,(gsub('\\..*','',colnames(sp))=='sd_r1')]
  r2=sp[,(gsub('\\..*','',colnames(sp))=='sd_r2')]
  temp_var_r_q_list_kl[[i]]<- cbind(q1,q2,r1,r2,v.r1=c(r1^2/(q1^2+r1^2)),v.r2=c(r2^2/(q2^2+r2^2)))
}

files_kw_par_2<- list.files(here('outputs','species parameter estimates','Key West','model 2'))
temp_var_r_q_list_kw<- list()
for(i in 1:nrow(key_west)){
  sp<- read.csv(here('outputs','species parameter estimates','Key West','model 2',files_kw_par_2[i]))
  q1=sp[,(gsub('\\..*','',colnames(sp))=='sd_q1')]
  q2=sp[,(gsub('\\..*','',colnames(sp))=='sd_q2')]
  r1=sp[,(gsub('\\..*','',colnames(sp))=='sd_r1')]
  r2=sp[,(gsub('\\..*','',colnames(sp))=='sd_r2')]
  temp_var_r_q_list_kw[[i]]<- cbind(q1,q2,r1,r2,v.r1=c(r1^2/(q1^2+r1^2)),v.r2=c(r2^2/(q2^2+r2^2)))
  temp_var_r_q_list_kw[[i]]<- cbind(temp_var_r_q_list_kw[[i]],diff.r1.r2=temp_var_r_q_list_kw[[i]][,5]-temp_var_r_q_list_kw[[i]][,6])
}

##Difference between RVC and REEF
sd_kl<- data.frame(sp=key_largo$SP,v.r1=NA,v.r1.l95=NA,v.r1.u95=NA,v.r2=NA,v.r2.l95=NA,v.r2.u95=NA,diff=NA,diff.l95=NA,diff.u95=NA,sd.v.r1=NA,sd.v.r2=NA)
for(i in 1:nrow(key_largo)){
  sd_kl[i,2]=median(temp_var_r_q_list_kl[[i]][,5])
  sd_kl[i,3]=quantile(temp_var_r_q_list_kl[[i]][,5],0.05)
  sd_kl[i,4]=quantile(temp_var_r_q_list_kl[[i]][,5],0.95)
  sd_kl[i,5]=median(temp_var_r_q_list_kl[[i]][,6])
  sd_kl[i,6]=quantile(temp_var_r_q_list_kl[[i]][,6],0.05)
  sd_kl[i,7]=quantile(temp_var_r_q_list_kl[[i]][,6],0.95)
  sd_kl[i,8]=median(sort(temp_var_r_q_list_kl[[i]][,5])-sort(temp_var_r_q_list_kl[[i]][,6]))
  sd_kl[i,9]=quantile(sort(temp_var_r_q_list_kl[[i]][,5])-sort(temp_var_r_q_list_kl[[i]][,6]),0.05)
  sd_kl[i,10]=quantile(sort(temp_var_r_q_list_kl[[i]][,5])-sort(temp_var_r_q_list_kl[[i]][,6]),0.95)
  sd_kl[i,11]=sd(temp_var_r_q_list_kl[[i]][,5])
  sd_kl[i,12]=sd(temp_var_r_q_list_kl[[i]][,6])
}
par(mfrow=c(2,2))
hist(sd_kl$v.r1,breaks=20)
hist(sd_kl$v.r2,breaks=20)
hist(sd_kl$diff,breaks=20)

mean_me_list_kl=list()
mean_me_kl=matrix(ncol=3,nrow=3000)
for(q in 1:nrow(temp_var_r_q_list_kl[[1]])){
  mean_me_list_kl[[q]]<- temp_var_r_q_list_kl[[1]][q,5:6]
  for(z in 2:length(mean_abund_list_kl)){
    mean_me_list_kl[[q]]<- rbind(mean_me_list_kl[[q]],temp_var_r_q_list_kl[[z]][q,5:6])
  }   
  mean_me_kl[q,1]<- mean(mean_me_list_kl[[q]][,1])
  mean_me_kl[q,2]<- mean(mean_me_list_kl[[q]][,2])
  mean_me_kl[q,3]<- mean(mean_me_list_kl[[q]][,1]-mean_me_list_kl[[q]][,2])
}
median(mean_me_kl[,1])
quantile(mean_me_kl[,1],0.025)
quantile(mean_me_kl[,1],0.975)
median(mean_me_kl[,2])
quantile(mean_me_kl[,2],0.025)
quantile(mean_me_kl[,2],0.975)
median(mean_me_kl[,3])
quantile(mean_me_kl[,3],0.025)
quantile(mean_me_kl[,3],0.975)

summary(sd_kl$diff)
summary(sd_kl$v.r1)
summary(sd_kl$v.r2)

nrow(subset(sd_kl,diff>0)) #Number of species with evidence of higher measurement error in RVC
nrow(subset(sd_kl,diff.l95>0)) #Number of species with significant evidence suggesting higher measurement error in RVC

sd_kw<- data.frame(sp=key_west$SP,v.r1=NA,v.r1.l95=NA,v.r1.u95=NA,v.r2=NA,v.r2.l95=NA,v.r2.u95=NA,diff=NA,diff.l95=NA,diff.u95=NA,sd.v.r1=NA,sd.v.r2=NA)
for(i in 1:nrow(key_west)){
  sd_kw[i,2]=median(temp_var_r_q_list_kw[[i]][,5])
  sd_kw[i,3]=quantile(temp_var_r_q_list_kw[[i]][,5],0.025)
  sd_kw[i,4]=quantile(temp_var_r_q_list_kw[[i]][,5],0.975)
  sd_kw[i,5]=median(temp_var_r_q_list_kw[[i]][,6])
  sd_kw[i,6]=quantile(temp_var_r_q_list_kw[[i]][,6],0.025)
  sd_kw[i,7]=quantile(temp_var_r_q_list_kw[[i]][,6],0.975)
  sd_kw[i,8]=median(sort(temp_var_r_q_list_kw[[i]][,5])-sort(temp_var_r_q_list_kw[[i]][,6]))
  sd_kw[i,9]=quantile(sort(temp_var_r_q_list_kw[[i]][,5])-sort(temp_var_r_q_list_kw[[i]][,6]),0.025)
  sd_kw[i,10]=quantile(sort(temp_var_r_q_list_kw[[i]][,5])-sort(temp_var_r_q_list_kw[[i]][,6]),0.975)
  sd_kw[i,11]=sd(temp_var_r_q_list_kw[[i]][,5])
  sd_kw[i,12]=sd(temp_var_r_q_list_kw[[i]][,6])
}
par(mfrow=c(2,2))
hist(sd_kw$v.r1,breaks=20)
hist(sd_kw$v.r2,breaks=20)
hist(sd_kw$diff,breaks=20)

summary(sd_kw$diff)
summary(sd_kw$v.r1)
median(sd_kw$v.r1)
summary(sd_kw$v.r2)
median(sd_kw$v.r2)

mean_me_list_kw=list()
mean_me_kw=matrix(ncol=3,nrow=3000)
for(q in 1:nrow(temp_var_r_q_list_kw[[1]])){
  mean_me_list_kw[[q]]<- temp_var_r_q_list_kw[[1]][q,5:6]
  for(z in 2:length(mean_abund_list_kw)){
    mean_me_list_kw[[q]]<- rbind(mean_me_list_kw[[q]],temp_var_r_q_list_kw[[z]][q,5:6])
  }   
  mean_me_kw[q,1]<- mean(mean_me_list_kw[[q]][,1])
  mean_me_kw[q,2]<- mean(mean_me_list_kw[[q]][,2])
  mean_me_kw[q,3]<- mean(mean_me_list_kw[[q]][,1]-mean_me_list_kw[[q]][,2])
}
median(mean_me_kw[,1])
quantile(mean_me_kw[,1],0.025)
quantile(mean_me_kw[,1],0.975)
median(mean_me_kw[,2])
quantile(mean_me_kw[,2],0.025)
quantile(mean_me_kw[,2],0.975)
median(mean_me_kw[,3])
quantile(mean_me_kw[,3],0.025)
quantile(mean_me_kw[,3],0.975)

nrow(subset(sd_kw,diff>0)) #Number of species with evidence of higher measurement error in RVC
nrow(subset(sd_kw,diff.l95>0)) #Number of species with strong differences in measurement error

sd_kl2<- cbind(sd_kl,key_largo)
sd_kl_drvc<- subset(sd_kl2,diff<0)
sd_kl_dreef<- subset(sd_kl2,diff>0)

sd_kw2<- cbind(sd_kw,key_west)
sd_kw_drvc<- subset(sd_kw2,diff<0)
sd_kw_dreef<- subset(sd_kw2,diff>0)

hist.rvc.kl<- hist(sd_kl$v.r1,breaks=20)
hist.reef.kl<- hist(sd_kl$v.r2,breaks=20)
hist.rvc.kw<- hist(sd_kw$v.r1,breaks=20)
hist.reef.kw<- hist(sd_kw$v.r2,breaks=20)

#Figure 3
par(mfrow=c(1,2))
plot(sd_kl$v.r1~rep(1,nrow(sd_kl)),ylim=c(0,1),xlim=c(1,1.2),bty='l',xaxt='n',ylab='Prop. Measurement Error',type='n',xlab='',data=sd_kl,main='Key Largo')
x.1.r<- jitter(rep(1.05,nrow(sd_kl_drvc)),amount=0.02)
x.1.rf<- jitter(rep(1.05,nrow(sd_kl_dreef)),amount=0.02)
x.2.r<- jitter(rep(1.15,nrow(sd_kl_drvc)),amount=0.02)
x.2.rf<- jitter(rep(1.15,nrow(sd_kl_dreef)),amount=0.02)
VerticalHist(x=1,hist=hist.rvc.kl,xscale=0.04,xwidth=0.01,rev=F)
VerticalHist(x=1.2,hist=hist.reef.kl,xscale=0.04,xwidth=0.01,rev=T)
for(i in 1:nrow(sd_kl_drvc)){
  lines(c(sd_kl_drvc$v.r1[i],sd_kl_drvc$v.r2[i])~c(x.1.r[i],x.2.r[i]),lwd=1,col=adjustcolor('deepskyblue4',alpha.f=0.4))
}
for(i in 1:nrow(sd_kl_dreef)){
  lines(c(sd_kl_dreef$v.r1[i],sd_kl_dreef$v.r2[i])~c(x.1.rf[i],x.2.rf[i]),lwd=1,col=adjustcolor('brown3',alpha.f=0.4))
}
points(sd_kl_drvc$v.r1~x.1.r,data=sd_kl_drvc,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(sd_kl_dreef$v.r1~x.1.rf,data=sd_kl_drvc,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(sd_kl_drvc$v.r2~x.2.r,data=sd_kl_dreef,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
points(sd_kl_dreef$v.r2~x.2.rf,data=sd_kl_dreef,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
mtext(side=1,at=c(1.05,1.15),c('RVC','REEF'),line=1)
text(par("usr")[1],par("usr")[4]*1.02,'A',cex=1,pos=3,xpd=T,font=2)

plot(sd_kw$v.r1~rep(1,nrow(sd_kw)),ylim=c(0,1),xlim=c(1,1.2),bty='l',xaxt='n',ylab='Prop. Measurement Error',type='n',xlab='',data=sd_kw,main='Key West')
x.1.r<- jitter(rep(1.05,nrow(sd_kw_drvc)),amount=0.02)
x.1.rf<- jitter(rep(1.05,nrow(sd_kw_dreef)),amount=0.02)
x.2.r<- jitter(rep(1.15,nrow(sd_kw_drvc)),amount=0.02)
x.2.rf<- jitter(rep(1.15,nrow(sd_kw_dreef)),amount=0.02)
VerticalHist(x=1,hist=hist.rvc.kw,xscale=0.04,xwidth=0.01,rev=F)
VerticalHist(x=1.2,hist=hist.reef.kw,xscale=0.04,xwidth=0.01,rev=T)
for(i in 1:nrow(sd_kw_drvc)){
  lines(c(sd_kw_drvc$v.r1[i],sd_kw_drvc$v.r2[i])~c(x.1.r[i],x.2.r[i]),lwd=1,col=adjustcolor('deepskyblue4',alpha.f=0.4))
}
for(i in 1:nrow(sd_kw_dreef)){
  lines(c(sd_kw_dreef$v.r1[i],sd_kw_dreef$v.r2[i])~c(x.1.rf[i],x.2.rf[i]),lwd=1,col=adjustcolor('brown3',alpha.f=0.4))
}
points(sd_kw_drvc$v.r1~x.1.r,data=sd_kw_drvc,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(sd_kw_dreef$v.r1~x.1.rf,data=sd_kw_drvc,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(sd_kw_drvc$v.r2~x.2.r,data=sd_kw_dreef,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
points(sd_kw_dreef$v.r2~x.2.rf,data=sd_kw_dreef,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
mtext(side=1,at=c(1.05,1.15),c('RVC','REEF'),line=1)
text(par("usr")[1],par("usr")[4]*1.02,'B',cex=1,pos=3,xpd=T,font=2)

dev.off()

###3 -Predictors in each survey#####
#By Family
sd_kl2$n.Fam<- NA
sd_kl2<- sd_kl2 %>% group_by(Family) %>% mutate(n=n())
sd_kl2_sub<- subset(sd_kl2,n>2) #removes singletons
fam_sum_kl<- sd_kl2_sub %>% group_by(Family) %>% summarize(n=n())

kl_meas_fam<- par_agg(x=files_kl_par_2,path=here('outputs','species parameter estimates','Key Largo','model 2'),group='Family',pars=c('sd_r1','sd_r2','sd_q1','sd_q2'),ref=key_largo)
kl_meas_fam$v.r1<- kl_meas_fam$par.sd_r1/(kl_meas_fam$par.sd_r1+kl_meas_fam$par.sd_q1) #variance due to measurement error - RVC
kl_meas_fam$v.r2<- kl_meas_fam$par.sd_r2/(kl_meas_fam$par.sd_r2+kl_meas_fam$par.sd_q2) #variance due to measurement error - RVC
kl_meas_fam$diff.r1.r2<- kl_meas_fam$v.r1-kl_meas_fam$v.r2

kl_meas_comp_rvc<- par_comp(x=kl_meas_fam,levels=fam_sum_kl$Family,par='v.r1',samps=3000)
kl_meas_comp_reef<- par_comp(x=kl_meas_fam,levels=fam_sum_kl$Family,par='v.r2',samps=3000)
kl_meas_diff<- par_comp(x=kl_meas_fam,levels=fam_sum_kl$Family,par='diff.r1.r2',samps=3000)


sd_kw2$n.Fam<- NA
sd_kw2<- sd_kw2 %>% group_by(Family) %>% mutate(n=n())
sd_kw2_sub<- subset(sd_kw2,n>2) #removes singletons
fam_sum_kw<- sd_kw2_sub %>% group_by(Family) %>% summarize(n=n())

kw_meas_fam<- par_agg(x=files_kw_par_2,path=here('outputs','species parameter estimates','Key West','model 2'),group='Family',pars=c('sd_r1','sd_r2','sd_q1','sd_q2'),ref=key_west)
kw_meas_fam$v.r1<- kw_meas_fam$par.sd_r1^2/(kw_meas_fam$par.sd_r1^2+kw_meas_fam$par.sd_q1^2) #variance due to measurement error - RVC
kw_meas_fam$v.r2<- kw_meas_fam$par.sd_r2^2/(kw_meas_fam$par.sd_r2^2+kw_meas_fam$par.sd_q2^2) #variance due to measurement error - RVC
kw_meas_fam$diff.r1.r2<- kw_meas_fam$v.r1-kw_meas_fam$v.r2

kw_meas_comp_rvc<- par_comp(x=kw_meas_fam,levels=fam_sum_kw$Family,par='v.r1',samps=3000)
kw_meas_comp_reef<- par_comp(x=kw_meas_fam,levels=fam_sum_kw$Family,par='v.r2',samps=3000)
kw_meas_diff<- par_comp(x=kw_meas_fam,levels=fam_sum_kw$Family,par='diff.r1.r2',samps=3000)


par(mfrow=c(1,2))
plot(kl_meas_comp_rvc$par.m~seq(1:nrow(kl_meas_comp_rvc)),ylim=c(0,1),xlim=c(0.5,12.5),bty='l',xaxt='n',ylab='Prop. Measurement Error',type='n',xlab='',data=sd_kl,main='Key Largo')
for(i in 1:nrow(kl_meas_comp_rvc)){
  lines(c(kl_meas_comp_rvc$par.l90[i],kl_meas_comp_rvc$par.u90[i])~rep(i-0.15,2),col=adjustcolor('darkblue',alpha.f=0.4))
}
for(i in 1:nrow(kl_meas_comp_rvc)){
  lines(c(kl_meas_comp_reef$par.l90[i],kl_meas_comp_reef$par.u90[i])~rep(i+0.15,2),col=adjustcolor('darkred',alpha.f=0.4))
}
points(kl_meas_comp_rvc$par.m~c(seq(1:12)-0.15),cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kl_meas_comp_reef$par.m~c(seq(1:12)+0.15),cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
text(c(seq(from=1,to=14,by=1)+0.1),-0.01,fam_sum_kl$n)
text(0.4,-0.01,'n =')
text(c(seq(from=0.5,to=11.5,by=1)),par("usr")[3] - 0.05,kw_meas_comp_reef$group,cex=0.8,srt = 45,pos=1,xpd=T)

plot(kw_meas_comp_rvc$par.m~seq(1:nrow(kw_meas_comp_rvc)),ylim=c(0,1),xlim=c(0.5,12.5),bty='l',xaxt='n',ylab='Prop. Measurement Error',type='n',xlab='',data=sd_kw,main='Key West')
for(i in 1:nrow(kw_meas_comp_rvc)){
  lines(c(kw_meas_comp_rvc$par.l90[i],kw_meas_comp_rvc$par.u90[i])~rep(i-0.15,2),col=adjustcolor('darkblue',alpha.f=0.4))
}
for(i in 1:nrow(kw_meas_comp_rvc)){
  lines(c(kw_meas_comp_reef$par.l90[i],kw_meas_comp_reef$par.u90[i])~rep(i+0.15,2),col=adjustcolor('darkred',alpha.f=0.4))
}
points(kw_meas_comp_rvc$par.m~c(seq(1:12)-0.15),cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kw_meas_comp_reef$par.m~c(seq(1:12)+0.15),cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
text(c(seq(from=1,to=13,by=1)+0.1),-0.01,fam_sum_kw$n)
text(0.4,-0.01,'n =')
text(c(seq(from=0.5,to=11.5,by=1)),par("usr")[3] - 0.05,kw_meas_comp_reef$group,cex=0.8,srt = 45,pos=1,xpd=T)


##3a Size###
plot(sd_kl$v.r1~log10(key_largo$size),ylim=c(0,1),bty='l',xaxt='n',ylab='Prop. Measurement Error',type='n',xlab='',data=sd_kl,main='Key Largo')
for(i in 1:nrow(sd_kl_drvc)){
  lines(c(sd_kl_drvc$v.r1[i],sd_kl_drvc$v.r2[i])~c(log10(sd_kl_drvc$size)[i],log10(sd_kl_drvc$size)[i]),lwd=1,col=adjustcolor('deepskyblue4',alpha.f=0.4))
}

for(i in 1:nrow(sd_kl_dreef)){
lines(c(sd_kl_dreef$v.r1[i],sd_kl_dreef$v.r2[i])~c(log10(sd_kl_dreef$size)[i],log10(sd_kl_dreef$size)[i]),lwd=1,col=adjustcolor('brown3',alpha.f=0.4))
}
points(sd_kl$v.r1~log10(key_largo$size),data=sd_kl,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(sd_kl$v.r2~log10(key_largo$size),data=sd_kl,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
axis(1, col="black", at=seq(0,3,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(1),expression(10),expression(100),expression(1000)))
pow <- 0:3
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
mtext(expression(paste("Total Length (cm)")),side=1,line = 3,cex=1)


plot(sd_kw$v.r1~log10(key_west$size),ylim=c(0,1),bty='l',xaxt='n',ylab='Prop. Measurement Error',type='n',xlab='',data=sd_kw,main='Key west')
for(i in 1:nrow(sd_kw_drvc)){
  lines(c(sd_kw_drvc$v.r1[i],sd_kw_drvc$v.r2[i])~c(log10(sd_kw_drvc$size)[i],log10(sd_kw_drvc$size)[i]),lwd=1,col=adjustcolor('deepskyblue4',alpha.f=0.4))
}

for(i in 1:nrow(sd_kw_dreef)){
  lines(c(sd_kw_dreef$v.r1[i],sd_kw_dreef$v.r2[i])~c(log10(sd_kw_dreef$size)[i],log10(sd_kw_dreef$size)[i]),lwd=1,col=adjustcolor('brown3',alpha.f=0.4))
}
points(sd_kw$v.r1~log10(key_west$size),data=sd_kw,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(sd_kw$v.r2~log10(key_west$size),data=sd_kw,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
axis(1, col="black", at=seq(0,3,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(1),expression(10),expression(100),expression(1000)))
pow <- 0:3
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
mtext(expression(paste("Total Length (cm)")),side=1,line = 3,cex=1)


cor_size_kw_rvc<- NA
cor_size_kw_reef<- NA
for(i in 1:nrow(temp_var_r_q_list_kw[[1]])){
  s<- temp_var_r_q_list_kw[[1]][i,5:6]
  for(q in 2:nrow(key_west)){
  s<- rbind(s,temp_var_r_q_list_kw[[q]][i,5:6])  
  }
  s<- cbind(s,key_west$size)
  cor_size_kw_rvc[i]<- cor(qlogis(s[,1]),log10(s[,3]))
  cor_size_kw_reef[i]<- cor(qlogis(s[,2]),log10(s[,3]))
}

median(cor_size_kw_rvc)
quantile(cor_size_kw_rvc,0.025)
quantile(cor_size_kw_rvc,0.975)

median(cor_size_kw_reef)
quantile(cor_size_kw_reef,0.025)
quantile(cor_size_kw_reef,0.975)

cor_size_kl_rvc<- NA
cor_size_kl_reef<- NA
for(i in 1:nrow(temp_var_r_q_list_kl[[1]])){
  s<- temp_var_r_q_list_kl[[1]][i,5:6]
  for(q in 2:nrow(key_largo)){
    s<- rbind(s,temp_var_r_q_list_kl[[q]][i,5:6])  
  }
  s<- cbind(s,key_largo$size)
  cor_size_kl_rvc[i]<- cor(qlogis(s[,1]),log10(s[,3]))
  cor_size_kl_reef[i]<- cor(qlogis(s[,2]),log10(s[,3]))
}

median(cor_size_kl_rvc)
quantile(cor_size_kl_rvc,0.025)
quantile(cor_size_kl_rvc,0.975)

median(cor_size_kl_reef)
quantile(cor_size_kl_reef,0.025)
quantile(cor_size_kl_reef,0.975)


###Mean Abundance###
plot(sd_kl$v.r1~log10(key_largo$mean.x.rvc),ylim=c(0,1),xlim=c(-2.5,2),bty='l',xaxt='n',ylab='Prop. Measurement Error',type='n',xlab='',data=sd_kl,main='Key Largo')
for(i in 1:nrow(sd_kl)){
  lines(c(sd_kl$v.r1[i],sd_kl$v.r1[i])~c(log10(key_largo$mean.x.rvc.l95)[i],log10(key_largo$mean.x.rvc.u95)[i]),lwd=1,col=adjustcolor('darkblue',alpha.f=0.4))
}
for(i in 1:nrow(sd_kl)){
  lines(c(sd_kl$v.r2[i],sd_kl$v.r2[i])~c(log10(key_largo$mean.x.reef.l95)[i],log10(key_largo$mean.x.reef.u95)[i]),lwd=1,col=adjustcolor('darkred',alpha.f=0.4))
}
points(sd_kl$v.r1~log10(key_largo$mean.x.rvc),data=sd_kl,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(sd_kl$v.r2~log10(key_largo$mean.x.reef),data=sd_kl,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
axis(1, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -2:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
mtext(expression(paste("Mean Abundance")),side=1,line = 3,cex=1)

plot(sd_kw$v.r1~log10(key_west$mean.x.rvc),ylim=c(0,1),xlim=c(-2.5,2),bty='l',xaxt='n',ylab='Prop. Measurement Error',type='n',xlab='',data=sd_kw,main='Key west')
for(i in 1:nrow(sd_kw)){
  lines(c(sd_kw$v.r1[i],sd_kw$v.r1[i])~c(log10(key_west$mean.x.rvc.l95)[i],log10(key_west$mean.x.rvc.u95)[i]),lwd=1,col=adjustcolor('darkblue',alpha.f=0.4))
}
for(i in 1:nrow(sd_kw)){
  lines(c(sd_kw$v.r2[i],sd_kw$v.r2[i])~c(log10(key_west$mean.x.reef.l95)[i],log10(key_west$mean.x.reef.u95)[i]),lwd=1,col=adjustcolor('darkred',alpha.f=0.4))
}
points(sd_kw$v.r1~log10(key_west$mean.x.rvc),data=sd_kw,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(sd_kw$v.r2~log10(key_west$mean.x.reef),data=sd_kw,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
axis(1, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -2:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
mtext(expression(paste("Mean Abundance")),side=1,line = 3,cex=1)

cor_abund_kw_rvc<- NA
cor_abund_kw_reef<- NA
for(i in 1:nrow(temp_var_r_q_list_kw[[1]])){
  s<- cbind(temp_var_r_q_list_kw[[1]][i,5],temp_var_r_q_list_kw[[1]][i,6],mean_abund_list_kw[[1]][i,1],mean_abund_list_kw[[1]][i,2])
  for(q in 2:nrow(key_west)){
    s<- rbind(s,cbind(temp_var_r_q_list_kw[[q]][i,5],temp_var_r_q_list_kw[[q]][i,6],mean_abund_list_kw[[q]][i,1],mean_abund_list_kw[[q]][i,2]))  
  }
  cor_abund_kw_rvc[i]<- cor(qlogis(s[,1]),log10(s[,3]))
  cor_abund_kw_reef[i]<- cor(qlogis(s[,2]),log10(s[,4]))
}

median(cor_abund_kw_rvc)
quantile(cor_abund_kw_rvc,0.025)
quantile(cor_abund_kw_rvc,0.975)

median(cor_abund_kw_reef)
quantile(cor_abund_kw_reef,0.025)
quantile(cor_abund_kw_reef,0.975)


cor_abund_kl_rvc<- NA
cor_abund_kl_reef<- NA
for(i in 1:nrow(temp_var_r_q_list_kl[[1]])){
  s<- cbind(temp_var_r_q_list_kl[[1]][i,5],temp_var_r_q_list_kl[[1]][i,6],mean_abund_list_kl[[1]][i,1],mean_abund_list_kl[[1]][i,2])
  for(q in 2:nrow(key_largo)){
    s<- rbind(s,cbind(temp_var_r_q_list_kl[[q]][i,5],temp_var_r_q_list_kl[[q]][i,6],mean_abund_list_kl[[q]][i,1],mean_abund_list_kl[[q]][i,2]))  
  }
  cor_abund_kl_rvc[i]<- cor(qlogis(s[,1]),log10(s[,3]))
  cor_abund_kl_reef[i]<- cor(qlogis(s[,2]),log10(s[,4]))
}

median(cor_abund_kl_rvc)
quantile(cor_abund_kl_rvc,0.025)
quantile(cor_abund_kl_rvc,0.975)

median(cor_abund_kl_reef)
quantile(cor_abund_kl_reef,0.025)
quantile(cor_abund_kl_reef,0.975)


#Schooling
kl_sd_scho0<- subset(sd_kl2,Schooling==0)
kl_sd_scho0rvc<- subset(kl_sd_scho0,diff<0)
kl_sd_scho0reef<- subset(kl_sd_scho0,diff>0)
kl_sd_scho1<- subset(sd_kl2,Schooling==1)
kl_sd_scho1rvc<- subset(kl_sd_scho1,diff<0)
kl_sd_scho1reef<- subset(kl_sd_scho1,diff>0)

kw_sd_scho0<- subset(sd_kw2,Schooling==0)
kw_sd_scho0rvc<- subset(kw_sd_scho0,diff<0)
kw_sd_scho0reef<- subset(kw_sd_scho0,diff>0)
kw_sd_scho1<- subset(sd_kw2,Schooling==1)
kw_sd_scho1rvc<- subset(kw_sd_scho1,diff<0)
kw_sd_scho1reef<- subset(kw_sd_scho1,diff>0)


kl_meas_schol<- par_agg(x=files_kl_par_2,path=here('outputs','species parameter estimates','Key Largo','model 2'),group='Schooling',pars=c('sd_r1','sd_r2','sd_q1','sd_q2'),ref=key_largo)
kl_meas_schol$v.r1<- kl_meas_schol$par.sd_r1/(kl_meas_fam$par.sd_r1+kl_meas_fam$par.sd_q1) #variance due to measurement error - RVC
kl_meas_schol$v.r2<- kl_meas_schol$par.sd_r2/(kl_meas_fam$par.sd_r2+kl_meas_fam$par.sd_q2) #variance due to measurement error - RVC

kl_meas_schol_rvc<- par_comp(x=kl_meas_schol,levels=c(0,1),par='v.r1',samps=3000)
kl_meas_comp_reef<- par_comp(x=kl_meas_schol,levels=c(0,1),par='v.r2',samps=3000)

kw_meas_schol<- par_agg(x=files_kw_par_2,path=here('outputs','species parameter estimates','Key West','model 2'),group='Schooling',pars=c('sd_r1','sd_r2','sd_q1','sd_q2'),ref=key_west)
kw_meas_schol$v.r1<- kw_meas_schol$par.sd_r1/(kw_meas_fam$par.sd_r1+kw_meas_fam$par.sd_q1) #variance due to measurement error - RVC
kw_meas_schol$v.r2<- kw_meas_schol$par.sd_r2/(kw_meas_fam$par.sd_r2+kw_meas_fam$par.sd_q2) #variance due to measurement error - RVC

kw_meas_schol_rvc<- par_comp(x=kw_meas_schol,levels=c(0,1),par='v.r1',samps=3000)
kw_meas_comp_reef<- par_comp(x=kw_meas_schol,levels=c(0,1),par='v.r2',samps=3000)


par(mfrow=c(1,2))
plot(sd_kl$v.r1~rep(1,nrow(sd_kl)),ylim=c(0,1),xlim=c(1,1.2),bty='l',xaxt='n',ylab='Prop. Measurement Error',type='n',xlab='',data=sd_kl,main='Key Largo')
x.1.r<- jitter(rep(1.05,nrow(kl_sd_scho0rvc)),amount=0.02)
x.1.rf<- jitter(rep(1.05,nrow(kl_sd_scho0reef)),amount=0.02)
x.2.r<- jitter(rep(1.15,nrow(kl_sd_scho1rvc)),amount=0.02)
x.2.rf<- jitter(rep(1.15,nrow(kl_sd_scho1reef)),amount=0.02)
for(i in 1:nrow(kl_sd_scho0rvc)){
  lines(c(kl_sd_scho0rvc$v.r1[i],kl_sd_scho0rvc$v.r2[i])~c(x.1.r[i],x.1.r[i]),lwd=1,col=adjustcolor('deepskyblue4',alpha.f=0.4))
}
for(i in 1:nrow(kl_sd_scho0reef)){
  lines(c(kl_sd_scho0reef$v.r1[i],kl_sd_scho0reef$v.r2[i])~c(x.1.rf[i],x.1.rf[i]),lwd=1,col=adjustcolor('brown3',alpha.f=0.4))
}
points(kl_sd_scho0rvc$v.r1~x.1.r,data=kl_sd_scho0rvc,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kl_sd_scho0reef$v.r1~x.1.rf,data=kl_sd_scho0reef,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kl_sd_scho0rvc$v.r2~x.1.r,data=kl_sd_scho0rvc,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
points(kl_sd_scho0reef$v.r2~x.1.rf,data=kl_sd_scho0reef,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)

for(i in 1:nrow(kl_sd_scho1rvc)){
  lines(c(kl_sd_scho1rvc$v.r1[i],kl_sd_scho1rvc$v.r2[i])~c(x.2.r[i],x.2.r[i]),lwd=1,col=adjustcolor('deepskyblue4',alpha.f=0.4))
}
for(i in 1:nrow(kl_sd_scho1reef)){
  lines(c(kl_sd_scho1reef$v.r1[i],kl_sd_scho1reef$v.r2[i])~c(x.2.rf[i],x.2.rf[i]),lwd=1,col=adjustcolor('brown3',alpha.f=0.4))
}
points(kl_sd_scho1rvc$v.r1~x.2.r,data=kl_sd_scho1rvc,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kl_sd_scho1reef$v.r1~x.2.rf,data=kl_sd_scho1reef,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kl_sd_scho1rvc$v.r2~x.2.r,data=kl_sd_scho1rvc,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
points(kl_sd_scho1reef$v.r2~x.2.rf,data=kl_sd_scho1reef,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
mtext(side=1,at=c(1.05,1.15),c('Non-Schooling','Schooling'),line=1)

plot(sd_kw$v.r1~rep(1,nrow(sd_kw)),ylim=c(0,1),xlim=c(1,1.2),bty='l',xaxt='n',ylab='Prop. Measurement Error',type='n',xlab='',data=sd_kw,main='Key West')
x.1.r<- jitter(rep(1.05,nrow(kw_sd_scho0rvc)),amount=0.02)
x.1.rf<- jitter(rep(1.05,nrow(kw_sd_scho0reef)),amount=0.02)
x.2.r<- jitter(rep(1.15,nrow(kw_sd_scho1rvc)),amount=0.02)
x.2.rf<- jitter(rep(1.15,nrow(kw_sd_scho1reef)),amount=0.02)
for(i in 1:nrow(kw_sd_scho0rvc)){
  lines(c(kw_sd_scho0rvc$v.r1[i],kw_sd_scho0rvc$v.r2[i])~c(x.1.r[i],x.1.r[i]),lwd=1,col=adjustcolor('deepskyblue4',alpha.f=0.4))
}
for(i in 1:nrow(kw_sd_scho0reef)){
  lines(c(kw_sd_scho0reef$v.r1[i],kw_sd_scho0reef$v.r2[i])~c(x.1.rf[i],x.1.rf[i]),lwd=1,col=adjustcolor('brown3',alpha.f=0.4))
}
points(kw_sd_scho0rvc$v.r1~x.1.r,data=kw_sd_scho0rvc,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kw_sd_scho0reef$v.r1~x.1.rf,data=kw_sd_scho0reef,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kw_sd_scho0rvc$v.r2~x.1.r,data=kw_sd_scho0rvc,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
points(kw_sd_scho0reef$v.r2~x.1.rf,data=kw_sd_scho0reef,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)

for(i in 1:nrow(kw_sd_scho1rvc)){
  lines(c(kw_sd_scho1rvc$v.r1[i],kw_sd_scho1rvc$v.r2[i])~c(x.2.r[i],x.2.r[i]),lwd=1,col=adjustcolor('deepskyblue4',alpha.f=0.4))
}
for(i in 1:nrow(kw_sd_scho1reef)){
  lines(c(kw_sd_scho1reef$v.r1[i],kw_sd_scho1reef$v.r2[i])~c(x.2.rf[i],x.2.rf[i]),lwd=1,col=adjustcolor('brown3',alpha.f=0.4))
}
points(kw_sd_scho1rvc$v.r1~x.2.r,data=kw_sd_scho1rvc,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kw_sd_scho1reef$v.r1~x.2.rf,data=kw_sd_scho1reef,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kw_sd_scho1rvc$v.r2~x.2.r,data=kw_sd_scho1rvc,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
points(kw_sd_scho1reef$v.r2~x.2.rf,data=kw_sd_scho1reef,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
mtext(side=1,at=c(1.05,1.15),c('Non-Schooling','Schooling'),line=1)

dev.off()
##Analysis


#Color
kl_sd_col0<- subset(sd_kl2,color=='drab')
kl_sd_col0rvc<- subset(kl_sd_col0,diff<0)
kl_sd_col0reef<- subset(kl_sd_col0,diff>0)
kl_sd_col1<- subset(sd_kl2,color=='colorful')
kl_sd_col1rvc<- subset(kl_sd_col1,diff<0)
kl_sd_col1reef<- subset(kl_sd_col1,diff>0)

kw_sd_col0<- subset(sd_kw2,color=='drab')
kw_sd_col0rvc<- subset(kw_sd_col0,diff<0)
kw_sd_col0reef<- subset(kw_sd_col0,diff>0)
kw_sd_col1<- subset(sd_kw2,color=='colorful')
kw_sd_col1rvc<- subset(kw_sd_col1,diff<0)
kw_sd_col1reef<- subset(kw_sd_col1,diff>0)

par(mfrow=c(1,2))
plot(sd_kl$v.r1~rep(1,nrow(sd_kl)),ylim=c(0,1),xlim=c(1,1.2),bty='l',xaxt='n',ylab='Prop. Measurement Error',type='n',xlab='',data=sd_kl,main='Key Largo')
x.1.r<- jitter(rep(1.05,nrow(kl_sd_col0rvc)),amount=0.02)
x.1.rf<- jitter(rep(1.05,nrow(kl_sd_col0reef)),amount=0.02)
x.2.r<- jitter(rep(1.15,nrow(kl_sd_col1rvc)),amount=0.02)
x.2.rf<- jitter(rep(1.15,nrow(kl_sd_col1reef)),amount=0.02)
for(i in 1:nrow(kl_sd_col0rvc)){
  lines(c(kl_sd_col0rvc$v.r1[i],kl_sd_col0rvc$v.r2[i])~c(x.1.r[i],x.1.r[i]),lwd=1,col=adjustcolor('deepskyblue4',alpha.f=0.4))
}
for(i in 1:nrow(kl_sd_col0reef)){
  lines(c(kl_sd_col0reef$v.r1[i],kl_sd_col0reef$v.r2[i])~c(x.1.rf[i],x.1.rf[i]),lwd=1,col=adjustcolor('brown3',alpha.f=0.4))
}
points(kl_sd_col0rvc$v.r1~x.1.r,data=kl_sd_col0rvc,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kl_sd_col0reef$v.r1~x.1.rf,data=kl_sd_col0reef,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kl_sd_col0rvc$v.r2~x.1.r,data=kl_sd_col0rvc,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
points(kl_sd_col0reef$v.r2~x.1.rf,data=kl_sd_col0reef,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)

for(i in 1:nrow(kl_sd_col1rvc)){
  lines(c(kl_sd_col1rvc$v.r1[i],kl_sd_col1rvc$v.r2[i])~c(x.2.r[i],x.2.r[i]),lwd=1,col=adjustcolor('deepskyblue4',alpha.f=0.4))
}
for(i in 1:nrow(kl_sd_col1reef)){
  lines(c(kl_sd_col1reef$v.r1[i],kl_sd_col1reef$v.r2[i])~c(x.2.rf[i],x.2.rf[i]),lwd=1,col=adjustcolor('brown3',alpha.f=0.4))
}
points(kl_sd_col1rvc$v.r1~x.2.r,data=kl_sd_col1rvc,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kl_sd_col1reef$v.r1~x.2.rf,data=kl_sd_col1reef,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kl_sd_col1rvc$v.r2~x.2.r,data=kl_sd_col1rvc,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
points(kl_sd_col1reef$v.r2~x.2.rf,data=kl_sd_col1reef,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
mtext(side=1,at=c(1.05,1.15),c('Drab Species','Colorful Species'),line=1)


plot(sd_kw$v.r1~rep(1,nrow(sd_kw)),ylim=c(0,1),xlim=c(1,1.2),bty='l',xaxt='n',ylab='Prop. Measurement Error',type='n',xlab='',data=sd_kw,main='Key West')
x.1.r<- jitter(rep(1.05,nrow(kw_sd_col0rvc)),amount=0.02)
x.1.rf<- jitter(rep(1.05,nrow(kw_sd_col0reef)),amount=0.02)
x.2.r<- jitter(rep(1.15,nrow(kw_sd_col1rvc)),amount=0.02)
x.2.rf<- jitter(rep(1.15,nrow(kw_sd_col1reef)),amount=0.02)
for(i in 1:nrow(kw_sd_col0rvc)){
  lines(c(kw_sd_col0rvc$v.r1[i],kw_sd_col0rvc$v.r2[i])~c(x.1.r[i],x.1.r[i]),lwd=1,col=adjustcolor('deepskyblue4',alpha.f=0.4))
}
for(i in 1:nrow(kw_sd_col0reef)){
  lines(c(kw_sd_col0reef$v.r1[i],kw_sd_col0reef$v.r2[i])~c(x.1.rf[i],x.1.rf[i]),lwd=1,col=adjustcolor('brown3',alpha.f=0.4))
}
points(kw_sd_col0rvc$v.r1~x.1.r,data=kw_sd_col0rvc,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kw_sd_col0reef$v.r1~x.1.rf,data=kw_sd_col0reef,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kw_sd_col0rvc$v.r2~x.1.r,data=kw_sd_col0rvc,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
points(kw_sd_col0reef$v.r2~x.1.rf,data=kw_sd_col0reef,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)

for(i in 1:nrow(kw_sd_col1rvc)){
  lines(c(kw_sd_col1rvc$v.r1[i],kw_sd_col1rvc$v.r2[i])~c(x.2.r[i],x.2.r[i]),lwd=1,col=adjustcolor('deepskyblue4',alpha.f=0.4))
}
for(i in 1:nrow(kw_sd_col1reef)){
  lines(c(kw_sd_col1reef$v.r1[i],kw_sd_col1reef$v.r2[i])~c(x.2.rf[i],x.2.rf[i]),lwd=1,col=adjustcolor('brown3',alpha.f=0.4))
}
points(kw_sd_col1rvc$v.r1~x.2.r,data=kw_sd_col1rvc,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kw_sd_col1reef$v.r1~x.2.rf,data=kw_sd_col1reef,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kw_sd_col1rvc$v.r2~x.2.r,data=kw_sd_col1rvc,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
points(kw_sd_col1reef$v.r2~x.2.rf,data=kw_sd_col1reef,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
mtext(side=1,at=c(1.05,1.15),c('Drab Species','Colorful Species'),line=1)

kl_meas_coul<- par_agg(x=files_kl_par_2,path=here('outputs','species parameter estimates','Key Largo','model 2'),group='color',pars=c('sd_r1','sd_r2','sd_q1','sd_q2'),ref=key_largo)
kl_meas_coul$v.r1<- kl_meas_coul$par.sd_r1/(kl_meas_fam$par.sd_r1+kl_meas_fam$par.sd_q1) #variance due to measurement error - RVC
kl_meas_coul$v.r2<- kl_meas_coul$par.sd_r2/(kl_meas_fam$par.sd_r2+kl_meas_fam$par.sd_q2) #variance due to measurement error - RVC

kl_meas_coul_rvc<- par_comp(x=kl_meas_coul,levels=levels(as.factor(key_largo$color)),par='v.r1',samps=2400)
kl_meas_coul_reef<- par_comp(x=kl_meas_coul,levels=levels(as.factor(key_largo$color)),par='v.r2',samps=2400)

kw_meas_coul<- par_agg(x=files_kw_par_2,path=here('outputs','species parameter estimates','Key West','model 2'),group='color',pars=c('sd_r1','sd_r2','sd_q1','sd_q2'),ref=key_west)
kw_meas_coul$v.r1<- kw_meas_coul$par.sd_r1/(kw_meas_fam$par.sd_r1+kw_meas_fam$par.sd_q1) #variance due to measurement error - RVC
kw_meas_coul$v.r2<- kw_meas_coul$par.sd_r2/(kw_meas_fam$par.sd_r2+kw_meas_fam$par.sd_q2) #variance due to measurement error - RVC

kw_meas_coul_rvc<- par_comp(x=kw_meas_coul,levels=levels(as.factor(key_west$color)),par='v.r1',samps=2400)
kw_meas_coul_reef<- par_comp(x=kw_meas_coul,levels=levels(as.factor(key_west$color)),par='v.r2',samps=2400)


#behaviour
kl_sd_cryp0<- subset(sd_kl2,behavior=='conspicuous')
kl_sd_cryp0rvc<- subset(kl_sd_cryp0,diff<0)
kl_sd_cryp0reef<- subset(kl_sd_cryp0,diff>0)
kl_sd_cryp1<- subset(sd_kl2,behavior=='cryptic')
kl_sd_cryp1rvc<- subset(kl_sd_cryp1,diff<0)
kl_sd_cryp1reef<- subset(kl_sd_cryp1,diff>0)

kw_sd_cryp0<- subset(sd_kw2,behavior=='conspicuous')
kw_sd_cryp0rvc<- subset(kw_sd_cryp0,diff<0)
kw_sd_cryp0reef<- subset(kw_sd_cryp0,diff>0)
kw_sd_cryp1<- subset(sd_kw2,behavior=='cryptic')
kw_sd_cryp1rvc<- subset(kw_sd_cryp1,diff<0)
kw_sd_cryp1reef<- subset(kw_sd_cryp1,diff>0)

par(mfrow=c(1,2))
plot(sd_kl$v.r1~rep(1,nrow(sd_kl)),ylim=c(0,1),xlim=c(1,1.2),bty='l',xaxt='n',ylab='Prop. Measurement Error',type='n',xlab='',data=sd_kl,main='Key Largo')
x.1.r<- jitter(rep(1.05,nrow(kl_sd_cryp0rvc)),amount=0.02)
x.1.rf<- jitter(rep(1.05,nrow(kl_sd_cryp0reef)),amount=0.02)
x.2.r<- jitter(rep(1.15,nrow(kl_sd_cryp1rvc)),amount=0.02)
x.2.rf<- jitter(rep(1.15,nrow(kl_sd_cryp1reef)),amount=0.02)
for(i in 1:nrow(kl_sd_cryp0rvc)){
  lines(c(kl_sd_cryp0rvc$v.r1[i],kl_sd_cryp0rvc$v.r2[i])~c(x.1.r[i],x.1.r[i]),lwd=1,col=adjustcolor('deepskyblue4',alpha.f=0.4))
}
for(i in 1:nrow(kl_sd_cryp0reef)){
  lines(c(kl_sd_cryp0reef$v.r1[i],kl_sd_cryp0reef$v.r2[i])~c(x.1.rf[i],x.1.rf[i]),lwd=1,col=adjustcolor('brown3',alpha.f=0.4))
}
points(kl_sd_cryp0rvc$v.r1~x.1.r,data=kl_sd_cryp0rvc,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kl_sd_cryp0reef$v.r1~x.1.rf,data=kl_sd_cryp0reef,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kl_sd_cryp0rvc$v.r2~x.1.r,data=kl_sd_cryp0rvc,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
points(kl_sd_cryp0reef$v.r2~x.1.rf,data=kl_sd_cryp0reef,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)

for(i in 1:nrow(kl_sd_cryp1rvc)){
  lines(c(kl_sd_cryp1rvc$v.r1[i],kl_sd_cryp1rvc$v.r2[i])~c(x.2.r[i],x.2.r[i]),lwd=1,col=adjustcolor('deepskyblue4',alpha.f=0.4))
}
for(i in 1:nrow(kl_sd_cryp1reef)){
  lines(c(kl_sd_cryp1reef$v.r1[i],kl_sd_cryp1reef$v.r2[i])~c(x.2.rf[i],x.2.rf[i]),lwd=1,col=adjustcolor('brown3',alpha.f=0.4))
}
points(kl_sd_cryp1rvc$v.r1~x.2.r,data=kl_sd_cryp1rvc,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kl_sd_cryp1reef$v.r1~x.2.rf,data=kl_sd_cryp1reef,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kl_sd_cryp1rvc$v.r2~x.2.r,data=kl_sd_cryp1rvc,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
points(kl_sd_cryp1reef$v.r2~x.2.rf,data=kl_sd_cryp1reef,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
mtext(side=1,at=c(1.05,1.15),c('Conspicuous Species','Cryptic Species'),line=1)


plot(sd_kw$v.r1~rep(1,nrow(sd_kw)),ylim=c(0,1),xlim=c(1,1.2),bty='l',xaxt='n',ylab='Prop. Measurement Error',type='n',xlab='',data=sd_kw,main='Key West')
x.1.r<- jitter(rep(1.05,nrow(kw_sd_cryp0rvc)),amount=0.02)
x.1.rf<- jitter(rep(1.05,nrow(kw_sd_cryp0reef)),amount=0.02)
x.2.r<- jitter(rep(1.15,nrow(kw_sd_cryp1rvc)),amount=0.02)
x.2.rf<- jitter(rep(1.15,nrow(kw_sd_cryp1reef)),amount=0.02)
for(i in 1:nrow(kw_sd_cryp0rvc)){
  lines(c(kw_sd_cryp0rvc$v.r1[i],kw_sd_cryp0rvc$v.r2[i])~c(x.1.r[i],x.1.r[i]),lwd=1,col=adjustcolor('deepskyblue4',alpha.f=0.4))
}
for(i in 1:nrow(kw_sd_cryp0reef)){
  lines(c(kw_sd_cryp0reef$v.r1[i],kw_sd_cryp0reef$v.r2[i])~c(x.1.rf[i],x.1.rf[i]),lwd=1,col=adjustcolor('brown3',alpha.f=0.4))
}
points(kw_sd_cryp0rvc$v.r1~x.1.r,data=kw_sd_cryp0rvc,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kw_sd_cryp0reef$v.r1~x.1.rf,data=kw_sd_cryp0reef,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kw_sd_cryp0rvc$v.r2~x.1.r,data=kw_sd_cryp0rvc,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
points(kw_sd_cryp0reef$v.r2~x.1.rf,data=kw_sd_cryp0reef,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)

for(i in 1:nrow(kw_sd_cryp1rvc)){
  lines(c(kw_sd_cryp1rvc$v.r1[i],kw_sd_cryp1rvc$v.r2[i])~c(x.2.r[i],x.2.r[i]),lwd=1,col=adjustcolor('deepskyblue4',alpha.f=0.4))
}
for(i in 1:nrow(kw_sd_cryp1reef)){
  lines(c(kw_sd_cryp1reef$v.r1[i],kw_sd_cryp1reef$v.r2[i])~c(x.2.rf[i],x.2.rf[i]),lwd=1,col=adjustcolor('brown3',alpha.f=0.4))
}
points(kw_sd_cryp1rvc$v.r1~x.2.r,data=kw_sd_cryp1rvc,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kw_sd_cryp1reef$v.r1~x.2.rf,data=kw_sd_cryp1reef,cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kw_sd_cryp1rvc$v.r2~x.2.r,data=kw_sd_cryp1rvc,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
points(kw_sd_cryp1reef$v.r2~x.2.rf,data=kw_sd_cryp1reef,cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
mtext(side=1,at=c(1.05,1.15),c('Conspicuous Species','Cryptic Species'),line=1)

kl_meas_crypsis<- par_agg(x=files_kl_par_2,path=here('outputs','species parameter estimates','Key Largo','model 2'),group='behavior',pars=c('sd_r1','sd_r2','sd_q1','sd_q2'),ref=key_largo)
kl_meas_crypsis$v.r1<- kl_meas_crypsis$par.sd_r1/(kl_meas_fam$par.sd_r1+kl_meas_fam$par.sd_q1) #variance due to measurement error - RVC
kl_meas_crypsis$v.r2<- kl_meas_crypsis$par.sd_r2/(kl_meas_fam$par.sd_r2+kl_meas_fam$par.sd_q2) #variance due to measurement error - RVC

kl_meas_crypsis_rvc<- par_comp(x=kl_meas_crypsis,levels=levels(as.factor(key_largo$behavior)),par='v.r1',samps=2400)
kl_meas_crypsis_reef<- par_comp(x=kl_meas_crypsis,levels=levels(as.factor(key_largo$behavior)),par='v.r2',samps=2400)

kw_meas_crypsis<- par_agg(x=files_kw_par_2,path=here('outputs','species parameter estimates','Key West','model 2'),group='behavior',pars=c('sd_r1','sd_r2','sd_q1','sd_q2'),ref=key_west)
kw_meas_crypsis$v.r1<- kw_meas_crypsis$par.sd_r1/(kw_meas_fam$par.sd_r1+kw_meas_fam$par.sd_q1) #variance due to measurement error - RVC
kw_meas_crypsis$v.r2<- kw_meas_crypsis$par.sd_r2/(kw_meas_fam$par.sd_r2+kw_meas_fam$par.sd_q2) #variance due to measurement error - RVC

kw_meas_crypsis_rvc<- par_comp(x=kw_meas_crypsis,levels=levels(as.factor(key_west$behavior)),par='v.r1',samps=2400)
kw_meas_crypsis_reef<- par_comp(x=kw_meas_crypsis,levels=levels(as.factor(key_west$behavior)),par='v.r2',samps=2400)


####4. Population trends #####
#Key Largo
files_kl_par_2<- list.files(here('outputs','species parameter estimates','Key Largo','model 2'))
files_kl_par_1<- list.files(here('outputs','species parameter estimates','Key Largo','model 1'))
r_mat_rvc<- list()
trend_vec_rvc_kl<- list()
r_mat_reef<- list()
trend_vec_reef_kl<- list()
r_mat_comb<- list()
trend_vec_comb_kl<- list()
trend_list_kl<- list()

trends_kl<- data.frame(sp=key_largo$SP,rvc.u.median=NA,rvc.u.l90=NA,rvc.u.u90=NA,rvc.u.l95=NA,rvc.u.u95=NA,reef.u.median=NA,reef.u.l90=NA,reef.u.u90=NA,reef.u.l95=NA,reef.u.u95=NA,comb.u.median=NA,comb.u.l95=NA,comb.u.u95=NA,comb.u.l90=NA,comb.u.u90=NA,diff.trend.m=NA,diff.trend.l95=NA,diff.trend.u95=NA)
for(i in 1:nrow(key_largo)){
  sp2<- read.csv(here('outputs','species parameter estimates','Key Largo','model 2',files_kl_par_2[i]))
  x1=sp2[,(gsub('\\..*','',colnames(sp2))=='x1')]
  x2=sp2[,(gsub('\\..*','',colnames(sp2))=='x2')]
  c=sp2[,2:5][(gsub('\\..*','',names(sp2[,2:5]))=='cut')]
  m_abund_rvc<- exp(x1)
  r_mat_rvc[[i]]<- matrix(data=NA,nrow=3000,ncol=ncol(x1)-1)
  trend_vec_rvc_kl[[i]]<- rep(0,3000)
  for(t in 1:nrow(m_abund_rvc)){
    for(z in 1:(ncol(x1)-1)){
      r_mat_rvc[[i]][t,z]=log(m_abund_rvc[t,z+1])-log(m_abund_rvc[t,z])
    }
    trend_vec_rvc_kl[[i]][t]<- (log(gm_mean(exp(r_mat_rvc[[i]][t,]))))
  }
  trends_kl[i,2]<- median(trend_vec_rvc_kl[[i]])
  trends_kl[i,3]<- quantile(trend_vec_rvc_kl[[i]],0.05)
  trends_kl[i,4]<- quantile(trend_vec_rvc_kl[[i]],0.95)
  trends_kl[i,5]<- quantile(trend_vec_rvc_kl[[i]],0.025)
  trends_kl[i,6]<- quantile(trend_vec_rvc_kl[[i]],0.975)
  
  m_abund_reef<- mean_ord_to_n(x=x2,c=c)
  r_mat_reef[[i]]<- matrix(data=NA,nrow=3000,ncol=ncol(x1)-1)
  trend_vec_reef_kl[[i]]<- rep(0,3000)
  for(t in 1:nrow(m_abund_reef)){
    for(z in 1:(ncol(x2)-1)){
      r_mat_reef[[i]][t,z]=log(m_abund_reef[t,z+1])-log(m_abund_reef[t,z])
    }
    trend_vec_reef_kl[[i]][t]<- log(gm_mean(exp(r_mat_reef[[i]][t,])))
  }
  trends_kl[i,7]<- median(trend_vec_reef_kl[[i]])
  trends_kl[i,8]<- quantile(trend_vec_reef_kl[[i]],0.05)
  trends_kl[i,9]<- quantile(trend_vec_reef_kl[[i]],0.95)
  trends_kl[i,10]<- quantile(trend_vec_reef_kl[[i]],0.025)
  trends_kl[i,11]<- quantile(trend_vec_reef_kl[[i]],0.975)
  
  sp1<- read.csv(here('outputs','species parameter estimates','Key Largo','model 1',files_kl_par_1[i]))
  x=sp1[,(gsub('\\..*','',colnames(sp1))=='x')]
  c=sp1[,(gsub('\\..*','',colnames(sp1))=='cut')]
  a=sp1[,(gsub('\\..*','',colnames(sp1))=='a')]
  
  m_abund_comb<- exp(x)
  r_mat_comb[[i]]<- matrix(data=NA,nrow=3000,ncol=ncol(x)-1)
  trend_vec_comb_kl[[i]]<- rep(0,3000)
  for(t in 1:nrow(m_abund_rvc)){
    for(z in 1:(ncol(x1)-1)){
      r_mat_comb[[i]][t,z]=log(m_abund_comb[t,z+1])-log(m_abund_comb[t,z])
    }
    trend_vec_comb_kl[[i]][t]<- log(gm_mean(exp(r_mat_comb[[i]][t,])))
  }
  trends_kl[i,12]<- median(trend_vec_comb_kl[[i]])
  trends_kl[i,13]<- quantile(trend_vec_comb_kl[[i]],0.025)
  trends_kl[i,14]<- quantile(trend_vec_comb_kl[[i]],0.975)
  trends_kl[i,15]<- quantile(trend_vec_comb_kl[[i]],0.05)
  trends_kl[i,16]<- quantile(trend_vec_comb_kl[[i]],0.95)
  trends_kl[i,17]<- median(trend_vec_rvc_kl[[i]]-trend_vec_reef_kl[[i]])
  trends_kl[i,18]<- quantile(trend_vec_rvc_kl[[i]]-trend_vec_reef_kl[[i]],0.025)
  trends_kl[i,19]<- quantile(trend_vec_rvc_kl[[i]]-trend_vec_reef_kl[[i]],0.975)
  
  trend_list_kl[[i]]<- data.frame(sp=rep(key_largo$SP[i],3000),t.rvc=as.numeric(trend_vec_rvc_kl[[i]]),t.reef=as.numeric(trend_vec_reef_kl[[i]]),t.comb=as.numeric(trend_vec_comb_kl[[i]]))
}

kl_dat<- cbind(key_largo,trends_kl)
nrow(subset(kl_dat,diff.trend.m>0)) #64 species have higher trends in RVC vs REEF

#Key West
files_kw_par_2<- list.files(here('outputs','species parameter estimates','Key West','model 2'))
files_kw_par_1<- list.files(here('outputs','species parameter estimates','Key West','model 1'))
r_mat_rvc<- list()
trend_vec_rvc_kw<- list()
r_mat_reef<- list()
trend_vec_reef_kw<- list()
r_mat_comb<- list()
trend_vec_comb_kw<- list()
trend_list_kw<- list()

trends_kw<- data.frame(sp=key_west$SP,rvc.u.median=NA,rvc.u.l90=NA,rvc.u.u90=NA,rvc.u.l95=NA,rvc.u.u95=NA,reef.u.median=NA,reef.u.l90=NA,reef.u.u90=NA,reef.u.l95=NA,reef.u.u95=NA,comb.u.median=NA,comb.u.l95=NA,comb.u.u95=NA,comb.u.l90=NA,comb.u.u90=NA,diff.trend.m=NA,diff.trend.l95=NA,diff.trend.u95=NA)
for(i in 1:nrow(key_west)){
  sp2<- read.csv(here('outputs','species parameter estimates','Key West','model 2',files_kw_par_2[i]))
  x1=sp2[,(gsub('\\..*','',colnames(sp2))=='x1')]
  x2=sp2[,(gsub('\\..*','',colnames(sp2))=='x2')]
  c=sp2[,2:5][(gsub('\\..*','',names(sp2[,2:5]))=='cut')]
  m_abund_rvc<- exp(x1)
  r_mat_rvc[[i]]<- matrix(data=NA,nrow=3000,ncol=ncol(x1)-1)
  trend_vec_rvc_kw[[i]]<- rep(0,3000)
  for(t in 1:nrow(m_abund_rvc)){
    for(z in 1:(ncol(x1)-1)){
      r_mat_rvc[[i]][t,z]=log(m_abund_rvc[t,z+1])-log(m_abund_rvc[t,z])
    }
    trend_vec_rvc_kw[[i]][t]<- log(gm_mean(exp(r_mat_rvc[[i]][t,])))
  }
  trends_kw[i,2]<- median(trend_vec_rvc_kw[[i]])
  trends_kw[i,3]<- quantile(trend_vec_rvc_kw[[i]],0.05)
  trends_kw[i,4]<- quantile(trend_vec_rvc_kw[[i]],0.95)
  trends_kw[i,5]<- quantile(trend_vec_rvc_kw[[i]],0.025)
  trends_kw[i,6]<- quantile(trend_vec_rvc_kw[[i]],0.975)
  
  m_abund_reef<- mean_ord_to_n(x=x2,c=c)
  r_mat_reef[[i]]<- matrix(data=NA,nrow=3000,ncol=ncol(x1)-1)
  trend_vec_reef_kw[[i]]<- rep(0,3000)
  for(t in 1:nrow(m_abund_reef)){
    for(z in 1:(ncol(x2)-1)){
      r_mat_reef[[i]][t,z]=log(m_abund_reef[t,z+1])-log(m_abund_reef[t,z])
    }
    trend_vec_reef_kw[[i]][t]<- log(gm_mean(exp(r_mat_reef[[i]][t,])))
  }
  trends_kw[i,7]<- median(trend_vec_reef_kw[[i]])
  trends_kw[i,8]<- quantile(trend_vec_reef_kw[[i]],0.05)
  trends_kw[i,9]<- quantile(trend_vec_reef_kw[[i]],0.95)
  trends_kw[i,10]<- quantile(trend_vec_reef_kw[[i]],0.025)
  trends_kw[i,11]<- quantile(trend_vec_reef_kw[[i]],0.975)
  
  sp1<- read.csv(here('outputs','species parameter estimates','Key West','model 1',files_kw_par_1[i]))
  x=sp1[,(gsub('\\..*','',colnames(sp1))=='x')]
  c=sp1[,(gsub('\\..*','',colnames(sp1))=='cut')]
  a=sp1[,(gsub('\\..*','',colnames(sp1))=='a')]
  
  m_abund_comb<- exp(x)
  r_mat_comb[[i]]<- matrix(data=NA,nrow=3000,ncol=ncol(x)-1)
  trend_vec_comb_kw[[i]]<- rep(0,3000)
  for(t in 1:nrow(m_abund_rvc)){
    for(z in 1:(ncol(x1)-1)){
      r_mat_comb[[i]][t,z]=log(m_abund_comb[t,z+1])-log(m_abund_comb[t,z])
    }
    trend_vec_comb_kw[[i]][t]<- log(gm_mean(exp(r_mat_comb[[i]][t,])))
  }
  trends_kw[i,12]<- median(trend_vec_comb_kw[[i]])
  trends_kw[i,13]<- quantile(trend_vec_comb_kw[[i]],0.025)
  trends_kw[i,14]<- quantile(trend_vec_comb_kw[[i]],0.975)
  trends_kw[i,15]<- quantile(trend_vec_comb_kw[[i]],0.05)
  trends_kw[i,16]<- quantile(trend_vec_comb_kw[[i]],0.95)
  trends_kw[i,17]<- median(trend_vec_rvc_kw[[i]]-trend_vec_reef_kw[[i]])
  trends_kw[i,18]<- quantile(trend_vec_rvc_kw[[i]]-trend_vec_reef_kw[[i]],0.025)
  trends_kw[i,19]<- quantile(trend_vec_rvc_kw[[i]]-trend_vec_reef_kw[[i]],0.975)
  
  trend_list_kw[[i]]<- data.frame(sp=rep(key_west$SP[i],3000),t.rvc=as.numeric(trend_vec_rvc_kw[[i]]),t.reef=as.numeric(trend_vec_reef_kw[[i]]),t.comb=as.numeric(trend_vec_comb_kw[[i]]))
  
}

kw_dat<- cbind(key_west,trends_kw)
nrow(subset(kw_dat,diff.trend.m>0)) #64 species have higher trends in RVC vs REEF
#Mean trend
#Key Largo
t_kl_rvc_m<- NA
t_kl_reef_m<- NA
mean_t_kl<- list()
for(q in 1:length(trend_vec_rvc_kl[[1]])){
  mean_t_kl[[q]]<- c(trend_vec_rvc_kl[[1]][q],trend_vec_reef_kl[[1]][q])
  for(z in 2:length(trend_vec_rvc_kl)){
    mean_t_kl[[q]]<- rbind(mean_t_kl[[q]],c(trend_vec_rvc_kl[[z]][q],trend_vec_reef_kl[[z]][q]))
  }   
  t_kl_rvc_m[q]<- mean(mean_t_kl[[q]][,1])
  t_kl_reef_m[q]<- mean(mean_t_kl[[q]][,2])
}

(exp(mean(t_kl_rvc_m))-1)*100 #Mean trend in RVC Converted into a % change per annum
(exp(mean(quantile(t_kl_rvc_m,0.025)))-1)*100 #lower 95%CI
(exp(mean(quantile(t_kl_rvc_m,0.975)))-1)*100 #upper 95% CI

(exp(mean(t_kl_reef_m))-1)*100 #Mean trend in REEF Converted into a % change per annum
(exp(mean(quantile(t_kl_reef_m,0.025)))-1)*100 #lower 95%CI
(exp(mean(quantile(t_kl_reef_m,0.975)))-1)*100 #upper 95% CI

t_kw_rvc_m<- NA
t_kw_reef_m<- NA
mean_t_kw<- list()
for(q in 1:length(trend_vec_rvc_kw[[1]])){
  mean_t_kw[[q]]<- c(trend_vec_rvc_kw[[1]][q],trend_vec_reef_kw[[1]][q])
  for(z in 2:length(trend_vec_rvc_kw)){
    mean_t_kw[[q]]<- rbind(mean_t_kw[[q]],c(trend_vec_rvc_kw[[z]][q],trend_vec_reef_kw[[z]][q]))
  }   
  t_kw_rvc_m[q]<- mean(mean_t_kw[[q]][,1])
  t_kw_reef_m[q]<- mean(mean_t_kw[[q]][,2])
}

(exp(mean(t_kw_rvc_m))-1)*100 #Mean trend in RVC Converted into a % change per annum
(exp(mean(quantile(t_kw_rvc_m,0.025)))-1)*100 #lower 95%CI
(exp(mean(quantile(t_kw_rvc_m,0.975)))-1)*100 #upper 95% CI

(exp(mean(t_kw_reef_m))-1)*100 #Mean trend in REEF Converted into a % change per annum
(exp(mean(quantile(t_kw_reef_m,0.025)))-1)*100 #lower 95%CI
(exp(mean(quantile(t_kw_reef_m,0.975)))-1)*100 #upper 95% CI


#Trends by Family####
kl_dat$n.Fam<- NA
kl_dat<- kl_dat %>% group_by(Family) %>% mutate(n=n())
kl_dat_sub<- subset(kl_dat,n>2) #removes singletons
fam_sum_kl<- kl_dat_sub %>% group_by(Family) %>% summarize(n=n())
#fam_sum_kl$Family.com=c('Surgeonfishes','Butterflyfishes','Groupers','Gobies','Grunts','Wrasses','Snappers','Filefishes','Angelfishes','Damselfishes','Parrotfishes','Hamlets')

kl_trend_fam<- do.call(rbind, lapply(trend_list_kl, data.frame, stringsAsFactors=FALSE))
kl_trend_fam$t.diff<- kl_trend_fam$t.rvc-kl_trend_fam$t.reef
m<- match(kl_trend_fam$sp,key_largo$SP)
kl_trend_fam$grp<- key_largo$Family[m]
kl_trend_comp_rvc<- par_comp(x=kl_trend_fam,levels=fam_sum_kl$Family,par='t.rvc',samps=2400)
kl_trend_comp_reef<- par_comp(x=kl_trend_fam,levels=fam_sum_kl$Family,par='t.reef',samps=2400)
#kl_trend_comp_comb<- par_comp(x=kl_trend_fam,levels=fam_sum_kl$Family,par='t.comb',samps=2400)
kl_trend_comp_diff<- par_comp(x=kl_trend_fam,levels=fam_sum_kl$Family,par='t.diff',samps=2400)
#kl_trend_comp_comb<- kl_trend_comp_comb[order(kl_trend_comp_comb$par.m),]
#kl_trend_comp_reef<- kl_trend_comp_reef[rownames(kl_trend_comp_comb),]
#kl_trend_comp_rvc<- kl_trend_comp_rvc[rownames(kl_trend_comp_comb),]
#fam_sum_kl<- fam_sum_kl[rownames(kl_trend_comp_comb),]
#diff<- data.frame(Family=kl_trend_comp_rvc$group,diff=kl_trend_comp_reef$par.m-kl_trend_comp_rvc$par.m)
#diff<- diff[order(diff$diff),]
#kl_trend_comp_comb<- kl_trend_comp_comb[rownames(diff),]
#kl_trend_comp_reef<- kl_trend_comp_reef[rownames(diff),]
#kl_trend_comp_rvc<- kl_trend_comp_rvc[rownames(diff),]
#kl_trend_comp_comb$Family.c=fam_sum_kl$Family.com[match(kl_trend_comp_comb$group,fam_sum_kl$Family)]
#kl_trend_comp_reef$Family.c=fam_sum_kl$Family.com[match(kl_trend_comp_reef$group,fam_sum_kl$Family)]
#kl_trend_comp_rvc$Family.c=fam_sum_kl$Family.com[match(kl_trend_comp_rvc$group,fam_sum_kl$Family)]


kw_dat$n.Fam<- NA
kw_dat<- kw_dat %>% group_by(Family) %>% mutate(n=n())
kw_dat_sub<- subset(kw_dat,n>2) #removes singletons
fam_sum_kw<- kw_dat_sub %>% group_by(Family) %>% summarize(n=n())
#fam_sum_kw$Family.com=c('Surgeonfishes','Butterflyfishes','Groupers','Gobies','Grunts','Wrasses','Snappers','Filefishes','Angelfishes','Damselfishes','Parrotfishes','Hamlets')

kw_trend_fam<- do.call(rbind, lapply(trend_list_kw, data.frame, stringsAsFactors=FALSE))
kw_trend_fam$t.diff<- kw_trend_fam$t.rvc-kw_trend_fam$t.reef
m<- match(kw_trend_fam$sp,key_west$SP)
kw_trend_fam$grp<- key_west$Family[m]
kw_trend_comp_rvc<- par_comp(x=kw_trend_fam,levels=fam_sum_kw$Family,par='t.rvc',samps=2400)
kw_trend_comp_reef<- par_comp(x=kw_trend_fam,levels=fam_sum_kw$Family,par='t.reef',samps=2400)
#kw_trend_comp_comb<- par_comp(x=kw_trend_fam,levels=fam_sum_kw$Family,par='t.comb',samps=2400)
kw_trend_comp_diff<- par_comp(x=kw_trend_fam,levels=fam_sum_kl$Family,par='t.diff',samps=2400)
#diff<- data.frame(Family=kw_trend_comp_rvc$group,diff=kw_trend_comp_reef$par.m-kw_trend_comp_rvc$par.m)
#diff<- diff[order(diff$diff),]
#kw_trend_comp_comb<- kw_trend_comp_comb[order(kw_trend_comp_comb$par.m),]
#kw_trend_comp_reef<- kw_trend_comp_reef[rownames(kw_trend_comp_comb),]
#kw_trend_comp_rvc<- kw_trend_comp_rvc[rownames(kw_trend_comp_comb),]
#fam_sum_kw<- fam_sum_kw[rownames(kw_trend_comp_comb),]
#kw_trend_comp_comb$Family.c=fam_sum_kw$Family.com[match(kw_trend_comp_comb$group,fam_sum_kw$Family)]
#kw_trend_comp_reef$Family.c=fam_sum_kw$Family.com[match(kw_trend_comp_reef$group,fam_sum_kw$Family)]
#kw_trend_comp_rvc$Family.c=fam_sum_kw$Family.com[match(kw_trend_comp_rvc$group,fam_sum_kw$Family)]


trends_kl$family<- key_largo$Family[match(trends_kl$sp,key_largo$SP)]
trends_kw$family<- key_west$Family[match(trends_kw$sp,key_west$SP)]
trends_kl<- trends_kl[order(trends_kl$diff.trend.m),]
trends_kw<- trends_kw[order(trends_kw$diff.trend.m),]

par(mfrow=c(2,2))
plot(trends_kl$rvc.u.median~seq(1:nrow(trends_kl)),ylim=c(-0.2,0.2),xlim=c(0.75,nrow(trends_kl)+0.25),bty='l',xaxt='n',ylab='Population Trend (1993-2018)',type='n',xlab='',main='Key Largo')
abline(h=0,lty=5)
for(i in 1:nrow(trends_kl)){
  lines(c(trends_kl[i,5],trends_kl[i,6])~rep(i-0.15,2),col=adjustcolor('darkblue',alpha.f=0.4))
}
for(i in 1:nrow(trends_kl)){
  lines(c(trends_kl[i,10],trends_kl[i,11])~rep(i+0.15,2),col=adjustcolor('darkred',alpha.f=0.4))
}
points(trends_kl[,2]~c(seq(1:nrow(trends_kl))-0.15),cex=1,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(trends_kl[,7]~c(seq(1:nrow(trends_kl))+0.15),cex=1,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
text(par("usr")[1],par("usr")[4]*1.1,'A',cex=1,pos=3,xpd=T,font=2)

plot(trends_kw$rvc.u.median~seq(1:nrow(trends_kw)),ylim=c(-0.3,0.3),xlim=c(0.75,nrow(trends_kw)+0.25),bty='l',xaxt='n',ylab='Population Trend (1998-2018)',type='n',xlab='',main='Key West')
abline(h=0,lty=5)
for(i in 1:nrow(trends_kw)){
  lines(c(trends_kw[i,5],trends_kw[i,6])~rep(i-0.15,2),col=adjustcolor('darkblue',alpha.f=0.4))
}
for(i in 1:nrow(trends_kw)){
  lines(c(trends_kw[i,10],trends_kw[i,11])~rep(i+0.15,2),col=adjustcolor('darkred',alpha.f=0.4))
}
points(trends_kw[,2]~c(seq(1:nrow(trends_kw))-0.15),cex=1,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(trends_kw[,7]~c(seq(1:nrow(trends_kw))+0.15),cex=1,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
text(par("usr")[1],par("usr")[4]*1.1,'B',cex=1,pos=3,xpd=T,font=2)

plot(kl_trend_comp_rvc$par.m~seq(1:nrow(kl_trend_comp_rvc)),ylim=c(-0.1,0.1),xlim=c(0.3,12.5),bty='l',xaxt='n',ylab='Mean Population Trend',type='n',xlab='',main='')
abline(h=0,lty=5)
for(i in 1:nrow(kl_trend_comp_rvc)){
  lines(c(kl_trend_comp_rvc[i,5],kl_trend_comp_rvc[i,6])~rep(i-0.15,2),col=adjustcolor('darkblue',alpha.f=0.4))
}
for(i in 1:nrow(kl_trend_comp_reef)){
  lines(c(kl_trend_comp_reef[i,5],kl_trend_comp_reef[i,6])~rep(i+0.15,2),col=adjustcolor('darkred',alpha.f=0.4))
}
points(kl_trend_comp_rvc[,2]~c(seq(1:12)-0.15),cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kl_trend_comp_reef[,2]~c(seq(1:12)+0.15),cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
text(c(seq(from=1,to=12,by=1)+0.1),-0.1,fam_sum_kl$n)
text(0.3,-0.1005,'n =')
text(c(seq(from=0.5,to=11.5,by=1)),par("usr")[3]*1.2,kl_trend_comp_rvc[,1],cex=0.8,srt = 45,pos=1,xpd=T)
text(par("usr")[1],par("usr")[4]*1.1,'C',cex=1,pos=3,xpd=T,font=2)
legend(10,0.1*1.23,c(as.expression(bquote(bold("RVC"))),as.expression(bquote(bold('REEF')))),text.col=c('darkblue','darkred'),bty='n',y.intersp=0.5,cex=0.8)

plot(kw_trend_comp_rvc$par.m~seq(1:nrow(kw_trend_comp_rvc)),ylim=c(-0.1,0.1),xlim=c(0.3,12.5),bty='l',xaxt='n',ylab='Mean Population Trend',type='n',xlab='',main='')
abline(h=0,lty=5)
for(i in 1:nrow(kw_trend_comp_rvc)){
  lines(c(kw_trend_comp_rvc[i,5],kw_trend_comp_rvc[i,6])~rep(i-0.15,2),col=adjustcolor('darkblue',alpha.f=0.4))
}
for(i in 1:nrow(kw_trend_comp_reef)){
  lines(c(kw_trend_comp_reef[i,5],kw_trend_comp_reef[i,6])~rep(i+0.15,2),col=adjustcolor('darkred',alpha.f=0.4))
}
points(kw_trend_comp_rvc[,2]~c(seq(1:12)-0.15),cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kw_trend_comp_reef[,2]~c(seq(1:12)+0.15),cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
text(c(seq(from=1,to=12,by=1)+0.1),-0.1,fam_sum_kw$n)
text(0.3,-0.1005,'n =')
text(c(seq(from=0.5,to=11.5,by=1)),par("usr")[3]*1.2,kw_trend_comp_rvc[,1],cex=0.8,srt = 45,pos=1,xpd=T)
text(par("usr")[1],par("usr")[4]*1.1,'D',cex=1,pos=3,xpd=T,font=2)
legend(10,0.1*1.23,c(as.expression(bquote(bold("RVC"))),as.expression(bquote(bold('REEF')))),text.col=c('darkblue','darkred'),bty='n',y.intersp=0.5,cex=0.8)


#
par(mfrow=c(1,2))
plot(kl_trend_comp_rvc$par.m~seq(1:12),ylim=c(-0.1,0.1),xlim=c(0.3,12.5),bty='l',xaxt='n',ylab='Mean Annual Trend',type='n',xlab='',main='Key Largo')
abline(h=0,lty=5)
for(i in 1:nrow(kl_trend_comp_rvc)){
  lines(c(kl_trend_comp_rvc[i,5],kl_trend_comp_rvc[i,6])~rep(i-0.25,2),col=adjustcolor('darkblue',alpha.f=0.4))
}
for(i in 1:nrow(kl_trend_comp_reef)){
  lines(c(kl_trend_comp_reef[i,5],kl_trend_comp_reef[i,6])~rep(i+0.25,2),col=adjustcolor('darkred',alpha.f=0.4))
}
points(kl_trend_comp_rvc[,2]~c(seq(1:12)-0.25),cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
#points(kl_trend_comp_comb[,2]~c(seq(1:12)),cex=2,col='white',bg=adjustcolor('darkslateblue',alpha.f=0.7),pch=21)
points(kl_trend_comp_reef[,2]~c(seq(1:12)+0.25),cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
text(c(seq(from=1,to=12,by=1)+0.1),-0.1,fam_sum_kl$n)
text(0.3,-0.1005,'n =')
text(c(seq(from=0.5,to=11.5,by=1)),par("usr")[3] - 0.01,kl_trend_comp_rvc[,1],cex=0.8,srt = 45,pos=1,xpd=T)

plot(kw_trend_comp_rvc$par.m~seq(1:11),ylim=c(-0.1,0.1),xlim=c(0.3,11.5),bty='l',xaxt='n',ylab='Mean Annual Trend',type='n',xlab='',main='Key West')
abline(h=0,lty=5)
for(i in 1:nrow(kw_trend_comp_rvc)){
  lines(c(kw_trend_comp_rvc[i,5],kw_trend_comp_rvc[i,6])~rep(i-0.25,2),col=adjustcolor('darkblue',alpha.f=0.4))
}
for(i in 1:nrow(kw_trend_comp_reef)){
  lines(c(kw_trend_comp_reef[i,5],kw_trend_comp_reef[i,6])~rep(i+0.25,2),col=adjustcolor('darkred',alpha.f=0.4))
}
points(kw_trend_comp_rvc[,2]~c(seq(1:12)-0.25),cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
#points(kw_trend_comp_comb[,2]~c(seq(1:12)),cex=2,col='white',bg=adjustcolor('goldenrod',alpha.f=0.7),pch=21)
points(kw_trend_comp_reef[,2]~c(seq(1:12)+0.25),cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
text(c(seq(from=1,to=11,by=1)+0.1),-0.1,fam_sum_kw$n)
text(0.3,-0.1005,'n =')
text(c(seq(from=0.5,to=10.5,by=1)),par("usr")[3] - 0.01,kw_trend_comp_rvc[,1],cex=0.8,srt = 45,pos=1,xpd=T)


#Trends by Size ####
par(mfrow=c(1,2))
plot(kl_dat$rvc.u.median~log10(key_largo$size),ylim=c(-0.2,0.2),bty='l',xaxt='n',ylab='Mean Annual Trend',type='n',xlab='',data=sd_kl,main='Key Largo')
abline(h=0,lty=5)
for(i in 1:nrow(kl_dat)){
  lines(c(kl_dat$rvc.u.l95[i],kl_dat$rvc.u.u95[i])~rep(log10(key_largo$size)[i],2),col=adjustcolor('darkblue',alpha.f=0.4))
}
#for(i in 1:nrow(kl_dat)){
#  lines(c(kl_dat$comb.u.l95[i],kl_dat$comb.u.u95[i])~rep(log10(key_largo$size)[i],2),col=adjustcolor('darkslateblue',alpha.f=0.4))
#}
for(i in 1:nrow(kl_dat)){
  lines(c(kl_dat$reef.u.l95[i],kl_dat$reef.u.u95[i])~rep(log10(key_largo$size)[i],2),col=adjustcolor('darkred',alpha.f=0.4))
}
points(kl_dat$rvc.u.median~log10(key_largo$size),cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kl_dat$reef.u.median~log10(key_largo$size),cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
axis(1, col="black", at=seq(0,3,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(1),expression(10),expression(100),expression(1000)))
pow <- 0:3
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
mtext(expression(paste("Total Length (cm)")),side=1,line = 3,cex=1)

plot(kw_dat$rvc.u.median~log10(key_west$size),ylim=c(-0.2,0.2),bty='l',xaxt='n',ylab='Mean Annual Trend',type='n',xlab='',data=sd_kw,main='Key West')
abline(h=0,lty=5)
for(i in 1:nrow(kw_dat)){
  lines(c(kw_dat$rvc.u.l95[i],kw_dat$rvc.u.u95[i])~rep(log10(key_west$size)[i],2),col=adjustcolor('darkblue',alpha.f=0.4))
}
#for(i in 1:nrow(kw_dat)){
#  lines(c(kw_dat$comb.u.l95[i],kw_dat$comb.u.u95[i])~rep(log10(key_largo$size)[i],2),col=adjustcolor('darkslateblue',alpha.f=0.4))
#}
for(i in 1:nrow(kw_dat)){
  lines(c(kw_dat$reef.u.l95[i],kw_dat$reef.u.u95[i])~rep(log10(key_west$size)[i],2),col=adjustcolor('darkred',alpha.f=0.4))
}
points(kw_dat$rvc.u.median~log10(key_west$size),cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kw_dat$reef.u.median~log10(key_west$size),cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
axis(1, col="black", at=seq(0,3,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(1),expression(10),expression(100),expression(1000)))
pow <- 0:3
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
mtext(expression(paste("Total Length (cm)")),side=1,line = 3,cex=1)

#Trends by abundance####
par(mfrow=c(1,2))
plot(kl_dat$rvc.u.median~log10(key_largo$mean.x.rvc),ylim=c(-0.2,0.2),xlim=c(-2.2,1.5),bty='l',xaxt='n',ylab='Mean Annual Trend',type='n',xlab='',data=sd_kl,main='Key Largo')
abline(h=0,lty=5)
for(i in 1:nrow(kl_dat)){
  lines(c(kl_dat$rvc.u.l95[i],kl_dat$rvc.u.u95[i])~rep(log10(key_largo$mean.x.rvc)[i],2),col=adjustcolor('darkblue',alpha.f=0.4))
}

#for(i in 1:nrow(kl_dat)){
#  lines(c(kl_dat$comb.u.l95[i],kl_dat$comb.u.u95[i])~rep(log10(key_largo$size)[i],2),col=adjustcolor('darkslateblue',alpha.f=0.4))
#}
for(i in 1:nrow(kl_dat)){
  lines(c(kl_dat$reef.u.l95[i],kl_dat$reef.u.u95[i])~rep(log10(key_largo$mean.x.reef)[i],2),col=adjustcolor('darkred',alpha.f=0.4))
}
points(kl_dat$rvc.u.median~log10(key_largo$mean.x.rvc),cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kl_dat$reef.u.median~log10(key_largo$mean.x.reef),cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
axis(1, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -2:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
mtext(expression(paste("Mean Abundance")),side=1,line = 3,cex=1)

plot(kw_dat$rvc.u.median~log10(key_west$mean.x.rvc),ylim=c(-0.2,0.2),xlim=c(-2.2,1.5),bty='l',xaxt='n',ylab='Mean Annual Trend',type='n',xlab='',main='Key West')
abline(h=0,lty=5)
for(i in 1:nrow(kw_dat)){
  lines(c(kw_dat$rvc.u.l95[i],kw_dat$rvc.u.u95[i])~rep(log10(key_west$mean.x.rvc)[i],2),col=adjustcolor('darkblue',alpha.f=0.4))
}
#for(i in 1:nrow(kw_dat)){
#  lines(c(kw_dat$comb.u.l95[i],kw_dat$comb.u.u95[i])~rep(log10(key_largo$size)[i],2),col=adjustcolor('darkslateblue',alpha.f=0.4))
#}
for(i in 1:nrow(kw_dat)){
  lines(c(kw_dat$reef.u.l95[i],kw_dat$reef.u.u95[i])~rep(log10(key_west$mean.x.reef)[i],2),col=adjustcolor('darkred',alpha.f=0.4))
}
points(kw_dat$rvc.u.median~log10(key_west$mean.x.rvc),cex=2,col='white',bg=adjustcolor('darkblue',alpha.f=0.7),pch=21)
points(kw_dat$reef.u.median~log10(key_west$mean.x.reef),cex=2,col='white',bg=adjustcolor('darkred',alpha.f=0.7),pch=21)
axis(1, col="black", at=seq(-2,2,by=1),   tcl=-0.45, cex.axis=1.1,
     labels=c(expression(0.01),expression(0.1),expression(1),expression(10),expression(100)))
pow <- -2:2
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
mtext(expression(paste("Mean Abundance")),side=1,line = 3,cex=1)


###Variance decomposition####
#For fixed effects (bottom time, depth, etc.) - need variance of effects from data for conditional variance
load(here('data','data_runs.RData'))

files_kl_par_2<- list.files(here('outputs','species parameter estimates','Key Largo','model 2'))
files_kl_par_1<- list.files(here('outputs','species parameter estimates','Key Largo','model 1'))

var_frame_reef_kl<- data.frame(key_largo$SP,prop.site=NA,prop.hab=NA,prop.strat=NA,prop.diver=NA,prop.d.cluster=NA,prop.m.cluster=NA,prop.process=NA,prop.obserror=NA,prop.btime=NA,prop.depth=NA,prop.vis=NA,prop.current=NA,prop.expert=NA)
var_frame_rvc_kl<- data.frame(key_largo$SP,prop.psu=NA,prop.hab=NA,prop.strat=NA,prop.process=NA,prop.obserror=NA,prop.depth=NA)
for(i in 1:nrow(key_largo)){
  sp<- read.csv(here('outputs','species parameter estimates','Key Largo','model 2',files_kl_par_2[i]))
  reef_occs<- reef_occs_3403[[i]]
  rvc_occs<- rvc_occs_3403[[i]]
  betas1<-sp[,(gsub('\\..*','',colnames(sp))=='beta1')]
  betas2<-sp[,(gsub('\\..*','',colnames(sp))=='beta2')]
  X1<- matrix(data=c(scale(as.numeric(rvc_occs$DEPTH))),ncol=1,nrow=nrow(rvc_occs))
  X2<- matrix(data=c(scale(as.numeric(reef_occs$btime)),scale(as.numeric(reef_occs$averagedepth)),scale(as.numeric(reef_occs$visibility)),scale(as.numeric(reef_occs$current)),reef_occs$exp_binary),ncol=5,nrow=nrow(reef_occs))
  var_beta1=NA
  for(z in 1:length(betas1)){ 
    var_beta1[z]<- var(betas1[z]*X1) #variance attributed to fixed effects - Nakagawa & Schielzeth 20
  }
  var_beta2<- data.frame(betas2)
  for(q in 1:ncol(betas2)){
      for(z in 1:nrow(betas2)){ 
        var_beta2[z,q]<- var(betas2[z,q]*X2[,q]) #variance attributed to fixed effects - Nakagawa & Schielzeth 20
      }
  }
  var_full1<- data.frame(sp$sd_psu^2,sp$sd_hab1^2,sp$sd_strat1^2,sp$sd_q1^2,sp$sd_r1^2,var_beta1)
  var_full2<- data.frame(sp$sd_site^2,sp$sd_hab2^2,sp$sd_strat2^2,sp$sd_dv^2,sp$sd_dmy^2,sp$sd_my^2,sp$sd_q2^2,sp$sd_r2^2,var_beta2)

  var_frame_rvc_kl[i,2]=median(var_full1[,1])  
  var_frame_rvc_kl[i,3]=median(var_full1[,2])  
  var_frame_rvc_kl[i,4]=median(var_full1[,3])  
  var_frame_rvc_kl[i,5]=median(var_full1[,4])  
  var_frame_rvc_kl[i,6]=median(var_full1[,5])  
  var_frame_rvc_kl[i,7]=median(var_full1[,6]) 
  
  var_frame_reef_kl[i,2]=median(var_full2[,1])  
  var_frame_reef_kl[i,3]=median(var_full2[,2])  
  var_frame_reef_kl[i,4]=median(var_full2[,3])  
  var_frame_reef_kl[i,5]=median(var_full2[,4])  
  var_frame_reef_kl[i,6]=median(var_full2[,5])  
  var_frame_reef_kl[i,7]=median(var_full2[,6]) 
  var_frame_reef_kl[i,8]=median(var_full2[,7]) 
  var_frame_reef_kl[i,9]=median(var_full2[,8]) 
  var_frame_reef_kl[i,10]=median(var_full2[,9])
  var_frame_reef_kl[i,11]=median(var_full2[,10])
  var_frame_reef_kl[i,12]=median(var_full2[,11])
  var_frame_reef_kl[i,13]=median(var_full2[,12])
  var_frame_reef_kl[i,14]=median(var_full2[,13])

}
write.csv(var_frame_rvc_kl,here('outputs','variance decomp','species_variances_rvc_kl.csv'))
write.csv(var_frame_reef_kl,here('outputs','variance decomp','species_variances_reef_kl.csv'))


files_kw_par_2<- list.files(here('outputs','species parameter estimates','Key West','model 2'))

var_frame_reef_kw<- data.frame(key_west$SP,prop.site=NA,prop.hab=NA,prop.strat=NA,prop.diver=NA,prop.d.cluster=NA,prop.m.cluster=NA,prop.process=NA,prop.obserror=NA,prop.btime=NA,prop.depth=NA,prop.vis=NA,prop.current=NA,prop.expert=NA)
var_frame_rvc_kw<- data.frame(key_west$SP,prop.psu=NA,prop.hab=NA,prop.strat=NA,prop.process=NA,prop.obserror=NA,prop.depth=NA)
for(i in 1:nrow(key_west)){
  sp<- read.csv(here('outputs','species parameter estimates','Key West','model 2',files_kw_par_2[i]))
  reef_occs<- reef_occs_3408[[i]]
  rvc_occs<- rvc_occs_3408[[i]]
  betas1<-sp[,(gsub('\\..*','',colnames(sp))=='beta1')]
  betas2<-sp[,(gsub('\\..*','',colnames(sp))=='beta2')]
  X1<- matrix(data=c(scale(as.numeric(rvc_occs$DEPTH))),ncol=1,nrow=nrow(rvc_occs))
  X2<- matrix(data=c(scale(as.numeric(reef_occs$btime)),scale(as.numeric(reef_occs$averagedepth)),scale(as.numeric(reef_occs$visibility)),scale(as.numeric(reef_occs$current)),reef_occs$exp_binary),ncol=5,nrow=nrow(reef_occs))
  var_beta1=NA
  for(z in 1:length(betas1)){ 
    var_beta1[z]<- var(betas1[z]*X1) #variance attributed to fixed effects - Nakagawa & Schielzeth 20
  }
  var_beta2<- data.frame(betas2)
  for(q in 1:ncol(betas2)){
    for(z in 1:nrow(betas2)){ 
      var_beta2[z,q]<- var(betas2[z,q]*X2[,q]) #variance attributed to fixed effects - Nakagawa & Schielzeth 20
    }
  }
  var_full1<- data.frame(sp$sd_psu^2,sp$sd_hab1^2,sp$sd_strat1^2,sp$sd_q1^2,sp$sd_r1^2,var_beta1)
  var_full2<- data.frame(sp$sd_site^2,sp$sd_hab2^2,sp$sd_strat2^2,sp$sd_dv^2,sp$sd_dmy^2,sp$sd_my^2,sp$sd_q2^2,sp$sd_r2^2,var_beta2)
  
  var_frame_rvc_kw[i,2]=median(var_full1[,1])  
  var_frame_rvc_kw[i,3]=median(var_full1[,2])  
  var_frame_rvc_kw[i,4]=median(var_full1[,3])  
  var_frame_rvc_kw[i,5]=median(var_full1[,4])  
  var_frame_rvc_kw[i,6]=median(var_full1[,5])  
  var_frame_rvc_kw[i,7]=median(var_full1[,6]) 
  
  var_frame_reef_kw[i,2]=median(var_full2[,1])  
  var_frame_reef_kw[i,3]=median(var_full2[,2])  
  var_frame_reef_kw[i,4]=median(var_full2[,3])  
  var_frame_reef_kw[i,5]=median(var_full2[,4])  
  var_frame_reef_kw[i,6]=median(var_full2[,5])  
  var_frame_reef_kw[i,7]=median(var_full2[,6]) 
  var_frame_reef_kw[i,8]=median(var_full2[,7]) 
  var_frame_reef_kw[i,9]=median(var_full2[,8]) 
  var_frame_reef_kw[i,10]=median(var_full2[,9])
  var_frame_reef_kw[i,11]=median(var_full2[,10])
  var_frame_reef_kw[i,12]=median(var_full2[,11])
  var_frame_reef_kw[i,13]=median(var_full2[,12])
  var_frame_reef_kw[i,14]=median(var_full2[,13])
  
}
write.csv(var_frame_rvc_kw,here('outputs','variance decomp','species_variances_rvc_kw.csv'))
write.csv(var_frame_reef_kw,here('outputs','variance decomp','species_variances_reef_kw.csv'))


####