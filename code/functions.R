`%notin%`<- Negate(`%in%`)

abund_tranfs<- function(probs){
  sum<- probs[1]*0+probs[2]*1+probs[3]*2+probs[4]*11+probs[5]*101
  return(sum)
}

mean_ord_to_n<- function(x,c){
  ord_yr=list()
  abund_x=matrix(ncol=ncol(x),nrow=nrow(c))
  if(ncol(c)==2){
    for(i in 1:ncol(x)){
      ord_yr[[i]]<- matrix(ncol=5,nrow=nrow(c))
      ord_yr[[i]][,1]<- plogis(c[,1]-x[,i])
      ord_yr[[i]][,2]<-plogis(c[,2]-x[,i])-plogis(c[,1]-x[,i])
      ord_yr[[i]][,3]<-1-plogis(c[,2]-x[,i])
      ord_yr[[i]][,4]<- 0
      ord_yr[[i]][,5]<- 0
      abund_x[,i]<- apply(ord_yr[[i]][,1:5],1,abund_tranfs)
    }
  }
  if(ncol(c)==3){
    for(i in 1:ncol(x)){
      ord_yr[[i]]<- matrix(ncol=5,nrow=nrow(c))
      ord_yr[[i]][,1]<- plogis(c[,1]-x[,i])
      ord_yr[[i]][,2]<-plogis(c[,2]-x[,i])-plogis(c[,1]-x[,i])
      ord_yr[[i]][,3]<-plogis(c[,3]-x[,i])-plogis(c[,2]-x[,i])
      ord_yr[[i]][,4]<- 1-plogis(c[,3]-x[,i])
      ord_yr[[i]][,5]<- 0
      abund_x[,i]<- apply(ord_yr[[i]][,1:5],1,abund_tranfs)
    }
  }
  if(ncol(c)==4){
    for(i in 1:ncol(x)){
      ord_yr[[i]]<- matrix(ncol=5,nrow=nrow(c))
      ord_yr[[i]][,1]<- plogis(c[,1]-x[,i])
      ord_yr[[i]][,2]<-plogis(c[,2]-x[,i])-plogis(c[,1]-x[,i])
      ord_yr[[i]][,3]<-plogis(c[,3]-x[,i])-plogis(c[,2]-x[,i])
      ord_yr[[i]][,4]<- plogis(c[,4]-x[,i])-plogis(c[,3]-x[,i])
      ord_yr[[i]][,5]<- 1-plogis(c[,4]-x[,i])
      abund_x[,i]<- apply(ord_yr[[i]][,1:5],1,abund_tranfs)
    }
  }
  return(abund_x)
}


ord_to_n<- function(x,c){
  abund_x<- numeric(length(x))
  p= matrix(ncol=5,nrow=length(x))
  if(ncol(c)==2){
    p[,1]=plogis(c[,1]-x)
    p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
    p[,3]=1-plogis(c[,2]-x)
    p[,4]=0
    p[,5]=0
    for(i in 1:length(x)){
      abund_x[i]=abund_tranfs(p[i,])  
    }
  }
  if(ncol(c)==3){
    p[,1]=plogis(c[,1]-x)
    p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
    p[,3]=plogis(c[,3]-x)-plogis(c[,2]-x)
    p[,4]=1-plogis(c[,3]-x)
    p[,5]=0
    for(i in 1:length(x)){
      abund_x[i]=abund_tranfs(p[i,])  
    }
  }
  if(ncol(c)==4){
    p[,1]=plogis(c[,1]-x)
    p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
    p[,3]=plogis(c[,3]-x)-plogis(c[,2]-x)
    p[,4]=plogis(c[,4]-x)-plogis(c[,3]-x)
    p[,5]=1-plogis(c[,4]-x)
    for(i in 1:length(x)){
      abund_x[i]=abund_tranfs(p[i,])  
    }
  }
  return(abund_x)
}


rvc_filter = function(x,GZ,sp,sites,year.start, year.end,strat.omit,hab.omit){ #Function to filter through RVC data
  x= x %>% subset(YEAR>=year.start) %>% subset(YEAR<=year.end) %>% subset(LAT_LON %in% sites$lat_lon) #Only keep designated years
  x$SSU_YEAR= paste(x$YEAR,x$PRIMARY_SAMPLE_UNIT,x$STATION_NR,sep='_')
  x1= x %>% subset(region.id %in% GZ) %>% select(SSU_YEAR,SPECIES_CD,everything())
  x2= complete(x1,SSU_YEAR,nesting(SPECIES_CD),fill=list(NUM=0)) #ensures all non-sightings are recorded
  zeros = anti_join(x2,x1)
  zeros[,3:29]= x1[match(zeros$SSU_YEAR,x1$SSU_YEAR),3:29]
  x3= rbind(x1,zeros)
  x3$HAB_CD2= gsub('\\_.*','',x3$HABITAT_CD)
  year_index= data.frame(yr=seq(year.start,year.end),y.ind=seq(1,year.end-year.start+1))
  x3$year_index=x3$y.ind[match(x3$YEAR,year_index$yr)]
  
  rvc_occs<- list()
  for(i in 1:nrow(sp)){
    x4= subset(x3,SPECIES_CD==sp$rvc_code[i])
    x5=  x4 %>% dplyr::group_by(SSU_YEAR) %>%
      dplyr::summarise(NUM.total=sum(NUM),occ=NA) %>% #Sums up the number of counts per SSU
      mutate(occ=ifelse(NUM.total>0,1,0)) %>% arrange(SSU_YEAR) #Also scores presence/absence at the SSU level
    x5[,4:32]<- x4[match(x5$SSU_YEAR,x4$SSU_YEAR),2:30]
    x5<- transform(x5,psu_id=match(LAT_LON,unique(LAT_LON)))
    x5$NUM.total2<- ceiling(x5$NUM.total)
    x5$STRAT=sites$stratum[match(x5$LAT_LON,sites$lat_lon)]
    x5$PSU_YEAR= paste(x5$YEAR,x5$PRIMARY_SAMPLE_UNIT,sep='_')
    x6<- subset(x5,STRAT %notin% strat.omit) %>% subset (HAB_CD2 %notin% hab.omit)
    
    rvc_occs[[i]]=x6
  }
  return(rvc_occs)
}

ts_rvc = function(x,miss){ #Takes the output from the previous function
  ts<- list()
  
  for(i in 1:length(x)){
    x1 = x[[i]]
    x2= x1 %>% group_by(YEAR) %>% summarize(n.occ=sum(occ),n.surv=n(),p.occ=n.occ/n.surv,sp=unique(SPECIES_CD))
    x3= x1 %>% group_by(YEAR,PRIMARY_SAMPLE_UNIT) %>% summarize(psu_abund=mean(NUM.total2)) %>% group_by(YEAR) %>% summarize(mean_abund=mean(psu_abund),sd_abund=sd(psu_abund))
    x4=left_join(x2,x3) %>% complete(YEAR=seq(min(x3$YEAR),max(x3$YEAR)))
    
    if(miss=='T')
      for(z in 1:nrow(x4)){
        if(is.na(x4$p.occ[z])==T){next}
        if(x4$p.occ[z]==0){
          x4$p.occ[z]=NA
          x4$mean_abund[z]=NA
        }
      }
    
    ts[[i]]=x4
  }
  return(ts)
}

reef_filter = function(R,GZ,sp,geog,year.start,year.end,strat.omit,hab.omit){ #function to trim dataframe and add in occurrence data by species
  occ_list<- list()
  TempDat<-R %>% subset(site4 %in% GZ) %>% subset(year>=year.start) %>% subset(year<=year.end) %>% select('formid','speciesid','abundance',everything())
  TempDat$hab_class<- geog$hab_class[match(TempDat$geogr,geog$geogid)]
  TempDat$stratum<- geog$stratum[match(TempDat$geogr,geog$geogid)]
  TempDat$lat_dd<- geog$lat_dd[match(TempDat$geogr,geog$geogid)]
  TempDat$lon_dd<- geog$lon_dd[match(TempDat$geogr,geog$geogid)]
  Zeros<- complete(TempDat,formid,nesting(speciesid),
                   fill=list(abundance=0)) %>%
    anti_join(TempDat) #Creates a dataset where each species is represented in every survey (adding in zero observations)
  m<- match(Zeros$formid,TempDat$formid) #Match the zero data frame to the rest of the data
  Zeros[,4:ncol(Zeros)]<- TempDat[m,4:ncol(TempDat)] #Replicate the survey-level data (except abundance)
  TempDat2<- rbind(TempDat,Zeros) %>% arrange(formid,speciesid) #Combine presence and absence data for all species
  TempDat2$occ=ifelse(TempDat2$abundance>0,1,0) #Code presence/absence based on abundance
  TempDat2<- TempDat2 %>% mutate(exp_binary=ifelse(exp=='E',1,0)) #Binary classification for expert or novice diver status
  TempDat2$abundance2<- TempDat2$abundance+1 #For modelling abundance categories in the ordinal regression, has to start at 1
  TempDat2$abund_trans <-NA #Transform abundance categories into the minimum counts
  TempDat2<- TempDat2%>%
    mutate(
      abund_trans=ifelse(abundance==0,0,abund_trans),
      abund_trans=ifelse(abundance==1,1,abund_trans),
      abund_trans=ifelse(abundance==2,2,abund_trans),
      abund_trans=ifelse(abundance==3,11,abund_trans),
      abund_trans=ifelse(abundance==4,101,abund_trans)
    )
  surveyors<- TempDat2 %>% group_by(fish_memberid) %>% summarize(n=length(unique(formid))) #Determine survey effort of each diver in the region
  surveyors_trim<- subset(surveyors,n>=5) #only keep surveys by members with 5 or more dives
  
  TempDat3<- subset(TempDat2, fish_memberid %in% surveyors_trim$fish_memberid) #Trim out surveys from divers with less than 5 surveys in this region 
  site_subsets<-TempDat3 %>% group_by(geogr) %>% summarize(n=n_distinct(formid),hab_class=unique(hab_class)) %>% subset(n>=5 & is.na(hab_class)==F) #Calculate surveys per site
  
  TempDat4<- TempDat3 %>% subset(geogr %in% site_subsets$geogr) #Only keep well surveyed sites (n>=5)
  TempDat4$hab_class2<- gsub('\\_.*','',TempDat4$hab_class)
  survs<- TempDat4 %>% distinct(formid,.keep_all = T) %>% group_by(geogr,day,month,year) %>% summarize(n=n()) %>% mutate(site_dmy=paste(geogr,day,month,year,sep='')) %>% subset(n>1) #Determine spatiotemporal clustering
  TempDat4<- TempDat4 %>% mutate(site_dmy=ifelse(paste(geogr,day,month,year,sep='') %in% survs$site_dmy,paste(geogr,day,month,year,sep=''),'baseline')) #For surveys where it was the only one that day on a site - put it in the global 'baseline' pool
  TempDat4$mth_cluster<- paste(TempDat4$year,TempDat4$month,sep='_')
  TempDat5<- subset(TempDat4, stratum %notin% strat.omit) %>% subset (hab_class2 %notin% hab.omit)
  
  for(i in 1:nrow(sp)){
    occ_list[[i]]<- subset(TempDat5,speciesid==sp$speciesid[i]) #Subset out each species in the provided dataframe
  }
  return(occ_list)
}

ts_reef = function(X,sp){
  ts_list<- list()
  for(i in 1:length(X)){
    sp_x<-X[[i]] #Subset out each species in the provided dataframe
    occ_by_year<- sp_x %>% group_by(year) %>% summarize(n.occ=sum(occ),n.surv=n(),p.occ=n.occ/n.surv,sp=sp$commonname[match(unique(sp_x$speciesid),sp$speciesid)])
    abun_by_year<- sp_x %>% group_by(year,geogr) %>% summarize(site_abun=mean(abund_trans),sd_abund=sd(abund_trans),n.surv=n()) %>% group_by(year) %>% summarize(mean_abund=mean(site_abun),n.survs=sum(n.surv),n.sites=n())
    total_sd<- sp_x %>% group_by(year) %>% summarize(sd=sd(abund_trans))
    
    comb<- left_join(occ_by_year,abun_by_year)
    comb2<- left_join(comb,total_sd)
    ts_list[[i]]<- comb2
  }
  return(ts_list)
}



###Time-series plot function###
scaled_timeseries_plot<- function(i,ts1,ts2,sp,GZ,mod,params1,params2,path,TT,TT.rvc,n.iter,yr.start,yr.end){
  pdf(file.path(path,paste(paste(i,sp,mod,GZ,sep='_'),'.pdf',sep='')),width=8,height=6)
  
  #Extract parameters
  cuts1=as.data.frame(params1)[grepl('cut',colnames(params1))]
  cuts2=as.data.frame(params2)[grepl('cut',colnames(params2))]
  x=as.data.frame(params1)[grepl('x',colnames(params1))]; x=x[,-1] #remove x0
  x1=as.data.frame(params2)[grepl('x1',colnames(params2))]; 
  x2=as.data.frame(params2)[grepl('x2',colnames(params2))];
  y1.1=as.data.frame(params1)[grepl('a_yr1',colnames(params1))]; 
  y1.2=as.data.frame(params1)[grepl('a_yr2',colnames(params1))]; 
  y2.1=as.data.frame(params2)[grepl('a_yr1',colnames(params2))]; 
  y2.2=as.data.frame(params2)[grepl('a_yr2',colnames(params2))]; 
  
  lambda_mat<- list()  
  for(i in 1:TT){
    reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,n.iter))
    
    if(ncol(cuts1)==2){
      if(mod=='model1'){
        reef_coef[,1]<- plogis(cuts1[,1]-y1.2[,i])
        reef_coef[,2]<-plogis(cuts1[,2]-y1.2[,i])-plogis(cuts1[,1]-y1.2[,i])
        reef_coef[,3]<-1-plogis(cuts1[,2]-y1.2[,i])
        reef_coef[,4]<- 0
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(cuts1[,1]-(x[,i]+params1$a))
        reef_coef[,7]<-plogis(cuts1[,2]-(x[,i]+params1$a))-plogis(cuts1[,1]-(x[,i]+params1$a))
        reef_coef[,8]<-1-plogis(cuts1[,2]-(x[,i]+params1$a))
        reef_coef[,9]<- 0
        reef_coef[,10]<- 0
      }
      if(mod=='model2'){
        reef_coef[,1]<- plogis(cuts2[,1]-y2.2[,i])
        reef_coef[,2]<-plogis(cuts2[,2]-y2.2[,i])-plogis(cuts2[,1]-y2.2[,i])
        reef_coef[,3]<-1-plogis(cuts2[,2]-y2.2[,i])
        reef_coef[,4]<- 0
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(cuts2[,1]-(x2[,i]))
        reef_coef[,7]<-plogis(cuts2[,2]-(x2[,i]))-plogis(cuts2[,1]-(x2[,i]))
        reef_coef[,8]<-1-plogis(cuts2[,2]-(x2[,i]))
        reef_coef[,9]<- 0
        reef_coef[,10]<- 0
      }
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(cuts1)==3){
      if(mod=='model1'){
        reef_coef[,1]<- plogis(cuts1[,1]-y1.2[,i])
        reef_coef[,2]<-plogis(cuts1[,2]-y1.2[,i])-plogis(cuts1[,1]-y1.2[,i])
        reef_coef[,3]<-plogis(cuts1[,3]-y1.2[,i])-plogis(cuts1[,2]-y1.2[,i])
        reef_coef[,4]<- 1-plogis(cuts1[,3]-y1.2[,i])
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(cuts1[,1]-(x[,i]+params1$a))
        reef_coef[,7]<-plogis(cuts1[,2]-(x[,i]+params1$a))-plogis(cuts1[,1]-(x[,i]+params1$a))
        reef_coef[,8]<-plogis(cuts1[,3]-(x[,i]+params1$a))-plogis(cuts1[,2]-(x[,i]+params1$a))
        reef_coef[,9]<- 1-plogis(cuts1[,3]-(x[,i]+params1$a))
        reef_coef[,10]<- 0
      }
      if(mod=='model2'){
        reef_coef[,1]<- plogis(cuts2[,1]-y2.2[,i])
        reef_coef[,2]<-plogis(cuts2[,2]-y2.2[,i])-plogis(cuts2[,1]-y2.2[,i])
        reef_coef[,3]<-plogis(cuts2[,3]-y2.2[,i])-plogis(cuts2[,2]-y2.2[,i])
        reef_coef[,4]<- 1-plogis(cuts2[,3]-y2.2[,i])
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(cuts2[,1]-(x2[,i]))
        reef_coef[,7]<-plogis(cuts2[,2]-(x2[,i]))-plogis(cuts2[,1]-(x2[,i]))
        reef_coef[,8]<-plogis(cuts2[,3]-(x2[,i]))-plogis(cuts2[,2]-(x2[,i]))
        reef_coef[,9]<- 1-plogis(cuts2[,3]-(x2[,i]))
        reef_coef[,10]<- 0
      }
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(params1[grepl('cut',colnames(params1))])==4){
      if(mod=='model1'){
        reef_coef[,1]<- plogis(cuts1[,1]-y1.2[,i])
        reef_coef[,2]<-plogis(cuts1[,2]-y1.2[,i])-plogis(cuts1[,1]-y1.2[,i])
        reef_coef[,3]<-plogis(cuts1[,3]-y1.2[,i])-plogis(cuts1[,2]-y1.2[,i])
        reef_coef[,4]<- plogis(cuts1[,4]-y1.2[,i])-plogis(cuts1[,3]-y1.2[,i])
        reef_coef[,5]<- 1- plogis(cuts1[,4]-y1.2[,i])
        reef_coef[,6]=plogis(cuts1[,1]-(x[,i]+params1$a))
        reef_coef[,7]<-plogis(cuts1[,2]-(x[,i]+params1$a))-plogis(cuts1[,1]-(x[,i]+params1$a))
        reef_coef[,8]<-plogis(cuts1[,3]-(x[,i]+params1$a))-plogis(cuts1[,2]-(x[,i]+params1$a))
        reef_coef[,9]<- plogis(cuts1[,4]-(x[,i]+params1$a))-plogis(cuts1[,3]-(x[,i]+params1$a))
        reef_coef[,10]<- 1-plogis(cuts1[,4]-(x[,i]+params1$a))
      }
      if(mod=='model2'){
        reef_coef[,1]<- plogis(cuts2[,1]-y2.2[,i])
        reef_coef[,2]<-plogis(cuts2[,2]-y2.2[,i])-plogis(cuts2[,1]-y2.2[,i])
        reef_coef[,3]<-plogis(cuts2[,3]-y2.2[,i])-plogis(cuts2[,2]-y2.2[,i])
        reef_coef[,4]<- plogis(cuts2[,4]-y2.2[,i])-plogis(cuts2[,3]-y2.2[,i])
        reef_coef[,5]<- 1- plogis(cuts2[,4]-y2.2[,i])
        reef_coef[,6]<-plogis(cuts2[,1]-(x2[,i]))
        reef_coef[,7]<-plogis(cuts2[,2]-(x2[,i]))-plogis(cuts2[,1]-(x2[,i]))
        reef_coef[,8]<-plogis(cuts2[,3]-(x2[,i]))-plogis(cuts2[,2]-(x2[,i]))
        reef_coef[,9]<- plogis(cuts2[,4]-(x2[,i]))-plogis(cuts2[,3]-(x2[,i]))
        reef_coef[,10]<- 1-plogis(cuts2[,4]-(x2[,i]))
      }
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    lambda_mat[[i]]=reef_coef
    
  }  
  median.x=median(do.call(rbind, lapply(lambda_mat, data.frame, stringsAsFactors=FALSE))$lambda.x)
  median.y=median(do.call(rbind, lapply(lambda_mat, data.frame, stringsAsFactors=FALSE))$lambda.y)
  
  
  x_mat<- data.frame(median.rvc=NA,l.95=NA,u.95=NA,median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  for(i in 1:TT){
    if(mod=='model1'){
      x_mat[i,1]=median(exp(x[,i])/median(as.vector(t(exp(x)))))
      x_mat[i,2]=quantile(exp(x[,i])/median(as.vector(t(exp(x)))),0.05)
      x_mat[i,3]=quantile(exp(x[,i])/median(as.vector(t(exp(x)))),0.95)
    }
    if(mod=='model2'){
      x_mat[i,1]=median(exp(x1[,i])/median(as.vector(t(exp(x1)))))
      x_mat[i,2]=quantile(exp(x1[,i])/median(as.vector(t(exp(x1)))),0.05)
      x_mat[i,3]=quantile(exp(x1[,i])/median(as.vector(t(exp(x1)))),0.95)
    }
    x_mat[i,4]=median(lambda_mat[[i]]$lambda.x/median.x)
    x_mat[i,5]=quantile(lambda_mat[[i]]$lambda.x/median.x,0.05)
    x_mat[i,6]=quantile(lambda_mat[[i]]$lambda.x/median.x,0.95)
  }
  
  y_mat_rvc<- data.frame(year=ts1$YEAR[complete.cases(ts1)],median.rvc=NA,l.95.rvc=NA,u.95.rvc=NA)
  y_mat_reef<- data.frame(year=seq(yr.start,yr.end),median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  for(i in 1:TT.rvc){
    if(mod=='model1'){
      y_mat_rvc[i,2]=median(exp(y1.1[,i])/median(as.vector(t(exp(y1.1)))))
      y_mat_rvc[i,3]=quantile(exp(y1.1[,i])/median(as.vector(t(exp(y1.1)))),0.05)
      y_mat_rvc[i,4]=quantile(exp(y1.1[,i])/median(as.vector(t(exp(y1.1)))),0.95)
    }
    if(mod=='model2'){
      y_mat_rvc[i,2]=median(exp(y2.1[,i])/median(as.vector(t(exp(y2.1)))))
      y_mat_rvc[i,3]=quantile(exp(y2.1[,i])/median(as.vector(t(exp(y2.1)))),0.05)
      y_mat_rvc[i,4]=quantile(exp(y2.1[,i])/median(as.vector(t(exp(y1.1)))),0.95)
    }
  }
  
  for(i in 1:TT){
    y_mat_reef[i,2]=median(lambda_mat[[i]]$lambda.y/median.y)
    y_mat_reef[i,3]=quantile(lambda_mat[[i]]$lambda.y/median.y,0.05)
    y_mat_reef[i,4]=quantile(lambda_mat[[i]]$lambda.y/median.y,0.95)
  }
  y_mat<- full_join(y_mat_reef,y_mat_rvc)
  
  par(xpd=T)
  plot(y_mat$median.rvc~y_mat$year,type='n',ylim=c(min(x_mat),max(c(max(y_mat[,2]),max(x_mat)))),col='darkblue',bty='l',ylab=expression('Scaled Relative Abundance'),xlab='Year',main=paste(sp,GZ,sep=' - '))
 
  x<- c(y_mat$year, rev(y_mat$year))
  y1<- c(x_mat[,2], rev(x_mat[,3]))
  y2<-  c(x_mat[,5], rev(x_mat[,6]))
  polygon(x, y1, col = adjustcolor('dodgerblue4', alpha = 0.2), border=NA) # Add uncertainty polygon
  polygon(x, y2, col = adjustcolor('firebrick4', alpha = 0.2), border=NA) # Add uncertainty polygon
  lines(x_mat[,1]~y_mat$year,lty=5,lwd=2,col='dodgerblue3')
  lines(x_mat[,4]~y_mat$year,lty=5,lwd=2,col='firebrick')
  
  lines(y_mat$median.rvc~y_mat$year,col='white',lwd=3)
  lines(y_mat$median.rvc~y_mat$year,col='dodgerblue4',lwd=2)
  points(y_mat$median.rvc~y_mat$year,col='white',pch=21,bg='dodgerblue4',cex=1.5,lwd=1)
  
  lines(y_mat$median.reef~y_mat$year,col='white',lwd=3)
  lines(y_mat$median.reef~c(seq(yr.start,yr.end)),col='firebrick4',lwd=2)
  points(y_mat$median.reef~c(seq(yr.start,yr.end)),col='white',pch=21,bg='firebrick4',cex=1.5,lwd=1)
  
  legend(yr.end-5,c(max(c(max(y_mat[,2]),max(x_mat)))*1.1),c('RVC surveys','REEF surveys'),text.col=c('dodgerblue4','firebrick4'),bty='n',y.intersp=1.4)
  dev.off()
}

timeseries_plot<- function(i,ts1,ts2,sp,GZ,mod,params1,params2,path,TT,TT.rvc,n.iter,yr.start,yr.end){
 # pdf(file.path(path,paste(paste(i,sp,mod,GZ,sep='_'),'.pdf',sep='')),width=8,height=6)
  
  lambda_mat<- list()  
  for(i in 1:TT){
    reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,n.iter))
    
    if(ncol(params_1[grepl('cut',colnames(params_1))])==2){
      if(mod=='model1'){
        reef_coef[,1]<- plogis(params1$cut[,1]-params1$a_yr2[,i])
        reef_coef[,2]<-plogis(params1$cut[,2]-params1$a_yr2[,i])-plogis(params1$cut[,1]-params1$a_yr2[,i])
        reef_coef[,3]<-1-plogis(params1$cut[,2]-params1$a_yr2[,i])
        reef_coef[,4]<- 0
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$cut[,1]-(params1$x[,i]+params1$a))
        reef_coef[,7]<-plogis(params1$cut[,2]-(params1$x[,i]+params1$a))-plogis(params1$cut[,1]-(params1$x[,i]+params1$a))
        reef_coef[,8]<-1-plogis(params1$cut[,2]-(params1$x[,i]+params1$a))
        reef_coef[,9]<- 0
        reef_coef[,10]<- 0
      }
      if(mod=='model2'){
        reef_coef[,1]<- plogis(params2$cut[,1]-params2$a_yr2[,i])
        reef_coef[,2]<-plogis(params2$cut[,2]-params2$a_yr2[,i])-plogis(params2$cut[,1]-params2$a_yr2[,i])
        reef_coef[,3]<-1-plogis(params2$cut[,2]-params2$a_yr2[,i])
        reef_coef[,4]<- 0
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params2$cut[,1]-params2$x2[,i])
        reef_coef[,7]<-plogis(params2$cut[,2]-params2$x2[,i])-plogis(params2$cut[,1]-params2$x2[,i])
        reef_coef[,8]<-1-plogis(params2$cut[,2]-params2$x2[,i])
        reef_coef[,9]<- 0
        reef_coef[,10]<- 0
      }
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(params_1[grepl('cut',colnames(params_1))])==3){
      if(mod=='model1'){
        reef_coef[,1]<- plogis(params1$cut[,1]-params1$a_yr2[,i])
        reef_coef[,2]<-plogis(params1$cut[,2]-params1$a_yr2[,i])-plogis(params1$cut[,1]-params1$a_yr2[,i])
        reef_coef[,3]<-plogis(params1$cut[,3]-params1$a_yr2[,i])-plogis(params1$cut[,2]-params1$a_yr2[,i])
        reef_coef[,4]<- 1-plogis(params1$cut[,3]-params1$a_yr2[,i])
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$cut[,1]-(params1$x[,i]+params1$a))
        reef_coef[,7]<-plogis(params1$cut[,2]-(params1$x[,i]+params1$a))-plogis(params1$cut[,1]-(params1$x[,i]+params1$a))
        reef_coef[,8]<-plogis(params1$cut[,3]-(params1$x[,i]+params1$a))-plogis(params1$cut[,2]-(params1$x[,i]+params1$a))
        reef_coef[,9]<- 1-plogis(params1$cut[,3]-(params1$x[,i]+params1$a))
        reef_coef[,10]<- 0
      }
      if(mod=='model2'){
        reef_coef[,1]<- plogis(params2$cut[,1]-params2$a_yr2[,i])
        reef_coef[,2]<-plogis(params2$cut[,2]-params2$a_yr2[,i])-plogis(params2$cut[,1]-params2$a_yr2[,i])
        reef_coef[,3]<-plogis(params2$cut[,3]-params2$a_yr2[,i])-plogis(params2$cut[,2]-params2$a_yr2[,i])
        reef_coef[,4]<- 1-plogis(params2$cut[,3]-params2$a_yr2[,i])
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params2$cut[,1]-params2$x2[,i])
        reef_coef[,7]<-plogis(params2$cut[,2]-params2$x2[,i])-plogis(params2$cut[,1]-params2$x2[,i])
        reef_coef[,8]<-plogis(params2$cut[,3]-params2$x2[,i])-plogis(params2$cut[,2]-params2$x2[,i])
        reef_coef[,9]<- 1-plogis(params2$cut[,3]-params2$x2[,i])
        reef_coef[,10]<- 0
      }
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(params_1[grepl('cut',colnames(params_1))])==4){
      if(mod=='model1'){
        reef_coef[,1]<- plogis(params1$cut[,1]-params1$a_yr2[,i])
        reef_coef[,2]<-plogis(params1$cut[,2]-params1$a_yr2[,i])-plogis(params1$cut[,1]-params1$a_yr2[,i])
        reef_coef[,3]<-plogis(params1$cut[,3]-params1$a_yr2[,i])-plogis(params1$cut[,2]-params1$a_yr2[,i])
        reef_coef[,4]<- plogis(params1$cut[,4]-params1$a_yr2[,i])-plogis(params1$cut[,3]-params1$a_yr2[,i])
        reef_coef[,5]<- 1- plogis(params1$cut[,4]-params1$a_yr2[,i])
        reef_coef[,6]=plogis(params1$cut[,1]-(params1$x[,i]+params1$a))
        reef_coef[,7]<-plogis(params1$cut[,2]-(params1$x[,i]+params1$a))-plogis(params1$cut[,1]-(params1$x[,i]+params1$a))
        reef_coef[,8]<-plogis(params1$cut[,3]-(params1$x[,i]+params1$a))-plogis(params1$cut[,2]-(params1$x[,i]+params1$a))
        reef_coef[,9]<- plogis(params1$cut[,4]-(params1$x[,i]+params1$a))-plogis(params1$cut[,3]-(params1$x[,i]+params1$a))
        reef_coef[,10]<- 1-plogis(params1$cut[,4]-(params1$x[,i]+params1$a))
      }
      if(mod=='model2'){
        reef_coef[,1]<- plogis(params2$cut[,1]-params2$a_yr2[,i])
        reef_coef[,2]<-plogis(params2$cut[,2]-params2$a_yr2[,i])-plogis(params2$cut[,1]-params2$a_yr2[,i])
        reef_coef[,3]<-plogis(params2$cut[,3]-params2$a_yr2[,i])-plogis(params2$cut[,2]-params2$a_yr2[,i])
        reef_coef[,4]<- plogis(params2$cut[,4]-params2$a_yr2[,i])-plogis(params2$cut[,3]-params2$a_yr2[,i])
        reef_coef[,5]<- 1- plogis(params2$cut[,4]-params2$a_yr2[,i])
        reef_coef[,6]<-plogis(params2$cut[,1]-params2$x2[,i])
        reef_coef[,7]<-plogis(params2$cut[,2]-params2$x2[,i])-plogis(params2$cut[,1]-params2$x2[,i])
        reef_coef[,8]<-plogis(params2$cut[,3]-params2$x2[,i])-plogis(params2$cut[,2]-params2$x2[,i])
        reef_coef[,9]<- plogis(params2$cut[,4]-params2$x2[,i])-plogis(params2$cut[,3]-params2$x2[,i])
        reef_coef[,10]<- 1-plogis(params2$cut[,4]-params2$x2[,i])
      }
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    lambda_mat[[i]]=reef_coef
    
  }  
  
  x_mat<- data.frame(median.rvc=NA,l.95=NA,u.95=NA,median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  for(i in 1:TT){
    if(mod=='model1'){
      x_mat[i,1]=median(exp(params1$x[,i]))
      x_mat[i,2]=quantile(exp(params1$x[,i]),0.05)
      x_mat[i,3]=quantile(exp(params1$x[,i]),0.95)
    }
    if(mod=='model2'){
      x_mat[i,1]=median(exp(params2$x1[,i]))
      x_mat[i,2]=quantile(exp(params2$x1[,i]),0.05)
      x_mat[i,3]=quantile(exp(params2$x1[,i]),0.95)
    }
    x_mat[i,4]=median(lambda_mat[[i]]$lambda.x)
    x_mat[i,5]=quantile(lambda_mat[[i]]$lambda.x,0.05)
    x_mat[i,6]=quantile(lambda_mat[[i]]$lambda.x,0.95)
  }
  
  y_mat_rvc<- data.frame(year=ts1$YEAR[complete.cases(ts1)],median.rvc=NA,l.95.rvc=NA,u.95.rvc=NA)
  y_mat_reef<- data.frame(year=seq(yr.start,yr.end),median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  for(i in 1:TT.rvc){
    if(mod=='model1'){
      y_mat_rvc[i,2]=median(exp(params1$a_yr1[,i]))
      y_mat_rvc[i,3]=quantile(exp(params1$a_yr1[,i]),0.05)
      y_mat_rvc[i,4]=quantile(exp(params1$a_yr1[,i]),0.95)
    }
    if(mod=='model2'){
      y_mat_rvc[i,2]=median(exp(params2$a_yr1[,i]))
      y_mat_rvc[i,3]=quantile(exp(params2$a_yr1[,i]),0.05)
      y_mat_rvc[i,4]=quantile(exp(params2$a_yr1[,i]),0.95)
    }
  }
  
  for(i in 1:TT){
    y_mat_reef[i,2]=median(lambda_mat[[i]]$lambda.y)
    y_mat_reef[i,3]=quantile(lambda_mat[[i]]$lambda.y,0.05)
    y_mat_reef[i,4]=quantile(lambda_mat[[i]]$lambda.y,0.95)
  }
  y_mat<- full_join(y_mat_reef,y_mat_rvc)
  
  par(xpd=T)
  plot(y_mat$median.rvc~y_mat$year,type='n',ylim=c(min(x_mat),max(c(max(na.omit(y_mat[,2:7])),max(x_mat)))),col='darkblue',bty='l',ylab=expression('Scaled Relative Abundance'),xlab='Year',main=paste(sp,GZ,sep=' - '))
  
  x<- c(y_mat$year, rev(y_mat$year))
  y1<- c(x_mat[,2], rev(x_mat[,3]))
  y2<-  c(x_mat[,5], rev(x_mat[,6]))
  polygon(x, y1, col = adjustcolor('dodgerblue4', alpha = 0.2), border=NA) # Add uncertainty polygon
  polygon(x, y2, col = adjustcolor('firebrick4', alpha = 0.2), border=NA) # Add uncertainty polygon
  lines(x_mat[,1]~y_mat$year,lty=5,lwd=2,col='dodgerblue3')
  lines(x_mat[,4]~y_mat$year,lty=5,lwd=2,col='firebrick')
  
  lines(y_mat$median.rvc~y_mat$year,col='white',lwd=3)
  lines(y_mat$median.rvc~y_mat$year,col='dodgerblue4',lwd=2)
  points(y_mat$median.rvc~y_mat$year,col='white',pch=21,bg='dodgerblue4',cex=1.5,lwd=1)
  
  lines(y_mat$median.reef~y_mat$year,col='white',lwd=3)
  lines(y_mat$median.reef~c(seq(yr.start,yr.end)),col='firebrick4',lwd=2)
  points(y_mat$median.reef~c(seq(yr.start,yr.end)),col='white',pch=21,bg='firebrick4',cex=1.5,lwd=1)
  
  legend(yr.end-5,c(max(c(max(y_mat[,2]),max(x_mat)))*1.1),c('RVC surveys','REEF surveys'),text.col=c('dodgerblue4','firebrick4'),bty='n',y.intersp=1.4)
 # dev.off()
  
}


`%notin%`<- Negate(`%in%`)

gm_mean<- function(x){
  prod(x)^(1/length(x))
}
sdlog10<- function(x) {sqrt(log10(1+var(x)/(mean(x))^2))
}

##Posterior parameter aggregation function##
par_agg<- function(x,path,pars,group,ref){ #x = list of files, pars = parameter of interest, group = grouping variable, ref = reference data-frame for grouping variables
  par_list<- list()
  for(i in 1:length(x)){
    dat<- read.csv(here(path,x[i]))
    pars_i<- dat[,gsub('\\..*','',colnames(dat)) %in% c(pars)]
    par_list[[i]]<- data.frame(par=pars_i,sp=rep(ref$SP[i],length(pars_i)),grp=rep(ref[i,colnames(ref)==group],length(pars_i)))
  }
  par_frame<- do.call(rbind, lapply(par_list, data.frame, stringsAsFactors=FALSE))
  return(par_frame)
}

##Posterior group comparison function##
par_comp<- function(x,levels,par,samps){ # x = parameters in aggregate from above, levels = levels of group to examine, pars = parameter of interest, samps = number of samples
  x_sub<- list()
  par_comp_dat<- data.frame(group=levels,par.m=NA,par.l90=NA,par.u90=NA,par.l95=NA,par.u95=NA)
  for(i in 1:length(levels)){
    x_sub[[i]]=subset(x, grp %in% levels[i])
    s=matrix(ncol=length(unique(x_sub[[i]]$sp)),nrow=samps)
    for(z in 1:length(unique(x_sub[[i]]$sp))){
      sp=subset(x_sub[[i]], sp %in% unique(x_sub[[i]]$sp)[z])
      s[,z]= sp[sample(nrow(sp),samps,replace=F),c(par)]
    }
    
    par_comp_dat[i,2]=mean(apply(s,1,mean))
    par_comp_dat[i,3]=quantile(apply(s,1,mean),0.05)
    par_comp_dat[i,4]=quantile(apply(s,1,mean),0.95)
    par_comp_dat[i,5]=quantile(apply(s,1,mean),0.025)
    par_comp_dat[i,6]=quantile(apply(s,1,mean),0.975)
  }
  return(par_comp_dat)
}

par_diff<- function(x,pars,samps){
  s=matrix(ncol=length(unique(x$sp)),nrow=samps)
  for(z in 1:length(unique(x_sub[[i]]$sp))){
    sp=subset(x, sp %in% unique(x$sp)[z])
    s[,z]= sp[sample(nrow(sp),samps,replace=F),colnames(x) %in% pars]
  }
}

VerticalHist <- function(x, xscale = NULL, xwidth, hist,
                         fillCol = adjustcolor("darkgray",alpha.f=0.4), lineCol ="white",rev){
  ## x (required) is the x position to draw the histogram
  ## xscale (optional) is the "height" of the tallest bar (horizontally),
  ##   it has sensible default behavior
  ## xwidth (required) is the horizontal spacing between histograms
  ## hist (required) is an object of type "histogram"
  ##    (or a list / df with $breaks and $density)
  ## fillCol and lineCol... exactly what you think.
  binWidth <- hist$breaks[2] - hist$breaks[1]
  if (is.null(xscale)) xscale <- xwidth * 0.90 / max(hist$density)
  n <- length(hist$density)
  x.l <- rep(x, n)
  x.r <- x.l + hist$density * xscale
  y.b <- hist$breaks[1:n]
  y.t <- hist$breaks[2:(n + 1)]
  if(rev==F){
    rect(xleft = x.l, ybottom = y.b, xright = x.r, ytop = y.t,
         col = fillCol, border = lineCol)
  }
  if(rev==T){
    x.r <- x.l - hist$density * xscale
    rect(xleft = x.r, ybottom = y.b, xright = x.l, ytop = y.t,
         col = fillCol, border = lineCol)
  }
}

ppd_plot_out=function(y,y_rep,file.path,mod,rvc){
  if(rvc==1){
    pdf(file.path(path,paste(paste(i,sp,mod,GZ,'RVC',sep='_'),'.pdf',sep='')),width=8,height=6)
    bayesplot::ppc_dens_overlay(y,y_rep)
    dev.off()
  }else{
    pdf(file.path(path,paste(paste(i,sp,mod,GZ,'REEF',sep='_'),'.pdf',sep='')),width=8,height=6)
    bayesplot::ppc_dens_overlay(y,y_rep)
    dev.off()
  }

}

