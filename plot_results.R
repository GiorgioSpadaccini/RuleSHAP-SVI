library(ggplot2)
library(dplyr)
nrep=20

#MSE
MSEMat=readRDS('output/MSEMat.Rda')
df=data.frame(MSE=c(MSEMat),
              Model=c('BART','RuleSHAP-GS','RuleSHAP-MS','RuleSHAP-SVI','OLS'),
              n=rep(c(1e3,3e3,5e3,1e4),each=nrep*nrow(MSEMat)),
              rep=rep(1:nrep,each=nrow(MSEMat)))

ggplot(df[df$Model!='OLS',], aes(x=as.factor(n), y=MSE, fill=Model)) +
  geom_boxplot()+
  theme(legend.position = "top", 
        legend.box = "horizontal",
        legend.direction = "horizontal")+
  guides(fill = guide_legend(nrow = 1))+
  xlab('Sample size')+ylab('MSE')










#Computational time
TimeMat=readRDS('output/TimeMat.Rda')

rowMeans(TimeMat[,1:nrep])
rowMeans(TimeMat[,(nrep+1):(2*nrep)])
rowMeans(TimeMat[,(2*nrep+1):(3*nrep)])
rowMeans(TimeMat[,(3*nrep+1):(4*nrep)])
rowMeans(TimeMat[,(4*nrep+1):(5*nrep)])




##### Local Inference
library(dplyr)
#Show results
nrep=20
ns=c(1e3,3e3,5e3,1e4)
lowsignal_quantile=0.10
SmallSigCounts=LargeSigCounts=NoiseCounts=list(matrix(0,ncol=nrep,nrow=4),
                                               matrix(0,ncol=nrep,nrow=4),
                                               matrix(0,ncol=nrep,nrow=4),
                                               matrix(0,ncol=nrep,nrow=4),
                                               matrix(0,ncol=nrep,nrow=4))
rownames(LargeSigCounts[[1]])=rownames(LargeSigCounts[[2]])=
  rownames(LargeSigCounts[[3]])=rownames(LargeSigCounts[[4]])=
  rownames(LargeSigCounts[[5]])=rownames(SmallSigCounts[[1]])=
  rownames(SmallSigCounts[[2]])=
  rownames(SmallSigCounts[[3]])=rownames(SmallSigCounts[[4]])=
  rownames(SmallSigCounts[[5]])=rownames(NoiseCounts[[1]])=
  rownames(NoiseCounts[[2]])=rownames(NoiseCounts[[3]])=
  rownames(NoiseCounts[[4]])=rownames(NoiseCounts[[5]])=
  c('BART','RuleSHAP-GS','RuleSHAP-MS','RuleSHAP-SVI')

names(SmallSigCounts)=names(LargeSigCounts)=
  names(NoiseCounts)=paste0('n=',ns)
n_vec=rep(1:4,each=nrep)

for(i in 1:length(n_vec)){
  
  n_index=n_vec[i]
  n=1e3
  
  bart_df=bart_dfs[[i]]
  if(n_index<=3){
    #If GS not available, keep the old one,
    #makes calculations convenient
    GSruleshap_df=GSruleshap_dfs[[i]]
  }
  MSruleshap_df=MSruleshap_dfs[[i]]
  SVIruleshap_df=SVIruleshap_dfs[[i]]
  
  #Assess local significance
  bart_df$sig=sign(bart_df$CIinf*bart_df$CIsup)==1
  GSruleshap_df$sig=sign(GSruleshap_df$CIinf*GSruleshap_df$CIsup)==1
  MSruleshap_df$sig=sign(MSruleshap_df$CIinf*MSruleshap_df$CIsup)==1
  SVIruleshap_df$sig=sign(SVIruleshap_df$CIinf*SVIruleshap_df$CIsup)==1
  
  #Asses size of signal
  target=bart_df[,c('x','value','predictor')]
  target$value=0
  target$value[1:n]=
    (10*sapply(target$x[1:n],\(x){integrate(\(z){sin(x*z*pi)},0,1)$value})+
       10*sin(target$x[1:n]*target$x[n+(1:n)]*pi))/2
  target$value[1:n]=target$value[1:n]-mean(target$value[1:n])
  target$value[n+(1:n)]=
    (10*sapply(target$x[n+(1:n)],\(x){integrate(\(z){sin(x*z*pi)},0,1)$value})+
       10*sin(target$x[1:n]*target$x[n+(1:n)]*pi))/2
  target$value[n+(1:n)]=target$value[n+(1:n)]-mean(target$value[n+(1:n)])
  target$value[2*n+(1:n)]=20*(target$x[2*n+(1:n)]-0.5)^2
  target$value[2*n+(1:n)]=target$value[2*n+(1:n)]-mean(target$value[2*n+(1:n)])
  target$value[3*n+(1:n)]=10*(target$x[3*n+(1:n)]-0.5)
  target$value[4*n+(1:n)]=5*(target$x[4*n+(1:n)]-0.5)
  
  target = target %>% 
    group_by(predictor) %>%
    mutate(thresh = quantile(abs(value),lowsignal_quantile)) %>%
    ungroup() %>%
    as.data.frame()
  
  small_signal=which(abs(target$value[1:(5*n)])<target$thresh[1:(5*n)])
  
  #Create variable determining signal/noise
  bart_df$signal=2*(bart_df$predictor %in% paste0('x.',1:5))
  GSruleshap_df$signal=2*(GSruleshap_df$predictor %in% paste0('x.',1:5))
  SVIruleshap_df$signal=2*(SVIruleshap_df$predictor %in% paste0('x.',1:5))
  MSruleshap_df$signal=2*(MSruleshap_df$predictor %in% paste0('x.',1:5))
  
  bart_df$signal[small_signal]=GSruleshap_df$signal[small_signal]=
    MSruleshap_df$signal[small_signal]=SVIruleshap_df$signal[small_signal]=1
  
  bart_counts=bart_df %>% group_by(signal) %>% summarize(npoints=sum(sig))
  GSruleshap_counts=GSruleshap_df %>% group_by(signal) %>% summarize(npoints=sum(sig))
  MSruleshap_counts=MSruleshap_df %>% group_by(signal) %>% summarize(npoints=sum(sig))
  SVIruleshap_counts=SVIruleshap_df %>% group_by(signal) %>% summarize(npoints=sum(sig))
  
  
  SmallSigCounts[[n_index]][1,i%%nrep+1]=
    sum(bart_counts$npoints*(bart_counts$signal==1))
  SmallSigCounts[[n_index]][2,i%%nrep+1]=
    sum(GSruleshap_counts$npoints*(GSruleshap_counts$signal==1))
  SmallSigCounts[[n_index]][3,i%%nrep+1]=
    sum(MSruleshap_counts$npoints*(MSruleshap_counts$signal==1))
  SmallSigCounts[[n_index]][4,i%%nrep+1]=
    sum(SVIruleshap_counts$npoints*(SVIruleshap_counts$signal==1))
  
  LargeSigCounts[[n_index]][1,i%%nrep+1]=
    sum(bart_counts$npoints*(bart_counts$signal==2))
  LargeSigCounts[[n_index]][2,i%%nrep+1]=
    sum(GSruleshap_counts$npoints*(GSruleshap_counts$signal==2))
  LargeSigCounts[[n_index]][3,i%%nrep+1]=
    sum(MSruleshap_counts$npoints*(MSruleshap_counts$signal==2))
  LargeSigCounts[[n_index]][4,i%%nrep+1]=
    sum(SVIruleshap_counts$npoints*(SVIruleshap_counts$signal==2))
  
  
  NoiseCounts[[n_index]][1,i%%nrep+1]=
    sum(bart_counts$npoints*(bart_counts$signal==0))
  NoiseCounts[[n_index]][2,i%%nrep+1]=
    sum(GSruleshap_counts$npoints*(GSruleshap_counts$signal==0))
  NoiseCounts[[n_index]][3,i%%nrep+1]=
    sum(MSruleshap_counts$npoints*(MSruleshap_counts$signal==0))
  NoiseCounts[[n_index]][4,i%%nrep+1]=
    sum(SVIruleshap_counts$npoints*(SVIruleshap_counts$signal==0))
}


#Ridgeline density plots for local rejection rates
library(ggplot2)
library(ggridges)
lines_df=data.frame(npoints=c(t(NoiseCounts[[1]])/5000,
                              t(NoiseCounts[[2]])/5000,
                              t(NoiseCounts[[3]])/5000,
                              t(NoiseCounts[[4]])/5000,
                              t(SmallSigCounts[[1]])/(lowsignal_quantile*5000),
                              t(SmallSigCounts[[2]])/(lowsignal_quantile*5000),
                              t(SmallSigCounts[[3]])/(lowsignal_quantile*5000),
                              t(SmallSigCounts[[4]])/(lowsignal_quantile*5000),
                              t(LargeSigCounts[[1]])/((1-lowsignal_quantile)*5000),
                              t(LargeSigCounts[[2]])/((1-lowsignal_quantile)*5000),
                              t(LargeSigCounts[[3]])/((1-lowsignal_quantile)*5000),
                              t(LargeSigCounts[[4]])/((1-lowsignal_quantile)*5000)),
                    n=rep(ns,each=nrep*4),
                    model=rep(rownames(NoiseCounts[[1]]),each=nrep),
                    order=rep(c(4,1,2,3),each=nrep),
                    signal=rep(c('Noise','Small signal','Large signal'),each=nrep*4*4))
lines_df$signal=factor(lines_df$signal, levels=c('Noise','Small signal','Large signal'))
lines_df = lines_df %>%
  group_by(n,model,signal) %>%
  mutate(avg_npoints=mean(npoints,na.rm = T)) %>%
  ungroup() %>%
  as.data.frame()


#For the sake of plotting visualization, set any NoiseCount larger than 6% to 6%
NoiseCounts[[1]]=pmin(NoiseCounts[[1]],0.06*5000)
NoiseCounts[[2]]=pmin(NoiseCounts[[2]],0.06*5000)
NoiseCounts[[3]]=pmin(NoiseCounts[[3]],0.06*5000)
NoiseCounts[[4]]=pmin(NoiseCounts[[4]],0.06*5000)

#And now redefine npoints (this way avg_npoints is correct)
lines_df$npoints=c(t(NoiseCounts[[1]])/5000,
                   t(NoiseCounts[[2]])/5000,
                   t(NoiseCounts[[3]])/5000,
                   t(NoiseCounts[[4]])/5000,
                   t(SmallSigCounts[[1]])/(lowsignal_quantile*5000),
                   t(SmallSigCounts[[2]])/(lowsignal_quantile*5000),
                   t(SmallSigCounts[[3]])/(lowsignal_quantile*5000),
                   t(SmallSigCounts[[4]])/(lowsignal_quantile*5000),
                   t(LargeSigCounts[[1]])/((1-lowsignal_quantile)*5000),
                   t(LargeSigCounts[[2]])/((1-lowsignal_quantile)*5000),
                   t(LargeSigCounts[[3]])/((1-lowsignal_quantile)*5000),
                   t(LargeSigCounts[[4]])/((1-lowsignal_quantile)*5000))

models=c('BART','RuleSHAP-GS','RuleSHAP-MS','RuleSHAP-SVI')

ggplot(lines_df[lines_df$n<=5e3 | lines_df$model != 'RuleSHAP-GS',], aes(x = npoints,
                     y= model,
                     fill = model,
                     color = model,
                     group = model,
                     height=after_stat(density))) +
  geom_density_ridges(alpha=0.45,
                      stat='density',
                      trim=TRUE,
                      scale=2) +
  geom_point(aes(x = avg_npoints,
                 y = model,
                 color = model,
                 #shape = factor(model,levels=models)
                ),
             inherit.aes = FALSE)+
  facet_grid(as.factor(n)~signal,scales='free')+
  theme_bw()+
  theme(strip.text.x = element_text(size = 10))+
  xlim(c(0,NA))+
  labs(x = "Fraction of observations with significant feature effect",
       y = "Frequency across replicates (density)")+
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 1),
    shape = guide_legend(order = 1)
  )+
  theme(legend.position = "top", 
        legend.box = "horizontal",
        legend.direction = "horizontal",
        panel.spacing = unit(0.5, "lines"))








#Plot distance of local effects

#Load data and compute MSE
n_vec=c(1e3,3e3,5e3,1e4)
p=10
models=c('BART','RuleSHAP-GS','RuleSHAP-MS','RuleSHAP-SVI')
library(dplyr)
library(ggplot2)

#Re-scale distances by theoretical variance of contribution
x1x2_avg=rmutil::int2(function(x,z){sin(x*z*pi)},a=c(0,0),b=c(1,1))
theo_var=c(400*integrate(function(x){((x-.5)^2-1/12)^2},0,1)$value, #x.3
           100/12,25/12, #x.4,x.5
           rep(1,p-5), #noise
           100*rmutil::int2(function(x,z){(sin(x*z*pi)-x1x2_avg)^2},a=c(0,0),b=c(1,1))) #x.1x.2

#Create dataframe that computes MSE of reconstruction, per repetition
final_df=NULL
for(i in 1:length(SVIruleshap_dfs)){
  #Retrieve n
  n=1e3
  #Select dataset and use it to compute target
  SVIruleshap_df=SVIruleshap_dfs[[i]]
  #Extract x1*x2 for the interactions
  x1x2=SVIruleshap_df$x[1:n]*SVIruleshap_df$x[(n+1):(2*n)]
  #Delete x1 and x2 individually cause you don't need them
  SVIruleshap_df=rbind(SVIruleshap_df,
                      data.frame(x=x1x2,
                                 value=SVIruleshap_df$value[1:n]+SVIruleshap_df$value[(n+1):(2*n)],
                                 CIinf=NA, CIsup=NA,
                                 predictor = 'x.1 & x.2'))
  SVIruleshap_df=SVIruleshap_df[!(SVIruleshap_df$predictor %in% c('x.1','x.2')),]
  #Add target
  target=c(20*(SVIruleshap_df$x[1:n]-.5)^2-5/3,10*(SVIruleshap_df$x[(n+1):(2*n)]-.5),5*(SVIruleshap_df$x[(2*n+1):(3*n)]-.5),
           rep(0,(p-5)*n),10*sin(x1x2*pi)-
             10*rmutil::int2(function(x,z){sin(x*z*pi)},a=c(0,0),b=c(1,1)))
  
  #Load other dataframes (and add joint effect for x.1 and x.2)
  if(i<=3*nrep){
  GSruleshap_df=rbind(GSruleshap_dfs[[i]],
                      data.frame(x=x1x2,
                                 value=GSruleshap_df$value[1:n]+GSruleshap_df$value[(n+1):(2*n)],
                                 CIinf=NA, CIsup=NA,
                                 predictor = 'x.1 & x.2'))
  GSruleshap_df=GSruleshap_df[!(GSruleshap_df$predictor %in% c('x.1','x.2')),]
  }
  MSruleshap_df=rbind(MSruleshap_dfs[[i]],
                      data.frame(x=x1x2,
                                 value=MSruleshap_df$value[1:n]+MSruleshap_df$value[(n+1):(2*n)],
                                 CIinf=NA, CIsup=NA,
                                 predictor = 'x.1 & x.2'))
  MSruleshap_df=MSruleshap_df[!(MSruleshap_df$predictor %in% c('x.1','x.2')),]
  bart_df=rbind(bart_dfs[[i]],
                data.frame(x=x1x2,
                           value=bart_df$value[1:n]+bart_df$value[(n+1):(2*n)],
                           CIinf=NA, CIsup=NA,
                           predictor = 'x.1 & x.2'))
  bart_df=bart_df[!(bart_df$predictor %in% c('x.1','x.2')),]
  
  
  df=rbind(cbind(model='RuleSHAP-SVI',SVIruleshap_df),
           cbind(model='RuleSHAP-GS',GSruleshap_df),
           cbind(model='RuleSHAP-MS',MSruleshap_df),
           cbind(model='BART',bart_df))
  df$target=target #Had to remove distance from joint x1x2 effects cause I forgot that
  df$theovar=rep(theo_var,each=n)
  df$predictor[df$predictor %in% paste0('x.',6:p)]='Noise'
  
  #Compute MSE
  final_df=rbind(final_df, cbind(
    rep=i,n=n_vec[ceiling(i/nrep)],df %>% group_by(predictor,model) %>% summarize(MSE=mean((value-target)^2/theovar))))
}


#Plot it
ggplot(final_df[(final_df$n<=5e3 | final_df$model != 'RuleSHAP-GS') & final_df$predictor != 'x.1 & x.2',], aes(x=predictor, y=MSE, fill=as.factor(model))) +
  geom_boxplot(width=0.6, position = position_dodge(width = .75))+
  stat_summary(fun = median, geom = "point", aes(group = as.factor(model),shape = as.factor(model)),
               size = 3, position = position_dodge(width = .75))+
  facet_grid(as.factor(n)~.,scales='free')+
  theme(legend.position = "top", 
        legend.box = "horizontal",
        legend.direction = "horizontal")+
  guides(fill = guide_legend(nrow = 1))+
  xlab('Predictor')+ylab('Mean Squared Distance from target (rescaled)')+
  scale_fill_manual(values=c(colorspace::rainbow_hcl(length(models))),
                    breaks=models,name='Model')+
  scale_shape_manual(values=1:(length(models)),
                     breaks=models, name='Model')


#Plot CPO to assess if it stabilized
niter_VB=length(c(readRDS('output/VB_CPO_n4_i5.Rda')[-1],readRDS('output/VB_CPOBurnIn_n4_i5.Rda')[-1]))
CPO_df=rbind(data.frame(CPO=c(
  readRDS('output/GS_CPO_n1_i1.Rda')[1:1e4],
  readRDS('output/GS_CPO_n1_i2.Rda')[1:1e4],
  readRDS('output/GS_CPO_n1_i3.Rda')[1:1e4],
  readRDS('output/GS_CPO_n1_i4.Rda')[1:1e4],
  readRDS('output/GS_CPO_n1_i5.Rda')[1:1e4],
  readRDS('output/GS_CPO_n2_i1.Rda')[1:1e4],
  readRDS('output/GS_CPO_n2_i2.Rda')[1:1e4],
  readRDS('output/GS_CPO_n2_i3.Rda')[1:1e4],
  readRDS('output/GS_CPO_n2_i4.Rda')[1:1e4],
  readRDS('output/GS_CPO_n2_i5.Rda')[1:1e4],
  readRDS('output/GS_CPO_n3_i1.Rda')[1:1e4],
  readRDS('output/GS_CPO_n3_i2.Rda')[1:1e4],
  readRDS('output/GS_CPO_n3_i3.Rda')[1:1e4],
  readRDS('output/GS_CPO_n3_i4.Rda')[1:1e4],
  readRDS('output/GS_CPO_n3_i5.Rda')[1:1e4],
  readRDS('output/GS_CPO_n3_i1.Rda')[1:1e4],
  readRDS('output/GS_CPO_n3_i2.Rda')[1:1e4],
  readRDS('output/GS_CPO_n3_i3.Rda')[1:1e4],
  readRDS('output/GS_CPO_n3_i4.Rda')[1:1e4],
  readRDS('output/GS_CPO_n3_i5.Rda')[1:1e4],
  readRDS('output/MS_CPO_n1_i1.Rda')[1:1e4],
  readRDS('output/MS_CPO_n1_i2.Rda')[1:1e4],
  readRDS('output/MS_CPO_n1_i3.Rda')[1:1e4],
  readRDS('output/MS_CPO_n1_i4.Rda')[1:1e4],
  readRDS('output/MS_CPO_n1_i5.Rda')[1:1e4],
  readRDS('output/MS_CPO_n2_i1.Rda')[1:1e4],
  readRDS('output/MS_CPO_n2_i2.Rda')[1:1e4],
  readRDS('output/MS_CPO_n2_i3.Rda')[1:1e4],
  readRDS('output/MS_CPO_n2_i4.Rda')[1:1e4],
  readRDS('output/MS_CPO_n2_i5.Rda')[1:1e4],
  readRDS('output/MS_CPO_n3_i1.Rda')[1:1e4],
  readRDS('output/MS_CPO_n3_i2.Rda')[1:1e4],
  readRDS('output/MS_CPO_n3_i3.Rda')[1:1e4],
  readRDS('output/MS_CPO_n3_i4.Rda')[1:1e4],
  readRDS('output/MS_CPO_n3_i5.Rda')[1:1e4],
  readRDS('output/MS_CPO_n4_i1.Rda')[1:1e4],
  readRDS('output/MS_CPO_n4_i2.Rda')[1:1e4],
  readRDS('output/MS_CPO_n4_i3.Rda')[1:1e4],
  readRDS('output/MS_CPO_n4_i4.Rda')[1:1e4],
  readRDS('output/MS_CPO_n4_i5.Rda'))[1:1e4],
  iterations=1:1e4,
  model=rep(c('RuleSHAP-GS','RuleSHAP-MS'),each=5*4*1e4),
  replicate=rep(1:5,each=1e4),
  thresholding=NA,
  n=rep(ns,each=5*1e4)
),
data.frame(CPO=c(
  readRDS('output/VB_CPOBurnIn_n1_i1.Rda')[-1],
  readRDS('output/VB_CPO_n1_i1.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n1_i2.Rda')[-1],
  readRDS('output/VB_CPO_n1_i2.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n1_i3.Rda')[-1],
  readRDS('output/VB_CPO_n1_i3.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n1_i4.Rda')[-1],
  readRDS('output/VB_CPO_n1_i4.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n1_i5.Rda')[-1],
  readRDS('output/VB_CPO_n1_i5.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n2_i1.Rda')[-1],
  readRDS('output/VB_CPO_n2_i1.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n2_i2.Rda')[-1],
  readRDS('output/VB_CPO_n2_i2.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n2_i3.Rda')[-1],
  readRDS('output/VB_CPO_n2_i3.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n2_i4.Rda')[-1],
  readRDS('output/VB_CPO_n2_i4.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n2_i5.Rda')[-1],
  readRDS('output/VB_CPO_n2_i5.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n3_i1.Rda')[-1],
  readRDS('output/VB_CPO_n3_i1.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n3_i2.Rda')[-1],
  readRDS('output/VB_CPO_n3_i2.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n3_i3.Rda')[-1],
  readRDS('output/VB_CPO_n3_i3.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n3_i4.Rda')[-1],
  readRDS('output/VB_CPO_n3_i4.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n3_i5.Rda')[-1],
  readRDS('output/VB_CPO_n3_i5.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n4_i1.Rda')[-1],
  readRDS('output/VB_CPO_n4_i1.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n4_i2.Rda')[-1],
  readRDS('output/VB_CPO_n4_i2.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n4_i3.Rda')[-1],
  readRDS('output/VB_CPO_n4_i3.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n4_i4.Rda')[-1],
  readRDS('output/VB_CPO_n4_i4.Rda')[-1],
  readRDS('output/VB_CPOBurnIn_n4_i5.Rda')[-1],
  readRDS('output/VB_CPO_n4_i5.Rda')[-1]),
  iterations=1:niter_VB,
  model='RuleSHAP-SVI',
  replicate=rep(1:5,each=niter_VB),
  thresholding=c(rep(TRUE,niter_VB-length(readRDS('output/VB_CPO_n4_i5.Rda')[-1])),rep(FALSE,length(readRDS('output/VB_CPO_n4_i5.Rda')[-1]))),
  n=rep(ns,each=5*niter_VB)
))


#Plot it (example for one case)
i_select=CPO_df$model=='RuleSHAP-MS' &
  CPO_df$n==1e4
ggplot(CPO_df[i_select,], aes(x=iterations, y=CPO)) +
  geom_line()+
  facet_grid(as.factor(replicate)~.,scales='free')+
  theme(legend.position = "top", 
        legend.box = "horizontal",
        legend.direction = "horizontal")+
  guides(fill = guide_legend(nrow = 1))+
  xlab('Iterations')+ylab('CPO (up to additive constant)')+
  scale_fill_manual(values=c(colorspace::rainbow_hcl(length(models))),
                    breaks=models,name='Model')+
  scale_shape_manual(values=1:(length(models)),
                     breaks=models, name='Model')
