burn.in=2e3
nmc=1e4
n_vec=c(1e3,3e3,5e3,1e4,3e4)
formula = y ~ .
n_max=max(n_vec)
p=10
nrep=20



## Load packages
library("pre")
library('dplyr')
library("partykit")
library("randomForest")
library("pre")
library("ggplot2")
library("ggthemes")
library('treeshap')
library('xgboost')
library('horseshoe')
library('horseshoenlm')
library('ggforce')
library('Matrix')


suppressMessages(source("ruleSHAP/R/DataGen.R"))
source("ruleSHAP/R/ModelFits.R")
source("ruleSHAP/R/ModelFitsMhorse.R")
source("ruleSHAP/R/VBFits.R")
source("ruleSHAP/R/RuleGen.R")
source("ruleSHAP/R/SHAP.R")
source("ruleSHAP/R/ruleSHAP-package.R")
source("ruleSHAP/R/Convergence.R")
Rcpp::sourceCpp('ruleSHAP/src/hs_VB.cpp')


#Shuffle samples and split into trains and test
set.seed(1)
friedman_trains=gendata.friedman1(n=nrep*n_max,p=p,sd=10,rho=0.3+diag(.7,nrow=p,ncol=p))
friedman_test=gendata.friedman1(n=1e4,p=p,sd=10,rho=0.3+diag(.7,nrow=p,ncol=p))

GSruleshap_dfs=MSruleshap_dfs=VBruleshap_dfs=bart_dfs=ols_dfs=list()
GSruleshapMat=MSruleshapMat=VBruleshapMat=OLSMat=matrix(nrow=p,ncol=nrep*length(n_vec))
TimeMat=MSEMat=matrix(ncol = nrep*length(n_vec), nrow=5)
rownames(TimeMat)=c('Baseline','BART','Horseshoe-GS','Mhorseshoe','Horseshoe-SVI')
rownames(MSEMat)=c('BART','Horseshoe-GS','Mhorseshoe','Horseshoe-SVI','OLS')
coeffs_list=list()

#Fit models
set.seed(42)
for(i in 1:nrep){
  for(n_index in 1:length(n_vec)){
    
    #Set repetition-specific quantities
    set.seed(length(n_vec)*(i-1)+n_index)
    n=n_vec[n_index]
    data = friedman_trains[(i-1)*n_max+(1:n),]
    
    #Verbose
    print(paste('Fitting on sample size n =',n))
    
    #Before starting, fit a small VB fit just to obtain the rules to use
    small_fit=ruleSHAP_VB(formula=formula, data=data,family='gaussian',verbose=T,burn.in = burn.in,nmc = nmc,
                          maxiter=2, maxiter_full = 2, batch_size=40) #!!! Check which arguments should be imput here
                          
    rule_objects=RuleMats(small_fit$rules,data[,1:p])
    rule_objects_test=RuleMats(small_fit$rules,data[1:min(1e3,n),1:p])
    WeightMatrices=ShapleyMats(data=data[,1:p],data_test=data[1:min(1e3,n),1:p],
                               Rs=rule_objects$Rs,Rs_test=rule_objects_test$Rs,
                               id_mat=rule_objects$RulePredMat,
                               interactions=F)
    
    
    
    
    
    #Now use a second small fit to determine the baseline time for the rest of the algorithm (without large Horseshoe model)
    base_time=Sys.time()
    set.seed(length(n_vec)*(i-1)+n_index)
    small_fit=ruleSHAP_VB(formula=formula, data=data,family='gaussian',verbose=T,burn.in = burn.in,nmc = nmc,
                          maxiter=2, maxiter_full=2,rules=small_fit$rules) #!!! Check which arguments should be imput here
    TimeMat[1,nrep*(n_index-1)+i]=difftime(Sys.time(),base_time,units='hours')
    
    #Start off with BART
    base_time=Sys.time()
    set.seed(length(n_vec)*(i-1)+n_index)
    BART_fit <- dbarts::bart(data[,1:p], data$y, keeptrees = TRUE, verbose = FALSE)
    TimeMat[2,nrep*(n_index-1)+i]=difftime(Sys.time(),base_time,units='hours')
    bart_dfs[[(n_index-1)*nrep+i]]=bartSHAP(BART_fit,data[1:1e3,1:p])
    
    #Now RuleSHAP with GS
    if(n_index<=3){
    base_time=Sys.time()
    set.seed(length(n_vec)*(i-1)+n_index)
    GS_fit=ruleSHAP(formula=formula, data=data,family='gaussian',verbose=T,burn.in = burn.in,nmc = nmc,
                    rules=small_fit$rules)
    TimeMat[3,nrep*(n_index-1)+i]=difftime(Sys.time(),base_time,units='hours')
    GSruleshap_dfs[[(n_index-1)*nrep+i]]=compute_SHAP(GS_fit,data[1:1e3,1:p],WeightMatrices=WeightMatrices)$marginal
    }
    
    #Now RuleSHAP with GS, using Mhorseshoe
    base_time=Sys.time()
    set.seed(length(n_vec)*(i-1)+n_index)
    MS_fit=ruleSHAP_Mhorse(formula=formula, data=data,family='gaussian',verbose=T,burn.in = burn.in,nmc = nmc,
                           rules=small_fit$rules)
    TimeMat[4,nrep*(n_index-1)+i]=difftime(Sys.time(),base_time,units='hours')
    MSruleshap_dfs[[(n_index-1)*nrep+i]]=compute_SHAP(MS_fit,data[1:1e3,1:p],WeightMatrices=WeightMatrices)$marginal
    
    #print(TimeMat[4,nrep*(n_index-1)+i])
    
    #Now RuleSHAP with VB
    base_time=Sys.time()
    set.seed(length(n_vec)*(i-1)+n_index)
    VB_fit=ruleSHAP_VB(formula=formula, data=data,family='gaussian',verbose=T,burn.in = burn.in,nmc = nmc,
                       maxiter=5e3,maxiter_full=2e3,
                       rules=MS_fit$rules,batch_size=200) #!!! Check which arguments should be imput here
    TimeMat[5,nrep*(n_index-1)+i]=difftime(Sys.time(),base_time,units='hours')
    VBruleshap_dfs[[(n_index-1)*nrep+i]]=compute_SHAP_VB(VB_fit,data[1:1e3,1:p],WeightMatrices=WeightMatrices)$marginal
    
    #print(TimeMat[5,nrep*(n_index-1)+i])
    
    #Now OLS
    OLS_fit=lm(formula=formula,data=data)
    ols_dfs[[(n_index-1)*nrep+i]]=data.frame(x = c(as.matrix(data[,1:p])),
      val = c(t(coef(OLS_fit)[-1]*t(as.matrix(scale(data[,1:p],scale=F))))),
      CIinf = c(t((coef(OLS_fit)[-1]+summary(OLS_fit)$coefficients[-1,2]*qnorm(.025))*t(as.matrix(scale(data[,1:p],scale=F))))),
      CIsup = c(t((coef(OLS_fit)[-1]+summary(OLS_fit)$coefficients[-1,2]*qnorm(.975))*t(as.matrix(scale(data[,1:p],scale=F))))),
      predictor = rep(names(data[,1:p]),each=n))
      
      
      
    #Store coeficients  
    VBruleshapMat[,nrep*(n_index-1)+i]=VB_fit$BetaHat[2:(p+1)]
    GSruleshapMat[,nrep*(n_index-1)+i]=GS_fit$BetaHat[2:(p+1)]
    MSruleshapMat[,nrep*(n_index-1)+i]=MS_fit$BetaHat[2:(p+1)]
    OLSMat[,nrep*(n_index-1)+i]=coef(OLS_fit)[-1]
    
    #coeffs_list[[nrep*(n_index-1)+i]]=list(VB=VB_fit$BetaHat, GS=GS_fit$BetaHat, MS=MS_fit$BetaHat)
    
    
    #Compute test MSE
    Ytr=data[, as.character(formula[[2]])]
    Xtr=data[, colnames(data) != as.character(formula[[2]])]
    my_modmat=cbind(1,pre:::get_modmat(formula = formula, data = friedman_test,
                               rules = small_fit$rules,
                               type = "both", x_names = colnames(Xtr),
                               winsfrac = 0, normalize = F, y_names = as.character(formula[[2]]))$x)
    
    MSEMat[,(n_index-1)*nrep+i]=colMeans((cbind(
          colMeans(predict(BART_fit,newdata=friedman_test)),
          c(my_modmat %*% GS_fit$BetaHat),
          c(my_modmat %*% MS_fit$BetaHat),
          c(my_modmat %*% VB_fit$BetaHat),
          predict(OLS_fit, newdata = friedman_test))-friedman_test$y)^2)
          
    #print(1-MSEMat[,(n_index-1)*nrep+i]/var(friedman_test$y))
    
    #Save results
    if(i<=5){
      saveRDS(GS_fit$LogCPO,paste0('output/GS_CPO_n',n_index,'_i',i,'.Rda'))
      saveRDS(MS_fit$LogCPO,paste0('output/MS_CPO_n',n_index,'_i',i,'.Rda'))
      saveRDS(VB_fit$CPO,paste0('output/VB_CPO_n',n_index,'_i',i,'.Rda'))
      saveRDS(VB_fit$CPOBurnIn,paste0('output/VB_CPOBurnIn_n',n_index,'_i',i,'.Rda'))
    }
    saveRDS(bart_dfs,'output/bart_dfs.Rda')
    saveRDS(VBruleshap_dfs,'output/VBruleshap_dfs.Rda')
    saveRDS(GSruleshap_dfs,'output/GSruleshap_dfs.Rda')
    saveRDS(MSruleshap_dfs,'output/MSruleshap_dfs.Rda')
    saveRDS(ols_dfs,'output/ols_dfs.Rda')
    saveRDS(VBruleshapMat,'output/VBruleshapMat.Rda')
    saveRDS(GSruleshapMat,'output/GSruleshapMat.Rda')
    saveRDS(MSruleshapMat,'output/MSruleshapMat.Rda')
    saveRDS(OLSMat,'output/OLSMat.Rda')
    saveRDS(TimeMat,'output/TimeMat.Rda')
    saveRDS(MSEMat,'output/MSEMat.Rda')
    
    }
}
