burn.in=2e3
nmc=2e4
maxdepth=3
rounding.digits=2
formula = chol ~ .
ntree=500
n=1e4
nb=2e2


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
source("ruleSHAP/R/Convergence.R")
source("ruleSHAP/R/ModelFitsMhorse.R")
source("ruleSHAP/R/VBFits.R")
source("ruleSHAP/R/RuleGen.R")
source("ruleSHAP/R/SHAP.R")
source("ruleSHAP/R/ruleSHAP-package.R")
Rcpp::sourceCpp('ruleSHAP/src/hs_VB.cpp')


## Load packages and data
fulldata=readRDS('helius.Rda')
print(names(fulldata))

#Compute number of predictors
p=ncol(fulldata)-2

#Select outcome
heliusdata=cbind(fulldata[,1:p],chol=fulldata$chol)

#Shuffle samples and split into trains and test
set.seed(1)
base_ind=sample(nrow(heliusdata))[1:n]
helius_train=heliusdata[base_ind,]

#Fit model
set.seed(42)
data = helius_train
VB_fit=ruleSHAP_VB(formula=formula, data=data,family='gaussian',verbose=T,burn.in = burn.in,nmc = nmc,
                   maxiter=2e4, maxiter_full=7e3, batch_size=nb,step_size=sqrt(nb/(n*log(3*n))),
                   step_size_sigma=sqrt(nb/(n*log(3*n))),
                   eps=1e-3,dvar=FALSE,invert_sigma=FALSE)

shapleys_df=compute_SHAP_VB(VB_fit,data[1:5e3,1:p],df=n,interactions=T)

saveRDS(shapleys_df,'output/shapleys_chol_longrun_n1e4.Rda')