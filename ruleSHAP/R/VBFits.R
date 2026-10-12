# Function to minimize to estimate scale
expscale2avg=function(x,avg){
  return(abs((x-1)/(x+1)/2/log(x)-avg))
}


#' Fit a RuleSHAP model using Variational Bayes
#'
#' This function fits a RuleSHAP model for continuous, binary and survival data.
#' It is based upon fitting the following tweaked horseshoe prior:
#' \deqn{y|X,a,b,\sigma^2 \sim \mathcal{N}(\beta_0+\sum_{j=1}^pb_jx_j+\sum_{k=1}^qa_kr_k,\sigma^2I),}
#' \deqn{b_j|\lambda^{(l)}_j,\tau,\sigma^2 \sim \mathcal{N}\Big(0,\big(\lambda^{(l)}_j\tau\tau_0\sigma\big)^2\Big),}
#' \deqn{a_k|\lambda^{(r)}_k,\tau,\sigma^2 \sim \mathcal{N}\Big(0,\big(\lambda^{(r)}_k\tau\tau_1\sigma\big)^2),}
#' \deqn{\lambda^{(l)}_j \sim \mathcal{C}^+(0,A^{(l)}),}
#' \deqn{\lambda^{(r)}_j \sim \mathcal{C}^+(0,A^{(r)}_j),}
#' \deqn{A^{(r)}_j=\frac{\big(2\cdot \min(1-\overline{r}_j,\overline{r}_j)\big)^\mu}{(m_j)^\eta},}
#' \deqn{\tau_0,\tau_1,\tau \sim \mathcal{C}^+(0,1),}
#' \deqn{\sigma^2 \sim \sigma^{-2}d\sigma^2,}
#' where \eqn{m_j} is the number of predictors involved in the rule \eqn{r_j}.
#'
#' @param formula an object of class "\link{formula}" which represents the
#' desired predicted outcome (Left-hand side) and the desired predictors
#' (Right-hand side). If intercept should not be included, do not actively
#' remove it from the formula. Instead, use the parameter \code{intercept}
#' @param data a data frame containing the data to be used to fit the model
#' @param intercept a logical entry determining whether an intercept should be
#' used (\code{TRUE}) or not used (\code{FALSE}). The former is recommended.
#' @param family A character string denoting the type of (Generalized) linear model
#' that is being fitted.
#' @param mu,eta hyparameters defining how strong the prior shrinkage of complex
#' rules should be. A higher value of \code{mu} induces stronger shrinkage on
#' rules with extreme support, a higher value of \code{eta} induces a stronger
#' shrinkage on the rules involving many predictors.
#' @param lin.relax the level of shrinkage of linear terms and intercept.
#' This is equal to the hyperparameter \eqn{A^{(l)}}.
#' @param method.tau.lin the type of prior distribution for the global shrinkage
#' \eqn{\tau_0} which is shared by the coefficients of the linear terms. Set to
#' \code{"halfCauchy"} by default, but may be set to \code{"fixed"} to have it
#' fixed to a value equal to \eqn{A^{(l)}} (see parameter \code{lin.relax}).
#' @param ntree number of trees that should be used when generating the rules
#' @param maxdepth maximum depth of the trees used to generate the rules. In particular,
#' this induces a maximum on the order of interactions between predictors.
#' @param disaggregate a logical value determining whether rules should be extracted
#' by also considering all possible subrules (\code{TRUE}) or if rules should only
#' be extracted in a top-down manner as in RuleFit and HorseRule (\code{FALSE}).
#' @param rounding.digits Number of decimals that the splitting points
#' are rounded to in the rule generation. This rounding simplifies readability
#' of rules and prevents rules from being very similar but not identical.
#' @param burn.in Number of burn-in MCMC iterations that are discarded before
#' the draws start to be recorded.
#' @param nmc Number of MCMC iterations after the initial burn-in period
#' @param thin Thinning parameter: if saving all MCMC draws is computationally
#' infeasible, thinning makes the model only record every \code{thin}-th iteration.
#' @param ping If the model is selected to be verbose, this denotes every how many
#' MCMC draws an update message should be printed
#' @param verbose logical value denoting whether the model should send printed
#' updates on how far the fitting procedure is
#' @param rules Either \code{NULL} or a character vector containing the rules that
#' should be used to fit the RuleSHAP model, if one wants to specify them in advance.
#' If left to its default \code{NULL}, rules will be generated using a
#' Parametric Random Forest.
#' @param wins.frac A numeric entry representing the fraction by which the dataset
#' is winsorized. It also applies to rules by inducing a minimum node size.
#'
#' @return rules a character vector containing the rules that were generated and used
#' for the model fitting
#' @return BetaSamples A matrix where every column is a MCMC sample of the coefficients.
#' The first row represents the intercept. The following \eqn{p} represent
#' the linear terms. The remaining rows are for the rules, in the same order as in
#' output \code{rules}
#' @return BetaHat A vector containing the posterior mean of the coefficients.
#' The first entry represents the intercept. The following \eqn{p} represent
#' the linear terms. The remaining entries are for the rules, in the same order as in
#' output \code{rules}
#' @export

ruleSHAP_VB=function(formula,data,intercept=TRUE,family=c('gaussian','binomial','survival'),
                  mu=0.5,eta=2,lin.relax=1,ntree=500,maxdepth=3,disaggregate=TRUE,rounding.digits=2,
                  maxiter=2e4,maxiter_full=1e4,eps=1e-3,burn.in=2e3,nmc=2e4,thin=1,
                  verbose=TRUE,rules=NULL,wins.frac=0.025,
                  batch_size=min(100,nrow(X)),step_size=sqrt(batch_size/(nrow(X)*log(3*nrow(X)))),decay=1-step_size,
                  step_size_sigma=-log(2*(1-1/nrow(X))^batch_size-(1-2/nrow(X))^batch_size)/log(5*nrow(X)*(nrow(X)-1)),decay_sigma=1-step_size_sigma,
                  ping=ceiling((10/step_size)/10^floor(log10((10/step_size))))*10^floor(log10((10/step_size))),
                  threshold.type=c("both","Mhorseshoe","fixed"),invert_sigma=FALSE,
                  threshold=0,threshold_update=1,adapt_p0=0,adapt_p1=0,dvar=FALSE, filter_r=TRUE){


  #Retrieve basic quantities (Survival is a bit special)
  family=match.arg(family)
  threshold.type=match.arg(threshold.type) #line specific to VB implementation
  X_names=labels(terms(formula, data=data))
  if(family == 'survival'){
    surv_name=as.character(formula[[2]])
    y_name=surv_name[2:3]
  } else{
    y_name=as.character(formula[[2]])
  }
  
  #Winsorize data and store winsorization points
  wins.points=matrix(nrow=2,ncol=ncol(data))
  which_numeric=which(sapply(data,is.numeric) & (names(data) %in% X_names))
  for(i in which_numeric){
    wins.points[,i]=quantile(data[,i],probs=c(wins.frac,1-wins.frac))
    data[data[,i]<wins.points[1,i],i]=wins.points[1,i]
    data[data[,i]>wins.points[2,i],i]=wins.points[2,i]
  }
  
  
  #Subset to relevant variables only
  y=data[, y_name]
  X=data[, X_names]
  data=data[,c(y_name,X_names)]
  if(family=='survival'){
    censored = which(y[,2] == 0)
  }
  
  
  #Retrieve linear terms (in case some variables are factors) and re-scale them
  if(family=='survival'){
    nosurv_formula=as.formula(paste("~", deparse(formula[[3]])))
    X_lin=model.matrix(nosurv_formula,data=data[,!(colnames(data) %in% y_name)])
  } else{
    X_lin=model.matrix(formula,data=data)
  }
  #Remove the intercept if requested
  if(!intercept){X_lin=X_lin[,-1]}
  #Standardize X_lin (but store info before doing so)
  muX_lin=apply(X_lin,MARGIN=2,FUN=mean)
  sdX_lin=apply(X_lin,MARGIN=2,FUN=sd)
  #Keep into account intercept in scaling
  if(intercept){
    X_lin[,-1]=scale(X_lin[,-1])
  } else{
    X_lin=scale(X_lin)
  }
  
  #Now actually for categorical entries do contrast-coding standardized as if balanced
  #Obtain number of categories and dummies used
  nlevs=unlist(lapply(X,nlevels))
  if(intercept){nlevs=c(0,nlevs)}
  ndummies=pmax(1,nlevs-1)
  #Compute which columns each variable uses
  end_i=cumsum(pmax(ndummies,1))
  start_i=c(1,end_i[-length(end_i)]+1)
  #For each categorical predictor, overwrite previous stdzation
  for(i in 1:length(nlevs)){
    if(nlevs[i]==0){ #Skip continuous columns (including intercept)
      next
    }
    else{ #Re-scale otherwise (use X_lin on RHS instead of X_mat)
      X_lin[,start_i[i]:end_i[i]]=(nlevs[i]*X_lin[,start_i[i]:end_i[i]]-1)/sqrt(nlevs[i]-1)
      #Update muX and sdX_lin acccordingly
      muX_lin[i]=1/nlevs[i]
      sdX_lin[i]=sqrt(nlevs[i]-1)/nlevs[i]
    }
  }
  
  #Get p (notice that p potentially includes the intercept) and n
  p=ncol(X_lin)
  n=nrow(data)
  
  
  #If gaussian, standardize Y. If survival, standardize log(Y)
  if(family=='gaussian'){
    MuY=mean(y)
    SDy=sd(y)
    if(intercept){
      y=(y-MuY)/SDy
    }
  } else if(family =='survival'){
    logY=log(y[,1])
    MuY=mean(logY)
    SDy=sd(logY)
    if(intercept){
      logY=(logY-MuY)/SDy
      y[,1]=exp(logY)
    }
  }
  
  #Fit model with just linear terms (used to generate rules and give starting point)
  HS_fit=hs(y=y, X=X_lin,family=family,
            method.tau.lin='halfCauchy', tau.lin=1,
            method.tau.rules = 'fixed', tau.rules = 1,
            method.tau = 'fixed', tau = 1,
            burn=burn.in, nmc=nmc,
            method.sigma="Jeffreys")
  
  
  #Round numeric predictors and create the dataset for rule generation
    which_numeric=sapply(X,is.numeric)
    rulegen_data=X
    rulegen_data[,which_numeric]=round(X[,which_numeric],digits=rounding.digits)
    if(family=='gaussian'){
      rulegen_data$res=y-c(X_lin%*%HS_fit$BetaHat)
    } else if(family=='binomial'){
      #For now I'll define the working residuals
      phat=1/(1+exp(-c(X_lin%*%HS_fit$BetaHat)))
      rulegen_data$res=(y-phat)/(phat*(1-phat))
    } else if(family == 'survival'){
      rulegen_data$res=logY-c(X_lin%*%HS_fit$BetaHat)
      #For censored observations, use E[r|r>residual], under normal distribution
      lowlim=rulegen_data$res[censored]/sqrt(HS_fit$Sigma2Hat)
      #rulegen_data$res[censored]=sqrt(HS_fit$Sigma2Hat)*dnorm(lowlim)/(1-pnorm(lowlim))
      rulegen_data$res[censored]=sqrt(HS_fit$Sigma2Hat)*
        exp(dnorm(lowlim,log = TRUE)-pnorm(lowlim,log.p = TRUE, lower.tail = FALSE))
      
      #actually, try to overwrite with working residuals, deal with censoring same way
      yhat=exp(c(X_lin%*%HS_fit$BetaHat))
      lowlim=(logY[censored]-log(yhat[censored]))/sqrt(HS_fit$Sigma2Hat)
      rulegen_data$res=(y[,1]-yhat)/yhat#=y/yhat-1=exp(log(Y)-Xbeta)-1
      rulegen_data$res[censored]=exp(HS_fit$Sigma2Hat/2
                                     +pnorm(lowlim-sqrt(HS_fit$Sigma2Hat),log.p = TRUE,
                                            lower.tail = FALSE)
                                     -pnorm(lowlim,log.p = TRUE, lower.tail = FALSE)
      )-1
    }
  
  #Fit parametric random forest on residuals if rules not pre-specified
  if(is.null(rules)){
    #Use the linear terms to compute residuals
    if(verbose){print('Generating rules')}
    
    rules=PRF(x=rulegen_data[,X_names],y=rulegen_data$res,
              ntree=ntree,rounding.digits=rounding.digits,
              maxdepth=maxdepth,disaggregate=disaggregate)
    
    #Filter out duplicate rules
    if(filter_r){
    rules=filter_rules(rules,X,filtering_tol=0,
                       priority=1/rowSums(RuleMats(rules,X)$RulePredMat))
    }
  }
  
  
  #Extend X_lin with the rules and standardize.
  #Since intercept is controlled externally, you can't just take the linear terms
  #from get_modmat, cause those change according to whether formula contains
  #the intercept or not. So, to be safe, attach rules to X_lin manually
  X_rules=pre:::get_modmat(formula = formula, data = data,
                           rules = rules, type = "both",
                           x_names = X_names, winsfrac = 0,
                           normalize = F, y_names = y_name)$x
  X_rules=X_rules[,(ncol(X_rules)-length(rules)+1):ncol(X_rules)]
  
  #Store info of X before rescaling
  muX=c(muX_lin,apply(X_rules,MARGIN=2,FUN=mean))
  X_mat=cbind(X_lin,scale(X_rules,scale=F))




  ####### Start the fitting process

  #Compute structured penalization for rules
  penaliz=(2*pmin(muX[-(1:p)],1-muX[-(1:p)]))^mu/
    sqrt(2*pmax(muX[-(1:p)],1-muX[-(1:p)]))/
    (rowSums(RuleMats(rules,X)$RulePredMat)^eta)
  prior_lambda=c(rep(lin.relax,p),penaliz)





  #Use CV Ridge regression to guesstimate a good MCMC starting point
  Ridge_fit=glmnet::cv.glmnet(x=t(t(X_rules)*penaliz),alpha=0,
                      y=rulegen_data$res,standardize=F,
                      family='gaussian',
                      weights=HS_fit$OmegaHat)
                      
  
  Beta_init=c(HS_fit$BetaHat,Ridge_fit$glmnet.fit$beta[,Ridge_fit$index[1]]*penaliz)
  TauLin_init=HS_fit$TauLinHat
  Tau_init=1
  TauRules_init=c(sqrt(sqrt(1/(n*Ridge_fit$lambda[Ridge_fit$index[1]]))))
  Lambda_init=c(HS_fit$LambdaHat,rep(TauRules_init,ncol(X_rules)))
  if(family == 'binomial'){
  Sigma2_init=1
  Omega_scale_init=rep(sd(y),n)
  for(i in 1:n){
  Omega_scale_init[i]=log(stats::optimise(expscale2avg,avg=HS_fit$OmegaMedian[i],
  lower=1,upper=1e5)$minimum)
  }
  TauRules_init=TauRules_init/sqrt(sqrt(mean(HS_fit$OmegaMedian))) #!!! is sqrt once because tau instead of tau2 and once because you split contribution between lambda and tau.
  Lambda_init[-(1:p)]=Lambda_init[-(1:p)]/sqrt(sqrt(mean(HS_fit$OmegaMedian)))
  } else{
  Sigma2_init=HS_fit$Sigma2Hat
  Omega_scale_init=rep(1,n)
  }

  if(verbose){print("Fitting model")}
  
  #Specific to VB: survival requires you to split outcomes
  if(family=='survival'){
    delta <- y[,2]
    y <- y[,1]
  } else{
    delta=rep(1,n)
  }

  #First fit with thresholding

  HS_fit=hs_VB(y=y, delta=delta, X=X_mat, family=substring(family,1,1),
               p_lin = p, max_iter=maxiter, eps=eps*10, #!!!
               prior_lambda = prior_lambda, ping=ping,
               Beta_init=Beta_init,Sigma2=Sigma2_init,Lambda_init=Lambda_init,
               tau_lin=TauLin_init*Tau_init,tau_rules=TauRules_init*Tau_init,
               Omega_scale_init=Omega_scale_init, verbose=verbose,
               threshold_type=substring(threshold.type,1,1),batch_size=batch_size,
               threshold=threshold,threshold_update=threshold_update,
               adapt_p0=adapt_p0,adapt_p1=adapt_p1,invert_sigma=invert_sigma,update_tau=TRUE,
               dvar=dvar)
  CPO_old=HS_fit$CPO
  
  #Then finish off without thresholding
  HS_fit=hs_VB(y=y, delta=delta, X=X_mat, family=substring(family,1,1),
               p_lin = p, max_iter=maxiter_full, eps=eps,
               prior_lambda = prior_lambda, ping=ping,
               Beta_init=HS_fit$BetaHat,Sigma2=1/(HS_fit$Sigma2InvHat),Lambda_init=1/sqrt(HS_fit$Lambda2InvHat),
               tau_lin=1/sqrt(HS_fit$TauL2InvHat),tau_rules=1/sqrt(HS_fit$TauR2InvHat),
               Omega_scale_init=Omega_scale_init, verbose=verbose, #!!! should have the new scale, not the current one
               threshold_type='f',batch_size=batch_size,
               threshold=1e15,threshold_update=1e3+1,
               adapt_p0=adapt_p0,adapt_p1=adapt_p1,update_tau=FALSE,
               dvar=dvar)
               HS_fit$CPOBurnIn=CPO_old
  

  #For gaussian and survival case, Y or logY was standardized.
  #Keep that into account, also for intercept
  if(family=='gaussian' || family=='survival'){
    HS_fit$BetaPrec=HS_fit$BetaPrec/(SDy^2)
    HS_fit$BetaHat=HS_fit$BetaHat*SDy
    if(intercept){
      HS_fit$BetaHat[1]=HS_fit$BetaHat[1]+MuY
    }
  }
  
  #Define output
  output=HS_fit
  output$rules=rules
  output$wins.points=wins.points
  
  #X was also centered and partially scaled.
  #Keep that into account for intercept
  if(intercept){
    #Update coefficients
    output$BetaHat[2:p]=output$BetaHat[2:p]/sdX_lin[-1]
    #Update intercept
    output$BetaHat[1]=
      output$BetaHat[1]-muX[-1]%*%output$BetaHat[-1]
    #Update precision matrix  
    output$BetaPrec[2:p,]=output$BetaPrec[2:p,]*sdX_lin[-1]
    output$BetaPrec[,2:p]=t(t(output$BetaPrec[,2:p])*sdX_lin[-1])
  } else{
    #Update coefficients
    output$BetaHat[1:p]=output$BetaHat[1:p]/sdX_lin
    #Update precision matrix  
    output$BetaPrec[1:p,]=output$BetaPrec[1:p,]*sdX_lin
    output$BetaPrec[,1:p]=t(t(output$BetaPrec[,1:p])*sdX_lin)
  }
  
  #Return fit model
  return(output)
}















compute_SHAP_VB=function(model,data,df=nrow(data),intercept=T,alpha=0.05,interactions=F,
                      combine=NULL,block_size=5e3,WeightMatrices=NULL){
  #Winsorize data (store unwinsorized version too)
  data_unwins=data
  for(i in 1:ncol(data)){
    if(!any(is.na(model$wins.points[,i]))){
      data[data[,i]<model$wins.points[1,i],i]=model$wins.points[1,i]
      data[data[,i]>model$wins.points[2,i],i]=model$wins.points[2,i]
    }
  }


  #Compute basic quantities
  p=ncol(data)
  n_test=nrow(data)
  if(is.null(WeightMatrices)){
    rule_objects=RuleMats(model$rules,data_unwins)
    WeightMatrices=ShapleyMats(data=data,Rs=rule_objects$Rs,
                             id_mat=rule_objects$RulePredMat,
                             interactions=interactions)
  }
  
  if(intercept){
    CoeffMat=model$BetaHat[-1]
    PrecMat=model$BetaPrec[-1,-1,drop=F] #should be transposed for later but it does not matter, it's symmetric
  }else{
    CoeffMat=model$BetaHat
    PrecMat=model$BetaPrec
  }

  #If combinations are requested, add them to WeightMatrix of marginal
  #To do that, do left-matrix-multiplication by a sparse matrix.
  if(!is.null(combine)){
    #Compute the matrix that groups up weights
    combi_lengths=sapply(combine,length)
    combi_i=rep(1:(n_test*length(combine)),times=rep(combi_lengths,each=n_test))
    combi_j=sapply(combine,function(combination){
      rep((which(names(data) %in% combination)-1)*n_test,times=n_test)+
        rep(1:n_test,each=length(combination))
    })
    combi_mat=sparseMatrix(i=combi_i,j=combi_j,x=1,
                           dims=c(length(combine)*n_test,nrow(WeightMatrices$marginal)))

    #Left multiply and append to WeightMatrices
    WeightMatrices$marginal=rbind(WeightMatrices$marginal,
                                  combi_mat%*%WeightMatrices$marginal)
  }


  #First do main effects
  vals=c(as.matrix(WeightMatrices$marginal%*%CoeffMat))
  if(is.na(alpha)){
    half_CI=NA
  } else{
    half_CI=qt(1-alpha/2,df=df)*sqrt(rowSums(t(Matrix:::solve(PrecMat,Matrix:::t(WeightMatrices$marginal)))*WeightMatrices$marginal))
  }

  shapleys_df=data.frame(x=c(unlist(lapply(data_unwins,as.numeric)),
                             rep(NA,nrow(data)*length(combine))),
                         value=vals,
                         CIinf=vals-half_CI,
                         CIsup=vals+half_CI,
                         predictor=rep(c(names(data),
                                         unlist(sapply(combine,paste,collapse=' & '))),
                                       each=nrow(data)))

  #Now do interactions, if requested
  shapleys_inter_df=NULL
  if(interactions){
    vals=c(as.matrix(WeightMatrices$interactions%*%CoeffMat))
    if(is.na(alpha)){
      half_CI=NA
    } else{
      half_CI=qt(1-alpha/2,df=df)*sqrt(rowSums(t(solve(PrecMat,Matrix:::t(WeightMatrices$interactions)))*WeightMatrices$interactions))
    }
    shapleys_inter_df=data.frame(x1=c(unlist(lapply(data_unwins,as.numeric))),
                                 x2=c(unlist(lapply(data_unwins,function(x){rep(as.numeric(x),times=p)}))),
                                 value=vals,
                                 CIinf=vals-half_CI,
                                 CIsup=vals+half_CI,
                                 predictor1=rep(names(data),each=n_test),
                                 predictor2=rep(names(data),each=n_test*p))
  }

  return(list(marginal=shapleys_df,interactions=shapleys_inter_df))
}
