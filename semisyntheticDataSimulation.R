library(ROCR)
library(glmnet)
library(parcor)
library(monomvn)
library(caret)                      ###Calculate PPV; tune alpha and lambda in elastic net#####
library(flare)                      ###Dantzig Selector######    
library(c060)
library(hdi)
library(ncvreg)


niter=64
metric_seq=c("pauc","rmse","tpr","ppv")
METHODS=c("lasso","henet","ridge","scad","stability","adalasso")   

load("tcgaExpression.RData")
TCGA=tcga=t(tcgaExpression)
N=nrow(tcga)    ############The number of samples in the TCGA dataset########


################################################
# select a scenario
################################################

scenarios <- read.table("semisyntheticDataScenarios.txt", header = TRUE,  stringsAsFactors = FALSE)
# Any scenario (row) in this data frame can be selected, or a new alternative scenario can be used
list2env(scenarios[1, ], .GlobalEnv)

if (s0B!=0 & pB==0 | s0B==0 & pB!=0) {
  stop()
}
if (s0B>0 & s0/s0B*pB>p) {
  stop()
}


##########################################



data_generation  <- function(n, p, s0, SNR, corDesign, pB, s0B) {

  ###############################################
  # subsample covariate data and generate outcome variable
  ###############################################
  
  genes=sample(1:ncol(TCGA),p,replace=FALSE)
  other.genes=setdiff(1:ncol(TCGA),genes)
  tcga=TCGA[,genes]
  XR=tcga
  
  if(s0B==0)
  {signal.positions=sample(1:p,size=s0,replace=FALSE)}
  
  if(s0B>0)
  {
    nblock=s0/s0B                  #######How many blocks have signals?#####
    XR=tcga 
    xx=XR 
    cor.xx=cor(xx);   
    
    num.large.cor=function(x) {sum(x>0.5)}
    rest_variables=1:p
    non_blocked_variables=1:p
    
    
    for(block.iter in 1:nblock)
    {
      cor.xx.tmp=abs(cor.xx); diag(cor.xx.tmp)=0
      #sub.all.potential.group=which(apply(cor.xx.tmp,1,max)>0.5)
      #sub.core=sample(sub.all.potential.group,size=1)
      sub.core=which.max(apply(cor.xx.tmp,1,max))
      
      
      sub.selected_variables=order(cor.xx[sub.core,],decreasing=TRUE)[1:pB]
      selected_variables=rest_variables[sub.selected_variables]
      assign(sprintf("block%s",block.iter),selected_variables)
      rest_variables=setdiff(rest_variables,selected_variables)
      non_blocked_variables=setdiff(non_blocked_variables,selected_variables)
      cor.xx=cor.xx[-sub.selected_variables,-sub.selected_variables]
    }
    signal.positions=c(); 
    for(block.id in 1:nblock)
    {block=get(sprintf("block%s",block.id)); signal.positions=c(signal.positions,block[1:s0B])}
  }   #######End of if(s0B>0)#########
  
  
  
  if(length(unique(signal.positions))!=s0) {stop("WRONG!")}
  row.idx=sample(1:nrow(XR),n,replace=FALSE); test.row.idx=setdiff(1:nrow(XR),row.idx); 
  Xr=XR[row.idx,]; rownames(Xr)=1:n; Xr.test=XR[test.row.idx,];n.test=length(test.row.idx)
  
  for(i in 1:ncol(Xr)){temp=Xr[,i];temp=(temp-mean(temp))/sd(temp);Xr[,i]=temp}
  for(i in 1:ncol(Xr.test)){temp=Xr.test[,i];temp=(temp-mean(temp))/sd(temp);Xr.test[,i]=temp}
  
  beta0=rep(0,p)
  beta0[signal.positions]=3
  Y0=Xr%*%beta0 ;sigmae=as.numeric(sqrt(t(beta0)%*%t(Xr)%*%Xr%*%beta0/(n*SNR^2)))
  {e.training=rnorm(n,0,sigmae); e.test=rnorm(n.test,0,sigmae)}
  
  
  Y=Y0+e.training; Y=Y-mean(Y)
  Y0.test=Xr.test%*%beta0
  Y.test=Y0.test+e.test; Y.test=Y.test-mean(Y.test)

  result=list()
  result$Xr=Xr; result$Y=Y; result$Xr.test=Xr.test; result$Y.test=Y.test; result$beta0=beta0; result$sigma=sigma
  result
}


pauc.lasso=c(); rmse.lasso=c(); tpr.lasso=c(); ppv.lasso=c()
pauc.henet=c(); rmse.henet=c(); tpr.henet=c(); ppv.henet=c()
pauc.scad=c(); rmse.scad=c(); tpr.scad=c(); ppv.scad=c()
pauc.ridge=c(); rmse.ridge=c(); tpr.ridge=c(); ppv.ridge=c()
pauc.stability=c(); rmse.stability=c(); tpr.stability=c(); ppv.stability=c()
pauc.adalasso=c(); rmse.adalasso=c(); tpr.adalasso=c(); ppv.adalasso=c()



for(iter in 1:niter)
{
  
  data = data_generation(n, p, s0, SNR, corDesign, pB, s0B)
  Xr = data$Xr; Y = data$Y; Xr.test = data$Xr.test; Y.test = data$Y.test; beta0 = data$beta0
  xlim=c(0,50/(p-s0))


  ###################################################
  # run methods
  ###################################################
  
  ##################################Lasso####################
  cv=cv.glmnet(x=Xr,y=Y,alpha=1,nlambda=100)
  lasso_run=glmnet(x=Xr,y=Y,alpha=1,lambda=cv$lambda)
  solution_path=Matrix(lasso_run$beta,sparse=TRUE)
  colnames(solution_path)=cv$lambda
  if(nrow(solution_path)!=p) {stop("Something wrong with solution path!")}
  lambda.min_idx=which.min(abs(cv$lambda-cv$lambda.min))
  beta.cv=solution_path[,lambda.min_idx]
  
  
  
  score=rep(0,p)
  for(var_idx in 1:p)
  {
    if(sum(solution_path[var_idx,]!=0)==0) {score[var_idx]=0}  ######Never entered solution path#####
    if(sum(solution_path[var_idx,]==0)==0) {score[var_idx]=cv$lambda[1]}
    
    if(sum(solution_path[var_idx,]!=0)!=0 & sum(solution_path[var_idx,]==0)!=0)
    {
      idx=max(which(solution_path[var_idx,]==0))
      if(idx!=length(cv$lambda)) {idx=idx+1; score[var_idx]=cv$lambda[idx]}
    }
  }    #####End of for(var_idx in 1:p)######
  pred=prediction(abs(score),ifelse(beta0!=0,1,0))
  auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
  
  
  pauc.lasso[iter]=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
  rmse.lasso[iter]=sqrt(mean((Y.test-Xr.test%*%beta.cv)^2))
  tp=sum(beta0!=0 & beta.cv!=0); fp=sum(beta0==0 & beta.cv!=0); tn=sum(beta0==0 & beta.cv==0); fn=sum(beta0!=0 & beta.cv==0)
  tpr.lasso[iter]=tp/(tp+fn)
  ppv.lasso[iter]=tp/(tp+fp)
  
  ############################################################
  
  
  cv=cv.glmnet(x=Xr,y=Y,alpha=0.3,nlambda=100)
  lasso_run=glmnet(x=Xr,y=Y,alpha=0.3,lambda=cv$lambda)
  solution_path=Matrix(lasso_run$beta,sparse=TRUE)
  colnames(solution_path)=cv$lambda
  if(nrow(solution_path)!=p) {stop("Something wrong with solution path!")}
  lambda.min_idx=which.min(abs(cv$lambda-cv$lambda.min))
  beta.cv=solution_path[,lambda.min_idx]
  
  
  score=rep(0,p)
  for(var_idx in 1:p)
  {
    if(sum(solution_path[var_idx,]!=0)==0) {score[var_idx]=0}  ######Never entered solution path#####
    if(sum(solution_path[var_idx,]==0)==0) {score[var_idx]=cv$lambda[1]}
    
    if(sum(solution_path[var_idx,]!=0)!=0 & sum(solution_path[var_idx,]==0)!=0)
    {
      idx=max(which(solution_path[var_idx,]==0))
      if(idx!=length(cv$lambda)) {idx=idx+1; score[var_idx]=cv$lambda[idx]}
    }
  }    #####End of for(var_idx in 1:p)######
  pred=prediction(abs(score),ifelse(beta0!=0,1,0))
  auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
  
  
  pauc.henet[iter]=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
  rmse.henet[iter]=sqrt(mean((Y.test-Xr.test%*%beta.cv)^2))
  tp=sum(beta0!=0 & beta.cv!=0); fp=sum(beta0==0 & beta.cv!=0); tn=sum(beta0==0 & beta.cv==0); fn=sum(beta0!=0 & beta.cv==0)
  tpr.henet[iter]=tp/(tp+fn)
  ppv.henet[iter]=tp/(tp+fp)
  
  
  ##################################Ridge####################
  cv=cv.glmnet(x=Xr,y=Y,alpha=0,nlambda=100)
  lasso_run=glmnet(x=Xr,y=Y,alpha=0,lambda=cv$lambda)
  solution_path=Matrix(lasso_run$beta,sparse=TRUE)
  colnames(solution_path)=cv$lambda
  if(nrow(solution_path)!=p) {stop("Something wrong with solution path!")}
  lambda.min_idx=which.min(abs(cv$lambda-cv$lambda.min))
  beta.cv=solution_path[,lambda.min_idx]
  score=abs(beta.cv)
  pred=prediction(abs(score),ifelse(beta0!=0,1,0))
  auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
  
  
  pauc.ridge[iter]=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
  rmse.ridge[iter]=sqrt(mean((Y.test-Xr.test%*%beta.cv)^2))
  
  
  #################################SCAD###########################
  cv=cv.ncvreg(X=Xr,y=Y,penalty="SCAD",max.iter=5000,nlambda=100)
  scad_run=ncvreg(X=Xr,y=Y,family="gaussian",penalty="SCAD",max.iter=5000,lambda=cv$lambda)
  solution_path=Matrix(scad_run$beta[-1,],sparse=TRUE)
  colnames(solution_path)=cv$lambda
  if(nrow(solution_path)!=p) {stop("Something wrong with solution path!")}
  lambda.min_idx=which.min(abs(cv$lambda-cv$lambda.min))
  beta.cv=solution_path[,lambda.min_idx]
  
  
  tp=sum(beta.cv!=0 & beta0!=0)
  fp=sum(beta.cv!=0 & beta0==0)
  tn=sum(beta.cv==0 & beta0==0)
  fn=sum(beta.cv==0 & beta0!=0)
  if(tp+fp+tn+fn!=p) {warning("wrong tp/fp/tn/fn!")}
  
  
  score=rep(0,p)
  for(var_idx in 1:p)
  {
    if(sum(solution_path[var_idx,]!=0)==0) {score[var_idx]=0}  ######Never entered solution path#####
    if(sum(solution_path[var_idx,]==0)==0) {score[var_idx]=cv$lambda[1]}
    
    if(sum(solution_path[var_idx,]!=0)!=0 & sum(solution_path[var_idx,]==0)!=0)
    {
      idx=max(which(solution_path[var_idx,]==0))
      if(idx!=length(cv$lambda)) {idx=idx+1; score[var_idx]=cv$lambda[idx]}
    }
  }    #####End of for(var_idx in 1:p)######
  pred=prediction(abs(score),ifelse(beta0!=0,1,0))
  auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
  
  
  
  
  pauc.scad[iter]=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
  rmse.scad[iter]=sqrt(mean((Y.test-Xr.test%*%beta.cv)^2))
  tp=sum(beta0!=0 & beta.cv!=0); fp=sum(beta0==0 & beta.cv!=0); tn=sum(beta0==0 & beta.cv==0); fn=sum(beta0!=0 & beta.cv==0)
  
  tpr.scad[iter]=tp/sum(beta0!=0)
  ppv.scad[iter]=tp/sum(beta.cv!=0)
  
  
  ########################################Stability selection######################################
  stability_run=stabpath(y=Y,x=Xr,size=0.632,steps=100,weakness=1,nlambda=100,mc.cores=6)
  solution_path=stability_run$x
  if(nrow(solution_path)!=p) {stop("Something wrong with solution path!")}
  beta.full=apply(solution_path,1,max)
  
  
  
  #beta.cv=beta.full
  #beta.cv[beta.cv<0.6]=0
  stab=stability(x=Xr,y=Y,EV=10,threshold=0.6,B=100,fraction=0.632,ncores=6);
  beta.cv=rep(0,p);beta.cv[stab$selected]=1
  
  
  
  model.size=sum(beta.cv!=0)
  tp=sum(beta.cv!=0 & beta0!=0)
  fp=sum(beta.cv!=0 & beta0==0)
  tn=sum(beta.cv==0 & beta0==0)
  fn=sum(beta.cv==0 & beta0!=0)
  tpr=tp/(tp+fn)
  ppv=tp/(tp+fp)
  if(tp+fp+tn+fn!=p) {warning("wrong tp/fp/tn/fn!")}
  if((tp+fp)==0) {mcc=0;ppv=0}
  if((tp+fp)!=0){mcc=((tp*tn-fp*fn)/p)/sqrt((tp+fp)*(tp+fn)/p^2*(tn+fp)*(tn+fn));ppv=tp/(tp+fp)}
  
  
  pred=prediction(abs(beta.full),ifelse(beta0!=0,1,0))
  auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
  
  
  
  pauc.stability[iter]=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
  tpr.stability[iter]=tpr
  ppv.stability[iter]=ppv
  
  
  
  ####################################Adalasso#######################################
  
  ridge1_cv=cv.glmnet(x = Xr, y = Y, alpha = 0)
  best_ridge_coef <- as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1]
  alasso1=glmnet(x= Xr, y = Y, alpha = 1,penalty.factor = 1 / abs(best_ridge_coef))
  cv=cv.glmnet(x = Xr, y = Y, alpha = 1, penalty.factor = 1 / abs(best_ridge_coef))
  beta.cv= as.numeric(coef(cv, s = cv$lambda.min))[-1]
  
  
  
  solution_path=Matrix(alasso1$beta,sparse=TRUE)
  #colnames(solution_path)=cv$lambda
  
  score=rep(0,p)
  for(var_idx in 1:p)
  {
    if(sum(solution_path[var_idx,]!=0)==0) {score[var_idx]=0}  ######Never entered solution path#####
    if(sum(solution_path[var_idx,]==0)==0) {score[var_idx]=cv$lambda[1]}
    
    if(sum(solution_path[var_idx,]!=0)!=0 & sum(solution_path[var_idx,]==0)!=0)
    {
      idx=max(which(solution_path[var_idx,]==0))
      if(idx!=length(cv$lambda)) {idx=idx+1; score[var_idx]=cv$lambda[idx]}
    }
  }    #####End of for(var_idx in 1:p)######
  pred=prediction(abs(score),ifelse(beta0!=0,1,0))
  auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
  
  
  pauc.adalasso[iter]=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
  rmse.adalasso[iter]=sqrt(mean((Y.test-Xr.test%*%beta.cv)^2))
  tp=sum(beta0!=0 & beta.cv!=0); fp=sum(beta0==0 & beta.cv!=0); tn=sum(beta0==0 & beta.cv==0); fn=sum(beta0!=0 & beta.cv==0)
  tpr.adalasso[iter]=tp/(tp+fn)
  ppv.adalasso[iter]=tp/(tp+fp)
  
  
  #######################################################################################################
} #######End of for(iter in 1:niter)######### 
