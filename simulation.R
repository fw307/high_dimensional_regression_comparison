Independence="Independence"; Pairwise="Pairwise"; Toeplitz="Toeplitz"
setwd(settings_directory) ########directory containing "SETTINGS.txt"##########
settings=as.vector(as.matrix(read.table("SETTINGS.txt",head=FALSE)))
setting=settings[1]

library(MASS)                       
library(glmnet)
library(parcor)
library(ncvreg)
library(monomvn)
library(flare)                  
library(c060)
library(Matrix)

alpha1=0.3; alpha2=0.6; beta_value=3
method_name=c("lasso","lenet","henet","ridge","dantzig","scad","stability")
for(name in method_name) {assign(sprintf("%s_list",name),list())}



setwd(script_directory)  #######directory containing "data_generation.R"##########
source("data_generation.R")
niter=1     







                                               for(iter in 1:niter)                                        
{
eval(parse(text=setting))
DATA=data_generation(setting)
X=DATA$X;Y=DATA$Y;X_test=DATA$X_test;Y_test=DATA$Y_test;
beta0=DATA$beta0; sigma=DATA$sigma

nfolds=5      
foldid=sample(rep(1:nfolds,n/nfolds),n)
xlim=c(0,50/(p-s0))





#############################################################Lasso###############################################
lasso=list()
cv=cv.glmnet(x=X,y=Y,alpha=1,foldid=foldid,nlambda=100)
ptm=proc.time()
lasso_run=glmnet(x=X,y=Y,alpha=1,lambda=cv$lambda)
run_time=proc.time()-ptm
solution_path=Matrix(lasso_run$beta,sparse=TRUE)
colnames(solution_path)=cv$lambda
lasso$run.time=run_time
lambda.min_idx=which.min(abs(cv$lambda-cv$lambda.min))
beta.cv=solution_path[,lambda.min_idx]
tp=sum(beta.cv!=0 & beta0!=0)
fp=sum(beta.cv!=0 & beta0==0)
tn=sum(beta.cv==0 & beta0==0)
fn=sum(beta.cv==0 & beta0!=0)
if((tp+fp)==0) {mcc=0;ppv=0}
if((tp+fp)!=0){mcc=((tp*tn-fp*fn)/p)/sqrt((tp+fp)*(tp+fn)/p^2*(tn+fp)*(tn+fn));ppv=tp/(tp+fp)}
tpr=tp/(tp+fn)
lasso$mcc=mcc
lasso$ppv=ppv
lasso$tpr=tpr


rmse=sqrt(mean((Y_test-X_test%*%beta.cv-mean(Y_test))^2))
lasso$rmse=rmse
score=rep(0,p)
for(var_idx in 1:p)
{
if(sum(solution_path[var_idx,]!=0)==0) {score[var_idx]=0}  ######Never entered solution path#####
if(sum(solution_path[var_idx,]==0)==0) {score[var_idx]=cv$lambda[1]}   ######Always in solution path, save the largest lambda####

        if(sum(solution_path[var_idx,]!=0)!=0 & sum(solution_path[var_idx,]==0)!=0)
    {
        idx=max(which(solution_path[var_idx,]==0))
        if(idx!=length(cv$lambda)) {idx=idx+1; score[var_idx]=cv$lambda[idx]}
    }
}    
pred=prediction(abs(score),ifelse(beta0!=0,1,0))
auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
pauc=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
lasso$pauc=pauc
lasso_list[[iter]]=lasso
















###############################################################Lenet#############################################
lenet=list()       
cv=cv.glmnet(x=X,y=Y,alpha=alpha2,foldid=foldid,nlambda=100)
ptm=proc.time()
lenet_run=glmnet(x=X,y=Y,alpha=alpha2,lambda=cv$lambda)
run_time=proc.time()-ptm
solution_path=Matrix(lenet_run$beta,sparse=TRUE)
colnames(solution_path)=cv$lambda
lenet$run.time=run_time
lambda.min_idx=which.min(abs(cv$lambda-cv$lambda.min))
beta.cv=solution_path[,lambda.min_idx]

tp=sum(beta.cv!=0 & beta0!=0)
fp=sum(beta.cv!=0 & beta0==0)
tn=sum(beta.cv==0 & beta0==0)
fn=sum(beta.cv==0 & beta0!=0)
if((tp+fp)==0) {mcc=0;ppv=0}
if((tp+fp)!=0){mcc=((tp*tn-fp*fn)/p)/sqrt((tp+fp)*(tp+fn)/p^2*(tn+fp)*(tn+fn));ppv=tp/(tp+fp)}
tpr=tp/(tp+fn)
lenet$mcc=mcc
lenet$ppv=ppv
lenet$tpr=tpr
rmse=sqrt(mean((Y_test-X_test%*%beta.cv-mean(Y_test))^2))
lenet$rmse=rmse

score=rep(0,p)
for(var_idx in 1:p)
{
    if(sum(solution_path[var_idx,]!=0)==0) {score[var_idx]=0}   
    if(sum(solution_path[var_idx,]==0)==0) {score[var_idx]=cv$lambda[1]}
    if(sum(solution_path[var_idx,]!=0)!=0 & sum(solution_path[var_idx,]==0)!=0)
    {
        idx=max(which(solution_path[var_idx,]==0))
        if(idx!=length(cv$lambda)) {idx=idx+1;score[var_idx]=cv$lambda[idx]}
    }
}     
pred=prediction(abs(score),ifelse(beta0!=0,1,0))
auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
pauc=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
lenet$pauc=pauc
lenet_list[[iter]]=lenet

















##############################################################Henet##################
henet=list()
cv=cv.glmnet(x=X,y=Y,alpha=alpha1,foldid=foldid,nlambda=100)
ptm=proc.time()
henet_run=glmnet(x=X,y=Y,alpha=alpha1,lambda=cv$lambda)
run_time=proc.time()-ptm
solution_path=Matrix(henet_run$beta,sparse=TRUE)
colnames(solution_path)=cv$lambda
henet$run.time=run_time

lambda.min_idx=which.min(abs(cv$lambda-cv$lambda.min))
beta.cv=solution_path[,lambda.min_idx]
tp=sum(beta.cv!=0 & beta0!=0)
fp=sum(beta.cv!=0 & beta0==0)
tn=sum(beta.cv==0 & beta0==0)
fn=sum(beta.cv==0 & beta0!=0)
if(tp+fp+tn+fn!=p) {warning("wrong tp/fp/tn/fn!")}
if((tp+fp)==0) {mcc=0;ppv=0}
if((tp+fp)!=0){mcc=((tp*tn-fp*fn)/p)/sqrt((tp+fp)*(tp+fn)/p^2*(tn+fp)*(tn+fn));ppv=tp/(tp+fp)}
tpr=tp/(tp+fn)
henet$mcc=mcc
henet$ppv=ppv
henet$tpr=tpr
rmse=sqrt(mean((Y_test-X_test%*%beta.cv-mean(Y_test))^2))
henet$rmse=rmse


score=rep(0,p)
for(var_idx in 1:p)
{
    if(sum(solution_path[var_idx,]!=0)==0) {score[var_idx]=0}   
    if(sum(solution_path[var_idx,]==0)==0) {score[var_idx]=cv$lambda[1]}
    if(sum(solution_path[var_idx,]!=0)!=0 & sum(solution_path[var_idx,]==0)!=0)
    {
        idx=max(which(solution_path[var_idx,]==0))
        if(idx!=length(cv$lambda)) {idx=idx+1; score[var_idx]=cv$lambda[idx]}
    }
}    
pred=prediction(abs(score),ifelse(beta0!=0,1,0))
auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
pauc=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
henet$pauc=pauc
henet_list[[iter]]=henet










                                               
                                            
                  



                                     

####################################################################Ridge Estimator#########################################################
ridge=list()
cv=cv.glmnet(x=X,y=Y,alpha=0,foldid=foldid,nlambda=100)
ptm=proc.time()
ridge_run=glmnet(x=X,y=Y,alpha=0,lambda=cv$lambda)
run_time=proc.time()-ptm

solution_path=Matrix(ridge_run$beta,sparse=TRUE)
colnames(solution_path)=cv$lambda
ridge$run.time=run_time

Beta=predict(ridge_run,type="coef",s=cv$lambda.min)
beta.cv=Beta[-1]
rmse=sqrt(mean((Y_test-X_test%*%beta.cv-mean(Y_test))^2))
ridge$rmse=rmse

pred=prediction(abs(beta.cv),ifelse(beta0!=0,1,0))
auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
pauc=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
ridge$pauc=pauc
ridge_list[[iter]]=ridge







####################################################################SCAD####################################################
scad=list()
cv=cv.ncvreg(X=X,y=Y,penalty="SCAD",nfolds=nfolds,max.iter=5000,nlambda=100)
ptm=proc.time()
scad_run=ncvreg(X=X,y=Y,family="gaussian",penalty="SCAD",max.iter=5000,lambda=cv$lambda)
run_time=proc.time()-ptm


solution_path=Matrix(scad_run$beta[-1,],sparse=TRUE)
colnames(solution_path)=cv$lambda
if(nrow(solution_path)!=p) {stop("Something wrong with solution path!")}
scad$run.time=run_time
lambda.min_idx=which.min(abs(cv$lambda-cv$lambda.min))
beta.cv=solution_path[,lambda.min_idx]

tp=sum(beta.cv!=0 & beta0!=0)
fp=sum(beta.cv!=0 & beta0==0)
tn=sum(beta.cv==0 & beta0==0)
fn=sum(beta.cv==0 & beta0!=0)
if(tp+fp+tn+fn!=p) {warning("wrong tp/fp/tn/fn!")}
if((tp+fp)==0) {mcc=0;ppv=0}
if((tp+fp)!=0){mcc=((tp*tn-fp*fn)/p)/sqrt((tp+fp)*(tp+fn)/p^2*(tn+fp)*(tn+fn));ppv=tp/(tp+fp)}
tpr=tp/(tp+fn)
scad$mcc=mcc
scad$ppv=ppv
scad$tpr=tpr
rmse=sqrt(mean((Y_test-X_test%*%beta.cv-mean(Y_test))^2))
scad$rmse=rmse

score=rep(0,p)
for(var_idx in 1:p)
{
    if(sum(solution_path[var_idx,]!=0)==0) {score[var_idx]=0}   
    if(sum(solution_path[var_idx,]==0)==0) {score[var_idx]=cv$lambda[1]}
    
    if(sum(solution_path[var_idx,]!=0)!=0 & sum(solution_path[var_idx,]==0)!=0)
    {
        idx=max(which(solution_path[var_idx,]==0))
        if(idx!=length(cv$lambda)) {idx=idx+1; score[var_idx]=cv$lambda[idx]}
    }
}     
pred=prediction(abs(score),ifelse(beta0!=0,1,0))
auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
pauc=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
scad$pauc=pauc
scad_list[[iter]]=scad








#######################################################Dantzig Selector#####################################
dantzig=list()
nlambda=10    ######number of lambda we try for finding the optimal one#####
ptm=proc.time()
dantzig_run=slim(X=X,Y=Y,method="dantzig",nlambda=nlambda,verbose=FALSE,lambda.min.ratio=0.01)
run_time=proc.time()-ptm

solution_path=Matrix(dantzig_run$beta,sparse=TRUE)
colnames(solution_path)=dantzig_run$lambda

#############k-fold cross-validation for each item of lambda_seq######
lambda_seq=dantzig_run$lambda
MSE_SEQ=c()   #####mean cross-validated MSE for each possible regularization level (lambda)######
MSE_SD_SEQ=c() #####Corresponding standard errors##########
k=nfolds  #####number of folds for cross-validation#####
samples=sample(1:n,n,replace=FALSE)

max=round(n/k)
x=seq_along(samples)
groups=split(samples,ceiling(x/max))   #####k groups(equal sizes) samples#####
for(i in 1:length(groups)){assign(sprintf("idx_%s",i),groups[[i]])}    

for(j in 1:length(lambda_seq))            #####for each possible regularization level#####
{
  lambda=lambda_seq[j]
  mse=c()
        for(k in 1:length(groups))
      {
         idx_test=get(sprintf("idx_%s",k))    
         idx_train=setdiff(1:n,idx_test)      
         X_test.part=X[idx_test,]
         Y_test.part=Y[idx_test]
         X_train.part=X[idx_train,]
         Y_train.part=Y[idx_train]
             
         dantzig_train=slim(X=X_train.part,Y=Y_train.part,method="dantzig",lambda=lambda,verbose=FALSE)
         beta_train=dantzig_train$beta
         Y_test_hat=X_test.part%*%beta_train+rep(dantzig_train$intercept,nrow(X_test.part))
         mse=c(mse,mean((Y_test.part-Y_test_hat)^2))
      }
   MSE_SEQ[j]=mean(mse)
   MSE_SD_SEQ[j]=sd(mse)
}


lambda.min=lambda_seq[match(min(MSE_SEQ),MSE_SEQ)]
cvm=MSE_SEQ
cvsd=MSE_SD_SEQ
dantzig$run.time=run_time

lambda.min_idx=which.min(abs(lambda_seq-lambda.min))
beta.cv=solution_path[,lambda.min_idx]
dantzig$beta.cv=beta.cv

tp=sum(beta.cv!=0 & beta0!=0)
fp=sum(beta.cv!=0 & beta0==0)
tn=sum(beta.cv==0 & beta0==0)
fn=sum(beta.cv==0 & beta0!=0)
if(tp+fp+tn+fn!=p) {warning("wrong tp/fp/tn/fn!")}
if((tp+fp)==0) {mcc=0;ppv=0}
if((tp+fp)!=0){mcc=((tp*tn-fp*fn)/p)/sqrt((tp+fp)*(tp+fn)/p^2*(tn+fp)*(tn+fn));ppv=tp/(tp+fp)}
tpr=tp/(tp+fn)
dantzig$mcc=mcc
dantzig$ppv=ppv
dantzig$tpr=tpr
rmse=sqrt(mean((Y_test-X_test%*%beta.cv-mean(Y_test))^2))
dantzig$rmse=rmse

score=rep(0,p)
for(var_idx in 1:p)
{
    if(sum(solution_path[var_idx,]!=0)==0) {score[var_idx]=0}  
    if(sum(solution_path[var_idx,]==0)==0) {score[var_idx]=lambda_seq[1]}
    if(sum(solution_path[var_idx,]!=0)!=0 & sum(solution_path[var_idx,]==0)!=0)
    {
        idx=max(which(solution_path[var_idx,]==0))
        if(idx!=length(lambda_seq)) {idx=idx+1; score[var_idx]=dantzig_run$lambda[idx]}
    }
}     
pred=prediction(abs(score),ifelse(beta0!=0,1,0))
auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
pauc=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
dantzig$pauc=pauc
dantzig_list[[iter]]=dantzig

























######################################################Stability Selection#################################
stability=list()
ptm=proc.time()
stability_run=stabpath(y=Y,x=X,size=0.632,steps=500,weakness=1,nlambda=100)
run_time=proc.time()-ptm
stability$run.time=run_time
solution_path=stability_run$x
beta.full=apply(solution_path,1,max)
beta.cv=beta.full
beta.cv[beta.cv<0.6]=0

tp=sum(beta.cv!=0 & beta0!=0)
fp=sum(beta.cv!=0 & beta0==0)
tn=sum(beta.cv==0 & beta0==0)
fn=sum(beta.cv==0 & beta0!=0)
if(tp+fp+tn+fn!=p) {warning("wrong tp/fp/tn/fn!")}
if((tp+fp)==0) {mcc=0;ppv=0}
if((tp+fp)!=0){mcc=((tp*tn-fp*fn)/p)/sqrt((tp+fp)*(tp+fn)/p^2*(tn+fp)*(tn+fn));ppv=tp/(tp+fp)}
tpr=tp/(tp+fn)
stability$mcc=mcc
stability$ppv=ppv
stability$tpr=tpr

pred=prediction(abs(beta.full),ifelse(beta0!=0,1,0))
auc.tmp=performance(pred,"auc",fpr.stop=xlim[2])
pauc=round(as.numeric(auc.tmp@y.values)/xlim[2],3)
stability$pauc=pauc
stability_list[[iter]]=stability
