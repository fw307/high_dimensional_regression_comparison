data_generation=function(setting)
{
   eval(parse(text = setting))
   n_training=n; n_test=500
   n=n_training+n_test  
   
   
   if(Design=="Independence")
   {
      X=matrix(rnorm(n*p),n,p)   
      gamma=rep(0,p); gamma[sample(1:p,size=s0,replace=FALSE)]=1
      beta0=rep(0,p); beta0[which(gamma==1)]=beta_value
   }
   if(Design=="Pairwise")
   {
      Block_rho=rho
      n_block=p/pB
      block=function()
      {
         signal_idx=sample(1:pB,size=s0_B,replace=FALSE)   ######Where you put the signals in each block#####
         e0=rep(1,pB)          #####A vector of 1's####
         Sigma=Block_rho*(e0%*%t(e0)-diag(pB))+diag(pB)
         result=mvrnorm(n=n,mu=rep(0,pB),Sigma=Sigma,tol=1e-12)    #####mvrnorm has realization in each ROW#######
         colnames(result)=rep(0,pB)
         colnames(result)[signal_idx]=1
         result
      }      
      
      X=block(); if(n_block>=2) {for(b in 1:(n_block-1)){X=cbind(X,block())}}
      gamma=as.numeric(colnames(X)); beta0=rep(0,p); beta0[which(gamma==1)]=beta_value
      
   }     
   
   if(Design=="Toeplitz")
   {
      rho_base=0.95
      n_block=p/pB
      cor_between_two_signals=rho
      no_signal_per_block=s0_B
      
      block=function()
      {
         Sigma=matrix(0,pB,pB)
         for(i in 1:ncol(Sigma))                 
         {
            for(j in 1:i)                            
            {Sigma[j,i]=rho_base^(abs(i-j))}
         }
         Sigma[lower.tri(Sigma)] = t(Sigma)[lower.tri(Sigma)]     
         result=mvrnorm(n=n,mu=rep(0,pB),Sigma=Sigma,tol=1e-12)    
         colnames(result)=rep(0,Block_size)
         if(no_signal_per_block!=2) {stop("Should be 2 signals per block for Toeplitz Design!")}
         dis=round(log(cor_between_two_signals,base=rho_base))
         if(1+dis>=pB) {stop("The block size is too small for Toeplitz Design!")}
         signal_idx=c(round(pB/2),round(pB/2)+dis)      
         colnames(result)[signal_idx]=1  
         result
      }
      X=block(); for(b in 1:(n_block-1)){X=cbind(X,block())}
      gamma=as.numeric(colnames(X)); beta0=rep(0,p); beta0[which(gamma==1)]=beta_value
   }     
   
   
   non_zero_idx=which(beta0!=0); if(length(non_zero_idx)>s0) {beta0[-non_zero_idx[1:s0]]=0} 
   #################for correlation designs, if not all blocks should have signals, those blocks besides first s0 signals are empty##########
   
   
   Xr=X; X_training=Xr[1:n_training,];X_test=Xr[(n_training+1):n,]
   X_training=scale(X_training); X_test=scale(X_test)
   sigmae=as.numeric(sqrt(t(beta0)%*%t(X_training)%*%X_training%*%beta0/(n_training*SNR^2)))
   Xr=rbind(X_training,X_test)              
   Y0=Xr%*%beta0
   
   Y0_training=Y0[1:n_training]
   e_training=rnorm(n_training,0,sigmae)
   Y_training=Y0_training+e_training
   Y_training=Y_training-mean(Y_training)
   Y0_test=Y0[(n_training+1):n]
   e_test=rnorm(n_test,0,sigmae)
   Y_test=Y0_test+e_test
   Y_test=Y_test-mean(Y_test)
   
   
   if(length(e_training)!=length(Y0_training) | length(e_test)!=length(Y0_test)) {stop("WRONG!")}
   
   
   result=list()
   result$X=X_training; result$Y=Y_training; result$X_test=X_test; result$Y_test=Y_test; result$beta0=beta0; result$sigma=sigma
   result
}                    