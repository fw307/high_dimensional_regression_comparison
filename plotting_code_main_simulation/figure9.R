par(mar=c(2.5,2.5,2,2.4))
par(oma=c(2,3.5,2,2.4))
par(mgp=c(2,1,0))  

setwd(sorted_address)        ######sorted_address is the folder for results of all settings######
all_folders=list.files()
rho_seq=c(0,0.5,0.7,0.9)
s0_B_seq=c(1,5)
pB_seq=c(10,100)
SNR_seq=c(4,1)


method_seq_all=c("lasso","lenet","henet","ridge","dantzig","scad","stability")
library(RColorBrewer);colors_all=(brewer.pal(length(method_seq_all),"Set1"));colors_all[1]="red";colors_all[5]="black";colors_all[6]="gold3";colors_all[7]="brown"
colors_all=colors_all[c(1,4,3,6,2,5,7)]
colors=colors_all

METHOD_SEQ=method_seq_all
COLORS=colors

method_seq.partial=METHOD_SEQ[c(1,3,6,7)]
colors.partial=COLORS[c(1,3,6,7)]

method_seq.all=METHOD_SEQ[c(1,3,5,6,7)]
colors.all=COLORS[c(1,3,5,6,7)]
lty_seq=c(1,2)
setting_seq=c("p=4000;s0=40","p=1000;s0=10")
n=300


metric1="tpr";metric2="ppv";metric3="tpr";metric4="ppv"
flag=0;XLIM=c(0,1)
pch_seq=c(0,1,2,3)
figure.count=1

par(mfrow=c(4,4))


 for(SNR.adjust in SNR_seq)
{
    for(pB in pB_seq)
 {
       for(s0_B in s0_B_seq) 
    {
          setting=setting_seq[1]
      

          eval(parse(text=setting)); SNR=SNR.adjust
          lty=lty_seq[match(setting,setting_seq)]
 if(p==1000){method_seq=method_seq.all;colors=colors.all}
 if(p!=1000){method_seq=method_seq.partial;colors=colors.partial}
 
           for(method in method_seq)
       {
           factorY1_seq=c()
           factorY2_seq=c()

                                                             for(rho in rho_seq)
                                                        {
                                             if(rho>0)                setting=sprintf("n=%s;p=%s;s0=%s;SNR=%s;Design=Pairwise;rho=%s;pB=%s;s0_B=%s",n,p,s0,SNR,rho,pB,s0_B)
                                             if(rho==0)                setting=sprintf("n=%s;p=%s;s0=%s;SNR=%s;Design=Independence;rho=0;pB=0;s0_B=0",n,p,s0,SNR)

                                                             setwd(sorted_address);setwd(setting)
                                                             factorY1=get(load(sprintf("%s_%s_%s.Rdata",metric1,method,setting)))
                                                             factorY2=get(load(sprintf("%s_%s_%s.Rdata",metric2,method,setting)))
                                                             col=colors[match(method,method_seq)];main=sprintf("%s: s0B=%s, pB=%s",LETTERS[figure.count],s0_B,pB) 

                                                             if(figure.count==1 | figure.count==5) main=sprintf("n=%s,SNR=%s \n%s: s0B=%s,pB=%s",n,SNR,LETTERS[figure.count],s0_B,pB) 
                                                             mean.factorY1=mean(factorY1); mean.factorY2=mean(factorY2);sd.factorY1=sd(factorY1); sd.factorY2=sd(factorY2)
                                                             factorY1_seq[match(rho,rho_seq)]=mean.factorY1; factorY2_seq[match(rho,rho_seq)]=mean.factorY2


xlab=toupper(metric1);ylab=toupper(metric2)
if(length(grep("pauc",metric1))>0){xlab="pAUC"}; if(length(grep("rmse",metric1))>0){xlab="RMSE"}; if(length(grep("smse",metric1))>0){xlab="SMSE"};if(length(grep("mcc",metric1))>0){xlab="MCC"};if(length(grep("tpr",metric1))>0){xlab="TPR"};if(length(grep("fscore",metric1))>0){xlab="F"} 
if(length(grep("pauc",metric2))>0){ylab="pAUC"}; if(length(grep("rmse",metric2))>0){ylab="RMSE"}; if(length(grep("smse",metric2))>0){ylab="SMSE"};if(length(grep("mcc",metric2))>0){ylab="MCC"};if(length(grep("tpr",metric2))>0){ylab="TPR"};if(length(grep("fscore",metric2))>0){ylab="F"};


                                                             plot(factorY1_seq,factorY2_seq,xlim=c(0,1),ylim=c(0,1),col=col,main="",xlab="",ylab="",pch=pch_seq,lty=lty,xaxt="n",yaxt="n")
                                                             lines(factorY1_seq,factorY2_seq,col=col,lty=lty)
                                                             axis(1,at=c(0,0.3,0.6,0.9),cex.axis=0.8)
                                                             axis(2,at=c(0,0.3,0.6,0.9),cex.axis=0.8)
mtext(side=1, text=xlab, line=2,cex=0.8)
mtext(side=2, text=ylab, line=2,cex=0.8)
mtext(text=LETTERS[figure.count],side=3,adj=0,cex=1)
                                                             par(new=TRUE)
                                                        }     

       }   
 
      
       if(figure.count %in% c(4,8)) {mtext(side=4,text=sprintf("SNR=%s",SNR),line=0.5,cex=1,font=2,las=1)}
       if(figure.count %in% 1:4) {mtext(side=3,text=bquote(bold('s'[0]^'B' == .(s0_B))),line=1,cex=1,font=2,las=1)}
       par(new=FALSE); figure.count=figure.count+1
    }
 }
}































 for(SNR.adjust in SNR_seq)
{

    for(pB in pB_seq)
 {
       for(s0_B in s0_B_seq) 
    {

           setting=setting_seq[2]
 
          eval(parse(text=setting)); SNR=SNR.adjust
          if(p==1000){method_seq=method_seq.all;colors=colors.all}
          if(p!=1000){method_seq=method_seq.partial;colors=colors.partial}
          lty=lty_seq[match(setting,setting_seq)]

           for(method in method_seq)
       {
           factorY1_seq=c()
           factorY2_seq=c()
   

                                                           for(rho in rho_seq)
                                                        {
                                              if(rho>0)               {setting=sprintf("n=%s;p=%s;s0=%s;SNR=%s;Design=Pairwise;rho=%s;pB=%s;s0_B=%s",n,p,s0,SNR,rho,pB,s0_B)}
                                              if(rho==0)               {setting=sprintf("n=%s;p=%s;s0=%s;SNR=%s;Design=Independence;rho=0;pB=0;s0_B=0",n,p,s0,SNR)}

                                                             setwd(sorted_address);setwd(setting)
                                                             factorY1=get(load(sprintf("%s_%s_%s.Rdata",metric3,method,setting)))
                                                             factorY2=get(load(sprintf("%s_%s_%s.Rdata",metric4,method,setting)))
                                                             col=colors[match(method,method_seq)];main=sprintf("%s: s0B=%s,pB=%s",LETTERS[figure.count],s0_B,pB) 

                                                             if(figure.count==9 | figure.count==13) main=sprintf("n=%s,SNR=%s  \n%s: s0B=%s,pB=%s",n,SNR,LETTERS[figure.count],s0_B,pB) 
                                                             mean.factorY1=mean(factorY1); mean.factorY2=mean(factorY2);sd.factorY1=sd(factorY1); sd.factorY2=sd(factorY2)
                                                             factorY1_seq[match(rho,rho_seq)]=mean.factorY1; factorY2_seq[match(rho,rho_seq)]=mean.factorY2

xlab="TPR"; ylab="PPV"
                                                             plot(factorY1_seq,factorY2_seq,xlim=c(0,1),ylim=c(0,ifelse(length(grep("rmse",metric4))>0,13,1)),col=col,main="",xlab="",ylab="",pch=pch_seq,xaxt="n",yaxt="n")
                                                             lines(factorY1_seq,factorY2_seq,col=col,lty=lty)
                                                             axis(1,at=c(0,0.3,0.6,0.9),cex.axis=0.8)
                                                             axis(2,at=c(0,0.3,0.6,0.9),cex.axis=0.8)
mtext(side=1, text=xlab, line=2,cex=0.8)
mtext(side=2, text=ylab, line=2,cex=0.8)
mtext(text=LETTERS[figure.count],side=3,adj=0,cex=1)
                                                             par(new=TRUE)
                                                        }      
        }    
  
       if(figure.count %in% c(12,16)) {mtext(side=4,text=sprintf("SNR=%s",SNR),line=0.5,cex=1,font=2,las=1)}
       par(new=FALSE);figure.count=figure.count+1                                    
    }
 }
}
mtext(text=bquote(bold('p'^'B' == 10)),at=0.25,side=3,outer=T,cex=1,font=2,line=0.2)  
mtext(text=bquote(bold('p'^'B' == 100)),at=0.75,side=3,outer=T,cex=1,font=2,line=0.2)  

mtext(gsub(";","\n",setting_seq[1]),at=0.75,side=2,outer=T,cex=1,font=2,las=2,line=-1.2)  
mtext(gsub(";","\n",setting_seq[2]),at=0.25,side=2,outer=T,cex=1,font=2,las=2,line=-1.2)  


all_methods_text=paste(method_seq,collapse=",") 
main=sprintf("Pairwise correlation,n=%s,p=%s,s0=%s,SNR=%s",n,p,s0,SNR) 
MAIN=main
par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE); plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend("bottom",method_seq.all,xpd=TRUE,horiz=TRUE,inset=c(0,0),lty=rep(1,4),bty="n",col=colors.all,cex=1,lwd=3)
legend("bottomleft",legend=paste("rho",rho_seq,sep="="),cex=0.9,pch=pch_seq,xpd=TRUE,horiz=FALSE,inset=c(0,0))
