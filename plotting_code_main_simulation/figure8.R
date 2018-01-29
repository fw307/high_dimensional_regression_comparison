par(mar=c(2.5,2.5,2.5,2.4))
par(oma=c(2,3.5,2,2.4))
par(mgp=c(2,1,0))  
setwd(sorted_address)        ######sorted_address is the folder for results of all settings######
all_folders=list.files()




n=300
setting_seq=c("p=4000;s0=40","p=1000;s0=10")
rho_seq=c(0,0.5,0.7,0.9);pB_seq=c(10,100);s0_B_seq=c(1,5)
SNR_seq=c(4,1)




method_seq_all=c("lasso","lenet","henet","ridge","dantzig","scad","stability")
library(RColorBrewer);colors_all=(brewer.pal(length(method_seq_all),"Set1"));colors_all[1]="red";colors_all[5]="black";colors_all[6]="gold3";colors_all[7]="brown"
colors_all=colors_all[c(1,4,3,6,2,5,7)]
colors=colors_all

method_seq.all=method_seq_all[-2]
colors.all=colors_all[-2]
method_seq.partial=method_seq_all[-c(2,5)]
colors.partial=colors_all[-c(2,5)]

metric="rmse"; Independence="Independence";YLIM=c(0,55)
pch_seq=c(0,2,5);lty_seq=c(5,3,1)
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
 if(p<=1000){method_seq=method_seq.all;colors=colors.all}
 if(p>1000){method_seq=method_seq.partial;colors=colors.partial}
 main=bquote(bold('s'[0]^'B' == .(s0_B) ~ ',' ~ 'p'^'B'== .(pB)))

 
           for(method in method_seq)
       {
                                                      mean.factorY_seq=c();sd.factorY_seq=c();col=colors[match(method,method_seq)]


                                                             for(rho in rho_seq)
                                                        {
                                           design=ifelse(rho==0,"Independence","Pairwise");
                                                      setting=sprintf("n=%s;p=%s;s0=%s;SNR=%s;Design=%s;rho=%s;pB=%s;s0_B=%s",n,p,s0,SNR,design,rho,ifelse(rho==0,0,pB),ifelse(rho==0,0,s0_B))
                                                      setwd(sorted_address);setwd(setting)
if(!(method=="dantzig" & p>1000))                                                      factorY=get(load(sprintf("%s_%s_%s.Rdata",metric,method,setting)))
if(!(method=="dantzig" & p>1000))                                                      mean.factorY=mean(factorY);sd.factorY=sd(factorY);mean.factorY_seq=c(mean.factorY_seq,mean.factorY);sd.factorY_seq=c(sd.factorY_seq,sd.factorY)

                                                        }    
 
           if(SNR==SNR_seq[1]){pch=pch_seq[1];lty=lty_seq[1]};if(SNR==SNR_seq[2]){pch=pch_seq[3];lty=lty_seq[3]}
           rho_seq.adjusted=rho_seq+(-0.02+0.005*match(method,method_seq))


           if(length(grep("pauc",metric))>0) {ylab="pAUC"}; if(length(grep("rmse",metric))>0) {ylab="RMSE"}; if(length(grep("smse",metric))>0) {ylab="SMSE"}
if(!(method=="dantzig" & p>1000))            {plot(rho_seq.adjusted,mean.factorY_seq,col=col,main="",ylim=YLIM,pch=pch,xlab="",ylab="",xaxt="n",xlim=range(rho_seq));lines(rho_seq.adjusted,mean.factorY_seq,col=col,lty=lty)}
if(!(method=="dantzig" & p>1000))            arrows(x0=rho_seq.adjusted,y0=pmax(mean.factorY_seq-sd.factorY_seq,0),x1=rho_seq.adjusted,y1=mean.factorY_seq+sd.factorY_seq,length=0.05,angle=90,code=3,col=col)
           mtext(side=2, text="RMSE", line=2,cex=1)
           mtext(side=1, text=expression(rho), line=2.5,cex=1)
           mtext(text=LETTERS[figure.count],side=3,adj=0,cex=1)

           par(new=TRUE)
       }   
       
       if(figure.count %in% c(4,8,12,16)) {mtext(side=4,text=sprintf("SNR=%s",SNR),line=0.5,cex=1,font=2,las=1)}
       if(figure.count %in% 1:4) {mtext(side=3,text=bquote(bold('s'[0]^'B' == .(s0_B))),line=1,cex=1,font=2,las=1)}
       par(new=FALSE);axis(las=2,at=rho_seq,side=1); figure.count=figure.count+1
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
          lty=lty_seq[match(setting,setting_seq)]
          if(p<=1000){method_seq=method_seq.all;colors=colors.all}
          if(p>1000){method_seq=method_seq.partial;colors=colors.partial}
          main=bquote(bold('s'[0]^'B' == .(s0_B) ~ ',' ~ 'p'^'B'== .(pB)))

 
           for(method in method_seq)
       {
                                                      mean.factorY_seq=c();sd.factorY_seq=c();col=colors[match(method,method_seq)]


                                                             for(rho in rho_seq)
                                                        {
                                                      design=ifelse(rho==0,"Independence","Pairwise");
                                                      setting=sprintf("n=%s;p=%s;s0=%s;SNR=%s;Design=%s;rho=%s;pB=%s;s0_B=%s",n,p,s0,SNR,design,rho,ifelse(rho==0,0,pB),ifelse(rho==0,0,s0_B))
                                                      setwd(sorted_address);setwd(setting)
if(!(method=="dantzig" & p>1000))                                                      factorY=get(load(sprintf("%s_%s_%s.Rdata",metric,method,setting)))
if(!(method=="dantzig" & p>1000))                                                      mean.factorY=mean(factorY);sd.factorY=sd(factorY);mean.factorY_seq=c(mean.factorY_seq,mean.factorY);sd.factorY_seq=c(sd.factorY_seq,sd.factorY)

                                                        }    
 
           if(SNR==SNR_seq[1]){pch=pch_seq[1];lty=lty_seq[1]};if(SNR==SNR_seq[2]){pch=pch_seq[3];lty=lty_seq[3]}
           rho_seq.adjusted=rho_seq+(-0.02+0.005*match(method,method_seq))


           if(length(grep("pauc",metric))>0) {ylab="pAUC"}; if(length(grep("rmse",metric))>0) {ylab="RMSE"}; if(length(grep("smse",metric))>0) {ylab="SMSE"}
if(!(method=="dantzig" & p>1000))            {plot(rho_seq.adjusted,mean.factorY_seq,col=col,main="",ylim=YLIM,pch=pch,xlab="",ylab="",xaxt="n",xlim=range(rho_seq));lines(rho_seq.adjusted,mean.factorY_seq,col=col,lty=lty)}
if(!(method=="dantzig" & p>1000))            arrows(x0=rho_seq.adjusted,y0=pmax(mean.factorY_seq-sd.factorY_seq,0),x1=rho_seq.adjusted,y1=mean.factorY_seq+sd.factorY_seq,length=0.05,angle=90,code=3,col=col)
           mtext(side=2, text="RMSE", line=2,cex=1)
           mtext(side=1, text=expression(rho), line=2.5,cex=1)
           mtext(text=LETTERS[figure.count],side=3,adj=0,cex=1)

           par(new=TRUE)
       }   
       if(figure.count %in% c(4,8,12,16)) {mtext(side=4,text=sprintf("SNR=%s",SNR),line=0.5,cex=1,font=2,las=1)}
       if(figure.count %in% 1:4) {mtext(side=3,text=bquote(bold('s'[0]^'B' == .(s0_B))),line=1,cex=1,font=2,las=1)}
       par(new=FALSE);axis(las=2,at=rho_seq,side=1); figure.count=figure.count+1
    }
 }
}




mtext(text=bquote(bold('p'^'B' == 10)),at=0.27,side=3,outer=T,cex=1,font=2,line=0.2)  
mtext(text=bquote(bold('p'^'B' == 100)),at=0.77,side=3,outer=T,cex=1,font=2,line=0.2)  

mtext(gsub(";","\n",setting_seq[1]),at=0.75,side=2,outer=T,cex=1,font=2,las=2,line=-1.2)  
mtext(gsub(";","\n",setting_seq[2]),at=0.25,side=2,outer=T,cex=1,font=2,las=2,line=-1.2)  


par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE); plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend("bottom",c("Lasso","HENet","Ridge","Dantzig","SCAD","Stability"),xpd=TRUE,horiz=TRUE,inset=c(0,0),lty=rep(1,4),bty="n",col=colors.all,cex=1,lwd=3)
