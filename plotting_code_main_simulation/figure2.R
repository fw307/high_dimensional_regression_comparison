par(mar=c(4,4,3,1))
par(oma=c(2,4,2,1))
setwd(sorted_address) ######sorted_address is the folder for results of all settings######
all_folders=list.files()


n_seq=c(300,200,100);s0_seq=c(10,40);SNR_seq=c(4,1)
p_seq=c(500,1000,2000,4000)
method_seq=c("lasso","lenet","henet","ridge","dantzig","scad","stability")
library(RColorBrewer);colors_all=(brewer.pal(length(method_seq),"Set1"));colors_all[1]="red";colors_all[5]="black";colors_all[6]="gold3";colors_all[7]="brown"
colors_all=colors_all[c(1,4,3,6,2,5,7)]
colors=colors_all
method_seq=method_seq[-c(2,5)]
colors=colors[-c(2,5)]
metric="pauc"; Independence="Independence";YLIM=c(0,1)
flag=0;
pch_seq=c(0,2,5);lty_seq=c(1,3,5)
figure.count=1




par(mfrow=c(3,4))

    for(n in n_seq)
 {
       for(SNR in SNR_seq) 
    {

           for(s0 in s0_seq)
       {
             main=sprintf("n=%s,SNR=%s,s0=%s",n,SNR,s0) 
                                                             for(method in method_seq)
                                                        {
                                                      mean.factorY_seq=c();sd.factorY_seq=c();col=colors[match(method,method_seq)] 
                                                                      for(p in p_seq) 
                                                                 { 
                                                      setting=sprintf("n=%s;p=%s;s0=%s;SNR=%s;Design=Independence;rho=%s;pB=%s;s0_B=%s",n,p,s0,SNR,0,0,0)
                                                      setwd(sorted_address);setwd(setting)
                                                      factorY=get(load(sprintf("%s_%s_%s.Rdata",metric,method,setting)))
                                                      mean.factorY=mean(factorY);sd.factorY=sd(factorY);mean.factorY_seq=c(mean.factorY_seq,mean.factorY);sd.factorY_seq=c(sd.factorY_seq,sd.factorY)
                                                                 }
           lty=lty_seq[1];pch=pch_seq[1]
           p_seq.adjusted=p_seq+(-30+10*match(method,method_seq))
           if(length(grep("rmse",metric))>0){ylab="RMSE"};  if(length(grep("smse",metric))>0){ylab="SMSE"}; if(length(grep("pauc",metric))>0){ylab="pAUC"}; 
           plot(p_seq.adjusted,mean.factorY_seq,col=col,main="",ylim=YLIM,pch=pch,xlab="",ylab="",xaxt="n",xlim=c(500,4000));lines(p_seq.adjusted,mean.factorY_seq,col=col,lty=lty)
           arrows(x0=p_seq.adjusted,y0=pmax(mean.factorY_seq-sd.factorY_seq,0),x1=p_seq.adjusted,y1=mean.factorY_seq+sd.factorY_seq,length=0.05,angle=90,code=3,col=col)

mtext(side=1, text="p", line=3,cex=0.8)
mtext(side=2, text=ylab, line=2,cex=1)
mtext(text=LETTERS[figure.count],side=3,adj=0,cex=1)
if(figure.count %in% c(1,5,9)) {mtext(side=2,text=sprintf("n=%s",n),line=4,cex=1,font=2,las=1)}
if(figure.count %in% 1:4) {mtext(side=3,text=sprintf("s0=%s",s0),line=1.5,cex=1,font=2,las=1)}

           par(new=TRUE)                                }
           par(new=FALSE);figure.count=figure.count+1  
           axis(las=2,at=p_seq,side=1)
       }
       par(new=FALSE)
          
    }     
 }  


mtext("SNR=4",at=0.27,side=3,outer=T,cex=1,font=2)  
mtext("SNR=1",at=0.77,side=3,outer=T,cex=1,font=2)  
all_methods_text=paste(method_seq,collapse=",") 
main=sprintf("Independence design") 
MAIN=main
par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE); plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend("bottom",method_seq,xpd=TRUE,horiz=TRUE,inset=c(0,0),lty=rep(1,4),bty="n",col=colors,cex=1,lwd=3)



