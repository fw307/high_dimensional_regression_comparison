par(mar=c(4,4,4,1))
par(oma=c(4,4,4,1))
setwd(results)        ######results is the folder for results of all settings######
all_folders=list.files()

n_seq=200
s0_seq=c(10,40)
SNR_seq=c(1,2,4)
p_seq=c(500,1000,4000)
figure.count=1



SNR_seq=c(4,1)
pch_seq=c(0:2,4)
method_seq_all=c("lasso","lenet","henet","ridge","dantzig","scad","stability")
library(RColorBrewer);colors_all=(brewer.pal(length(method_seq_all),"Set1"));colors_all[1]="red";colors_all[5]="black";colors_all[6]="gold3";colors_all[7]="brown"
colors_all=colors_all[c(1,4,3,6,2,5,7)]
colors=colors_all


method_seq=method_seq_all[-c(2,4)]
colors=colors[-c(2,4)]


metric1="tpr";metric2="ppv"
flag=0
Independence="Independence"


par(mfrow=c(length(n_seq),4))

    for(n in n_seq)
 {
       for(SNR in SNR_seq) 
    {

                    
           for(s0 in s0_seq)
       {
           lty=1
         
                                                             for(method in method_seq)
                                                        {
                                                             col=colors[match(method,method_seq)];main=sprintf("%s: n=%s,s0=%s",LETTERS[figure.count],n,s0) 
                                                             factorY1.p=c(); factorY2.p=c()               
                                                                    for(p in p_seq)
                                                                {   
                                                                 pch=pch_seq[match(p,p_seq)]
                                                                 setting=sprintf("n=%s;p=%s;s0=%s;SNR=%s;Design=Independence;rho=%s;pB=%s;s0_B=%s",n,p,s0,SNR,0,0,0)
                                                                 setwd(sorted_address);setwd(setting)
                                                                 new.factorY1=get(load(sprintf("%s_%s_%s.Rdata",metric1,method,setting)))
                                                                 new.factorY2=get(load(sprintf("%s_%s_%s.Rdata",metric2,method,setting)))
                                                                 factorY1.p=c(factorY1.p,mean(new.factorY1)); 
                                                                 factorY2.p=c(factorY2.p,mean(new.factorY2))
                                                                 xlab=toupper(metric1); ylab=toupper(metric2)
                                                                 plot(mean(new.factorY1),mean(new.factorY2),col=col,pch=pch,xlab="",ylab="",main="",xlim=c(0,1),ylim=c(0,1))
mtext(side=1, text=xlab, line=3,cex=0.8)
mtext(side=2, text=ylab, line=2,cex=1)
mtext(text=LETTERS[figure.count],side=3,adj=0,cex=1)
if(figure.count %in% c(1,5,9)) {mtext(side=2,text=sprintf("n=%s",n),line=4,cex=1,font=2,las=1)}
if(figure.count %in% 1:4) {mtext(side=3,text=sprintf("s0=%s",s0),line=1.5,cex=1,font=2,las=1)}
                                                                 par(new=TRUE)
                                                                }
                                                                 lines(factorY1.p,factorY2.p,col=col,lty=lty)
                                                                 par(new=TRUE)
                                                        }      
par(new=FALSE); figure.count=figure.count+1
           
       }                            
       min.r=round(n/(s0*max(log(p_seq))),2); max.r=round(n/(s0*min(log(p_seq))),2)
       par(new=FALSE)
     
    }
 }
mtext("SNR=4",at=0.27,side=3,outer=T,cex=1,font=2)  
mtext("SNR=1",at=0.77,side=3,outer=T,cex=1,font=2)  
all_methods_text=paste(method_seq,collapse=",") 
main=sprintf("Independence design"); #MAIN=paste(main,all_methods_text,sep="\n")
MAIN=main
par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE); plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend("bottom",method_seq,xpd=TRUE,horiz=TRUE,inset=c(0,0),lty=rep(1,4),bty="n",col=colors,cex=1,lwd=3)
legend("bottomleft",legend=paste("p",p_seq,sep="="),xpd=TRUE,horiz=FALSE,inset=c(0,0),cex=1,pch=pch_seq)
