par(mar=c(3,4,2.5,1))
par(oma=c(2.5,4,2,1))
par(mgp=c(2,1,0))  
setwd(sorted_address)        ######sorted_address is the folder for results of all settings######
all_folders=list.files()


setting_seq=c("pB=10","pB=100")
n_seq=c(300,100);s0_seq=c(10,40);SNR_seq=c(4,1)
p_seq=c(500,1000,2000,4000)
rho=0.7; s0_B=2
SNR_seq=c(4,1)


method_seq_all=c("lasso","lenet","henet","ridge","dantzig","scad","stability")
library(RColorBrewer);colors_all=(brewer.pal(length(method_seq_all),"Set1"));colors_all[1]="red";colors_all[5]="black";colors_all[6]="gold3";colors_all[7]="brown"
colors_all=colors_all[c(1,4,3,6,2,5,7)]
colors=colors_all

method_seq.all=method_seq_all[-2]
colors.all=colors_all[-2]
method_seq.partial=method_seq_all[-c(2,5)]
colors.partial=colors_all[-c(2,5)]

metric="pauc"; Independence="Independence";YLIM=c(0,1)
pch_seq=c(0,2,5);lty_seq=c(5,3,1)
figure.count=1
par(mfrow=c(4,4))


method_seq=method_seq.all
colors=colors.all

 for(n in n_seq)
{
    for(SNR in SNR_seq)
 {
       for(s0 in s0_seq) 
    {
          setting=setting_seq[1]
          eval(parse(text=setting))
          main=bquote(bold('s'[0]^'B' == .(s0_B) ~ ',' ~ 'p'^'B'== .(pB)))

 
           for(method in method_seq)
       {
                                                      mean.factorY_seq=c();sd.factorY_seq=c();col=colors[match(method,method_seq)]


                                                             for(p in p_seq)
                                                        {
                                                      design=ifelse(rho==0,"Independence","Pairwise");
                                                      setting=sprintf("n=%s;p=%s;s0=%s;SNR=%s;Design=%s;rho=%s;pB=%s;s0_B=%s",n,p,s0,SNR,design,rho,ifelse(rho==0,0,pB),ifelse(rho==0,0,s0_B))
                                                      setwd(sorted_address);setwd(setting)
if(!(method=="dantzig" & p>1000))                     factorY=get(load(sprintf("%s_%s_%s.Rdata",metric,method,setting)))
if(method=="dantzig" & p>1000)                        factorY=rep(NA,64)
                                                      mean.factorY=mean(factorY);sd.factorY=sd(factorY);mean.factorY_seq=c(mean.factorY_seq,mean.factorY);sd.factorY_seq=c(sd.factorY_seq,sd.factorY)

                                                        }    
 
           if(SNR==SNR_seq[1]){pch=pch_seq[1];lty=lty_seq[1]};if(SNR==SNR_seq[2]){pch=pch_seq[3];lty=lty_seq[3]}
           p_seq.adjusted=p_seq+(-30+10*match(method,method_seq))

plot(p_seq.adjusted,mean.factorY_seq,col=col,main="",ylim=YLIM,pch=pch,xlab="",ylab="",xaxt="n",xlim=range(p_seq));lines(p_seq.adjusted,mean.factorY_seq,col=col)
arrows(x0=p_seq.adjusted,y0=pmax(mean.factorY_seq-sd.factorY_seq,0),x1=p_seq.adjusted,y1=mean.factorY_seq+sd.factorY_seq,length=0.05,angle=90,code=3,col=col)


           mtext(side=2, text="pAUC", line=2,cex=1)
           mtext(side=1, text='p', line=3,cex=1)
           mtext(text=LETTERS[figure.count],side=3,adj=0,cex=1)

           par(new=TRUE)
       }   
if(figure.count %in% c(1,5,9,13)) {mtext(side=2,text=sprintf("n=%s",n),line=4,cex=1,font=2,las=1)}
if(figure.count %in% 1:4) {mtext(side=3,text=sprintf("s0=%s",s0),line=1.5,cex=1,font=2,las=1)}
par(new=FALSE);axis(las=2,at=p_seq,side=1); figure.count=figure.count+1
    }
 }
}





























 for(n in n_seq)
{
    for(SNR in SNR_seq)
 {
       for(s0 in s0_seq) 
    {
          setting=setting_seq[2]
          eval(parse(text=setting)) 
          main=bquote(bold('s'[0]^'B' == .(s0_B) ~ ',' ~ 'p'^'B'== .(pB)))

 
           for(method in method_seq)
       {
                                                      mean.factorY_seq=c();sd.factorY_seq=c();col=colors[match(method,method_seq)]


                                                             for(p in p_seq)
                                                        {
                                                      design=ifelse(rho==0,"Independence","Pairwise");
                                                      setting=sprintf("n=%s;p=%s;s0=%s;SNR=%s;Design=%s;rho=%s;pB=%s;s0_B=%s",n,p,s0,SNR,design,rho,ifelse(rho==0,0,pB),ifelse(rho==0,0,s0_B))
      if(!is.na(match(setting,all_folders)))
  {
setwd(sorted_address);setwd(setting)
if(!(method=="dantzig" & p>1000))                     factorY=get(load(sprintf("%s_%s_%s.Rdata",metric,method,setting)))
if(method=="dantzig" & p>1000)                        factorY=rep(NA,64)
                                                      mean.factorY=mean(factorY);sd.factorY=sd(factorY);mean.factorY_seq=c(mean.factorY_seq,mean.factorY);sd.factorY_seq=c(sd.factorY_seq,sd.factorY)
  }
      if(is.na(match(setting,all_folders)))
  {mean.factorY_seq=c(mean.factorY_seq,NA);sd.factorY_seq=c(sd.factorY_seq,NA)
}

                                                        }    
 
           if(SNR==SNR_seq[1]){pch=pch_seq[1];lty=lty_seq[1]};if(SNR==SNR_seq[2]){pch=pch_seq[3];lty=lty_seq[3]}
           p_seq.adjusted=p_seq+(-30+10*match(method,method_seq))
plot(p_seq.adjusted,mean.factorY_seq,col=col,main="",ylim=YLIM,pch=pch,xlab="",ylab="",xaxt="n",xlim=range(p_seq));lines(p_seq.adjusted,mean.factorY_seq,col=col)
arrows(x0=p_seq.adjusted,y0=pmax(mean.factorY_seq-sd.factorY_seq,0),x1=p_seq.adjusted,y1=mean.factorY_seq+sd.factorY_seq,length=0.05,angle=90,code=3,col=col)

           mtext(side=2, text="pAUC", line=2,cex=1)
           mtext(side=1, text="p", line=3,cex=1)
           mtext(text=LETTERS[figure.count],side=3,adj=0,cex=1)
           par(new=TRUE)
       }   

if(figure.count %in% c(1,5,9,13)) {mtext(side=2,text=sprintf("n=%s",n),line=4,cex=1,font=2,las=1)}
if(figure.count %in% 1:4) {mtext(side=3,text=sprintf("s0=%s",s0),line=1.5,cex=1,font=2,las=1)}
par(new=FALSE);axis(las=2,at=p_seq,side=1); figure.count=figure.count+1
    }
 }
}

mtext(bquote(bold('p'^'B'== 10)),at=0.75,side=2,outer=T,cex=1,font=2,las=2,line=-0.5)  
mtext(bquote(bold('p'^'B'== 100)),at=0.25,side=2,outer=T,cex=1,font=2,las=2,line=-1.2)  
mtext("SNR=4",at=0.27,side=3,outer=T,cex=1,font=2)  
mtext("SNR=1",at=0.77,side=3,outer=T,cex=1,font=2)  


all_methods_text=paste(method_seq,collapse=",") 
main=sprintf("Pairwise correlation,n=%s,p=%s,s0=%s,SNR=%s",n,p,s0,SNR) 
MAIN=main
par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE); plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend("bottom",c("Lasso","HENet","Ridge","Dantzig","SCAD","Stability"),xpd=TRUE,horiz=TRUE,inset=c(0,0),lty=rep(1,4),bty="n",col=colors.all,cex=1,lwd=3)
