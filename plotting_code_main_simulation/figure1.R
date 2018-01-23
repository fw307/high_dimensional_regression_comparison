method_seq_all=c("lasso","lenet","henet","ridge","dantzig","scad","stability")
library(RColorBrewer);colors_all=(brewer.pal(length(method_seq_all),"Set1"));colors_all[1]="red";colors_all[5]="black";colors_all[6]="gold3";colors_all[7]="brown"
colors_all=colors_all[c(1,4,3,6,2,5,7)]
figure.count=1

dev.off()
par(mar=c(3.5,3.5,1,1))
par(oma=c(2,1,1,1))
setwd(results)        ######results is the folder for results of all settings######
all_folders=list.files();Independence="Independence";target.snr=2
all_independence_folders=all_folders[grep("Independence",all_folders)]
independence_settings=all_independence_folders[grep(sprintf("SNR=%s",target.snr),all_independence_folders)]
r_seq=rep(0,length(independence_settings))
SNR_seq=rep(0,length(independence_settings))
for(i in 1:length(independence_settings)){setting=independence_settings[i];eval(parse(text=setting));r=n/(s0*log(p-s0));r_seq[i]=r;SNR_seq[i]=SNR}
unique.r.values=sort(unique(r_seq))

par(mfrow=c(2,2))
main=sprintf("Independence design,SNR=%s",target.snr);MAIN=main
setwd(sorted_address)
metric_seq=c("pauc","rmse","tpr","ppv")


     for(metric in metric_seq)
{
     if(length(grep("tpr",metric))>0){method_seq=c("lasso","lenet","henet","dantzig","scad","stability");colors=colors_all[-match("ridge",method_seq_all)];M=1;pch=1}
     if(length(grep("pauc",metric))>0) {method_seq=c("lasso","lenet","henet","ridge","dantzig","scad","stability");colors=colors_all;M=1;pch=1}
     if(length(grep("mse",metric))>0) {method_seq=c("lasso","lenet","henet","ridge","dantzig","scad");colors=colors_all[-match("stability",method_seq_all)];M=22;pch=1}
   
       for(method in method_seq)
  {
       metric_vec=rep(0,length(unique.r.values))
       col=colors[match(method,method_seq)]

        for(r in unique.r.values)
      { 
            individual.metric=c();settings=independence_settings[which(r_seq==r)]
            for(setting in settings)
          {
             setwd(setting)      
             individual.metric=c(individual.metric,mean(get(load(sprintf("%s_%s_%s.Rdata",metric,method,setting)))))
             setwd(sorted_address)
          }
       metric_vec[which(unique.r.values==r)]=mean(individual.metric)
     }
       order=order(unique.r.values)
       ylab=toupper(metric)
       if(length(grep("pauc",metric))>0) {ylab="pAUC"}; if(length(grep("rmse",metric))>0) {ylab="RMSE"}; if(length(grep("smse",metric))>0) {ylab="SMSE"}; if(length(grep("tpr",metric))>0) {ylab="TPR"}
      
       #plot(unique.r.values,metric_vec,xlim=c(0,5),ylim=c(0,M),xlab="r",ylab=ylab,main=LETTERS[figure.count],col=col,pch=pch)
plot(unique.r.values,metric_vec,xlim=c(0,5),ylim=c(0,M),xlab="",ylab="",main="",col=col,pch=pch)
mtext(side=1, text="r", line=2)
mtext(side=2, text=ylab, line=2)
mtext(text=LETTERS[figure.count],side=3,adj=0,cex=1.5)



       lines(unique.r.values[order],metric_vec[order],col=col)
       par(new=TRUE)
  }
       par(new=FALSE); figure.count=figure.count+1
}
method_all=c("lasso","lenet","henet","ridge","dantzig","scad","stability") 
colors=colors_all
par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE); plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend("bottom",method_all,xpd=TRUE,horiz=TRUE,inset=c(0,0),lty=rep(1,4),bty="n",col=colors,cex=0.8,lwd=3)
par(new=FALSE)
print(match(metric,metric_seq)==2)


