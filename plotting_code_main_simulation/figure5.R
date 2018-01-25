setwd(sorted_address)        ######sorted_address is the folder for results of all settings######
all_folders=list.files(); Independence="Independence"; Pairwise="Pairwise"; Toeplitz="Toeplitz"
ALL_METHODS=c("lasso","henet","ridge","scad","stability")

target.snr=2
LIM.rmse=c(5,23)


independence.list=c()
pairwise.list=c()
toeplitz.list=c()





                          for(folder in all_folders) 
{ 
  eval(parse(text=folder))
  if(s0_B==0 & SNR==target.snr) {independence.list=c(independence.list,folder)}
  if(Design=="Pairwise" & pB==100 & s0_B==2 & rho==0.7 & SNR==target.snr) {pairwise.list=c(pairwise.list,folder)}
  if(Design=="Toeplitz" & pB==100 & s0_B==2 & rho==0.7 & SNR==target.snr) {toeplitz.list=c(toeplitz.list,folder)}
}
metric_seq=c("pauc","rmse","tpr","ppv","mcc"); 






r_seq.independence=c()

                         for(folder in independence.list)
{
setwd(sorted_address); setwd(folder); 
eval(parse(text=folder)); r_seq.independence=c(r_seq.independence,n/(s0*log(p-s0))); 
            for(metric in metric_seq)
      {
            if(metric=="pauc") {all_methods=ALL_METHODS} 
            if(metric=="rmse") {all_methods=ALL_METHODS[-5]}
            for(method in ALL_METHODS) {assign(sprintf("%s.%s_n=%s,p=%s,s0=%s,Design=Independence",method,metric,n,p,s0),mean(get(load(sprintf("%s_%s_%s.Rdata",metric,method,folder)))))}
      }

}


r_seq.toeplitz=c()

                         for(folder in toeplitz.list)
{
setwd(sorted_address); setwd(folder); 
eval(parse(text=folder)); r_seq.toeplitz=c(r_seq.toeplitz,n/(s0*log(p-s0))); 
            for(metric in metric_seq)
      {
            if(metric=="pauc") {all_methods=ALL_METHODS} 
            if(metric=="rmse") {all_methods=ALL_METHODS[-5]}
            for(method in ALL_METHODS) {assign(sprintf("%s.%s_n=%s,p=%s,s0=%s,Design=Toeplitz",method,metric,n,p,s0),mean(get(load(sprintf("%s_%s_%s.Rdata",metric,method,folder)))))}
      }

}






r_seq.pairwise=c()
factor.array=c();

    for(folder in pairwise.list)
{
setwd(sorted_address); setwd(folder); 
eval(parse(text=folder))
factor.array=c(factor.array,sprintf("n=%s;p=%s;s0=%s",n,p,s0))
}
identity.array=unique(factor.array)                  


    for(sub_setting in identity.array)
  {
    eval(parse(text=sub_setting)); folders=pairwise.list[grep(sprintf("n=%s;p=%s;s0=%s",n,p,s0),pairwise.list)]
    r_seq.pairwise=c(r_seq.pairwise,n/(s0*log(p-s0))); 
    

            for(metric in metric_seq)
      {
            if(metric=="pauc") {all_methods=ALL_METHODS} 
            if(metric=="rmse") {all_methods=ALL_METHODS[-5]}
            for(method in ALL_METHODS) 
        {
            result.vec=c()
            for(folder in folders) {setwd(sorted_address); setwd(folder); result.vec=c(result.vec,get(load(sprintf("%s_%s_%s.Rdata",metric,method,folder))))}
            assign(sprintf("%s.%s_n=%s,p=%s,s0=%s,Design=Pairwise",method,metric,n,p,s0),mean(result.vec))

        }
      }
  }

identity.array.pairwise=identity.array 







method_seq_all=c("lasso","lenet","henet","ridge","dantzig","scad","stability")
library(RColorBrewer);colors_all=(brewer.pal(length(method_seq_all),"Set1"));colors_all[1]="red";colors_all[5]="black";colors_all[6]="gold3";colors_all[7]="brown"
colors_all=colors_all[c(1,4,3,6,2,5,7)]
colors=colors_all
method_seq_all=method_seq_all[-2]
all_methods=method_seq_all
colors=colors[-2]
all_designs=c("Independence","Pairwise","Toeplitz")
par(mfrow=c(2,2))
metric_seq=c("pauc","rmse","tpr","ppv")





par(mar=c(3,3,1,1))
par(oma=c(3.6,1,1,1))

                              for(indicator in 1:4)                     
{
design1="Independence"; design2="Pairwise"
metric=metric_seq[indicator]; all_methods=ALL_METHODS
xlim=ylim=c(0,1); if(metric=="rmse") {all_methods=ALL_METHODS[-5]; xlim=ylim=LIM.rmse}
main.text=toupper(metric); if(metric=="pauc") {main.text="pAUC"}; if(metric=="rmse") {main.text="RMSE"}
main=main.text
list.settings=identity.array.pairwise; if(design1=="Toeplitz" | design2=="Toeplitz") {list.settings=toeplitz.list}

                            par(mgp=c(2,1,0))  
                            for(method in all_methods) 
                              {
                                 vec.design1=c(); vec.design2=c()
                               for(folder in list.settings) 
                                {
                                 eval(parse(text=folder))
                                 vec.design1=c(vec.design1,get(sprintf("%s.%s_n=%s,p=%s,s0=%s,Design=%s",method,metric,n,p,s0,design1)))
                                 vec.design2=c(vec.design2,get(sprintf("%s.%s_n=%s,p=%s,s0=%s,Design=%s",method,metric,n,p,s0,design2)))
                                }
                          plot(pch=ceiling(ifelse(r_seq.pairwise<=3,r_seq.pairwise,4)),vec.design2,vec.design1,col=colors[match(method,method_seq_all)],ylab="",xlab="",cex=0.75,xlim=xlim,ylim=ylim,main=main);abline(0,1)
xlab=design2
ylab=design1                         
mtext(side=1, text=xlab, line=2,cex=1)
mtext(side=2, text=ylab, line=2,cex=1)
mtext(text=LETTERS[indicator],side=3,adj=0,cex=1)

                          par(new=TRUE)
                              }



                        method="dantzig"; r.dantzig=c(); vec.design1=c(); vec.design2=c()
                        for(folder in list.settings) 
                       {
                            eval(parse(text=folder));
                            if(p<=1000)          
                          {
                            r.dantzig=c(r.dantzig,n/(s0*log(p-s0)));                       
                            folder.name=sprintf("n=%s;p=%s;s0=%s;SNR=%s;Design=%s;rho=%s;pB=%s;s0_B=%s",n,p,s0,target.snr,design1,ifelse(design1=="Independence",0,0.7),ifelse(design1=="Independence",0,100),ifelse(design1=="Independence",0,2))
                            setwd(sorted_address); setwd(folder.name)
                            vec.design1=c(vec.design1,mean(get(load(sprintf("%s_dantzig_%s.Rdata",metric,folder.name)))))

                            folder.name=sprintf("n=%s;p=%s;s0=%s;SNR=%s;Design=%s;rho=%s;pB=%s;s0_B=%s",n,p,s0,target.snr,design2,ifelse(design2=="Independence",0,0.7),ifelse(design2=="Independence",0,100),ifelse(design2=="Independence",0,2))
                            setwd(sorted_address); setwd(folder.name)
                            vec.design2=c(vec.design2,mean(get(load(sprintf("%s_dantzig_%s.Rdata",metric,folder.name)))))
                          }
                      }
                          plot(pch=ceiling(ifelse(r.dantzig<=3,r.dantzig,4)),vec.design2,vec.design1,col=colors[match(method,method_seq_all)],xlab="",ylab="",cex=0.75, xlim=xlim,ylim=ylim,main=main);abline(0,1) 
                          par(new=FALSE) 

}

par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE); plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend(legend=c("0<r<=1", "1<r<=2","2<r<=3","r>3"),xpd=TRUE,horiz=TRUE,inset=c(0,0),pch=1:4,bty="n",cex=1.5,x=-0.55,y=-0.85)
legend(legend=c("Lasso","HENet","Ridge","Dantzig","SCAD","Stability"),xpd=TRUE,horiz=TRUE,inset=c(0,0),col=colors,pch=15,bty="n",cex=1.5,x=-1,y=-0.95)
