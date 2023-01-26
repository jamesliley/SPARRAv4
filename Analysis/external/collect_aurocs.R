## Collect AUROCs and AUPRCs from text files into a table
## Run this from the directory containing subfolder 'Analysis'

## Text files with AUROCs 
roc_list=list.files("Analysis/full_model/Performance/constituent_predictors/",
  full.names=TRUE,pattern="^roc_pred.*txt$")
roc_final="Analysis/full_model/Performance/v4_final/roc_super.txt"

## Text files with AUPRCs 
prc_list=list.files("Analysis/full_model/Performance/constituent_predictors/",
  full.names=TRUE,pattern="^prc_pred.*txt$")
prc_final="Analysis/full_model/Performance/v4_final/prc_super.txt"


## AUROCs
aurocs=c(); roc_se=c()
for (i in 1:length(roc_list)) {
  rx=readLines(roc_list[i])
  w0=which(rx=="***************")  
  com=paste(rx[(w0+1):length(rx)],collapse="\n")
  eval(parse(text=com))
  aurocs=rbind(aurocs,xroc$auc); roc_se=rbind(roc_se,xroc$se)
}  
rxf=readLines(roc_final)
w0=which(rxf=="***************")  
eval(parse(text=paste(rxf[(w0+1):length(rxf)],collapse="\n")))
aurocs=rbind(aurocs,xroc$auc); roc_se=rbind(roc_se,xroc$se)

rnames=c(gsub("roc_pred.","",gsub(".txt","",basename(roc_list))),"final")
cnames=c("fold1","fold2","fold3","avg")
rownames(aurocs)=rnames; rownames(roc_se)=rnames
colnames(aurocs)=cnames; colnames(roc_se)=cnames
aurocs=data.frame(aurocs)


## AUPRCs
auprcs=c(); prc_se=c()
for (i in 1:length(prc_list)) {
  rx=readLines(prc_list[i])
  w0=which(rx=="***************")  
  com=paste(rx[(w0+1):length(rx)],collapse="\n")
  eval(parse(text=com))
  auprcs=rbind(auprcs,xprc$auc); prc_se=rbind(prc_se,xprc$se)
}  
rxf=readLines(prc_final)
w0=which(rxf=="***************")  
eval(parse(text=paste(rxf[(w0+1):length(rxf)],collapse="\n")))
auprcs=rbind(auprcs,xprc$auc); prc_se=rbind(prc_se,xprc$se)

rnames=c(gsub("prc_pred.","",gsub(".txt","",basename(prc_list))),"final")
cnames=c("fold1","fold2","fold3","avg")
rownames(auprcs)=rnames; rownames(prc_se)=rnames
colnames(auprcs)=cnames; colnames(prc_se)=cnames
auprcs=data.frame(auprcs)


## Add coefficient in ensemble
efile="Analysis/full_model/Description/misc.txt"
rxx=readLines(efile)
w123=which(rxx=="Lambda_min")
coef=c(); 
for (i in 1:3) 
  coef=rbind(coef,as.numeric(unlist(strsplit(rxx[w123[i]+2]," "))))
mods=unlist(strsplit(rxx[w123[1]+1]," ")); mods=mods[which(!(mods==""))]
colnames(coef)=mods
rownames(coef)=c("fold1","fold2","fold3")
coef=data.frame(t(coef))

## Format table
ftab=cbind(
  aurocs$fold1,auprcs$fold1,c(coef$fold1,NA),
  aurocs$fold2,auprcs$fold2,c(coef$fold2,NA),
  aurocs$fold3,auprcs$fold3,c(coef$fold3,NA),
  aurocs$avg,auprcs$avg
)
colnames(ftab)=c(rep(c("AUROC","AUPRC","Coef."),3),
  c("AUROC","AUPRC"))
signif(ftab,digits=4)
rownames(ftab)=c("ANN","GLM",  "Naive Bayes",  "RF, max. 20",  
  "RF, max. 40",  "SPARRA v3",  "XG-boost, max. 4",  
  "XG-boost, max. 8",  "XG-boost, max. 3",  "Ensemble")


## Write table (csv)
write.csv(signif(ftab,digits=4),file="Analysis/external/auc.txt",quote=F)

## Maximum standard error for AUROCs
max(roc_se)

## Maximum standard error for AUORCs
max(prc_se)
