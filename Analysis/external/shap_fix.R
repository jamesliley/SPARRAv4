### Fix Shapley value plots: currently no title
## Run this from the directory containing subfolder 'Analysis'

library(dplyr)
source("SPARRAv4/auxiliary.R")


old_dir="Analysis/full_model/Shapley_values/All_old/"
new_dir="Analysis/full_model/Shapley_values/All/"

shap_names=paste0(new_dir,gsub(".pdf","",list.files(old_dir,pattern="*.pdf",full.names=FALSE)))
shap_var=gsub("_large_range","",gsub("_small_range","",gsub("shapley_","",basename(shap_names))))
sub=which(!is.na(longvarnames(shap_var)))
shap_names=shap_names[sub]

# Manual removals for some difficult inclusions
shap_names=setdiff(shap_names,c(
  "Analysis/full_model/Shapley_values/All/shapley_days_since_last_SMR04_large_range",
  "Analysis/full_model/Shapley_values/All/shapley_elective_bed_days_large_range",
  "Analysis/full_model/Shapley_values/All/shapley_emergency_bed_days_small_range",
  "Analysis/full_model/Shapley_values/All/shapley_num_ae2_attendances_small_range"))


large_range=c(-0.015,0.15)
small_range=c(-0.022,0.065)


# Get labels, because I forgot to do this
library(pdftools)
labs=list()
for (i in 1:length(shap_names)) {
  p1=pdf_text(paste0(gsub("All","All_old",shap_names[i]),".pdf"))
  p2=unlist(strsplit(p1,"\n")); p2=p2[which(!(p2==""))]; 
  p2=trimws(gsub("●","",p2)); p2=setdiff(p2,c("","Shapley value","SIMD decile"))
  if (grepl("small",basename(shap_names[i]))) 
    yv=c("0.06","0.04","0.02","0.00","−0.02") else
      yv=c("0.00","0.05","0.10","0.15")
  p2=setdiff(p2,yv)
  if (any(c("No","Yes") %in% p2)) p2=c("No","Yes") else p2=p2[which(nchar(p2)<15)]
  p2=gsub("−","-",p2)
  labs[[i]]=p2
  names(labs)[length(labs)]=basename(shap_names[i])
}

for (i in 1:length(shap_names)) {
  r1=readLines(paste0(gsub("All","All_old",shap_names[i]),".txt"))
  w=which(r1=="***************")
  r1s=paste0(r1[(w+1):length(r1)], collapse="\n")
  eval(parse(text=r1s))
  
  if (grepl("small",basename(shap_names[i]))) yvar=small_range else yvar=large_range
  varname=gsub("_small_range","",
               gsub("_large_range","",
                    gsub("shapley_","",basename(shap_names[i]))))
  pdf(paste0(shap_names[i],".pdf"),width=7,height=6)
  par(mar=c(7,4.1,4.1,2.1))
  xmain=longvarnames(varname)
  if (nchar(xmain)>50) {
    gsp=which(unlist(strsplit(xmain,""))==" "); gbr=gsp[round(length(gsp)/2)]
    xmain=paste0(substring(xmain,1,gbr),"\n",substring(xmain,gbr+1,nchar(xmain)))
  }
  plot(0,type="n",xlim=range(x0),ylim=yvar,main=xmain,xaxt="n",
       xlab="",ylab="Shapley value") # note fixed ylim
  lines(x0,mx,lty=2,col="red")
  segments(x0,mx-sdx,x0,mx+sdx,col="red")
  points(x0,mx,pch=16,col="red")
  axis(1,at=x0,labels=labs[[i]],las=2)
  abline(h=0,lty=2)
  dev.off()
 
  # Copy over .txt
  com=paste0("cp ",gsub("All","All_old",shap_names[i]),".txt ",shap_names[i],".txt")
  system(com)
}


