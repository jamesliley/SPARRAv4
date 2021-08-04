## Read in topic breakdown and output as table
## Run this from the directory containing subfolder 'Analysis'

## File with topic breakdown
source("Analysis/full_model/Shapley_values/topic_breakdown.txt")

## Sink to output file
sink("Analysis/external/topic_table.txt")


# Write as table
for (i in 1:30) {
  it=topic_summary1[[i]]$terms
  wb=it[which(substring(it,1,5)=="(BNF)")]
  wx=it[which(substring(it,1,7)=="(ICD10)")]
  wb=gsub("(BNF)  ","",wb,fixed=TRUE)
  wx=gsub("(ICD10) ","",wx,fixed=TRUE)
  wx=tolower(wx); 
  if (length(wx)>0) for (i in 1:length(wx)) 
    wx[i]=paste0(toupper(substring(wx[i],1,1)),substring(wx[i],2,nchar(wx[i])))
  lb=length(wb); lx=length(wx)
  wb=gsub("&","\\&",wb,fixed=TRUE)
  wx=gsub("&","\\&",wx,fixed=TRUE)
  if (lb>0) for (j in 1:lb) {
    outstr=paste0(wb[j]," & \\\\")
    cat(outstr)
    cat("\n")
  }
  if (lx>0) for (j in 1:lx) {
    outstr=paste0("\\emph{",wx[j],"} & \\\\")
    cat(outstr)
    cat("\n")
  }
  cat("\n & \\\\ \n")
}

sink()