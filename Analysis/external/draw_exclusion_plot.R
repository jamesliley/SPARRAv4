## Draw Venn diagram of exclusions
# Should be in directory with folders 'Analytics', 'Description' etc. Usually called 'full_model'

library(VennDiagram)

ex=readLines("Description/exclusions.txt")

# A, B, C, D
wA=which(ex=="A: Number of records for which individual died before time cutoff")
all=c(as.numeric(ex[wA+3*(0:13)+1]),0)
names(all)=c("A","B","C","D","AB","AC","AD","BC","BD","CD","ABC","ABD","ACD","BCD","ABCD")
for (i in 1:length(all)) assign(names(all)[i],all[i])

pdf("../../Diagrams/exclusions.pdf",width=4,height=4)
grid.newpage()
draw.quad.venn(area1=A, area2=B, area3=C, 
               area4 =D, n12=AB, n23=BC, n13=AC, 
               n14= AD,n24=BD, n34=CD, n123=ABC, 
               n124=ABD, n234=BCD, n134=ACD, n1234=ABCD, 
               category=c("Died","No v3","No SIMD","Unmatched"),
               col="Green",fill=c("Red","Pink","Blue","Orange"),lty="dashed")
dev.off()