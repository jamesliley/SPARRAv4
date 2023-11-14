# Scripts
source("Figures/util.R") # for import_sparra_expr()
source("SPARRAv4/auxiliary.R") # for roc_2panel() etc

## Age/SIMD comparison by Shapley values ####

eval(import_sparra_expr("Analysis/full_model/Shapley_values/equivalent_years_older_shapley.txt"))

pdf("Figures/pdfs/Unsorted/equivalent_years_older_shapley.pdf",width=5,height=5)

plot(0,type="n",xlim=range(aget),ylim=range(age_d-aget),
     xlab="For individuals of this age",ylab="SIMD1 vs 10 equiv. to this many years older")
lines(aget,age_d-aget)
dev.off()


## Equivalent ages ####

eval(import_sparra_expr("Analysis/full_model/Shapley_values/age_simd_equivalent.txt"))
xx=seq(2,80);
xcol=colorRampPalette(c("blue","gray","red"))(10)
pdf("Figures/pdfs/Unsorted/age_simd_equivalent.pdf",width=5,height=5)
plot(0,type="n",xlab="Chron. age",ylab="Effective age",xlim=c(15,max(xx)),ylim=c(20,max(xx)))
abline(0,1,col="black",lty=2,lwd=2)
for (i in 1:10) lines(xx,yeq[,i],col=xcol[i],lwd=2)
legend("bottomright",c("1","...","10"),title="SIMD",lty=1,col=c(xcol[1],xcol[5],xcol[10]),lwd=2)
dev.off()



