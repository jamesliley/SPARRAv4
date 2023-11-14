library("tidyverse")
library("patchwork")
library("glue")

source("Figures/util.R") # for import_sparra_expr()
source("SPARRAv4/auxiliary.R") # for roc_2panel()



eval(import_sparra_expr("Analysis/full_model/Performance/v3_vs_v4/v3_vs_v4_num_in_topN.txt"))



# Old Figure 2(e)
par(mar=c(4,4,0.5,2))
plot(0,0,type="n",ylim=c(0,max(tp_super-tp_v3)),xlim=range(topn),
     xlab="Number of patients",ylab="Difference")
points(topn,tp_super-tp_v3,pch=16,col="black")
lines(loess.smooth(topn,abs(topn*meancut_v3-tp_v3)), col="blue")
lines(topn,abs(topn*meancut_super-tp_super), col="red")
if(FALSE) segments(topn,tp_super-tp_v3,topn,0*topn,lty=2,col="black")

legend("topleft",c("v3","v4"),lty=1,col=c("blue","red"))



# New Figure 2(e)
# NB: dropping the lines since they are to indicate calibration, but we have
#     calibration plots separately already
p <- ggplot(data.frame(x = topn, y = tp_super-tp_v3)) +
  geom_point(aes(x = x, y = y)) +
  xlab("Number of Patients") + ylab("Difference") +
  theme_minimal(base_size = 8)

print(p)

ggsave("Figures/pdfs/Unsorted/v3v4_difference_topN.pdf",
       width = 7.5, height = 4.2, units = "cm",
       device = cairo_pdf)
