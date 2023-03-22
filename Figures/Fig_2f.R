library("tidyverse")
library("patchwork")
library("glue")

source("Figures/util.R") # for import_sparra_expr()
source("SPARRAv4/auxiliary.R") # for roc_2panel()



eval(import_sparra_expr("Analysis/full_model/Performance/v3_vs_v4/v3_vs_v4_ntreat.txt"))



# Old Figure 2(f)
plot(nadm,s_v3-s_super,col="red", #main="Number of patients needed to treat to target avoidable admissions",
     xlab="To target this many avoidable admissions",ylab="Treat this many fewer patients using v4",pch=16)



# New Figure 2(f)
p <- ggplot(data.frame(x = nadm, y = s_v3-s_super)) +
  geom_point(aes(x = x, y = y)) +
  xlab("To target this many avoidable admissions") + ylab("Treat this many fewer\npatients using v4") +
  theme_minimal(base_size = 8)

print(p)

ggsave("Figures/pdfs/Fig_2f.pdf",
       width = 7.5, height = 4.2, units = "cm",
       device = cairo_pdf)
