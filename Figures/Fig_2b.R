library("tidyverse")
library("patchwork")

source("Figures/util.R") # for import_sparra_expr()
source("SPARRAv4/auxiliary.R") # for roc_2panel()



eval(import_sparra_expr("Analysis/full_model/Performance/v4_final/prc_max.txt"))



# Old Figure 2(b)
auc3=signif(xprc3$auc[4],digits=3); se3=signif(xprc3$se[4],digits=2);
auc4=signif(xprc4$auc[4],digits=3); se4=signif(xprc4$se[4],digits=2);
aucm=signif(xprcm$auc[4],digits=3); sem=signif(xprcm$se[4],digits=2);
labs=c(paste0("v3: ",auc3," (",se3,")"),paste0("v4: ",auc4," (",se4,")"),
       paste0("Max: ",aucm," (",sem,")"))
prc_2panel(list(xprc3b,xprc4b,xprcmb),labels=labs,col=c("blue","red","black"),
           title="AUPRC (SE)",title.col="black",text.col=c("blue","red","black"))



# New Figure 2(b)
df <- rbind(
  cbind(Model = "v3", data.frame(sens = xprc3b$sens[1,],
                                 ppv = xprc3b$ppv[1,])),
  cbind(Model = "v4", data.frame(sens = xprc4b$sens[1,],
                                 ppv = xprc4b$ppv[1,])),
  cbind(Model = "Max", data.frame(sens = xprcmb$sens[1,],
                                  ppv = xprcmb$ppv[1,]))
) |> mutate(Model = fct_relevel(Model, "v3", "v4", "Max"))
df2 <- build_diff(df, sens, ppv)

p1 <- ggplot(df) +
  geom_line(aes(x = sens, y = ppv, col = Model), linewidth = 0.4) +
  xlim(0, 1) + ylim(0, 1) +
  xlab("") + ylab("Precision") +
  theme_minimal(base_size = 8) + theme(legend.justification = c(1,1),
                                       legend.position = c(1,1),
                                       legend.spacing = unit(0, "npc"),
                                       legend.margin = unit(0, "npc"),
                                       legend.background = element_rect(fill = "white", size = 0, colour = "white"))

p2 <- ggplot(df2 |>
               mutate(v4 = v4 - v3,
                      Max = Max - v3,
                      v3 = v3 - v3) |>
               pivot_longer(v3:Max, names_to = "Model", values_to = "delta_sens") |>
               mutate(Model = fct_relevel(Model, "v3", "v4", "Max"))) +
  geom_line(aes(x = sens, y = delta_sens, col = Model), linewidth = 0.4) +
  xlim(0, 1) +
  xlab("Recall") + ylab("Δ Precision") +
  theme_minimal(base_size = 8) + theme(legend.position = "none")

p <- p1 / p2 + plot_layout(heights = c(3,1))

print(p)

ggsave("Figures/pdfs/Fig_2b.pdf",
       width = 7.5, height = 7.5, units = "cm",
       device = cairo_pdf)
