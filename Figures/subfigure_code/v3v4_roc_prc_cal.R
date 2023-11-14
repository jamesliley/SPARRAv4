library("tidyverse")
library("patchwork")

source("Figures/util.R") # for import_sparra_expr()
source("SPARRAv4/auxiliary.R") # for roc_2panel()

eval(import_sparra_expr("Analysis/full_model/Performance/v4_final/roc_max.txt"))

# ROC ####
auc3=signif(xroc3$auc[4],digits=3); se3=signif(xroc3$se[4],digits=2);
auc4=signif(xroc4$auc[4],digits=3); se4=signif(xroc4$se[4],digits=2);
aucm=signif(xrocm$auc[4],digits=3); sem=signif(xrocm$se[4],digits=2);
labs=c(paste0("v3: ",auc3," (",se3,")"),paste0("v4: ",auc4," (",se4,")"),
       paste0("Max: ",aucm," (",sem,")"))
roc_2panel(list(xroc3b,xroc4b,xrocmb),labels=labs,
           col=c("blue","black","red"),xy_lty=2,xy_col="gray",
           text.col=c("blue","red","black"),title="AUROC (SE)",
           title.col="black")



# New Figure 2(a)
df <- rbind(
  cbind(Model = "v3", data.frame(sens = xroc3b$sens[1,],
                                 spec = xroc3b$spec[1,])),
  cbind(Model = "v4", data.frame(sens = xroc4b$sens[1,],
                                 spec = xroc4b$spec[1,])),
  cbind(Model = "Max", data.frame(sens = xrocmb$sens[1,],
                                  spec = xrocmb$spec[1,]))
) |> mutate(Model = fct_relevel(Model, "v3", "v4", "Max"))
df2 <- build_diff(df, spec, sens)

p1 <- ggplot(df |>
               mutate(spec = 1-spec) |>
               filter(Model != "Max")) +
  geom_line(aes(x = spec, y = sens, col = Model), linewidth = 0.4) +
  xlim(0, 1) + ylim(0, 1) +
  xlab("") + ylab("Sensitivity") +
  theme_minimal(base_size = 8) + theme(legend.justification = c(1,0),
                                       legend.position = c(1,0),
                                       legend.spacing = unit(0, "npc"),
                                       legend.margin = unit(0, "npc"),
                                       legend.background = element_rect(fill = "white", size = 0, colour = "white"))

p2 <- ggplot(df2 |>
               mutate(spec = 1-spec,
                      v4 = v4 - v3,
                      Max = Max - v3,
                      v3 = v3 - v3) |>
               pivot_longer(v3:Max, names_to = "Model", values_to = "delta_sens") |>
               mutate(Model = fct_relevel(Model, "v3", "v4", "Max")) |>
               filter(Model != "Max")) +
  geom_line(aes(x = spec, y = delta_sens, col = Model), linewidth = 0.4) +
  xlim(0, 1) +
  xlab("Recall") + ylab("Δ Sensitivity") +
  theme_minimal(base_size = 8) + theme(legend.position = "none")

p <- p1 / p2 + plot_layout(heights = c(3,1))

print(p)

ggsave("Figures/pdfs/Unsorted/roc_v3v4.pdf",
       width = 7.5, height = 7.25, units = "cm",
       device = cairo_pdf)






####### PRC #####
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

p1 <- ggplot(df |>
               filter(Model != "Max")) +
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
               mutate(Model = fct_relevel(Model, "v3", "v4", "Max")) |>
               filter(Model != "Max")) +
  geom_line(aes(x = sens, y = delta_sens, col = Model), linewidth = 0.4) +
  xlim(0, 1) +
  xlab("Recall") + ylab("Δ Precision") +
  theme_minimal(base_size = 8) + theme(legend.position = "none")

p <- p1 / p2 + plot_layout(heights = c(3,1))

print(p)

ggsave("Figures/pdfs/Unsorted/prc_v3v4.pdf",
       width = 7.5, height = 7.25, units = "cm",
       device = cairo_pdf)





### CAL #####

eval(import_sparra_expr("Analysis/full_model/Performance/v4_final/cal_max.txt"))



# Old Figure 2(c)
labs=c("v3", "v4","Max")
cic=c(rgb(0,0,1,alpha=0.5),rgb(1,0,1,alpha=0.5),rgb(0.5,0.5,0.5,alpha=0.5))
cal_2panel(list(xcal3,xcal4,xcalm),labels=labs,col=c("blue","red","black"),
           text.col=c("blue","red","black"),ci_col=cic,xy_col="gray",xy_lty=2)



# New Figure 2(c)
df <- rbind(
  cbind(Model = "v3", data.frame(pred = xcal3$x,
                                 obs  = xcal3$y[,1],
                                 del  = xcal3$y[,1] - xcal3$x,
                                 dell = xcal3$lower[,1] - xcal3$x,
                                 delu = xcal3$upper[,1] - xcal3$x)),
  cbind(Model = "v4", data.frame(pred = xcal4$x,
                                 obs  = xcal4$y[,1],
                                 del  = xcal4$y[,1] - xcal4$x,
                                 dell = xcal4$lower[,1] - xcal4$x,
                                 delu = xcal4$upper[,1] - xcal4$x)),
  cbind(Model = "Max", data.frame(pred = xcalm$x,
                                  obs  = xcalm$y[,1],
                                  del  = xcalm$y[,1] - xcalm$x,
                                  dell = xcalm$lower[,1] - xcalm$x,
                                  delu = xcalm$upper[,1] - xcalm$x))
) |> mutate(Model = fct_relevel(Model, "v3", "v4", "Max"))

p1 <- ggplot(df |>
               filter(Model != "Max")) +
  geom_line(aes(x = pred, y = obs, col = Model), linewidth = 0.4) +
  xlim(0, 1) + ylim(0, 1) +
  xlab("") + ylab("Observed") +
  theme_minimal(base_size = 8) + theme(legend.justification = c(0,1),
                                       legend.position = c(0,1),
                                       legend.spacing = unit(0, "npc"),
                                       legend.margin = unit(0, "npc"),
                                       legend.background = element_rect(fill = "white", size = 0, colour = "white"))

p2 <- ggplot(df |>
               filter(Model != "Max")) +
  geom_ribbon(aes(x = pred, ymin = dell, ymax = delu, fill = Model), alpha = 0.5) +
  geom_line(aes(x = pred, y = del, col = Model), linewidth = 0.4) +
  xlim(0, 1) +
  xlab("Predicted") + ylab("Δ from Calibrated") +
  theme_minimal(base_size = 8) + theme(legend.position = "none")

p <- p1 / p2 + plot_layout(heights = c(3,1))

print(p)

ggsave("Figures/pdfs/Unsorted/cal_v3v4.pdf",
       width = 7.5, height = 7.25, units = "cm",
       device = cairo_pdf)
