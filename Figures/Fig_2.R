library("tidyverse")
library("patchwork")

source("Figures/util.R") # for import_sparra_expr()
source("SPARRAv4/auxiliary.R") # for roc_2panel()

eval(import_sparra_expr("Analysis/full_model/Performance/v4_final/roc_max.txt"))

# Old Figure 2(a)
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
df2 <- data.frame(spec = xroc3b$spec[1,],
                  v3 = xroc3b$sens[1,],
                  v4 = xroc4b$sens[1,],
                  Max = xrocmb$sens[1,])

p1 <- ggplot(df |>
               mutate(spec = 1-spec)) +
  geom_line(aes(x = spec, y = sens, col = Model), linewidth = 0.75) +
  xlab("") + ylab("Sensitivity") + theme_minimal() + theme(legend.position = c(0.91,0.25))

p2 <- ggplot(df2 |>
               mutate(spec = 1-spec,
                      v4 = v4 - v3,
                      Max = Max - v3,
                      v3 = v3 - v3) |>
               pivot_longer(v3:Max, names_to = "Model", values_to = "delta_sens") |>
               mutate(Model = fct_relevel(Model, "v3", "v4", "Max"))) +
  geom_line(aes(x = spec, y = delta_sens, col = Model), linewidth = 0.75) +
  xlab("1 - Specificity") + ylab("Δ Sensitivity") + theme_minimal() + theme(legend.position = "none")

p1 / p2 + plot_layout(heights = c(3,1))


