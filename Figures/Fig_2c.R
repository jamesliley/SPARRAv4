library("tidyverse")
library("patchwork")

source("Figures/util.R") # for import_sparra_expr()
source("SPARRAv4/auxiliary.R") # for roc_2panel()



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

p1 <- ggplot(df) +
  geom_line(aes(x = pred, y = obs, col = Model), linewidth = 0.4) +
  xlim(0, 1) + ylim(0, 1) +
  xlab("") + ylab("Observed") +
  theme_minimal(base_size = 8) + theme(legend.justification = c(0,1),
                                       legend.position = c(0,1),
                                       legend.spacing = unit(0, "npc"),
                                       legend.margin = unit(0, "npc"),
                                       legend.background = element_rect(fill = "white", size = 0, colour = "white"))

p2 <- ggplot(df) +
  geom_ribbon(aes(x = pred, ymin = dell, ymax = delu, fill = Model), alpha = 0.5) +
  geom_line(aes(x = pred, y = del, col = Model), linewidth = 0.4) +
  xlim(0, 1) +
  xlab("Predicted") + ylab("Δ from Calibrated") +
  theme_minimal(base_size = 8) + theme(legend.position = "none")

p <- p1 / p2 + plot_layout(heights = c(3,1))

print(p)

ggsave("Figures/pdfs/Fig_2c.pdf",
       width = 7.5, height = 7.5, units = "cm",
       device = cairo_pdf)
