library("tidyverse")
library("patchwork")
library("glue")

source("Figures/util.R") # for import_sparra_expr()
source("SPARRAv4/auxiliary.R") # for roc_2panel()



eval(import_sparra_expr("Analysis/full_model/Performance/v3_vs_v4/v3_vs_v4_calibration_differential.txt"))



# Old Figure 2(d)
del=0.1
cic=c(rgb(1,0,0,alpha=0.5),rgb(1,0,0,alpha=0.5),rgb(0.5,0.5,0.5,alpha=0.5),rgb(0.5,0.5,0.5,alpha=0.5))
cal_2panel(list(cv3h3,cv3h4,cv4h3,cv4h4),labels=c("v3, v3>v4","v3, v4>v3","v4, v3>v4", "v4, v4>v3"),
           title=paste0("|v3-v4| > ",del),xy_col="gray",xy_lty=2,col=c("red","red","black","black"),
           lty=c(1,2,1,2),ci_col=cic)



# New Figure 2(d)
df <- rbind(
  cbind(Model = "v3",
        Order = "v3>v4",
        data.frame(pred = cv3h3$x,
                   obs  = cv3h3$y[,1],
                   del  = cv3h3$y[,1] - cv3h3$x,
                   dell = cv3h3$lower[,1] - cv3h3$x,
                   delu = cv3h3$upper[,1] - cv3h3$x)),
  cbind(Model = "v3",
        Order = "v4>v3",
        data.frame(pred = cv3h4$x,
                   obs  = cv3h4$y[,1],
                   del  = cv3h4$y[,1] - cv3h4$x,
                   dell = cv3h4$lower[,1] - cv3h4$x,
                   delu = cv3h4$upper[,1] - cv3h4$x)),
  cbind(Model = "v4",
        Order = "v3>v4",
        data.frame(pred = cv4h3$x,
                   obs  = cv4h3$y[,1],
                   del  = cv4h3$y[,1] - cv4h3$x,
                   dell = cv4h3$lower[,1] - cv4h3$x,
                   delu = cv4h3$upper[,1] - cv4h3$x)),
  cbind(Model = "v4",
        Order = "v4>v3",
        data.frame(pred = cv4h4$x,
                   obs  = cv4h4$y[,1],
                   del  = cv4h4$y[,1] - cv4h4$x,
                   dell = cv4h4$lower[,1] - cv4h4$x,
                   delu = cv4h4$upper[,1] - cv4h4$x))
) |> mutate(Model = fct_relevel(Model, "v3", "v4"),
            Order = fct_relevel(Order, "v3>v4", "v4>v3"))

p1 <- ggplot(df) +
  geom_line(aes(x = pred, y = obs, linetype = Model, col = Order), linewidth = 0.4) +
  xlim(0, 1) + ylim(0, 1) +
  xlab("") + ylab("Observed") +
  guides(linetype = guide_legend(title = glue("|v3 - v4| > {del}\n\nModel"), order = 1), col = guide_legend(order = 2)) +
  theme_minimal(base_size = 8) + theme(legend.justification = c(0,1),
                                       legend.position = c(0,1),
                                       legend.background = element_rect(fill = "white", size = 0))

p2 <- ggplot(df) +
  geom_ribbon(aes(x = pred, ymin = dell, ymax = delu, linetype = Model, fill = Order), alpha = 0.5) +
  geom_line(aes(x = pred, y = del, linetype = Model, col = Order), linewidth = 0.4) +
  xlim(0, 1) +
  xlab("Predicted") + ylab("Δ from Calibrated") +
  theme_minimal(base_size = 8) + theme(legend.position = "none")

p <- p1 / p2 + plot_layout(heights = c(3,1))

print(p)

ggsave("Figures/pdfs/Fig_2d.pdf",
       width = 7.5, height = 7.5, units = "cm",
       device = cairo_pdf)
