library("tidyverse")
library("patchwork")

source("Figures/util.R") # for import_sparra_expr()


# SIMD panels
simd_plot <- function(x, ylim.max = 0.2) {
  # Load the data
  xsimd <- eval(import_sparra_expr(x[1]))
  df <- data.frame(
    simd <- as.numeric(names(xsimd)),
    freq <- as.numeric(xsimd)
  )
  ggplot(data = df, aes(x = factor(simd), y = freq)) +
    geom_bar(stat="identity", width = 0.3) +
    xlab("SIMD decile") + ylab("Frequency") +
    ylim(0, ylim.max) +
    ggtitle(x[2]) +
    theme_minimal(base_size = 12)
}

simd <- list.files("Analysis/full_model/Description",
                   pattern="simd_.*\\.txt", full.names = TRUE)
simd.names <- gsub("Analysis/full_model/Description/simd_|.txt","", simd)
simd.names <- recode(simd.names,
                     "AE2_A_E" = "A&E",
                     "all" = "All",
                     "PIS_Prescr" = "Prescriptions",
                     "SMR00_Outpt" = "Outpatients",
                     "SMR01_Inp_day" = "Inpatients",
                     "SMR01E_SystemWatch_Ger_systwatch" = "Other",
                     "SMR04_MH_inp_day" = "Mental health inpatients")
names(simd) <- simd.names
simd <- cbind(simd, simd.names)
simd.plots <- apply(simd, 1, simd_plot)

# Place-holder to see how the plots look like
(simd.plots[[1]] + simd.plots[[2]]) /
  (simd.plots[[3]] + simd.plots[[4]]) /
  (simd.plots[[5]] + simd.plots[[6]]) /
  (simd.plots[[7]] + plot_spacer())





### Age
eval(import_sparra_expr_v2("Analysis/full_model/Description/age_all.txt"))


