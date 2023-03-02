library("tidyverse")
library("patchwork")

source("Figures/util.R") # for import_sparra_expr()

# SIMD panels

## Function used to generate SIMD plots per data source
simd_plot <- function(x, ylim.max = 0.2, title = TRUE,
                      base_size = 9) {
  # Load the data
  xsimd <- eval(import_sparra_expr(x[1]))
  df <- data.frame(
    simd <- as.numeric(names(xsimd)),
    freq <- as.numeric(xsimd)
  )
  p <- ggplot(data = df, aes(x = factor(simd), y = freq)) +
    geom_bar(stat="identity", width = 0.3) +
    xlab("SIMD decile") + ylab("Frequency") +
    ylim(0, ylim.max) +
    theme_minimal(base_size = base_size)
  if(title == TRUE)
    p <- p + ggtitle("") #ggtitle(x[2])
  return(p)
}

## Load the list of files extracted from the DSH
simd <- list.files("Analysis/full_model/Description",
                   pattern="simd_.*\\.txt", full.names = TRUE)

## Set the labels to be used in the figure
simd.names <- gsub("Analysis/full_model/Description/simd_|.txt","", simd)
simd.names <- recode(simd.names,
                     "AE2_A_E" = "Accidents and emergency",
                     "all" = "All",
                     "PIS_Prescr" = "GP prescribing",
                     "SMR00_Outpt" = "Outpatients",
                     "SMR01_Inp_day" = "Acute inpatients and day cases",
                     "SMR01E_SystemWatch_Ger_systwatch" = "Other",
                     "SMR04_MH_inp_day" = "Mental health inpatient and day cases")
names(simd) <- simd.names
simd <- cbind(simd, simd.names)

## Create the plots
simd.plots <- apply(simd, 1, simd_plot, title = TRUE)

## Place-holder to see how the plots look like
(simd.plots[[1]] + simd.plots[[2]]) /
  (simd.plots[[3]] + simd.plots[[4]]) /
  (simd.plots[[5]] + simd.plots[[6]]) /
  (simd.plots[[7]] + plot_spacer())


# Age/sex panels

## Function used to generate age/sex plots per data source
## DSH extract had to be modified due to missing object
age_plot <- function(x, ylim.max = 0.0225,
                     title = TRUE, withlegend = TRUE,
                     base_size = 9) {
  # Load the data
  eval(import_sparra_expr(x[1]))
  mydf <- age <- data.frame(age.f = df$x, dens.f = df$y,
                            age.m = dm$x, dens.m = dm$y)
  p <- ggplot(age) +
    geom_area(aes(x = age.f, y = dens.f, fill = "coral1"), alpha = 0.5) +
    geom_area(aes(x = age.m, y = dens.m, fill = "cornflowerblue"), alpha = 0.5) +
    scale_fill_discrete(name = 'Sex',
                        labels = c('F','M')) +
    xlab("Age") + ylab("Density") +
    ylim(0, ylim.max) +
    theme_minimal(base_size = base_size)
  if(title == TRUE)
    p <- p + ggtitle(x[2])
  if(withlegend == TRUE) {
    p <- p + theme(legend.position = c(0.9, 0.9),
                   legend.key.size = unit(0.2, 'cm'),
                   legend.background = element_rect(fill="white",
                                                    size=0.5, linetype="solid",
                                                    colour ="darkgrey"))
  } else {
    p <- p + theme(legend.position="none")
  }
  return(p)
}

## Load the list of files extracted from the DSH
age <- list.files("Analysis/full_model/Description",
                  pattern="age_.*\\_fix.txt", full.names = TRUE)

## Set the labels to be used in the figure
age.names <- gsub("Analysis/full_model/Description/age_|_fix.txt","", age)
age.names <- recode(age.names,
                    "AE2_A_E" = "Accidents and emergency",
                    "all" = "All",
                    "PIS_Prescr" = "GP prescribing",
                    "SMR00_Outpt" = "Outpatients",
                    "SMR01_Inp_day" = "Acute inpatients and day cases",
                    "SMR01E_SystemWatch_Ger_systwatch" = "Other",
                    "SMR04_MH_inp_day" = "Mental health inpatient and day cases")
names(age) <- age.names
age <- cbind(age, age.names)

## Create the plots
age.plots <- apply(age, 1, age_plot, title = TRUE, withlegend = TRUE)

# Place-holder to see how the plots look like
(age.plots[[1]] + age.plots[[2]]) /
  (age.plots[[3]] + age.plots[[4]]) /
  (age.plots[[5]] + age.plots[[6]]) /
  (age.plots[[7]] + plot_spacer())

# Combine plots based on the source table
persource.plot <- function(datasource, age.plots, simd.plots) {
  age.plots[[datasource]] + simd.plots[[datasource]] +
    plot_annotation(title = datasource,
                    theme = theme(plot.title = element_text(hjust = 0.5)))
}

all.plots <- lapply(names(age.plots), persource.plot,
                    age.plots = age.plots, simd.plots = simd.plots)
names(all.plots) <- gsub("&| ","", names(age.plots))

#for(i in seq_len(length(all.plots))) {
#  ggsave(paste0("Figures/pdfs/Fig_1b_", names(all.plots)[i] ,".pdf"),
#         plot = all.plots[[i]],
#         width = 10, height = 5, units = "cm",
#         device = cairo_pdf)
#}

# Final plots

## Figure 1b
age_plot(age["All",], title = FALSE, base_size = 10) +
  simd_plot(simd["All",], title = FALSE, base_size = 10)

ggsave("Figures/pdfs/Fig_1b.pdf",
       width = 12, height = 6, units = "cm",
       device = cairo_pdf)

## Supplementary
all.plots$All <- NULL
wrap_plots(all.plots, nrow = 3)

ggsave("Figures/pdfs/SupFig_1.pdf",
       width = 24, height = 18, units = "cm",
       device = cairo_pdf)

# Table S2

PaperTotal <- 478243585

AE2 <- c(7539454, 3031773)
deaths <- c(283554, 283554)
PIS <- c(393573549, 5589772)
SMR00 <- c(27463987, 3753240)
SMR01 <- c(26326889, 2205606)
SMR01E <- c(137340, 25015)
SMR04 <- c(111487, 51635)
SPARRALTC <- c(3286987, 1978171)
SystemWatch <- c(12890591, 1750009)
SMR01M <- c(6629747, 2311727)

AE2 + PIS + SMR00 + SMR01 + SMR01E + SMR04 + SystemWatch

AE2 + PIS + SMR00 + SMR01 + SMR01E + SMR04 + SystemWatch + SPARRALTC

AE2 + PIS + SMR00 + SMR01 + SMR01E + SMR04 + SystemWatch + SPARRALTC + SMR01M

AE2 + PIS + SMR00 + SMR01 + SMR01E + SMR04 + SystemWatch + SPARRALTC + SMR01M + deaths

# LTCs reconstructed
