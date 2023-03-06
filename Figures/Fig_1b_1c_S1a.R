################################################################
# Code used to generate Fig 1b, 1c (Venn diagram only) and S1a #
################################################################

library("tidyverse")
library("patchwork")

source("Figures/util.R") # for import_sparra_expr()

###############################################################
# Load and plot the distribution of records across EHR tables #
###############################################################

# SIMD

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

# Age/sex

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

# Combine plots based on the source table
persource.plot <- function(datasource, age.plots, simd.plots) {
  age.plots[[datasource]] + simd.plots[[datasource]] +
    plot_annotation(title = datasource,
                    theme = theme(plot.title = element_text(hjust = 0.5)))
}

all.plots <- lapply(names(age.plots), persource.plot,
                    age.plots = age.plots, simd.plots = simd.plots)
names(all.plots) <- gsub("&| ","", names(age.plots))

###############################################################
########################## Figure 1b ##########################
###############################################################

age_plot(age["All",], title = FALSE, base_size = 10) +
  simd_plot(simd["All",], title = FALSE, base_size = 10)

ggsave("Figures/pdfs/Fig_1b.pdf",
       width = 12, height = 6, units = "cm",
       device = cairo_pdf)

###############################################################
########################## Figure 1c ##########################
###############################################################

## Read the DSH extract
ex <- readLines("Analysis/full_model/Description/exclusions.txt")

### Total individuals
as.numeric(ex[which(ex == "Total patients recorded anywhere in data tables")+1])

### Total samples (individual-time pairs)
as.numeric(ex[which(ex == "Total patient-time pairs under consideration")+1])

## VennDiagram summarising exclusions
## Code source: Analysis/external/draw_exclusions_plot.R

library(VennDiagram)

wA=which(ex=="A: Number of records for which individual died before time cutoff")
all=c(as.numeric(ex[wA+3*(0:13)+1]),0)
names(all)=c("A","B","C","D","AB","AC","AD","BC","BD","CD","ABC","ABD","ACD","BCD","ABCD")
for (i in 1:length(all)) assign(names(all)[i],all[i])

pdf("Figures/pdfs/Fig_1c_venn.pdf",width=4,height=4)
grid.newpage()
draw.quad.venn(area1=A, area2=B, area3=C,
               area4 =D, n12=AB, n23=BC, n13=AC,
               n14= AD,n24=BD, n34=CD, n123=ABC,
               n124=ABD, n234=BCD, n134=ACD, n1234=ABCD,
               category=c("Died","No v3","No SIMD","Unmatched"),
               col="Green",fill=c("Red","Pink","Blue","Orange"),lty="dashed")
dev.off()

### Total samples (individual-time pairs) after exclusions
### Matches total pre-exclusions minus all numbers in Venn Diagram
as.numeric(ex[which(ex == "Total records excl. individuals with no v3 score")+1])

### Total individuals after exclusions
as.numeric(ex[which(ex == "Total patients in study")+1])

##############################################################
######################### Figure S1a #########################
##############################################################

all.plots$All <- NULL
wrap_plots(all.plots, nrow = 3)

ggsave("Figures/pdfs/SupFig_1a.pdf",
       width = 24, height = 18, units = "cm",
       device = cairo_pdf)

##############################################################
######################### Figure S1b #########################
##############################################################

# Outcome decomposition (EA vs death) by age - cumulative
eval(import_sparra_expr("Analysis/full_model/Analytics/admission_type_by_age_cumulative.txt"))

# The line above loads 3 objects detailed as follows:
# n_tot: total number of samples with an event per age and SIMD group
# n_d: total number of samples with a deaths without a prior EA per age and SIMD group
# n_b: total number of samples with a death & EA (both) per age and SIMD group

# Total number of samples with an event (composite outcome)
sum(n_tot)
# Total number of samples with EA (with our without death)
sum(n_tot) - sum(n_d)
# Total number of samples with a death without a prior EA
sum(n_d)
# Total number of samples with with EA and death
sum(n_b)

# Code source: Pipelines/main_pipeline.R

# Age cutoffs: aggregate groups
age_cuts=c(0,5,20 + 5*(0:14))
nt=rowSums(n_tot)[1:amax]
nd=rowSums(n_d)[1:amax]
nb=rowSums(n_b)[1:amax]
x=(c(age_cuts,100)+c(0,age_cuts))[2:(amax+1)]/2

# General plot
pdf("Figures/pdfs/SupFig_1b.pdf",width=5,height=5)
plot(0,type="n",xlim=range(x),ylim=c(0,1),xlab="Age",
     ylab="Cumulative proportion of events",xaxs="i",yaxs="i")
polygon(c(min(x),min(x),max(x),max(x)),c(0,1,1,0),col="gray",border=NA)
polygon(c(x,max(x),min(x)),c((nd+nb)/nt,0,0),col="red",border=NA)
polygon(c(x,max(x),min(x)),c((nd/nt),0,0),col="black",border=NA)
legend("topleft",
       c("EA","Death (without prior EA)","EA and subsequent death"),
       pch = 16, col = c("gray","black","red"), bg = "white")
dev.off()


