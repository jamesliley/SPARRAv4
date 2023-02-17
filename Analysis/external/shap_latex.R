### Output LaTeX code to draw Shapley Value plots
## Run this from the directory containing subfolder 'Analysis'

# Get names of plots
shap_names=gsub(".pdf","",list.files("Analysis/full_model/Shapley_values/All/",pattern="*.pdf"))

# Order alphabetically
ox=order(gsub("shapley_","",shap_names))
shap_names=shap_names[ox]

# Basic subfigure code
shap1=paste0("\\begin{subfigure}{0.3\\textwidth}\n",
             "\\includegraphics[width=\\textwidth]{{",shap_names,"}.pdf}\n",
             "\\caption{}\n",
             "\\label{supp_fig:",shap_names,"}\n",
             "\\end{subfigure}")


# 3x3 on each line
shap3=c()
ind=1
for (i in 1:floor(length(shap_names)/3)) {
  sub=paste0(
    "\\begin{subfigure}{\\textwidth}\n",
    shap1[ind],"\n",shap1[ind+1],"\n",shap1[ind+2],"\n",
    "\\end{subfigure}"
  )
  ind=ind+3
  shap3=c(shap3,sub)
}
subx=paste0("\\begin{subfigure}{\\textwidth}\n")
for (j in ind:(length(shap_names))) subx=paste0(subx,shap1[j],"\n")
subx=paste0(subx,"\\end{subfigure}")
shap3=c(shap3,subx)


# 9x9 on each page
shap9=c()
ind=1
for (i in 1:floor(length(shap_names)/9)) {
  sub=paste0(
    "\\begin{figure}[b]\n",
    shap3[ind],"\n",shap3[ind+1],"\n",shap3[ind+2],"\n",
    "\\caption{Influence profile on risk scores by variable and value, ",
    "in decreasing order of mean absolute Shapley value. Vertical lines ",
    "show standard deviations. Y axes are identical in all plots. Data ",
    "are anonymised to have $>10$ individuals for each x axis mark, so ",
    "some X axis spacings are irregular. Plots are not shown if data are ",
    "not sufficiently anonymous }\n",
    "\\label{supp_fig:shapley_all",i,"}\n",
    "\\end{figure}\n"
  )
  ind=ind+3
  shap9=c(shap9,sub)
}
subx=paste0("\\begin{figure}[b]\n")
for (j in ind:(length(shap3))) subx=paste0(subx,shap3[j],"\n")
subx=paste0(subx,    "\\caption{Influence profile on risk scores by variable and value, ",
            "in decreasing order of mean absolute Shapley value. Vertical lines ",
            "show standard deviations. Y axes are identical in all plots. Data ",
            "are anonymised to have $>10$ individuals for each x axis mark, so ",
            "some X axis spacings are irregular. Plots are not shown if data are ",
            "not sufficiently anonymous }\n",
            "\\label{supp_fig:shapley_all",ceiling(length(shap_names)/9),"}\n",
            "\\end{figure}\n"
)
shap9=c(shap9,subx)

shap_all=paste(shap9,collapse="\n\n")

sink("Analysis/external/shapley_latex.txt")
cat(shap_all)
sink()