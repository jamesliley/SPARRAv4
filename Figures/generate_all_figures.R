## Draw all figures used in paper

## Find names of R scripts/files
scriptnames=list.files("Figures",pattern="*.R",full.names=TRUE)
scriptnames=grep("util.R",scriptnames,val=T,invert = TRUE)
scriptnames=grep("generate_all_figures.R",scriptnames,val=T,invert = TRUE)
figurenames=gsub(".R","",basename(scriptnames),fixed=TRUE)

# Clear destination folders
for (i in 1:length(figurenames)) {
  flist=list.files(paste0("Figures/pdfs/",figurenames[i]),full.names=TRUE)
  if (length(flist)>0) {
    for (j in 1:length(flist)) {
     file.remove(flist[j]) 
    }
  }
}

# Run R scripts
for (iiii in 1:length(scriptnames)) { # long loop name is because we keep redefining 'i'
  source(scriptnames[iiii])
  print(paste0("Completed figure: ",scriptnames[iiii]))
}