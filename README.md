## SPARRAv4

# Overview

This is a general repository for SPARRAv4 pipelines, functions and figures. 

# Organisation

This repository contains folders:

1. **Analysis**: containing figures, summaries and co-ordinates
2. **Pipelines**: containing code to fit models and draw figures, and to demonstrate calibration estimator
3. **SPARRAv4**: containing functions
4. **Diagrams**: containing diagrams and drawings
5. **Figures**: containing code to draw figures used in paper.

# Use

To generate all figure panels used in the manuscript, run **Figures/generate_all_figures.R** with the working directory at the main directory for this repository. Figure panels are saved into the relevant subfolder in **Figures/pdfs**.


# Notes

A small number of analyses are performed by code in
 - **Analysis/external/**

To find the original code to generate a plot, search for its basename in these files.

Generally, for each plot  **X.pdf** there is an associated file **X.txt**. This file contains all the data that was plotted. Each file firstly contains a nominally human-readable form of the data, then a line '***************', and then some lines which can be read into R directly to recreate the variables used (using dput). 

Figures were drawn on the windows machine in the NSH, for which package versions are listed as part of the R script **SPARRAv4/ehr\_mockup.R** (object _ip\_windows_).