library(dplyr)
library(tidyr)
library(ggplot2)

source("rscript/functions.R")

# Read data files and select variably expressed tTREs.
source("rscript/fig1S.make.expression_variation_pvalue.R")

# Make exaple plots in SLFN and BCL genes
source("rscript/fig1.examples.R")

# Read data files and calculate background correlation distributions.
source("rscript/fig1S.correlation_background.R")

# Annotate mRNA TSS based on transcription activity
source("rscript/figS2.RNA_rennotation.R")

# Generate co-expression coefficient and p-value 
source("rscript/fig2.distance_coexpression.R")

# TF binding and correlation
source("rscript/fig3.TF_dependency_boxplots.R")
source("rscript/fig3.TF_dependency_AUC.R")
source("rscript/fig3.TF_dependency_FDR.R")
source("rscript/fig3S.make.auc.pval.R")

# Hi-C raw data analysis
source("rscript/fig4.HiC.R")

# Close range convergent transcription
source("rscript/fig5.close_range.R")

