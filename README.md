# GEDDplot

Author: Long H. Do (long.h.do@gmail.com)

This script takes as input .diff files produced by CuffDiff to produce plots that search for GEDDs
Requires R (runs Rscript version) and the preprocesCore package from R Bioconductor for the quantile normalizaton

To install preprocessCore:
source("http://bioconductor.org/biocLite.R")
biocLite("preprocessCore")

usage: ./thisscript cuffDiffSample1.diff cuffDiffSample2.diff

output: spearman.pdf file in the current directory
