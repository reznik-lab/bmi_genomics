# A quick script to import the MSK-IMPACT FACETS data and create a samplexgene matrix that encodes losses, gains, etc.

rm(list = ls())
setwd('BMI/impact_bmi/')
library(data.table)
library(ggplot2)
library(reshape2)

# read in the impact data and then reshape to sample*gene either before/after filtering for VUS/silent
facets = fread('msk_impact_facets_annotated.gene_level.txt.gz',header = TRUE,stringsAsFactors = FALSE,data.table = FALSE)
facetscohort = fread('msk_impact_facets_annotated.cohort.txt.gz',header = TRUE,stringsAsFactors = FALSE,data.table = FALSE)
rownames(facetscohort) = facetscohort$tumor_sample
facets$tumorsample = sapply(facets$sample,function(x){strsplit(x,'\\_')[[1]][1]})

facets2use = rownames(facetscohort)[which(facetscohort$facets_qc)]
facets = facets[which(facets$tumorsample %in% facets2use),]

facetslcn = dcast(data = facets,formula = tumorsample~gene, value.var = 'lcn')
facetstcn = dcast(data = facets,formula = tumorsample~gene, value.var = 'tcn')
rownames(facetslcn) = facetslcn[,1]
rownames(facetstcn) = facetstcn[,1]
facetslcn = facetslcn[,-1]
facetstcn = facetstcn[,-1]

g2use = read.csv('impact341_gene_panel.txt',header = FALSE,sep = '\t')
g2use = as.character(g2use[4,2:dim(g2use)[2]])
facetslcn = facetslcn[,intersect(g2use,colnames(facetslcn))]
facetstcn = facetstcn[,intersect(g2use,colnames(facetstcn))]

write.csv(facetslcn,'facetslcn.csv')
write.csv(facetstcn,'facetstcn.csv')
