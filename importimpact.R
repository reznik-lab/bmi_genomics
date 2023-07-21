# A quick script to imporat the MSK-IMPACT data and create a samplexgene matrix that encodes oncogenic mutations and optionally VUSs

# BUGS TO FIX: 
# 1: THE NUMBER OF SAMPLES IN SILENT, NS, AND ONCOGENIC TABLES IS NOT EQUAL

rm(list = ls())
setwd('/BMI/impact_bmi/')
library(data.table)
library(ggplot2)
library(reshape2)

# read in the impact data and then reshape to sample*gene either before/after filtering for VUS/silent
muts = fread('data_mutations_extended_somatic.oncokb.txt.gz',stringsAsFactors = FALSE,header = TRUE,data.table = FALSE)
clin = fread('data_clinical_sample_somatic.oncokb.txt',stringsAsFactors = FALSE,header = TRUE,data.table = FALSE)
mutssilent = fread('data_nonsignedout_mutations.txt',stringsAsFactors = FALSE,header = TRUE,data.table = FALSE)
mutssilent = mutssilent[which(mutssilent$Variant_Classification == 'Silent'),]
mutssilent$ID = 1

mutsog = muts[which(muts$ONCOGENIC %in% c('Likely Oncogenic','Oncogenic','Resistance')),]
mutsns = muts[which(muts$ONCOGENIC %in% c('Oncogenic','Likely Oncogenic','Resistance','Inconclusive','Unknown')),]

# for each of the above, make a table of samples x patients and annotate with cancer type and detailed cancer type
ogtable = dcast(data=mutsog,formula = Tumor_Sample_Barcode~Hugo_Symbol,value.var = 'ONCOGENIC')
rownames(ogtable) = ogtable[,1]
ogtable = ogtable[,-1]

nstable = dcast(data=mutsns,formula = Tumor_Sample_Barcode~Hugo_Symbol,value.var = 'ONCOGENIC')
rownames(nstable) = nstable[,1]
nstable = nstable[,-1]

silenttable = dcast(data = mutssilent,formula = Tumor_Sample_Barcode~Hugo_Symbol, value.var = 'ID')
rownames(silenttable) = silenttable[,1]
silenttable = silenttable[,-1]

# Restrict the silent table to just the genes in the oncogenic mutations table
silenttable = silenttable[,intersect(colnames(ogtable),colnames(silenttable))]

# We need to add back cases where we have fully wild-type samples with on apparent mutations
clin = clin[which(clin$GENE_PANEL%in%c('IMPACT341','IMPACT410','IMPACT468','IMPACT505')),]
diffog = setdiff(clin$SAMPLE_ID,rownames(ogtable))
diffns = setdiff(clin$SAMPLE_ID,rownames(nstable))
diffsilent = setdiff(clin$SAMPLE_ID,rownames(silenttable))

# Add these rows as zeros
ogtable[diffog,] = 0
nstable[diffns,] = 0
silenttable[diffsilent,] = 0

# FIX THIS: LET'S MAKE EXACTLY THE SAME PATIENTS FOR ALL 3 TABLES
ixall = intersect(rownames(ogtable),intersect(rownames(nstable),rownames(silenttable)))
ogtable = ogtable[ixall,]
nstable = nstable[ixall,]
silenttable = silenttable[ixall,]

# Write out tables
write.csv(ogtable,'IMPACT_Oncogenic_Table.csv')
write.csv(nstable,'IMPACT_AllMuts_Table.csv')
write.csv(silenttable,'IMPACT_SilentMuts_Table.csv')
