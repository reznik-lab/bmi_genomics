# Investigate the specific alleles in KRAS and EGFR in lung adenocarcinoma

setwd('/BMI/impact_bmi/')
library(data.table)
library(ggplot2)
library(sjPlot)
library(ComplexHeatmap)
####################################################################################
# Initialize the key variables to investigate
################################################################################

# Genes to investigate
g2use = c('KRAS','EGFR')

# Get the sample IDs from lung with EGFR/KRAS alterations and their obesity status
maf = fread('data_mutations_extended_somatic.oncokb.txt.gz',header = TRUE,stringsAsFactors = FALSE,data.table = FALSE)
muts = fread('IMPACT_Oncogenic_Table.csv',stringsAsFactors = FALSE,sep = ',',data.table = FALSE) # decide on oncogenic or all NS
rownames(muts) = muts[,1] 
muts = muts[,-1]

# get ancestry as well
ancestry = fread('cd-22-0312_supplementary_tables_s1_-_s6_suppst1_ancestry.csv',header = TRUE,stringsAsFactors = FALSE,sep = ',',data.table = FALSE)
ancestry$Patient = sapply(ancestry[,1],function(x){strsplit(x,'MSK-')[[1]][2]})
ancestry = ancestry[-which(duplicated(ancestry$Patient) | is.na(ancestry$Patient)),]
rownames(ancestry) = ancestry$Patient

# read in the clinical table and add BMI
clin = fread('data_clinical_sample_somatic.oncokb.txt',header = TRUE,stringsAsFactors = FALSE,data.table = FALSE)
rownames(clin) = clin[,1]
bmi = fread('bmi_for_reznik_3_21_23.csv',header = TRUE,sep = ',',data.table = FALSE,stringsAsFactors = FALSE)
bmi = bmi[-which(duplicated(bmi[,1])),]
rownames(bmi) = bmi[,1]

comorbcols2use = c('BMI','AGE','GENDER','Smoking_History_Prediction_Aggregated',
                   'CCI_DIABETES','CCI_HIV_AND_AIDS','CCI_RENAL_DISEASE','CCI_RHEUMATOLOGIC_DISEASE')
names2use = c('bmi','age','sex','smoking','diabetes','hiv','renal','rheum')

clin$ancestry = ancestry[clin$PATIENT_ID,2]
clin = clin[-which(is.na(clin$ancestry)),]
clin[,names2use] = bmi[clin$PATIENT_ID,comorbcols2use]

# Drop low BMIs, data is wrong
clin = clin[-which(clin$bmi < 1),]

# 1 sample per patient --> earliest primary if there are multiple
clin <- clin[
  with(clin, order(PATIENT_ID, DATE_ADDED, rev(SAMPLE_TYPE))),
]

clin <- clin[!duplicated(clin$PATIENT_ID),]

# Clinical variable to investigate
clinvar2use = 'bmi' # this indicates the clinical variable of interest


# Intersect all of our clinical data with mutation data
ixpats = intersect(rownames(clin),rownames(muts))
muts = muts[ixpats,]
clin = clin[ixpats,]

# Choose the lung adenocarcinoma patients
ixluad = rownames(clin)[which(clin$CANCER_TYPE_DETAILED == 'Lung Adenocarcinoma')]
muts = muts[ixluad,]
clin = clin[ixluad,]
clin$EGFR = muts[rownames(clin),'EGFR']
clin$KRAS = muts[rownames(clin),'KRAS']
clin[which(clin$bmi < 18.5),'bmicat'] = 'Underweight'
clin[which(clin$bmi >= 18.5 & clin$bmi < 25),'bmicat'] = 'Healthy'
clin[which(clin$bmi >= 25 & clin$bmi < 30),'bmicat'] = 'Overweight'
clin[which(clin$bmi >= 30),'bmicat'] = 'Obese'

# convert multiples to 1s
clin[which(clin$KRAS > 1),'KRAS'] = 1
clin[which(clin$EGFR > 1),'EGFR'] = 1
clin$bmicat = factor(clin$bmicat,levels =  c('Underweight','Healthy','Overweight','Obese'))
  
write.csv(clin,'LUADclin.csv') # use this file to make the key oncoprints in cbioportal

# Make the key forest plot
luadglm = glm(data = clin,formula = KRAS ~ bmicat + age + sex + smoking + ancestry,family = 'binomial')
p = plot_model(luadglm) + theme(text = element_text(size = 8)) + 
  font_size(axis_title.x = 8,title = 8, labels.x = 8,labels.y = 8)
pdf('LUAD KRAS.pdf',height = 3.29,width = 3.5)
print(p)
dev.off()

luadglmegfr = glm(data = clin,formula = EGFR ~ bmicat + age + sex + smoking + ancestry,family = 'binomial')










