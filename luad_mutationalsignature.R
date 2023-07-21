# DETERMINE WHETHER THE MUTATIONAL SIGNATURE FOR SMOKING (SBS4) HAS AN EFFECT ON OUR RESULTS IN LUNG ADENOCARCINOMA PATIENTS
rm(list = ls())
setwd('/BMI/impact_bmi/')
library(data.table)
library(ggplot2)
library(sjPlot)
library(ComplexHeatmap)
####################################################################################
# Initialize the key variables to investigate
################################################################################

# get 341 gene panel
g2use = read.csv('impact341_gene_panel.txt',header = FALSE,sep = '\t')
g2use = as.character(g2use[4,2:dim(g2use)[2]])

mutsonc = fread('IMPACT_Oncogenic_Table.csv',stringsAsFactors = FALSE,sep = ',',data.table = FALSE)
rownames(mutsonc) = mutsonc[,1]
mutsonc=mutsonc[,-1]
mutsns = fread('IMPACT_AllMuts_Table.csv',sep = ',',stringsAsFactors = FALSE,data.table = FALSE)
rownames(mutsns) = mutsns[,1]
mutsns = mutsns[,-1]
mutssilent = fread('IMPACT_SilentMuts_Table.csv',sep = ',',stringsAsFactors = FALSE,data.table = FALSE)
rownames(mutssilent) = mutssilent[,1]
mutssilent = mutssilent[,-1]

# convert all multiples to 1's
mutsonc[which(mutsonc > 1,arr.ind = TRUE)] = 1
mutsns[which(mutsns > 1,arr.ind = TRUE)] = 1
mutssilent[which(mutssilent > 1,arr.ind = TRUE)] = 1

# get ancestry as well
ancestry = fread('cd-22-0312_supplementary_tables_s1_-_s6_suppst1_ancestry.csv',header = TRUE,stringsAsFactors = FALSE,sep = ',',data.table = FALSE)
ancestry$Patient = sapply(ancestry[,1],function(x){strsplit(x,'MSK-')[[1]][2]})
ancestry = ancestry[-which(duplicated(ancestry$Patient) | is.na(ancestry$Patient)),]
rownames(ancestry) = ancestry$Patient

# read in the clinical table, add BMI
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
nrow(clin) #3486

# Intersect all of our clinical data with mutation data
ixpats = intersect(rownames(clin),rownames(mutsonc))
mutsonc = mutsonc[ixpats,]
mutsns = mutsns[ixpats,]
mutssilent = mutssilent[ixpats,]
clin = clin[ixpats,]

# Restrict analysis to just the cancer types where we have signal
ccounts = table(clin$CANCER_TYPE_DETAILED)
uqctypes = names(ccounts)[which(ccounts > 50)]


g2use = intersect(g2use,colnames(mutsonc)) # in case gene names change
clinvar2use = 'bmi' # this indicates the clinical variable of interest

#extract lung adenocarinoma patients 
la <- rownames(clin)[which(clin$CANCER_TYPE_DETAILED == 'Lung Adenocarcinoma')]

#add smoking signature to clinical data
ss <- read.csv("~/Documents/BMI/impact_bmi/nlp_impact_samples_AM.csv")

#add smoking signature to clinical data
clin$smoking_signature <-ss$Dominant_Signature[match(clin$SAMPLE_ID,ss$SAMPLE_ID)]
clin <- clin[-which(is.na(clin$smoking_signature)),] #remove where there aren't matches in ss dataframe
clin <- clin[-which(clin$smoking_signature == "Not Included in the previous study"),] #remove where signatures weren't calculated
clin$smoking_signature <- ifelse(clin$smoking_signature == "Smoking (SBS4)","yes","no")

################################################################################
# End initialization
################################################################################
# Make the data frames to store the relevant results
resonc = data.frame()
resns = data.frame()

resoncun = data.frame() # univariate results
resnsun = data.frame() # univariate results
ressilentun = data.frame() # univariate results

cc <- c("lung_adeno")



ii = intersect(la,clin$SAMPLE_ID)
tempclin = clin[ii,]

temponc = mutsonc[ii,]
tempns = mutsns[ii,]
tempsilent = mutssilent[ii,]


temponc[,c('clinvar','anc','age','sex','smoking_signature')] = clin[rownames(temponc),c(clinvar2use,'ancestry','age','sex','smoking_signature')]
tempns[,c('clinvar','anc','age','sex','smoking_signature')] = clin[rownames(tempns),c(clinvar2use,'ancestry','age','sex',"smoking_signature")]   
tempsilent[,c('clinvar','anc','age','sex','smoking_signature')] = clin[rownames(tempsilent),c(clinvar2use,'ancestry','age','sex','smoking_signature')]  

for (gg in g2use){
  
  temponcgg = temponc[,c(gg,'clinvar','anc','age','sex','smoking_signature')]
  if (sum(temponcgg[,gg]) < 10){next}
  colnames(temponcgg) = c('gene','clinvar','anc','age','sex','smoking_signature')
  tempnsgg = tempns[,c(gg,'clinvar','anc','age','sex','smoking_signature')]
  colnames(tempnsgg) = c('gene','clinvar','anc','age','sex','smoking_signature')
  tempsilentgg = tempsilent[,c(gg,'clinvar','anc','age','sex','smoking_signature')]
  colnames(tempsilentgg) = c('gene','clinvar','anc','age','sex','smoking_signature')
  
  #### Oncogenic mutations, all covariates
  if (length(unique(temponcgg$sex[!is.na(temponcgg$sex)])) > 1){ # if there are multiple genders for this disease
    tempglm = glm(data = temponcgg,formula = gene ~ clinvar + age + sex + anc + smoking_signature,family = 'binomial')
  }else{
    tempglm = glm(data = temponcgg,formula = gene ~ clinvar + age + anc + smoking_signature,family = 'binomial')
  }
  
  tempsummary = summary(tempglm)
  resonc[paste0(cc,':',gg),c('cancer','gene','estimate','p')] = c(cc,gg,tempsummary$coefficients['clinvar','Estimate'],tempsummary$coefficients['clinvar','Pr(>|z|)'])
  
  
  #### All non-synonymous mutations, all covariates
  if (length(unique(tempnsgg$sex[!is.na(temponcgg$sex)])) > 1){ # if there are multiple genders for this disease
    tempglmns = glm(data = tempnsgg,formula = gene ~ clinvar + age + sex + smoking_signature,family = 'binomial')
  }else{
    tempglmns = glm(data = tempnsgg,formula = gene ~ clinvar + age + smoking_signature,family = 'binomial')
  }
  tempsummaryns = summary(tempglmns)
  resns[paste0(cc,':',gg),c('cancer','gene','estimate','p')] = c(cc,gg,tempsummaryns$coefficients['clinvar','Estimate'],tempsummaryns$coefficients['clinvar','Pr(>|z|)'])
  
  #### Oncogenic mutations, univariate
  tempglmun = glm(data = temponcgg,formula = gene ~ clinvar ,family = 'binomial')
  tempsummaryun = summary(tempglmun)
  resoncun[paste0(cc,':',gg),c('cancer','gene','estimate','p')] = c(cc,gg,tempsummaryun$coefficients['clinvar','Estimate'],tempsummaryun$coefficients['clinvar','Pr(>|z|)'])
  
  #### All non-synonymous mutations, univariate
  tempglmunns = glm(data = tempnsgg,formula = gene ~ clinvar ,family = 'binomial')
  tempsummaryunns = summary(tempglmunns)
  resnsun[paste0(cc,':',gg),c('cancer','gene','estimate','p')] = c(cc,gg,tempsummaryunns$coefficients['clinvar','Estimate'],tempsummaryunns$coefficients['clinvar','Pr(>|z|)'])
  
  #### All silent mutations, univariate
  tempglmsilent = glm(data = tempsilentgg,formula = gene ~ clinvar ,family = 'binomial')
  tempsummarysilent = summary(tempglmsilent)
  ressilentun[paste0(cc,':',gg),c('cancer','gene','estimate','p')] = c(cc,gg,tempsummarysilent$coefficients['clinvar','Estimate'],tempsummarysilent$coefficients['clinvar','Pr(>|z|)'])
  
}




resonc$p = as.numeric(resonc$p)
resonc$padj = p.adjust(resonc$p,method = 'BH')
resonc = resonc[order(resonc$p,decreasing = FALSE),]

resoncun$p = as.numeric(resoncun$p)
resoncun$padj = p.adjust(resoncun$p,method = 'BH')
resoncun = resoncun[order(resoncun$p,decreasing = FALSE),]

resns$p = as.numeric(resns$p)
resns$padj = p.adjust(resns$p,method = 'BH')
resns = resns[order(resns$p,decreasing = FALSE),]

resnsun$p = as.numeric(resnsun$p)
resnsun$padj = p.adjust(resnsun$p,method = 'BH')
resnsun = resnsun[order(resnsun$p,decreasing = FALSE),]

ressilentun$p = as.numeric(ressilentun$p)
ressilentun$padj = p.adjust(ressilentun$p,method = 'BH')
ressilentun = ressilentun[order(ressilentun$p,decreasing = FALSE),]

write.csv(resonc,paste0('/results/ss/',clinvar2use,'_results_oncogenic_allcovariates.csv'))
write.csv(resns,paste0('/results/ss/',clinvar2use,'_results_allmuts_allcovariates.csv'))

write.csv(resoncun,paste0('/results/ss/',clinvar2use,'_results_oncogenic_univariate.csv'))
write.csv(resnsun,paste0('/results/ss/',clinvar2use,'_results_allmuts_univariate.csv'))
write.csv(ressilentun,paste0('/results/ss/',clinvar2use,'_results_silent_univariate.csv'))









