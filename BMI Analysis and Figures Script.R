library(dplyr)
library(lubridate)
library(data.table)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(sjPlot)

load("BMI Manuscript Data.RData") #load data downloaded from zenodo

#########Start analyses############
#restrict analysis to just the cancer types where we have signal
ccounts = table(clin$CANCER_TYPE_DETAILED)
uqctypes = names(ccounts)[which(ccounts > 50)]

g2use = intersect(g2use,colnames(mutsonc)) # in case gene names change
clinvar2use = 'bmi' # this indicates the clinical variable of interest

#make the data frames to store the relevant results
resonc = data.frame()
resns = data.frame()

resoncun = data.frame() # univariate results
resnsun = data.frame() # univariate results
ressilentun = data.frame() # univariate results

#run logistic regressions on individual gene-cancer type pairs
for (cc in uqctypes){
  print(cc)
  
  ii = rownames(clin)[which(clin$CANCER_TYPE_DETAILED == cc)]
  tempclin = clin[ii,]
  
  temponc = mutsonc[ii,]
  tempns = mutsns[ii,]
  tempsilent = mutssilent[ii,]
  
  #confirm there are at least 10 mutations
  if (sum(clin[rownames(temponc),clinvar2use],na.rm = TRUE) < 10){next}
  
  temponc[,c('clinvar','anc','age','sex','tmb')] = clin[rownames(temponc),c(clinvar2use,'ancestry','age','sex','tmb')]
  tempns[,c('clinvar','anc','age','sex','tmb')] = clin[rownames(tempns),c(clinvar2use,'ancestry','age','sex','tmb')]   
  tempsilent[,c('clinvar','anc','age','sex','tmb')] = clin[rownames(tempsilent),c(clinvar2use,'ancestry','age','sex','tmb')]   
  
  for (gg in g2use){
    
    temponcgg = temponc[,c(gg,'clinvar','anc','age','sex','tmb')]
    if (sum(temponcgg[,gg]) < 10){next}
    colnames(temponcgg) = c('gene','clinvar','anc','age','sex','tmb')
    tempnsgg = tempns[,c(gg,'clinvar','anc','age','sex','tmb')]
    colnames(tempnsgg) = c('gene','clinvar','anc','age','sex','tmb')
    tempsilentgg = tempsilent[,c(gg,'clinvar','anc','age','sex','tmb')]
    colnames(tempsilentgg) = c('gene','clinvar','anc','age','sex','tmb')
    
    ####oncogenic mutations, all covariates
    if (length(unique(temponcgg$sex[!is.na(temponcgg$sex)])) > 1){ # if there are multiple genders for this disease
      tempglm = glm(data = temponcgg,formula = gene ~ clinvar + age + sex + anc + tmb,family = 'binomial')
    }else{
      tempglm = glm(data = temponcgg,formula = gene ~ clinvar + age + anc + tmb,family = 'binomial')
    }
    
    tempsummary = summary(tempglm)
    resonc[paste0(cc,':',gg),c('cancer','gene','estimate','p')] = c(cc,gg,tempsummary$coefficients['clinvar','Estimate'],tempsummary$coefficients['clinvar','Pr(>|z|)'])
     
    ####all non-synonymous mutations, all covariates
    if (length(unique(tempnsgg$sex[!is.na(temponcgg$sex)])) > 1){ # if there are multiple genders for this disease
      tempglmns = glm(data = tempnsgg,formula = gene ~ clinvar + age + sex + tmb,family = 'binomial')
    }else{
      tempglmns = glm(data = tempnsgg,formula = gene ~ clinvar + age + tmb,family = 'binomial')
    }
    tempsummaryns = summary(tempglmns)
    resns[paste0(cc,':',gg),c('cancer','gene','estimate','p')] = c(cc,gg,tempsummaryns$coefficients['clinvar','Estimate'],tempsummaryns$coefficients['clinvar','Pr(>|z|)'])
    
    ####oncogenic mutations, univariate
    tempglmun = glm(data = temponcgg,formula = gene ~ clinvar ,family = 'binomial')
    tempsummaryun = summary(tempglmun)
    resoncun[paste0(cc,':',gg),c('cancer','gene','estimate','p')] = c(cc,gg,tempsummaryun$coefficients['clinvar','Estimate'],tempsummaryun$coefficients['clinvar','Pr(>|z|)'])
    
    ####all non-synonymous mutations, univariate
    tempglmunns = glm(data = tempnsgg,formula = gene ~ clinvar ,family = 'binomial')
    tempsummaryunns = summary(tempglmunns)
    resnsun[paste0(cc,':',gg),c('cancer','gene','estimate','p')] = c(cc,gg,tempsummaryunns$coefficients['clinvar','Estimate'],tempsummaryunns$coefficients['clinvar','Pr(>|z|)'])
    
    ####all silent mutations, univariate
    tempglmsilent = glm(data = tempsilentgg,formula = gene ~ clinvar ,family = 'binomial')
    tempsummarysilent = summary(tempglmsilent)
    ressilentun[paste0(cc,':',gg),c('cancer','gene','estimate','p')] = c(cc,gg,tempsummarysilent$coefficients['clinvar','Estimate'],tempsummarysilent$coefficients['clinvar','Pr(>|z|)'])
    
  }
  
}

# Write out the results
resonc$p = as.numeric(resonc$p) #oncogenic mutations multivariate
resonc$padj = p.adjust(resonc$p,method = 'BH')
resonc = resonc[order(resonc$p,decreasing = FALSE),]

resoncun$p = as.numeric(resoncun$p) #oncogenic mutations univariate
resoncun$padj = p.adjust(resoncun$p,method = 'BH')
resoncun = resoncun[order(resoncun$p,decreasing = FALSE),]

resns$p = as.numeric(resns$p) #all non-synonymous mutations multivariate
resns$padj = p.adjust(resns$p,method = 'BH')
resns = resns[order(resns$p,decreasing = FALSE),]

resnsun$p = as.numeric(resnsun$p) #all non-synonymous mutations univariate
resnsun$padj = p.adjust(resnsun$p,method = 'BH')
resnsun = resnsun[order(resnsun$p,decreasing = FALSE),]

ressilentun$p = as.numeric(ressilentun$p) #silent mutations univariate
ressilentun$padj = p.adjust(ressilentun$p,method = 'BH')
ressilentun = ressilentun[order(ressilentun$p,decreasing = FALSE),]

#########Start making figures############
#main figure 1, panel a
resoncun$label <- rownames(resoncun)
resoncun$estimate <- as.numeric(resoncun$estimate)
resoncun$sig <- ifelse(resoncun$padj < 0.05, "significant","not significant")

ggplot(resoncun,aes(x=estimate,y=-log10(p))) + geom_point(aes(color = sig)) + xlim(-0.25,0.25) + 
  geom_text_repel(data=subset(resoncun, sig == "significant"), aes(label=label),size=2,color="black") +
  geom_vline(xintercept = 0, color = "grey90", linetype = "dashed") +
  scale_color_manual(values=c("#D9D9D9", "#FF0000")) +
  theme_minimal() +
  xlab("Estimate") + ylab("-Log10 P-Value") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        text=element_text(size=6.5))


#main figure, panel b
ixluad = rownames(clin)[which(clin$CANCER_TYPE_DETAILED == 'Lung Adenocarcinoma')]
clin_la <- clin[ixluad,]
clin_la$KRAS <- mutsonc[ixluad,"KRAS"]
clin_la$EGFR <- mutsonc[ixluad,"EGFR"]
clin_la$bmi_cat <- NA
clin_la[which(clin_la$bmi < 18.5),'bmi_cat'] = 'Underweight'
clin_la[which(clin_la$bmi >= 18.5 & clin_la$bmi < 25),'bmi_cat'] = 'Healthy'
clin_la[which(clin_la$bmi >= 25 & clin_la$bmi < 30),'bmi_cat'] = 'Overweight'
clin_la[which(clin_la$bmi >= 30),'bmi_cat'] = 'Obese'
clin_la$bmi_cat = factor(clin_la$bmi_cat,levels =  c('Healthy','Underweight','Overweight','Obese'))

luadglm = glm(data = clin_la,formula = KRAS ~ bmi_cat + age + sex + smoking + ancestry + tmb,family = 'binomial')

plot_model(luadglm) + theme(plot.title = element_text(hjust = 0.5), text=element_text(size=6.5)) + 
  font_size(axis_title.x = 8,title = 8, labels.x = 8,labels.y = 8) + theme_minimal()



#main figure, panel c
sigres = resoncun[which(resoncun$padj < .05),]
sigres$estimate = as.numeric(sigres$estimate)
sigres$score = -log10(sigres$padj) * sign(sigres$estimate)
sigresun = resoncun[which(resoncun$padj < .05),]

maxscore = 3
sigres[which(sigres$score > maxscore),'score'] = maxscore
sigres[which(sigres$score < -maxscore),'score'] = -maxscore

rr <- rownames(sigres[2,]) #this is the row of lung adenocarcinoma:kras
gg = sigres[rr,'gene']
cc = sigres[rr,'cancer']
ii = rownames(clin)[which(clin$CANCER_TYPE_DETAILED == cc)]
temponc = mutsonc[ii,gg,drop = FALSE]
temponc$bmi = clin[rownames(temponc),clinvar2use]
colnames(temponc) = c('mutant','bmi')
temponc$anc = clin[rownames(temponc),'ancestry']
temponc$age = clin[rownames(temponc),'age']
temponc$sex = clin[rownames(temponc),'gender']
temponc$mutant = factor(temponc$mutant)
temponc$mutant <- as.integer(temponc$mutant)
temponc$mutant <- temponc$mutant -1

temponc$bmicat = NA
temponc[which(temponc$bmi < 18.5),'bmicat'] = 'Underweight'
temponc[which(temponc$bmi >= 18.5 & temponc$bmi < 25),'bmicat'] = 'Healthy'
temponc[which(temponc$bmi >= 25 & temponc$bmi < 30),'bmicat'] = 'Overweight'
temponc[which(temponc$bmi >= 30),'bmicat'] = 'Obese'
temponc$bmicat = factor(temponc$bmicat,levels = c('Underweight','Healthy','Overweight','Obese'))

mutation_fraction <- temponc %>%
  group_by(anc, bmicat) %>%
  summarise(
    Mean_Mutant_Fraction = mean(mutant),
    SampleCount = n()
  )
mutation_fraction <- as.data.frame(mutation_fraction)
mutation_fraction$SE <- sapply(1:nrow(mutation_fraction),function(i){sqrt(mutation_fraction[i,3]*(1-mutation_fraction[i,3])/mutation_fraction[i,4])})
use_anc <- c("ASJ","EAS","EUR")
mutation_fraction <- mutation_fraction[mutation_fraction$anc %in% use_anc,]

ggplot(mutation_fraction, aes(x = bmicat, y = Mean_Mutant_Fraction)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(x=bmicat, ymin=(Mean_Mutant_Fraction - SE),ymax=Mean_Mutant_Fraction),width=NA,color='white') +
  geom_errorbar(aes(x=bmicat,ymin=Mean_Mutant_Fraction,ymax=(Mean_Mutant_Fraction + SE)),width=NA,color='black') +
  facet_grid(anc~.,scales = 'fixed')+
  labs(title = "Mutation Fraction by Ancestry and BMI Category", y = "Mean Mutant Fraction") +
  theme_minimal()  + theme(text=element_text(size=6.5),plot.title = element_text(hjust = 0.5))


#main figure, panel d
prop_barplot <- function(rr = row){
  print(rr)
  gg = sigres[rr,'gene']
  cc = sigres[rr,'cancer']
  ii = rownames(clin)[which(clin$CANCER_TYPE_DETAILED == cc)]
  temponc = mutsonc[ii,gg,drop = FALSE]
  temponc$bmi = clin[rownames(temponc),clinvar2use]
  colnames(temponc) = c('mutant','bmi')
  temponc$anc = clin[rownames(temponc),'ancestry']
  temponc$age = clin[rownames(temponc),'age']
  temponc$sex = clin[rownames(temponc),'gender']
  temponc$mutant = factor(temponc$mutant)
  ggplot(temponc,aes(mutant,bmi)) + geom_violin() + geom_jitter() + theme_minimal() + ggtitle(paste(gg,cc))
  
  # make categories
  temponc$bmicat = NA
  temponc[which(temponc$bmi < 18.5),'bmicat'] = 'Underweight'
  temponc[which(temponc$bmi >= 18.5 & temponc$bmi < 25),'bmicat'] = 'Healthy'
  temponc[which(temponc$bmi >= 25 & temponc$bmi < 30),'bmicat'] = 'Overweight'
  temponc[which(temponc$bmi >= 30),'bmicat'] = 'Obese'
  temponc$bmicat = factor(temponc$bmicat,levels = c('Underweight','Healthy','Overweight','Obese'))
  propvals = prop.table(table(temponc$mutant,temponc$bmicat),margin = 2)
  bmi_count <- table(temponc$bmicat)
  prop_table <- sapply(1:4,function(i){sqrt(propvals[2,i]*(1-propvals[2,i])/bmi_count[i])})
  prop_table <- rbind(prop_table,propvals[2,])
  prop_table <- as.data.frame(t(prop_table))
  colnames(prop_table) <- c("se","prop")
  prop_table$cat <- rownames(prop_table)
  maxy = max(propvals[2,])
  prop_table$cat <- factor(prop_table$cat, levels = rownames(prop_table))
  
  tempcounts = table(temponc$bmicat)
  xlabels = paste0(names(tempcounts),'(',tempcounts,')')
  
  # for the purposes of plotting, drop NAs
  if (length(which(is.na(temponc$bmicat))) > 0){
    temponc = temponc[-which(is.na(temponc$bmicat)),]
  }
  
  sigres[rr,'chisqp'] = chisq.test(table(temponc$bmicat,temponc$mutant))$p.value
  
  plot <- ggplot(prop_table) + 
    geom_bar(aes(x=cat,y=prop),stat="identity") + 
    geom_errorbar(aes(x=cat, ymin=(prop-se),ymax=prop),width=NA,color='white') +
    geom_errorbar(aes(x=cat,ymin=prop,ymax=(prop+se)),width=NA,color='black') +
    ggtitle(paste0(gg,':',cc)) + 
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ylim(0,maxy+.05) +
    ylab('Proportion') + 
    scale_x_discrete(labels = xlabels) +
    scale_fill_manual(values = c('0' = 'white','1' = 'black')) + 
    xlab('') + ylab('Fraction Mutant Patients') + theme(legend.position = 'none',
                                                        text=element_text(size=6.5))
  return(plot)
}

prop_barplot(rownames(sigres)[1]) #Lung Adenocarcinoma:EGFR
prop_barplot(rownames(sigres)[2]) #Lung Adenocarcinoma:KRAS

#main figure, panel e
bmi_dfci$BMI <- c("Underweight(64)","Healthy(763)","Overweight(611)","Obese(394)")
bmi_dfci$Count <- c(64,763,611,394)
bmi_dfci$EGFR_prop <- bmi_dfci$EGFR.MUT/bmi_dfci$Count
bmi_dfci$KRAS_prop <- bmi_dfci$KRAS.MUT/bmi_dfci$Count
bmi_dfci$BMI <- factor(bmi_dfci$BMI, levels = c("Underweight(64)","Healthy(763)","Overweight(611)","Obese(394)"))
bmi_dfci$EGFR_se <- sapply(1:4,function(i){sqrt(bmi_dfci[i,7]*(1-bmi_dfci[i,7])/bmi_dfci[i,6])})
bmi_dfci$KRAS_se <- sapply(1:4,function(i){sqrt(bmi_dfci[i,8]*(1-bmi_dfci[i,8])/bmi_dfci[i,6])})

ggplot(bmi_dfci) + 
  geom_bar(aes(x=BMI,y=EGFR_prop),stat="identity") + 
  geom_errorbar(aes(x=BMI, ymin=EGFR_prop-EGFR_se,max=EGFR_prop),width=NA,color='white') +
  geom_errorbar(aes(x=BMI, ymin=EGFR_prop,max=EGFR_prop+EGFR_se),width=NA,color='black') +
  ggtitle("DFCI Lung Adenocarcinoma:EGFR") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), plot.title = element_text(hjust = 0.5)) +
  ylab('Proportion') + 
  xlab('') + ylab('Fraction Mutant Patients') + theme(legend.position = 'none', text=element_text(size=6.5))

ggplot(bmi_dfci) + 
  geom_bar(aes(x=BMI,y=KRAS_prop),stat="identity") + 
  geom_errorbar(aes(x=BMI, ymin=KRAS_prop-KRAS_se,max=KRAS_prop),width=NA,color='white') +
  geom_errorbar(aes(x=BMI, ymin=KRAS_prop,max=KRAS_prop+KRAS_se),width=NA,color='black') +
  ggtitle("DFCI Lung Adenocarcinoma:KRAS") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), plot.title = element_text(hjust = 0.5)) +
  ylab('Proportion') + 
  xlab('') + ylab('Fraction Mutant Patients') + theme(legend.position = 'none', text=element_text(size=6.5))




#supplemental figure 2
qqplot = function(d,name){
  
  # QQ plot function, expects ordered set with increasing p-values
  
  # Grab the p-values, get their quantiles, and then plot the negative logs against each other
  d$quantile = (1:dim(d)[1])/dim(d)[1]
  d$logpadj = -log10(d$p)
  d$logquant = -log10(d$quantile)
  
  
  
  d$significant = FALSE
  d[which(d$padj < 5e-2),'significant'] = TRUE
  d$label = rownames(d)
  d[which(d$significant == FALSE),'label'] = NA
  
  p = ggplot(d,aes(logquant,logpadj,color = significant,label = label)) + geom_point()  + geom_abline() +
    theme_minimal(base_size = 4) + xlim(0,max(c(d$logpadj,d$logquant))) +
    ylim(0,max(c(d$logpadj,d$logquant))) +
    ggtitle(name) +
    xlab('-Log10 Expected P Value') +
    ylab('-Log10 Observed P Value') +
    geom_text_repel(size = 2,force = 20) +
    scale_color_manual(values = c('FALSE' = 'gray','TRUE' = 'red')) +
    theme(legend.position = 'none')
  
    return(p)
}
qqplot(ressilentun,'Silent Mutations, Univariate')







