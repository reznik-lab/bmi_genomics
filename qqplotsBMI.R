# Analyze the results of the association analysis by making qqplots (panels 1a and 1b)

rm(list = ls())
setwd('BMI/impact_bmi/')
library(data.table)
library(ggplot2)
library(ggrepel)

# Read in the relevant data
resonc = fread('../results/bmi_results_oncogenic_allcovariates.csv',stringsAsFactors = FALSE,data.table = FALSE,header = TRUE)
resns = fread('../results/bmi_results_allmuts_allcovariates.csv',stringsAsFactors = FALSE,data.table = FALSE,header = TRUE)
resoncun = fread('../results/bmi_results_oncogenic_univariate.csv',stringsAsFactors = FALSE,data.table = FALSE,header = TRUE)
resnsun = fread('../results/bmi_results_allmuts_univariate.csv',stringsAsFactors = FALSE,data.table = FALSE,header = TRUE)
ressilentun = fread('../results/bmi_results_silent_univariate.csv',stringsAsFactors = FALSE,data.table = FALSE,header = TRUE)

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

  # pdf(paste0('../results/qq/',name,'.pdf'),height = 2,width = 2)
  # print(p)
  # dev.off()
   
   return(p)
    
  
}

qqplot(resonc,'Oncogenic Mutations, All Covariates')
qqplot(resns,'All Nonsynonymous Mutations, All Covariates')
qqplot(resoncun,'Oncogenic Mutations, Univariate')
qqplot(resnsun,'All Nonsynonymous Mutations, Univariate')
qqplot(ressilentun,'Silent Mutations, Univariate')


fig1a <- qqplot(resnsun,'All Nonsynonymous Mutations, Univariate')
ggsave("All Nonsynonymous Mutations, Univariate Scatterplot.pdf", width = 2.5, height = 2.5, units = "in")
fig1b <- qqplot(ressilentun,'Silent Mutations, Univariate')
ggsave("Silent Mutations, Univariate Scatterplot.pdf", width = 2.5, height = 2.5, units = "in")
