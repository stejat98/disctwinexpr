## Install & Load necesssary packages
#source("http://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")
require(gdata)
library(Biobase)
library(GEOquery)
require(data.table)
library(rmeta)

## get data pertaining to GSE16059 
gse16059 <- getGEO("GSE16059", GSEMatrix = TRUE,AnnotGPL = TRUE)

## used to extract phenotype/disease diagnosis corresponding to each twin sample
phenodiagnosis_GSE16059 <- gse16059$ GSE16059_series_matrix.txt.gz@phenoData@data$characteristics_ch1.2

## used to extract twin pair # and sex info corresponding to each twin sample
twp_num <- gse16059$GSE16059_series_matrix.txt.gz@phenoData@data$characteristics_ch1
sex_group <- gse16059$GSE16059_series_matrix.txt.gz@phenoData@data$characteristics_ch1.1 

## create vectors of 0s and 1s to indicate phenotype status
Healthy0.CFS1 <- vector(mode="integer", length=88)
Healthy0.CFS1[grep('unaffected', phenodiagnosis_GSE16059)] <- 0
Healthy0.CFS1[c(grep('CFS', phenodiagnosis_GSE16059),grep('ICF', phenodiagnosis_GSE16059))] <- 1

## extract expression table
gsm_tables_GSE16059 <- exprs(gse16059$'GSE16059_series_matrix.txt.gz')
gsm_tables_GSE16059 <- as.data.frame(gsm_tables_GSE16059)
gsm16059_names <- gse16059$'GSE16059_series_matrix.txt.gz'@phenoData@data$geo_accession


## Table with information (sex, age, twin pair number, etc.) for each sample
Sample_Info <- cbind(TWP_NUM = gsub('twin pair:','',twp_num), Healthy0.CFS1,Sex_Group = gsub('sex: ','',sex_group),GSM = as.character(gsm16059_names))
Sample_Info <- as.data.frame(Sample_Info)

## separate Sample_Info table into two tables: one only for affected samples and one for the unaffected
healthy_samples_male <- subset(Sample_Info, Healthy0.CFS1 ==0 & Sex_Group == 'male')
diseased_samples_male <- subset(Sample_Info, Healthy0.CFS1 == 1 & Sex_Group == 'male')


## affected group of samples
tone <- gsm_tables_GSE16059[,match(diseased_samples_male$GSM,colnames(gsm_tables_GSE16059))]

## unaffected group of samples (control)
ttwo <- gsm_tables_GSE16059[,match(healthy_samples_male$GSM,colnames(gsm_tables_GSE16059))]


## get annotation information for GPL570 platform
gpl570info <- getGEO('GPL570',AnnotGPL = TRUE)
str(gpl570info)

## Extract annotation table for GPL570 platform
gpl570infodf <- gpl570info@dataTable@table



## Extract microarray probe ids and assign them to a vector
probe_id_gse16059 <- as.character(t(gpl570infodf$ID))
## Extract Entrez Gene IDs and assign them to a vector
gene_id_gse16059 <- as.character(t(gpl570infodf$"Gene ID"))

###  transform data into format compatible with the (paired) t.test function
toneparameter <- t(tone)
ttwoparameter <- t(ttwo)

## function that assigns output values (mean difference, p-value, standard error) from t.test function into a data frame format
tloop <- function(a,b) { 
  
  ttest <- t.test(a,b, paired = TRUE)
  
  ttestoutput <- data.frame( meandiff = ttest$estimate, pvalue = ttest$p.value, (stanerror = (ttest$conf.int[2] - ttest$estimate)/1.96))
  
  return(ttestoutput)
  
}
## runs tloop() -- performs t-test --  over all the columns in both twin groups
twinttestoutputdf <- sapply(c(1: ncol(toneparameter)) , function(a) tloop(toneparameter[,a], ttwoparameter[,a]))

meandiff16059 <- twinttestoutputdf[1,] # vector of mean differences for gse16059 study twins data
pvalue16059 <- twinttestoutputdf[2,]    # vector of pvalue
stanerror16059 <- twinttestoutputdf[3,]     # vector of standard error values

## extract probe labels
gse16059probelabels <- rownames(tone)

## column bind probe labels with t-test results
seconditergpl16059df <- cbind(gse16059probelabels,meandiff16059,pvalue16059,stanerror16059)

## transform data into compatible format for meta-analysis of probe-level values and bind entrez geneids (from annotation table) to t-test results data frame
thirditeringpl570probeorderdf <- seconditergpl16059df[match(probe_id_gse16059,seconditergpl16059df[,1]),]
thirditeringpl570probeorderdf <- as.data.frame(thirditeringpl570probeorderdf)
thirditeringpl570probeorderdf$gse16059probelabels <- as.character(thirditeringpl570probeorderdf$gse16059probelabels)
thirditeringpl570probeorderdf$meandiff16059 <- as.numeric(as.character(thirditeringpl570probeorderdf$meandiff16059))
thirditeringpl570probeorderdf$pvalue16059 <- as.numeric(as.character(thirditeringpl570probeorderdf$pvalue16059))
thirditeringpl570probeorderdf$stanerror16059 <- as.numeric(as.character(thirditeringpl570probeorderdf$stanerror16059))
gene_id_gse16059 <- as.data.frame(gene_id_gse16059)
forthiterplusgeneiddf <- cbind(thirditeringpl570probeorderdf[,1],gene_id_gse16059, thirditeringpl570probeorderdf[,2], thirditeringpl570probeorderdf[,3],thirditeringpl570probeorderdf[,4])
colnames(forthiterplusgeneiddf) <- c("probeid","geneid","meandiff","pval","stanerror")

## convert data frame into data table format
dttemp <- data.table(forthiterplusgeneiddf)
## split rows with multiple geneids corresponding to one probe id into separate rows
dttemp.out <- dttemp[, list(probeid=probeid , geneid = unlist(strsplit(as.character(geneid), "///")), meandiff =meandiff,stanerror=stanerror,pval=pval), by=1:nrow(dttemp)]
## convert into data frame format
secondtemp16059df <- data.frame(lapply(dttemp.out, as.character), stringsAsFactors=FALSE)

## get data frame with frequencies of geneids
testinputloop <- as.data.frame(table(secondtemp16059df$geneid))
## get unique geneids
testinputloop$Var1


## define meta-analytic function
getgenelevelidvalfunct <- function(a)
{
  secondtemp16059df.sub1 <- subset(secondtemp16059df, geneid == as.character(a) , select = c(meandiff, stanerror,geneid))
  metasumtest1 <- meta.summaries(as.numeric(secondtemp16059df.sub1$meandiff), as.numeric(secondtemp16059df.sub1$stanerror), method = "fixed", logscale = FALSE)
  assign(paste("tempgenevector",a,sep=""), c(metasumtest1$summary,metasumtest1$se.summary,pnorm(abs(metasumtest1$test[1]), lower.tail=F)*2))
  
}

## apply over all geneids 
genelevelvalgse16059_male <- lapply(testinputloop$Var1,getgenelevelidvalfunct)
## convert to data frame
genelevelvalgse16059df_male <- as.data.frame(genelevelvalgse16059_male)
## transpose data frame
genelevelvaluegse16059df_male <- t(genelevelvalgse16059df_male)
colnames(genelevelvaluegse16059df_male) <- c("meandiff", "stanerror","pval")
genelevelvaluegse16059df_male <- as.data.frame(genelevelvaluegse16059df_male)

## FDR-correction of p-values for each gene
FDR_GSE16059_male <- p.adjust(genelevelvaluegse16059df_male$pval, method="fdr")
## bind FDR values column to gene-level value results data frame
genelevelvaluegse16059df_male <- cbind(GENEID=testinputloop$Var1,genelevelvaluegse16059df_male,FDR_GSE16059_male )

## save data frame to your desired sub-directory
#saveRDS(genelevelvaluegse16059df_male, '/home/st320/second_iter_gse_genelevel_dfs/Male_genelevelvalue_plain_t_test_gse16059df.RDS')

