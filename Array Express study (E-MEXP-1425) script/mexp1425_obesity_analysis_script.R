## Install & Load necesssary packages
#source("http://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")
require(gdata)
library(Biobase)
library(GEOquery)
library('ArrayExpress')
library(rmeta)
require(data.table)


mexp1425 = getAE("E-MEXP-1425", type = "full")
## Build a an ExpressionSet from the raw data
MEXP1425raw = ae2bioc(mageFiles = mexp1425)
## Build a an ExpressionSet from the processed data
cnames = getcolproc(mexp1425)
MEXP1425proc = procset(mexp1425, cnames[2])

## used to extract phenotype/disease diagnosis corresponding to each twin sample
phenodiagnosis_MEXP1425 <- MEXP1425proc@phenoData@data$Characteristics..disease.

## used to find which two samples should use from triplets -- not using "Fat twin 9 slimmer"
add_pheno <- MEXP1425proc@phenoData@data$Sample.Name 

twp_num <- MEXP1425proc@phenoData@data$Characteristics..family.relationship.
sex_group <- MEXP1425proc@phenoData@data$Characteristics..sex. 
Age <- MEXP1425proc@phenoData@data$Characteristics..age.

gsub('twin pair:','',twp_num)

gsub('sex: ','',sex_group)


grep('normal', phenodiagnosis_MEXP1425)
grep('obesity', phenodiagnosis_MEXP1425)

## create vector of 0s and 1st indicating phenotype status
Healthy0.Obesity1 <- vector(mode="integer", length=27)
Healthy0.Obesity1[grep('normal', phenodiagnosis_MEXP1425)] <- 0
Healthy0.Obesity1[grep('obesity', phenodiagnosis_MEXP1425)] <- 1

## extract expression tables
gsm_tables_MEXP1425 <- exprs(MEXP1425proc)
gsm_tables_MEXP1425 <- as.data.frame(gsm_tables_MEXP1425)

## Table with information (sex, age, twin pair number, etc.) for each sample
Sample_Info <- cbind(Colnames_eset= colnames(gsm_tables_MEXP1425),Healthy0.Obesity1, add_pheno,Phenodiagnosis=phenodiagnosis_MEXP1425,TWP_NUM = gsub('monozygotic twin pair ','',twp_num),Sex_Group = sex_group,Age)
Sample_Info <- as.data.frame(Sample_Info)

# remove third twin in TWP 9 as it consist of a triplet -- only want to look at pairs
Sample_Info <- Sample_Info[-14,]

## separate Sample_Info table into two tables: one only for affected samples and one for the unaffected
healthy_samples <- subset(Sample_Info, Healthy0.Obesity1 ==0)
diseased_samples <- subset(Sample_Info, Healthy0.Obesity1 == 1)

diseased_samples <- diseased_samples[match(healthy_samples$TWP_NUM,diseased_samples$TWP_NUM),]

## affected group of samples
tone <- gsm_tables_MEXP1425[,match(diseased_samples$Colnames_eset,colnames(gsm_tables_MEXP1425))]

## unaffected group of samples (control)
ttwo <- gsm_tables_MEXP1425[,match(healthy_samples$Colnames_eset,colnames(gsm_tables_MEXP1425))]



## get annotation information for GPL570 platform
gpl570info <- getGEO('GPL570',AnnotGPL = TRUE)
str(gpl570info)

## Extract annotation table for GPL570 platform
gpl570infodf <- gpl570info@dataTable@table



## Extract microarray probe ids and assign them to a vector
probe_id_MEXP1425 <- as.character(t(gpl570infodf$ID))
## Extract Entrez Gene IDs and assign them to a vector
gene_id_MEXP1425 <- as.character(t(gpl570infodf$"Gene ID"))

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

meandiffMEXP1425 <- twinttestoutputdf[1,] # vector of mean differences for MEXP1425 study twins data
pvalueMEXP1425 <- twinttestoutputdf[2,]    # vector of pvalue
stanerrorMEXP1425 <- twinttestoutputdf[3,]     # vector of standard error values

## extract probe labels
MEXP1425probelabels <- rownames(tone)

## column bind probe labels with t-test results
seconditergplMEXP1425df <- cbind(MEXP1425probelabels,meandiffMEXP1425,pvalueMEXP1425,stanerrorMEXP1425)

## transform data into compatible format for meta-analysis of probe-level values and bind entrez geneids (from annotation table) to t-test results data frame
thirditeringpl570probeorderdf <- seconditergplMEXP1425df[match(probe_id_MEXP1425,seconditergplMEXP1425df[,1]),]
thirditeringpl570probeorderdf <- as.data.frame(thirditeringpl570probeorderdf)
thirditeringpl570probeorderdf$MEXP1425probelabels <- as.character(thirditeringpl570probeorderdf$MEXP1425probelabels)
thirditeringpl570probeorderdf$meandiffMEXP1425 <- as.numeric(as.character(thirditeringpl570probeorderdf$meandiffMEXP1425))
thirditeringpl570probeorderdf$pvalueMEXP1425 <- as.numeric(as.character(thirditeringpl570probeorderdf$pvalueMEXP1425))
thirditeringpl570probeorderdf$stanerrorMEXP1425 <- as.numeric(as.character(thirditeringpl570probeorderdf$stanerrorMEXP1425))
gene_id_MEXP1425 <- as.data.frame(gene_id_MEXP1425)
forthiterplusgeneiddf <- cbind(thirditeringpl570probeorderdf[,1],gene_id_MEXP1425, thirditeringpl570probeorderdf[,2], thirditeringpl570probeorderdf[,3],thirditeringpl570probeorderdf[,4])
colnames(forthiterplusgeneiddf) <- c("probeid","geneid","meandiff","pval","stanerror")

## convert data frame into data table format
dttemp <- data.table(forthiterplusgeneiddf)
## split rows with multiple geneids corresponding to one probe id into separate rows
dttemp.out <- dttemp[, list(probeid=probeid , geneid = unlist(strsplit(as.character(geneid), "///")), meandiff =meandiff,stanerror=stanerror,pval=pval), by=1:nrow(dttemp)]
## convert into data frame format
secondtempMEXP1425df <- data.frame(lapply(dttemp.out, as.character), stringsAsFactors=FALSE)

## get data frame with frequencies of geneids
testinputloop <- as.data.frame(table(secondtempMEXP1425df$geneid))
## get unique geneids
testinputloop$Var1

## define meta-analytic function
getgenelevelidvalfunct <- function(a)
{
  secondtempMEXP1425df.sub1 <- subset(secondtempMEXP1425df, geneid == as.character(a) , select = c(meandiff, stanerror,geneid))
  metasumtest1 <- meta.summaries(as.numeric(secondtempMEXP1425df.sub1$meandiff), as.numeric(secondtempMEXP1425df.sub1$stanerror), method = "fixed", logscale = FALSE)
  assign(paste("tempgenevector",a,sep=""), c(metasumtest1$summary,metasumtest1$se.summary,pnorm(abs(metasumtest1$test[1]), lower.tail=F)*2))
  
}

## apply over all geneids 
genelevelvalMEXP1425 <- lapply(testinputloop$Var1,getgenelevelidvalfunct)
## convert to data frame
genelevelvalMEXP1425df <- as.data.frame(genelevelvalMEXP1425)
## transpose data frame
genelevelvalueMEXP1425df <- t(genelevelvalMEXP1425df)
colnames(genelevelvalueMEXP1425df) <- c("meandiff", "stanerror","pval")
genelevelvalueMEXP1425df <- as.data.frame(genelevelvalueMEXP1425df)

## FDR-correction of p-values for each gene
FDR_MEXP1425 <- p.adjust(genelevelvalueMEXP1425df$pval, method="fdr")

## bind FDR values column to gene-level value results data frame
genelevelvalueMEXP1425df <- cbind(GENEID=testinputloop$Var1,genelevelvalueMEXP1425df,FDR_MEXP1425)

## save data frame to your desired sub-directory 
#saveRDS(genelevelvalueMEXP1425df, '/home/st320/second_iter_gse_genelevel_dfs/genelevelvalue_plain_t_test_MEXP1425df.RDS')

