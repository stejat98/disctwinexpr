## Install & Load necesssary packages
#source("http://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")
require(gdata)
library(Biobase)
library(GEOquery)
library(rmeta)
require(data.table)



## get data pertaining to GSE22619 
gse22619 <- getGEO("GSE22619", GSEMatrix = TRUE,AnnotGPL = TRUE)


## extract expression data table
gsm_tables_GSE22619 <- exprs(gse22619$'GSE22619_series_matrix.txt.gz')
gsm_tables_GSE22619 <- as.data.frame(gsm_tables_GSE22619)
#gsm22619_names <- gse22619$'GSE22619_series_matrix.txt.gz'@phenoData@data$geo_accession

## Diseased vs. Healthy GSMs
#1 GSM560969   | GSM560961
#2 GSM1012178  | GSM1012177
#3 GSM1012180  | GSM1012179
#4 GSM560970   | GSM560962
#6 GSM560971   | GSM560963
#7 GSM560972   | GSM560964
#8 GSM560973   | GSM560965
#9 GSM560974   | GSM560966
#10GSM560975   | GSM560967
#12GSM560976   | GSM560968


## list of GSMs (vectors) corresponding to diseasesd and healthy respectively
GSM_diseased_GSE22619 <- c('GSM560969','GSM1012178','GSM1012180','GSM560970','GSM560971','GSM560972','GSM560973','GSM560974','GSM560975','GSM560976')
GSM_healthy_GSE22619 <- c('GSM560961','GSM1012177','GSM1012179','GSM560962','GSM560963','GSM560964','GSM560965','GSM560966','GSM560967','GSM560968')

## affected group of samples
tone <- gsm_tables_GSE22619[,match(GSM_diseased_GSE22619,colnames(gsm_tables_GSE22619))]
## unaffected group of samples (control)
ttwo <- gsm_tables_GSE22619[,match(GSM_healthy_GSE22619,colnames(gsm_tables_GSE22619))]

## get annotation information for GPL570 platform
gpl570info <- getGEO('GPL570',AnnotGPL = TRUE)
str(gpl570info)

## Extract annotation table for GPL570 platform
gpl570infodf <- gpl570info@dataTable@table

## Extract microarray probe ids and assign them to a vector
probe_id_gse22619 <- as.character(t(gpl570infodf$ID))
## Extract Entrez Gene IDs and assign them to a vector
gene_id_gse22619 <- as.character(t(gpl570infodf$"Gene ID"))

###  transform data into format compatible with the (paired) t.test function
toneparameter <- t(tone)
ttwoparameter <- t(ttwo)

## function that assigns output values (mean difference, p-value, standard error) from t.test function into a data frame format
tloop <- function(a,b) { 
  
  ttest <- t.test(a,b, paired = TRUE)
  
  ttestoutput <- data.frame( meandiff = ttest$estimate, pvalue = ttest$p.value, (stanerror = (ttest$conf.int[2] - ttest$estimate)/1.96))
  
  return(ttestoutput)
  
}
## runs tloop() over all the columns in both twin groups
twinttestoutputdf <- sapply(c(1: ncol(toneparameter)) , function(a) tloop(toneparameter[,a], ttwoparameter[,a]))

meandiff22619 <- twinttestoutputdf[1,] # vector of mean differences for gse22619 study twins data
pvalue22619 <- twinttestoutputdf[2,]    # vector of pvalue
stanerror22619 <- twinttestoutputdf[3,]     # vector of standard error values

## extract probe labels
gse22619probelabels <- rownames(tone)

## column bind probe labels with t-test results
seconditergpl22619df <- cbind(gse22619probelabels,meandiff22619,pvalue22619,stanerror22619)

## transform data into compatible format for meta-analysis of probe-level values and bind entrez geneids (from annotation table) to t-test results data frame
thirditeringpl570probeorderdf <- seconditergpl22619df[match(probe_id_gse22619,seconditergpl22619df[,1]),]
thirditeringpl570probeorderdf <- as.data.frame(thirditeringpl570probeorderdf)
thirditeringpl570probeorderdf$gse22619probelabels <- as.character(thirditeringpl570probeorderdf$gse22619probelabels)
thirditeringpl570probeorderdf$meandiff22619 <- as.numeric(as.character(thirditeringpl570probeorderdf$meandiff22619))
thirditeringpl570probeorderdf$pvalue22619 <- as.numeric(as.character(thirditeringpl570probeorderdf$pvalue22619))
thirditeringpl570probeorderdf$stanerror22619 <- as.numeric(as.character(thirditeringpl570probeorderdf$stanerror22619))
gene_id_gse22619 <- as.data.frame(gene_id_gse22619)
forthiterplusgeneiddf <- cbind(thirditeringpl570probeorderdf[,1],gene_id_gse22619, thirditeringpl570probeorderdf[,2], thirditeringpl570probeorderdf[,3],thirditeringpl570probeorderdf[,4])
colnames(forthiterplusgeneiddf) <- c("probeid","geneid","meandiff","pval","stanerror")

## convert data frame into data table format
dttemp <- data.table(forthiterplusgeneiddf)
## split rows with multiple geneids corresponding to one probe id into separate rows
dttemp.out <- dttemp[, list(probeid=probeid , geneid = unlist(strsplit(as.character(geneid), "///")), meandiff =meandiff,stanerror=stanerror,pval=pval), by=1:nrow(dttemp)]
## convert into data frame format
secondtemp22619df <- data.frame(lapply(dttemp.out, as.character), stringsAsFactors=FALSE)

## get data frame with frequencies of geneids
testinputloop <- as.data.frame(table(secondtemp22619df$geneid))
## get unique geneids
testinputloop$Var1

## define meta-analytic function
getgenelevelidvalfunct <- function(a)
{
  secondtemp22619df.sub1 <- subset(secondtemp22619df, geneid == as.character(a) , select = c(meandiff, stanerror,geneid))
  metasumtest1 <- meta.summaries(as.numeric(secondtemp22619df.sub1$meandiff), as.numeric(secondtemp22619df.sub1$stanerror), method = "fixed", logscale = FALSE)
  assign(paste("tempgenevector",a,sep=""), c(metasumtest1$summary,metasumtest1$se.summary,pnorm(abs(metasumtest1$test[1]), lower.tail=F)*2))
  
}

## apply over all geneids 
genelevelvalgse22619 <- lapply(testinputloop$Var1,getgenelevelidvalfunct)
## convert to data frame
genelevelvalgse22619df <- as.data.frame(genelevelvalgse22619)
## transpose data frame
genelevelvaluegse22619df <- t(genelevelvalgse22619df)
colnames(genelevelvaluegse22619df) <- c("meandiff", "stanerror","pval")
genelevelvaluegse22619df <- as.data.frame(genelevelvaluegse22619df)


## FDR-correction of p-values for each gene
FDR_GSE22619 <- p.adjust(genelevelvaluegse22619df$pval, method="fdr")
## bind FDR values column to gene-level value results data frame
genelevelvaluegse22619df <- cbind(GENEID=testinputloop$Var1,genelevelvaluegse22619df,FDR_GSE22619)

## save data frame to your desired sub-directory
#saveRDS(genelevelvaluegse22619df, '/home/st320/second_iter_gse_genelevel_dfs/genelevelvalue_plain_t_test_gse22619df.RDS')

