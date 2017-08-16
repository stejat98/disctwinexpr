## Install & Load necesssary packages
#source("http://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")
require(gdata)
library(Biobase)
library(GEOquery)
library(rmeta)
require(data.table)

# No Activity      Active
#1 GSM509171   | GSM509172
#2 GSM509174   | GSM509173
#3 GSM509175   | GSM509176
#4 GSM509177   | GSM509178
#5 GSM509179   | GSM509180
#6 GSM509181   | GSM509182
#7 GSM509184   | GSM509183
#8 GSM509185   | GSM509186
#9 GSM509188   | GSM509187
#10GSM509189   | GSM509190

## list of GSMs (vectors) corresponding to no activity and activity respectively
GSM_diseased_GSE20319 <- c('GSM509171','GSM509174','GSM509175','GSM509177','GSM509179','GSM509181','GSM509184','GSM509185','GSM509188','GSM509189')
GSM_healthy_GSE20319 <- c('GSM509172','GSM509173','GSM509176','GSM509178','GSM509180','GSM509182','GSM509183','GSM509186','GSM509187','GSM509190')


## get data pertaining to GSE20319 
gse20319 <- getGEO("GSE20319", GSEMatrix = TRUE,AnnotGPL = TRUE)

## extract expression data table
gsm_tables_GSE20319 <- exprs(gse20319$'GSE20319_series_matrix.txt.gz')
gsm_tables_GSE20319 <- as.data.frame(gsm_tables_GSE20319)
#gsm20319_names <- gse20319$'GSE20319_series_matrix.txt.gz'@phenoData@data$geo_accession

## (no activity) group of samples
tone <- gsm_tables_GSE20319[,match(GSM_diseased_GSE20319,colnames(gsm_tables_GSE20319))]
## (activity) group of samples 
ttwo <- gsm_tables_GSE20319[,match(GSM_healthy_GSE20319,colnames(gsm_tables_GSE20319))]

#gsm20319names <- gse20319$GSE20319_series_matrix.txt.gz@phenoData@data$geo_accession

## get annotation information for GPL6884 platform
gpl6884info <- getGEO('GPL6884',AnnotGPL = TRUE)
str(gpl6884info)

## Extract annotation table for GPL6884 platform
gpl6884infodf <- gpl6884info@dataTable@table

## Extract microarray probe ids and assign them to a vector
probe_id_gse20319 <- as.character(t(gpl6884infodf$ID))
## Extract Entrez Gene IDs and assign them to a vector
gene_id_gse20319 <- as.character(t(gpl6884infodf$"Gene ID"))

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

meandiff20319 <- twinttestoutputdf[1,] # vector of mean differences for gse20319 study twins data
pvalue20319 <- twinttestoutputdf[2,]    # vector of pvalue
stanerror20319 <- twinttestoutputdf[3,]     # vector of standard error values

## extract probe labels
gse20319probelabels <- rownames(tone)

## column bind probe labels with t-test results
seconditergpl20319df <- cbind(gse20319probelabels,meandiff20319,pvalue20319,stanerror20319)

## transform data into compatible format for meta-analysis of probe-level values and bind entrez geneids (from annotation table) to t-test results data frame
thirditeringpl6884probeorderdf <- seconditergpl20319df[match(probe_id_gse20319,seconditergpl20319df[,1]),]
thirditeringpl6884probeorderdf <- as.data.frame(thirditeringpl6884probeorderdf)
thirditeringpl6884probeorderdf$gse20319probelabels <- as.character(thirditeringpl6884probeorderdf$gse20319probelabels)
thirditeringpl6884probeorderdf$meandiff20319 <- as.numeric(as.character(thirditeringpl6884probeorderdf$meandiff20319))
thirditeringpl6884probeorderdf$pvalue20319 <- as.numeric(as.character(thirditeringpl6884probeorderdf$pvalue20319))
thirditeringpl6884probeorderdf$stanerror20319 <- as.numeric(as.character(thirditeringpl6884probeorderdf$stanerror20319))
gene_id_gse20319 <- as.data.frame(gene_id_gse20319)
forthiterplusgeneiddf <- cbind(thirditeringpl6884probeorderdf[,1],gene_id_gse20319, thirditeringpl6884probeorderdf[,2], thirditeringpl6884probeorderdf[,3],thirditeringpl6884probeorderdf[,4])
colnames(forthiterplusgeneiddf) <- c("probeid","geneid","meandiff","pval","stanerror")

## convert data frame into data table format
dttemp <- data.table(forthiterplusgeneiddf)
## split rows with multiple geneids corresponding to one probe id into separate rows
dttemp.out <- dttemp[, list(probeid=probeid , geneid = unlist(strsplit(as.character(geneid), "///")), meandiff =meandiff,stanerror=stanerror,pval=pval), by=1:nrow(dttemp)]
## convert into data frame format
secondtemp20319df <- data.frame(lapply(dttemp.out, as.character), stringsAsFactors=FALSE)

## get data frame with frequencies of geneids
testinputloop <- as.data.frame(table(secondtemp20319df$geneid))
## get unique geneids
testinputloop$Var1

## define meta-analytic function
getgenelevelidvalfunct <- function(a)
{
  secondtemp20319df.sub1 <- subset(secondtemp20319df, geneid == as.character(a) , select = c(meandiff, stanerror,geneid))
  metasumtest1 <- meta.summaries(as.numeric(secondtemp20319df.sub1$meandiff), as.numeric(secondtemp20319df.sub1$stanerror), method = "fixed", logscale = FALSE)
  assign(paste("tempgenevector",a,sep=""), c(metasumtest1$summary,metasumtest1$se.summary,pnorm(abs(metasumtest1$test[1]), lower.tail=F)*2))
  
}

## apply over all geneids
genelevelvalgse20319 <- lapply(testinputloop$Var1,getgenelevelidvalfunct)
## convert to data frame
genelevelvalgse20319df <- as.data.frame(genelevelvalgse20319)
## transpose data frame
genelevelvaluegse20319df <- t(genelevelvalgse20319df)
colnames(genelevelvaluegse20319df) <- c("meandiff", "stanerror","pval")
genelevelvaluegse20319df <- as.data.frame(genelevelvaluegse20319df)

## FDR-correction of p-values for each gene
FDR_GSE20319 <- p.adjust(genelevelvaluegse20319df$pval, method="fdr")
## bind FDR values column to gene-level value results data frame
genelevelvaluegse20319df <- cbind(GENEID=testinputloop$Var1,genelevelvaluegse20319df,FDR_GSE20319)

## save data frame to your desired sub-directory
#saveRDS(genelevelvaluegse20319df, '/home/st320/second_iter_gse_genelevel_dfs/genelevelvalue_plain_t_test_gse20319df.RDS')
