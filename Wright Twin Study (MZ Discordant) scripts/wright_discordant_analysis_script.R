## Install & Load necesssary packages
#source("http://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")
require(gdata)
library(Biobase)
library(GEOquery)
library(stringr)
require(data.table)
library(rmeta)


## read table with sample info and expression values for each twin sample
twinsExpression <- readRDS('/groups/patel/dbgap/twins_depression/PhenoGenotypeFiles/ExpressionFiles/twins_expression.RDS')

## subset data frame to get only MZ pairs
MZ_twinsExpression_df <- subset(twinsExpression, TWIN > 0)

## subset to only get twins with MDD
MZ_MDD_twinsExpression_df <- subset(MZ_twinsExpression_df, MDD == 1)

## subset to only get twins without MDD/ controls
MZ_Healthy_twinsExpression_df <- subset(MZ_twinsExpression_df, MDD == 0)

## get list of discordant twins
discordant_twin_pair_ids <-  intersect(MZ_MDD_twinsExpression_df$TWIN,MZ_Healthy_twinsExpression_df$TWIN)

## get annotation information for GPL13667 platform
gpl13667info <- getGEO('GPL13667',AnnotGPL = TRUE)
str(gpl13667info)

## Extract annotation table for GPL13667 platform
gpl13667infodf <- gpl13667info@dataTable@table



## Extract microarray probe ids and assign them to a vector
probe_id_gsewright_discordant_general <- as.character(t(gpl13667infodf$ID))

probe_id_gsewright_discordant_general <- unlist(lapply(as.character(t(gpl13667infodf$ID)), function(a){paste('x',a,sep="")}))
## Extract Entrez Gene IDs and assign them to a vector
gene_id_gsewright_discordant_general <- as.character(t(gpl13667infodf$"Entrez Gene"))



## subset to only include discordant twin samples
discordant_MZ_MDD_twinsExpression_df <- MZ_MDD_twinsExpression_df[match(discordant_twin_pair_ids,MZ_MDD_twinsExpression_df$TWIN),]

discordant_MZ_Healthy_twinsExpression_df <- MZ_Healthy_twinsExpression_df[match(discordant_twin_pair_ids,MZ_Healthy_twinsExpression_df$TWIN),]

## include only probe exression value columns
toneparameter <- discordant_MZ_MDD_twinsExpression_df[,grep('^x', colnames(twinsExpression))]

ttwoparameter <- discordant_MZ_Healthy_twinsExpression_df[,grep('^x', colnames(twinsExpression))]

## function that assigns output values (mean difference, p-value, standard error) from t.test function into a data frame format
tloop <- function(a,b) { 
  
  ttest <- t.test(a,b, paired = TRUE)
  
  ttestoutput <- data.frame( meandiff = ttest$estimate, pvalue = ttest$p.value, (stanerror = (ttest$conf.int[2] - ttest$estimate)/1.96))
  
  return(ttestoutput)
  
}
## runs tloop() over all the columns in both twin groups
twinttestoutputdf <- sapply(c(1: ncol(toneparameter)) , function(a) tloop(toneparameter[,a], ttwoparameter[,a]))

meandiffwright_discordant_general <- twinttestoutputdf[1,] # vector of mean differences for gsewright_discordant_general study twins data
pvaluewright_discordant_general <- twinttestoutputdf[2,]    # vector of pvalue
stanerrorwright_discordant_general <- twinttestoutputdf[3,]     # vector of standard error values

## get probe-level ids
gsewright_discordant_generalprobelabels <- colnames(toneparameter)
## create dataframe with probe-level ids and corresponding mean difference, p-value, and standard error values
seconditergplwright_discordant_generaldf <- cbind(gsewright_discordant_generalprobelabels,meandiffwright_discordant_general,pvaluewright_discordant_general,stanerrorwright_discordant_general)


## format problematic probe-level ids properly

temp_1_df <- gsub("xAFFX-HUMGAPDH/M33197_3_","xAFFX-HUMGAPDH/M33197_3_at",seconditergplwright_discordant_generaldf[,1])
temp_1_df <- gsub("xAFFX-HUMGAPDH/M33197_5_","xAFFX-HUMGAPDH/M33197_5_at",temp_1_df)
temp_1_df <- gsub("xAFFX-HUMGAPDH/M33197_M_","xAFFX-HUMGAPDH/M33197_M_at",temp_1_df)
temp_1_df <- gsub("xAFFX-HUMISGF3A/M97935_3","xAFFX-HUMISGF3A/M97935_3_at",temp_1_df)
temp_1_df <- gsub("xAFFX-HUMISGF3A/M97935_5","xAFFX-HUMISGF3A/M97935_5_at",temp_1_df)
temp_1_df <- gsub("^xAFFX-HUMISGF3A/M97935_M$","xAFFX-HUMISGF3A/M97935_MA_at",temp_1_df)
temp_1_df <- gsub("^xAFFX-HUMISGF3A/M97935_M.1$","xAFFX-HUMISGF3A/M97935_MB_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC03_a","xAFFX-Nonspecific-GC03_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC04_a","xAFFX-Nonspecific-GC04_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC05_a","xAFFX-Nonspecific-GC05_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC06_a","xAFFX-Nonspecific-GC06_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC07_a","xAFFX-Nonspecific-GC07_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC08_a","xAFFX-Nonspecific-GC08_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC09_a","xAFFX-Nonspecific-GC09_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC10_a","xAFFX-Nonspecific-GC10_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC11_a","xAFFX-Nonspecific-GC11_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC12_a","xAFFX-Nonspecific-GC12_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC13_a","xAFFX-Nonspecific-GC13_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC14_a","xAFFX-Nonspecific-GC14_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC15_a","xAFFX-Nonspecific-GC15_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC16_a","xAFFX-Nonspecific-GC16_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC17_a","xAFFX-Nonspecific-GC17_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC18_a","xAFFX-Nonspecific-GC18_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC19_a","xAFFX-Nonspecific-GC19_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC20_a","xAFFX-Nonspecific-GC20_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC21_a","xAFFX-Nonspecific-GC21_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC22_a","xAFFX-Nonspecific-GC22_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC23_a","xAFFX-Nonspecific-GC23_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC24_a","xAFFX-Nonspecific-GC24_at",temp_1_df)
temp_1_df <- gsub("xAFFX-Nonspecific-GC25_a","xAFFX-Nonspecific-GC25_at",temp_1_df)

## reassign formatted probe-level ids 
seconditergplwright_discordant_generaldf[,1] <- temp_1_df

## transform data into compatible format for meta-analysis of probe-level values and bind entrez geneids (from annotation table) to t-test results data frame
thirditeringpl13667probeorderdf <- seconditergplwright_discordant_generaldf[match(probe_id_gsewright_discordant_general,seconditergplwright_discordant_generaldf[,1]),]
thirditeringpl13667probeorderdf <- as.data.frame(thirditeringpl13667probeorderdf)
thirditeringpl13667probeorderdf$gsewright_discordant_generalprobelabels <- as.character(thirditeringpl13667probeorderdf$gsewright_discordant_generalprobelabels)
thirditeringpl13667probeorderdf$meandiffwright_discordant_general <- as.numeric(as.character(thirditeringpl13667probeorderdf$meandiffwright_discordant_general))
thirditeringpl13667probeorderdf$pvaluewright_discordant_general <- as.numeric(as.character(thirditeringpl13667probeorderdf$pvaluewright_discordant_general))
thirditeringpl13667probeorderdf$stanerrorwright_discordant_general <- as.numeric(as.character(thirditeringpl13667probeorderdf$stanerrorwright_discordant_general))
gene_id_gsewright_discordant_general <- as.data.frame(gene_id_gsewright_discordant_general)
forthiterplusgeneiddf <- cbind(thirditeringpl13667probeorderdf[,1],gene_id_gsewright_discordant_general, thirditeringpl13667probeorderdf[,2], thirditeringpl13667probeorderdf[,3],thirditeringpl13667probeorderdf[,4])
colnames(forthiterplusgeneiddf) <- c("probeid","geneid","meandiff","pval","stanerror")
## remove the 1708 probes the study didn't use (or for which there are no values recorded)
fifthiterplusgeneiddf <- forthiterplusgeneiddf[-which(is.na(forthiterplusgeneiddf$meandiff)),]

## convert data frame into data table format
dttemp <- data.table(fifthiterplusgeneiddf)
## split rows with multiple geneids corresponding to one probe id into separate rows
dttemp.out <- dttemp[, list(probeid=probeid , geneid = unlist(strsplit(as.character(geneid), "///")), meandiff =meandiff,stanerror=stanerror,pval=pval), by=1:nrow(dttemp)]
## convert into data frame format
secondtempwright_discordant_generaldf <- data.frame(lapply(dttemp.out, as.character), stringsAsFactors=FALSE)

secondtempwright_discordant_generaldf$geneid <- str_trim(secondtempwright_discordant_generaldf$geneid)
## get data frame with frequencies of geneids
testinputloop <- as.data.frame(table(secondtempwright_discordant_generaldf$geneid))
## get unique geneids
testinputloop$Var1


## define function to meta-analyze over probe-level values corresponding to each gene to obtain gene-level values
getgenelevelidvalfunct <- function(a)
{
  secondtempwright_discordant_generaldf.sub1 <- subset(secondtempwright_discordant_generaldf, geneid == as.character(a) , select = c(meandiff, stanerror,geneid))
  metasumtest1 <- meta.summaries(as.numeric(secondtempwright_discordant_generaldf.sub1$meandiff), as.numeric(secondtempwright_discordant_generaldf.sub1$stanerror), method = "fixed", logscale = FALSE)
  assign(paste("tempgenevector",a,sep=""), c(metasumtest1$summary,metasumtest1$se.summary,pnorm(abs(metasumtest1$test[1]), lower.tail=F)*2))
  
}

## apply meta-analytic function over each gene id
genelevelvalgsewright_discordant_general <- lapply(testinputloop$Var1,getgenelevelidvalfunct)

## convert to data frame
genelevelvalgsewright_discordant_generaldf <- as.data.frame(genelevelvalgsewright_discordant_general)
## transpose data frame
genelevelvaluegsewright_discordant_generaldf <- t(genelevelvalgsewright_discordant_generaldf)
colnames(genelevelvaluegsewright_discordant_generaldf) <- c("meandiff", "stanerror","pval")

## convert to data frame
genelevelvaluegsewright_discordant_generaldf <- as.data.frame(genelevelvaluegsewright_discordant_generaldf)
## FDR-correct p-values for all genes
FDR_GSEwright_discordant_general <- p.adjust(genelevelvaluegsewright_discordant_generaldf$pval, method="fdr")

## column bind the FDR values vector and geneid vector to gene-level values dataframe
genelevelvaluegsewright_discordant_generaldf <- cbind(GENEID=testinputloop$Var1,genelevelvaluegsewright_discordant_generaldf,FDR_GSEwright_discordant_general )
## remove non-entrez gene ids
genelevelvaluegsewright_discordant_generaldf <- genelevelvaluegsewright_discordant_generaldf[-c(1,19286:19311),]

## find number of FDR significant genes
sum(genelevelvaluegsewright_discordant_generaldf[,5] < 0.05)

## save dataframe to sub-directory
#saveRDS(genelevelvaluegsewright_discordant_generaldf, '/home/st320/second_iter_gse_genelevel_dfs/genelevelvalue_plain_t_test_gsewright_discordant_generaldf.RDS')

