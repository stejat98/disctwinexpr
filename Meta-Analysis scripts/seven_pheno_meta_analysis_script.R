library(metafor)
library(ggplot2)
library(mygene)
library(dplyr)

## read data frame with gene-level differences for each of the 7 studies (phenotypes)
seven_pheno_studies <- readRDS('/home/st320/second_iter_gse_genelevel_dfs/seven_pheno_studies.RDS')

## define meta-analytic function to meta-analyze gene-level values across phenotypes
getgenelevelidvalfunctmetafor <- function(a)
{
  df.sub1 <- subset(seven_pheno_studies, GENEID == as.character(a) , select = c(meandiff, stanerror,GENEID,Pheno_Disease))
  metasumtest1 <- rma.uni(yi=as.numeric(df.sub1$meandiff), sei=as.numeric(df.sub1$stanerror), method = "DL")
  assign(paste("tempgenevector",a,sep=""), c(nrow(df.sub1),metasumtest1$b,metasumtest1$se,metasumtest1$pval,metasumtest1$I2,metasumtest1$QE,metasumtest1$QEp))
}

## get vector of unique geneids
uniquegenesidsacrossallstudies <- unique(seven_pheno_studies$GENEID)
## apply meta-analytic function for each gene
genelevelvalacrossallstudies <- lapply(uniquegenesidsacrossallstudies,getgenelevelidvalfunctmetafor)
## convert to data frame
genelevelvalueacrossallstudies <- as.data.frame(genelevelvalacrossallstudies)
## transpose data frame
genelevelvalacrossallstudiesdf <- t(genelevelvalueacrossallstudies)
colnames(genelevelvalacrossallstudiesdf) <- c("NumPhenotypes","meandiff","stanerror","pval", "I2","QE","QEp")

## column bind geneids vector to meta-analytic gene-level values data frame
genelevelvalacrossallstudiesdfmod <- cbind(as.numeric(as.character(uniquegenesidsacrossallstudies)), genelevelvalacrossallstudiesdf)
colnames(genelevelvalacrossallstudiesdfmod) <- c("GENEID","NumPhenotypes","meandiff","stanerror","pval", "I2","QE","QEp")

## convert to data frame
genelevelvalacrossallstudiesdfmod <- as.data.frame(genelevelvalacrossallstudiesdfmod)
rownames(genelevelvalacrossallstudiesdfmod) <- NULL

## get gene symbols
temphgncsymbol_all_studies_meta<- queryMany(unique(genelevelvalacrossallstudiesdfmod$GENEID),scopes="entrezgene",fields="symbol",species="human")
head(temphgncsymbol_all_studies_meta)
colnames(temphgncsymbol_all_studies_meta)[1] <- 'id'
genelevelvalacrossallstudiesdfmod$GENEID <- as.character(genelevelvalacrossallstudiesdfmod$GENEID)
genelevelvaluegseall_studies_meta_generaldf_plus_symbols <- left_join(genelevelvalacrossallstudiesdfmod,temphgncsymbol_all_studies_meta,copy= TRUE,by=c('GENEID'='id'))


## FDR correction for both overall p-values and QEp values

FDR <- p.adjust(genelevelvaluegseall_studies_meta_generaldf_plus_symbols$pval, method='fdr')
FDR_QEp <- p.adjust(genelevelvaluegseall_studies_meta_generaldf_plus_symbols$QEp, method='fdr')

## Column bind FDR and FDR_QEp vectors to meta-analytic results data frame
genelevelvaluegseall_studies_meta_generaldf_plus_symbols_updated <- cbind(genelevelvaluegseall_studies_meta_generaldf_plus_symbols, FDR= FDR, FDR_QEp =  FDR_QEp)

## remove unnecessary columns
genelevelvaluegseall_studies_meta_generaldf_plus_symbols_updated <-  genelevelvaluegseall_studies_meta_generaldf_plus_symbols_updated[,-c(9,11,12)]

## save data frame to your desired sub-directory
#saveRDS(genelevelvaluegseall_studies_meta_generaldf_plus_symbols_updated, '/home/st320/genelevelvaluegse_meta_7_phenotypes_plus_symbols_updated_08_07_17.RDS')

## read data frame
#genelevelvaluegseall_studies_meta_generaldf_plus_symbols_updated <- readRDS('/home/st320/genelevelvaluegse_meta_7_phenotypes_plus_symbols_updated_08_07_17.RDS')

## subset to find FDR (of mean difference) significant genes
metasig_genes_df <- subset(genelevelvaluegseall_studies_meta_generaldf_plus_symbols_updated, FDR < 0.05 & NumPhenotypes > 1)

## subset to find FDR(of QEp -- heterogeneity) significant genes
meta_genes_df <- subset(genelevelvaluegseall_studies_meta_generaldf_plus_symbols_updated, FDR_QEp < 0.05 & NumPhenotypes > 1)

## order genes based on increasing FDR values
ndx_metasig <- order(metasig_genes_df$FDR, decreasing = F)
temp_df_metasig <- metasig_genes_df[ndx_metasig,]

## find number of overall significant genes
nrow(subset(metasig_genes_df,!is.na(symbol)))

nrow(subset(meta_genes_df,!is.na(symbol)))

## read data frame
seven_pheno_studies <- readRDS('/home/st320/second_iter_gse_genelevel_dfs/seven_pheno_studies.RDS')

## check if meta-analytic sig genes appear signficant in individual phenotypes as well
subset(seven_pheno_studies, symbol == 'RPSA')
subset(seven_pheno_studies, symbol == 'CCT3')
subset(seven_pheno_studies, symbol == 'SDHD')

## get data frame of all overall significant genes found significant in within studies(pheno)
all_df <- data.frame()
for(i in 1:nrow(temp_df_metasig))
{
  temp <- subset(seven_pheno_studies, GENEID == temp_df_metasig$GENEID[i]) 
  temp_2 <- subset(temp, FDR < 0.05)
  all_df <- rbind(all_df,temp_2)
  
}

## get data frame of all overall significant genes found in phenotypes other than IQ
subset(all_df, Pheno_Disease != 'IQ')




## find most heterogeneous genes
max(temp_df_metasig$I2)

nrow(subset(temp_df_metasig, !is.na(symbol) & QEp < 0.05))

## observe least heterogeneous genes
i2_zero_sig_genes <- subset(temp_df_metasig, I2 == 0)

## get -log10 of FDR_QEp values (vector)
genelevelvaluegseall_studies_meta_generaldf_plus_symbols_updated$neglog_FDR_QEp <- -log10(genelevelvaluegseall_studies_meta_generaldf_plus_symbols_updated$FDR_QEp)

## plot -log10(FDR_QEp) versus. I2
ggplot(genelevelvaluegseall_studies_meta_generaldf_plus_symbols_updated, aes(x= I2, y= neglog_FDR_QEp))+geom_point() + ylab(expression('-log'[10] ~ '(FDR'[QEp] ~ ')')) + xlab(expression('I'^2))

## save plot to your desired sub-directory
#savePlot('/home/st320/FDR_QEp_I2_plot_08_07_17.png')

## Empirical CDF plot 08-07-17

ecdf(genelevelvaluegseall_studies_meta_generaldf_plus_symbols_updated$I2)

## plot ECDF for I2
ggplot(genelevelvaluegseall_studies_meta_generaldf_plus_symbols_updated, aes(I2)) + stat_ecdf(geom = "step", pad = "FALSE") + xlab(expression('I'^2)) + ylab('Cumulative Distribution Function')

## save plot to your desired sub-directory
#savePlot('/home/st320/Empirical_CDF_plot_08_07_17.png')
