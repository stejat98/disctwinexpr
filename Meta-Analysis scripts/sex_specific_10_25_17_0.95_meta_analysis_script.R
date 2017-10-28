## load necessary (R) packages
library(metafor)
library(mygene)
library(dplyr)




## male-specific meta-analysis

male_genelevelvaluegsewright_discordant_generaldf <- readRDS('/home/st320/second_iter_gse_genelevel_dfs/male_genelevelvalue_plain_t_test_gsewright_discordant_generaldf.RDS')
genelevelvaluegse16059df_male <- readRDS('/home/st320/second_iter_gse_genelevel_dfs/Male_genelevelvalue_plain_t_test_gse16059df.RDS')
genelevelvalueMEXP1425df_male <- readRDS('/home/st320/second_iter_gse_genelevel_dfs/male_genelevelvalue_plain_t_test_MEXP1425df.RDS')


colnames(male_genelevelvaluegsewright_discordant_generaldf)[5] <- 'FDR'
colnames(genelevelvaluegse16059df_male)[5] <- 'FDR'
colnames(genelevelvalueMEXP1425df_male)[5] <- 'FDR'

## create phenotype label vector
Phenotype_MDD_twins <- rep('MDD',each = nrow(male_genelevelvaluegsewright_discordant_generaldf))
Phenotype_CFS_twins <- rep('CFS', each = nrow(genelevelvaluegse16059df_male))
Phenotype_OB_twins <- rep('OB', each = nrow(genelevelvalueMEXP1425df_male))
Phenotype <- c(Phenotype_MDD_twins,Phenotype_CFS_twins,Phenotype_OB_twins)

## row bind all the male-specific gene-level differences for each of the phenotypes 
all_3_male_studies_estimates_comb_df <- rbind(male_genelevelvaluegsewright_discordant_generaldf, genelevelvaluegse16059df_male, genelevelvalueMEXP1425df_male)
## column bind phenotype label vector to the gene-level differences data frame
all_3_male_studies_estimates_comb_df <- cbind(all_3_male_studies_estimates_comb_df, Phenotype)

## meta-analytic function to meta-analyze the male-specific values for each gene across all 3 phenotypes
getgenelevelidvalfunctmetafor <- function(a)
{
  df.sub1 <- subset(all_3_male_studies_estimates_comb_df, GENEID == as.character(a) , select = c(meandiff, stanerror,GENEID,Phenotype))
  metasumtest1 <- rma.uni(yi=as.numeric(df.sub1$meandiff), sei=as.numeric(df.sub1$stanerror), method = "DL")
  assign(paste("tempgenevector",a,sep=""), c(nrow(df.sub1),metasumtest1$b,metasumtest1$se,metasumtest1$pval,metasumtest1$I2,metasumtest1$QE,metasumtest1$QEp))
}

## get vector of all unique gene ids
male_uniquegenesidsacrossallstudies <- unique(all_3_male_studies_estimates_comb_df$GENEID)
## run meta-analysis for each gene
male_genelevelvalacrossallstudies <- lapply(male_uniquegenesidsacrossallstudies,getgenelevelidvalfunctmetafor)
## convert results to data frame 
male_genelevelvalueacrossallstudies <- as.data.frame(male_genelevelvalacrossallstudies)
## transpose data frame
male_genelevelvalacrossallstudiesdf <- t(male_genelevelvalueacrossallstudies)
colnames(male_genelevelvalacrossallstudiesdf) <- c("NumPhenotypes","meandiff","stanerror","pval", "I2","QE","QEp")

## column bind vector of geneids to results data
male_genelevelvalacrossallstudiesdfmod <- cbind(as.numeric(as.character(male_uniquegenesidsacrossallstudies)), male_genelevelvalacrossallstudiesdf)
colnames(male_genelevelvalacrossallstudiesdfmod) <- c("GENEID","NumPhenotypes","meandiff","stanerror","pval", "I2","QE","QEp")
## convert to data frame
male_genelevelvalacrossallstudiesdfmod <- as.data.frame(male_genelevelvalacrossallstudiesdfmod)

## save data frame to desired sub-directory
#saveRDS(male_genelevelvalacrossallstudiesdfmod, '/home/st320/sex_specific_meta_analysis/updated_08_08_17_male_genelevel_diffs_meta.RDS')




## female-specific meta-analysis


female_genelevelvaluegsewright_discordant_generaldf <- readRDS('/home/st320/second_iter_gse_genelevel_dfs/female_genelevelvalue_plain_t_test_gsewright_discordant_generaldf.RDS')
genelevelvaluegse16059df_female <- readRDS('/home/st320/second_iter_gse_genelevel_dfs/Female_genelevelvalue_plain_t_test_gse16059df.RDS')
genelevelvalueMEXP1425df_female <- readRDS('/home/st320/second_iter_gse_genelevel_dfs/female_genelevelvalue_plain_t_test_MEXP1425df.RDS')
colnames(female_genelevelvaluegsewright_discordant_generaldf)[5] <- 'FDR'
colnames(genelevelvaluegse16059df_female)[5] <- 'FDR'
colnames(genelevelvalueMEXP1425df_female)[5] <- 'FDR'

## create phenotype vector 
Phenotype_MDD_twins_female <- rep('MDD',each = nrow(female_genelevelvaluegsewright_discordant_generaldf))
Phenotype_CFS_twins_female <- rep('CFS', each = nrow(genelevelvaluegse16059df_female))
Phenotype_OB_twins_female <- rep('OB', each = nrow(genelevelvalueMEXP1425df_female))
Phenotype_female <- c(Phenotype_MDD_twins_female,Phenotype_CFS_twins_female,Phenotype_OB_twins_female)

## row bind all of the female-specific gene-level differences for each of the phenotypes
all_3_female_studies_estimates_comb_df <- rbind(female_genelevelvaluegsewright_discordant_generaldf, genelevelvaluegse16059df_female, genelevelvalueMEXP1425df_female)
## column bind phenotype label vector to the gene-level differences data frame
all_3_female_studies_estimates_comb_df <- cbind(all_3_female_studies_estimates_comb_df, Phenotype = Phenotype_female)

## meta-analytic function to meta-analyze the female-specific values for each gene across all 3 phenotypes
getgenelevelidvalfunctmetafor <- function(a)
{
  df.sub1 <- subset(all_3_female_studies_estimates_comb_df, GENEID == as.character(a) , select = c(meandiff, stanerror,GENEID,Phenotype))
  metasumtest1 <- rma.uni(yi=as.numeric(df.sub1$meandiff), sei=as.numeric(df.sub1$stanerror), method = "DL")
  assign(paste("tempgenevector",a,sep=""), c(nrow(df.sub1),metasumtest1$b,metasumtest1$se,metasumtest1$pval,metasumtest1$I2,metasumtest1$QE,metasumtest1$QEp))
}

## get vector of all unique geneids
female_uniquegenesidsacrossallstudies <- unique(all_3_female_studies_estimates_comb_df$GENEID)
## run meta-analysis for each gene
female_genelevelvalacrossallstudies <- lapply(female_uniquegenesidsacrossallstudies,getgenelevelidvalfunctmetafor)
## convert to data frame
female_genelevelvalueacrossallstudies <- as.data.frame(female_genelevelvalacrossallstudies)
## transpose data frame
female_genelevelvalacrossallstudiesdf <- t(female_genelevelvalueacrossallstudies)
colnames(female_genelevelvalacrossallstudiesdf) <- c("NumPhenotypes","meandiff","stanerror","pval", "I2","QE","QEp")

## column bind vector of geneids to results data
female_genelevelvalacrossallstudiesdfmod <- cbind(as.numeric(as.character(female_uniquegenesidsacrossallstudies)), female_genelevelvalacrossallstudiesdf)
colnames(female_genelevelvalacrossallstudiesdfmod) <- c("GENEID","NumPhenotypes","meandiff","stanerror","pval", "I2","QE","QEp")

## convert to data frame 
female_genelevelvalacrossallstudiesdfmod <- as.data.frame(female_genelevelvalacrossallstudiesdfmod)

## save data frame to desired sub-directory
#saveRDS(female_genelevelvalacrossallstudiesdfmod, '/home/st320/sex_specific_meta_analysis/updated_08_08_17_female_genelevel_diffs_meta.RDS')

## read meta-analytic RDS files for female and male groups
female_genelevelvalacrossallstudiesdfmod <- readRDS('/home/st320/sex_specific_meta_analysis/updated_08_08_17_female_genelevel_diffs_meta.RDS')
male_genelevelvalacrossallstudiesdfmod <-  readRDS('/home/st320/sex_specific_meta_analysis/updated_08_08_17_male_genelevel_diffs_meta.RDS')


## female (get gene symbols)

temphgncsymbol_all_studies_meta_female<- queryMany(unique(female_genelevelvalacrossallstudiesdfmod$GENEID),scopes="entrezgene",fields="symbol",species="human")
head(temphgncsymbol_all_studies_meta_female)
colnames(temphgncsymbol_all_studies_meta_female)[1] <- 'id'
female_genelevelvalacrossallstudiesdfmod$GENEID <- as.character(female_genelevelvalacrossallstudiesdfmod$GENEID)
female_genelevelvalues_meta_generaldf_plus_symbols <- left_join(female_genelevelvalacrossallstudiesdfmod,temphgncsymbol_all_studies_meta_female,copy= TRUE,by=c('GENEID'='id'))
female_genelevelvalues_meta_generaldf_plus_symbols <-  female_genelevelvalues_meta_generaldf_plus_symbols[,-c(9,11,12)]
colnames(female_genelevelvalues_meta_generaldf_plus_symbols)[2] <- 'NumPhenotypes'


## male (get gene symbols)

temphgncsymbol_all_studies_meta_male<- queryMany(unique(male_genelevelvalacrossallstudiesdfmod$GENEID),scopes="entrezgene",fields="symbol",species="human")
head(temphgncsymbol_all_studies_meta_male)
colnames(temphgncsymbol_all_studies_meta_male)[1] <- 'id'
male_genelevelvalacrossallstudiesdfmod$GENEID <- as.character(male_genelevelvalacrossallstudiesdfmod$GENEID)
male_genelevelvalues_meta_generaldf_plus_symbols <- left_join(male_genelevelvalacrossallstudiesdfmod,temphgncsymbol_all_studies_meta_male,copy= TRUE,by=c('GENEID'='id'))
male_genelevelvalues_meta_generaldf_plus_symbols <-  male_genelevelvalues_meta_generaldf_plus_symbols[,-c(9,11,12)]
colnames(male_genelevelvalues_meta_generaldf_plus_symbols)[2] <- 'NumPhenotypes'

## FDR-correction for p-values of mean differences for females
FDR_female <- p.adjust(female_genelevelvalues_meta_generaldf_plus_symbols$pval, method='fdr')
## FDR-correction for QEp values (significance of I2 -- heterogeneity) for females
FDR_QEp_female <- p.adjust(female_genelevelvalues_meta_generaldf_plus_symbols$QEp, method='fdr')

## FDR-correction for p-values of mean differences for males
FDR_male <- p.adjust(male_genelevelvalues_meta_generaldf_plus_symbols$pval, method='fdr')
## FDR-correction for QEp values (significance of I2 -- heterogeneity) for males
FDR_QEp_male <- p.adjust(male_genelevelvalues_meta_generaldf_plus_symbols$QEp, method='fdr')

## column bind FDR_female and FDR_QEp_female columns to female meta-analytic data frame
female_genelevelvalues_meta_generaldf_plus_symbols_updated <- cbind(female_genelevelvalues_meta_generaldf_plus_symbols, FDR= FDR_female, FDR_QEp = FDR_QEp_female)


## find number of overall (across multiple phenotypes) FDR (of mean difference) significant genes 
sum(female_genelevelvalues_meta_generaldf_plus_symbols_updated$FDR<0.05)

## save (female) data frame to your desired sub-directory
#saveRDS(female_genelevelvalues_meta_generaldf_plus_symbols_updated,'/home/st320/sex_specific_meta_analysis/female_genelevelvalues_meta_generaldf_plus_symbols_08_08_17.RDS')

## column bind FDR_male and FDR_QEp_male columns to male meta-anlytic data frame
male_genelevelvalues_meta_generaldf_plus_symbols_updated <- cbind(male_genelevelvalues_meta_generaldf_plus_symbols, FDR= FDR_male, FDR_QEp = FDR_QEp_male)

#sum(male_genelevelvalues_meta_generaldf_plus_symbols_updated$pval<0.05)

## find number of overall (across multiple phenotypes) FDR (of mean difference) significant genes
sum(male_genelevelvalues_meta_generaldf_plus_symbols_updated$FDR<0.05)

## save (male) data frame to your desired sub-directory
#saveRDS(male_genelevelvalues_meta_generaldf_plus_symbols_updated,'/home/st320/sex_specific_meta_analysis/male_genelevelvalues_meta_generaldf_plus_symbols_08_08_17.RDS')


## read updated meta-analytic RDS files for female and male groups
female_genelevelvalues_meta_generaldf_plus_symbols_updated <- readRDS('/home/st320/sex_specific_meta_analysis/female_genelevelvalues_meta_generaldf_plus_symbols_08_08_17.RDS')
male_genelevelvalues_meta_generaldf_plus_symbols_updated <- readRDS('/home/st320/sex_specific_meta_analysis/male_genelevelvalues_meta_generaldf_plus_symbols_08_08_17.RDS')

## male info -- find significant genes (in terms of differentially expressed, most heterogeneous, etc.)

male_metasig_genes_df_FDR <- subset(male_genelevelvalues_meta_generaldf_plus_symbols_updated, FDR < 0.05 & NumPhenotypes > 1 & abs(meandiff) > quantile(abs(male_genelevelvalues_meta_generaldf_plus_symbols_updated$meandiff),0.95))
male_ndx_metasig_FDR <- order(male_metasig_genes_df_FDR$FDR, decreasing = F)
male_temp_df_metasig_FDR <- male_metasig_genes_df_FDR[male_ndx_metasig_FDR,] 
length(which(!is.na(male_temp_df_metasig_FDR$symbol))) ## 9 genes without na symbols

#max(male_temp_df_metasig_FDR$I2)  ## max I2 value = 0

#nrow(subset(male_temp_df_metasig_FDR, I2 > 99.9))

i2_zero_sig_genes <- subset(male_temp_df_metasig_FDR, I2 == 0)

nrow(subset(male_temp_df_metasig_FDR, I2== 0))

## find significant genes that also have I2 of 0
i2_zero_sig_genes_male <- subset(male_temp_df_metasig_FDR, I2 == 0 & !is.na(symbol)) ## 9



## female info -- find significant genes (in terms of differentially expressed, most heterogeneous, etc.)
#female_genelevelvalues_meta_generaldf_plus_symbols_updated <- readRDS('/home/st320/sex_specific_meta_analysis/female_genelevelvalues_meta_generaldf_plus_symbols_08_08_17.RDS')
female_metasig_genes_df_FDR <- subset(female_genelevelvalues_meta_generaldf_plus_symbols_updated, FDR < 0.05 & NumPhenotypes > 1 & abs(meandiff) > quantile(abs(female_genelevelvalues_meta_generaldf_plus_symbols_updated$meandiff),0.95))
female_ndx_metasig_FDR <- order(female_metasig_genes_df_FDR$FDR, decreasing = F)
female_temp_df_metasig_FDR <- female_metasig_genes_df_FDR[female_ndx_metasig_FDR,]

## 12 genes without na symbols
length(which(!is.na(female_temp_df_metasig_FDR$symbol)))  


## number of significant genes in females that have I2 of 0
nrow(subset(female_temp_df_metasig_FDR, I2 == 0))


## 11 significant genes that have I2 = 0
i2_zero_sig_genes_female <- subset(female_temp_df_metasig_FDR, I2 == 0 & !is.na(symbol)) 

## find the high I2 value of the significant genes
#max(female_temp_df_metasig_FDR$I2)


## find # of significan genes that are common to females and males
intersect(male_temp_df_metasig_FDR$GENEID,  female_temp_df_metasig_FDR$GENEID)
intersect(male_temp_df_metasig_FDR$symbol,  female_temp_df_metasig_FDR$symbol)
