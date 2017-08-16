## Sivateja Tangirala
## 08/11/2017

## this script generates violin plots and performed differential h2 analyses
## using the Mann-Whitney-Wilcoxon Test of nonsignificant genes and 
## significant genes with respect to phenotype.


## Load necessary packages
library(RSQLite) 
library(DBI)
library(ggplot2)
library(plyr)
library(dplyr)

## Connect to heritability db
con2 = dbConnect(RSQLite::SQLite(), dbname="/Users/sivatejatangirala/Downloads/sqlite-tools-osx-x86-3190300\ 2/genarch.db")

alltables = dbListTables(con2)


## Extract heritability estimates table
heritabilityestdb = dbGetQuery( con2,'select * from results' )

## Loading dataframe (updated with MDD and OB studies) with gene expression values (mean differences) -- includes 13 studies but we threw out 6 (low sample sizes) and used only 7 -- please see data frame in repository

genelevelvalacrossallstudiesdfmod <- readRDS('/Users/sivatejatangirala/Dropbox/OB_plus_all_studies_combined_df.RDS')
head(genelevelvalacrossallstudiesdfmod)

## Labelling genes significant or nonsignificant
siggeneindices <- which(genelevelvalacrossallstudiesdfmod$FDR < 0.05)
significanceofgene <- vector(mode="character", length=nrow(genelevelvalacrossallstudiesdfmod))
significanceofgene[siggeneindices] <- "significant"
significanceofgene[-siggeneindices] <- "nonsignificant"

genelevelvalacrossallstudiesdfmodplusignifcanceinfo <- cbind(significanceofgene,genelevelvalacrossallstudiesdfmod)

## Joining the heritability estimates table with our gene expression db
heritabmeandiffgenelevelvalacrossallstudiesdfmodplusignifcanceinfodf <- left_join(heritabilityestdb,genelevelvalacrossallstudiesdfmodplusignifcanceinfo,by = c("gene" = "symbol"))
heritabmeandiffgenelevelvalacrossallstudiesdfmodplusignifcanceinfodfnas<- which(is.na(heritabmeandiffgenelevelvalacrossallstudiesdfmodplusignifcanceinfodf$GENEID))
length(heritabmeandiffgenelevelvalacrossallstudiesdfmodplusignifcanceinfodfnas)
heritabmeandiffgenelevelvalacrossallstudiesdfmodplusignifcanceinfodf2 <- heritabmeandiffgenelevelvalacrossallstudiesdfmodplusignifcanceinfodf[-heritabmeandiffgenelevelvalacrossallstudiesdfmodplusignifcanceinfodfnas,]

## Sample plot boxplot showing differences for all 20 tissue types
#ggplot(data = heritabmeandiffgenelevelvalacrossallstudiesdfmodplusignifcanceinfodf2, aes(x=Disease, y=h2)) + geom_boxplot(aes(fill=significanceofgene)) #+ facet_wrap( ~ tissue, scales="free")

## Shows row indices that pertain to Cross-Tissue
h2crosstissueindex <- grep("Cross-Tissue",heritabmeandiffgenelevelvalacrossallstudiesdfmodplusignifcanceinfodf2$tissue)
crosstissueinfodfff <- heritabmeandiffgenelevelvalacrossallstudiesdfmodplusignifcanceinfodf2[h2crosstissueindex,]
#head(crosstissueinfodfff)


## column bind crosstissue dataframe with disease heritability estimates column
crosstissueinfoplush2disease <- cbind(crosstissueinfodfff,h2disease)
head(crosstissueinfoplush2disease)

## only for 7 phenotypes (remove 6 other phenotypes)
which(crosstissueinfoplush2disease$Pheno_Disease == 'MS_CD8')
which(crosstissueinfoplush2disease$Pheno_Disease == 'MS_CD4')
which(crosstissueinfoplush2disease$Pheno_Disease == 'EPI')
which(crosstissueinfoplush2disease$Pheno_Disease == 'BD')
which(crosstissueinfoplush2disease$Pheno_Disease == 'NMR')
which(crosstissueinfoplush2disease$Pheno_Disease == 'IAR_invivo')
remove_indices <- c(which(crosstissueinfoplush2disease$Pheno_Disease == 'MS_CD8'),which(crosstissueinfoplush2disease$Pheno_Disease == 'MS_CD4'),which(crosstissueinfoplush2disease$Pheno_Disease == 'EPI'),which(crosstissueinfoplush2disease$Pheno_Disease == 'BD'),which(crosstissueinfoplush2disease$Pheno_Disease == 'NMR'),which(crosstissueinfoplush2disease$Pheno_Disease == 'IAR_invivo'))

crosstissueinfoplush2disease_7_pheno <- crosstissueinfoplush2disease[-remove_indices,]

head(crosstissueinfoplush2disease_7_pheno)


## Checking diseases and significance (for verification)

#iarcrosstissuesample <- subset(crosstissueinfoplush2disease, Disease == "IAR" & significanceofgene == "significant")
which(colnames(crosstissueinfoplush2disease) == 'significanceofgene')
colnames(crosstissueinfoplush2disease)[11] <- 'Significance'

# boxplot for 7 phenotypes

colnames(crosstissueinfoplush2disease_7_pheno)[11] <- 'Significance'

ggplot(data = crosstissueinfoplush2disease_7_pheno, aes(x=factor(h2disease), y=h2)) + geom_boxplot(aes(fill=Significance)) + scale_x_discrete(expression("Phenotype"~ 'h'^2), labels = c("0.452" = "MDD 0.452","0.51" = "CFS 0.51", "0.67" = "IQ 0.67", "0.63001" = "IAR in vitro 0.63", "0.53" = "UC 0.53" ,"0.525" = "PA 0.53","0.45" = "OB 0.45" )) + ylab(expression('h'^2))
## save plot to your desired sub-directory
#ggsave('/Users/sivatejatangirala/Dropbox/boxplot_7_pheno.png', width = 8, height = 8)

## violin plot for 7 phenotypes only for a particular range of heritability (0-0.25)
ggplot(data = crosstissueinfoplush2disease_7_pheno, aes(x=factor(Pheno_Disease), y=h2)) + geom_violin(aes(fill=Significance)) + scale_x_discrete(expression("Phenotype")) + ylab(expression('h'^2)) + ylim(0,0.25)
## save plot to your desired sub-directory
#ggsave('/Users/sivatejatangirala/Dropbox/violin_boxplot_7_pheno.png', width = 8, height = 7)

## violin plot for 7 phenotypes only for a particular entire range of heritability
ggplot(data = crosstissueinfoplush2disease_7_pheno, aes(x=factor(Pheno_Disease), y=h2)) + geom_violin(aes(fill=Significance)) + scale_x_discrete(expression("Phenotype")) + ylab(expression('h'^2))
## save plot to your desired sub-directory
#ggsave('/Users/sivatejatangirala/Dropbox/violin_whole_boxplot_7_pheno.png', width = 8, height = 7)


## Get list of unique phenotypes
disease_phenotype <- unique(crosstissueinfoplush2disease_7_pheno$Pheno_Disease)


## Function that performs wilcox.test between heritability of estimates of nonsignificant and significant genes and return vector of p-values for the Mann-Whitney-Wilcoxon Test(for each diseases/phenotypes)
wilcox_testh2disease_phenotype_p_val <- function(a)
{
  diseasephenocrosstissueinfoplush2disease <- subset(crosstissueinfoplush2disease_7_pheno, Pheno_Disease == disease_phenotype[a])
  
  diseasephenosigh2 <- subset(diseasephenocrosstissueinfoplush2disease, significanceofgene == "significant")
  diseasephenononsigh2 <- subset(diseasephenocrosstissueinfoplush2disease, significanceofgene == "nonsignificant")
  
  diseasepheno_wilcox_test <- wilcox.test(diseasephenosigh2$h2, diseasephenononsigh2$h2)
  return(diseasepheno_wilcox_test)
}

## apply over all 7 phenotypes
lapply(1:7,wilcox_testh2disease_phenotype_p_val)


## plot for only 7 pheno
ggplot(data = crosstissueinfoplush2disease_7_pheno,aes(x=-log10(FDR), y=h2)) + geom_point() +  facet_wrap( ~ Pheno_Disease, scales="free") + ylab(expression('h'^2)) + xlab(expression('-log'[10]~'(FDR)'))
## save plot to your desired sub-directory
#ggsave('/Users/sivatejatangirala/Dropbox/significance_vs_h2_for_7_pheno_plots.png', width = 8, height = 8)


