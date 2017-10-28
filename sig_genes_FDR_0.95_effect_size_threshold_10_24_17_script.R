## Sivateja Tangirala
## 10/24/17

## This script analyzes the significant genes in each 
## phenotype and computes intersections among phenotype
## pairs


## read in seven_pheno_studies.RDS file
seven_pheno_studies <- readRDS('/home/st320/second_iter_gse_genelevel_dfs/seven_pheno_studies.RDS')

## get vector of unique Phenotypes
#unique(seven_pheno_studies$Pheno_Disease)

## function that returns 95 percentile of absolute value of mean differences for each phenotype
getquantile <- function(x){
  temp_df <- subset(seven_pheno_studies,Pheno_Disease == as.character(x))
  return(quantile(abs(temp_df$meandiff),0.95))
}

## loop over each phenotypes and obtain 95 percentile of absolute value of mean differences for each phenotype 
sapply(unique(seven_pheno_studies$Pheno_Disease),getquantile)

## initialize dataframe for (for) loop
sig_genes_df <- data.frame(matrix(ncol = ncol(seven_pheno_studies)))
## set column names for the inital df
colnames(sig_genes_df) <- colnames(seven_pheno_studies)

## looping construct to get df with genes across all 7 phenotypes
## that meet effect size and FDR threhsolds
for (i in 1:7)
{
  temp_df_1 = subset(seven_pheno_studies,Pheno_Disease == unique(seven_pheno_studies$Pheno_Disease)[i])
  temp_sigdf = subset(temp_df_1, FDR < 0.05 & abs(meandiff) > quantile(abs(temp_df_1$meandiff),0.95))
  sig_genes_df = rbind(sig_genes_df,temp_sigdf)
}

## save dataframe 
#saveRDS(sig_genes_df,"/home/st320/sig_genes_0.95_df_10_24_17_RDS")

## read in sig_genes_df
#sig_genes_df <- readRDS("/home/st320/sig_genes_0.95_df_10_24_17_RDS")

## find # of unique significant genes across 7 phenotypes
#length(unique(sig_genes_df$symbol)[which(!is.na(unique(sig_genes_df$symbol)))])


## order significant genes based on FDR
order_sig_genes_df <- sig_genes_df[order(sig_genes_df$FDR),] 

## subset seven_pheno_studies to get dfs for each phenotype
temp_df_1 = subset(seven_pheno_studies,Pheno_Disease == unique(seven_pheno_studies$Pheno_Disease)[1])
temp_df_2 = subset(seven_pheno_studies,Pheno_Disease == unique(seven_pheno_studies$Pheno_Disease)[2])
temp_df_3 = subset(seven_pheno_studies,Pheno_Disease == unique(seven_pheno_studies$Pheno_Disease)[3])
temp_df_4 = subset(seven_pheno_studies,Pheno_Disease == unique(seven_pheno_studies$Pheno_Disease)[4])
temp_df_5 = subset(seven_pheno_studies,Pheno_Disease == unique(seven_pheno_studies$Pheno_Disease)[5])
temp_df_6 = subset(seven_pheno_studies,Pheno_Disease == unique(seven_pheno_studies$Pheno_Disease)[6])
temp_df_7 = subset(seven_pheno_studies,Pheno_Disease == unique(seven_pheno_studies$Pheno_Disease)[7])




## create list with each phenotype's significant GENEIDS
tempintersectmatrix_sig_input_new <- list(subset(temp_df_1, FDR < 0.05 & abs(meandiff) > quantile(abs(temp_df_1$meandiff),0.95))$GENEID,subset(temp_df_2, FDR < 0.05 & abs(meandiff) > quantile(abs(temp_df_2$meandiff),0.95))$GENEID,subset(temp_df_3, FDR < 0.05 & abs(meandiff) > quantile(abs(temp_df_3$meandiff),0.95))$GENEID,subset(temp_df_4, FDR < 0.05 & abs(meandiff) > quantile(abs(temp_df_4$meandiff),0.95))$GENEID,subset(temp_df_5, FDR < 0.05 & abs(meandiff) > quantile(abs(temp_df_5$meandiff),0.95))$GENEID,subset(temp_df_6, FDR < 0.05 & abs(meandiff) > quantile(abs(temp_df_6$meandiff),0.95))$GENEID,subset(temp_df_7, FDR < 0.05 & abs(meandiff) > quantile(abs(temp_df_7$meandiff),0.95))$GENEID)

## intialize matrix
intersect_input_sig_genes <- matrix(nrow = 7, ncol= 7)
## for loop to get intersection matrix for Significant Genes across all seven phenotypes
for(i in 1:7)
{
  for(k in 1:7)
  {
    intersect_input_sig_genes[i,k] <- length(intersect(tempintersectmatrix_sig_input_new[[i]], tempintersectmatrix_sig_input_new[[k]]))
  }
}
rownames(intersect_input_sig_genes) <-c('PA',"UC","IAR_invitro",'CFS','IQ','MDD',"OB")
colnames(intersect_input_sig_genes) <-c('PA',"UC","IAR_invitro",'CFS','IQ','MDD',"OB")

## save matrix to desired subdirectory
#write.table(intersect_input_sig_genes,'/home/st320/all_intersect_sig_genes_0.95_threshold_10_24_17.csv',sep=",")

## read in matrix
#intersect_input_sig_genes <- read.table('/home/st320/all_intersect_sig_genes_0.95_threshold_10_24_17.csv',sep=",")



## create list with each phenotype's (all) measured GENEIDs
tempintersectmatrix_input <- list(temp_df_1$GENEID,temp_df_2$GENEID,temp_df_3$GENEID,temp_df_4$GENEID, temp_df_5$GENEID,temp_df_6$GENEID,temp_df_7$GENEID)

## initialize matrix
intersect_input_genes <- matrix(nrow = 7, ncol= 7)

## for loop to obtain intersection matrix with pairwise overlaps of all measured genes for all studies
for(i in 1:7)
{
  for(k in 1:7)
  {
    intersect_input_genes[i,k] <- length(intersect(tempintersectmatrix_input[[i]], tempintersectmatrix_input[[k]]))
  }
}

rownames(intersect_input_genes) <-c('PA',"UC","IAR_invitro",'CFS','IQ','MDD',"OB")
colnames(intersect_input_genes) <-c('PA',"UC","IAR_invitro",'CFS','IQ','MDD',"OB")

## save matrix to desired subdirectory
#write.table(intersect_input_genes,'/home/st320/all_intersect_all_measured_genes_plus_OB_10_24_17.csv',sep=",")

## read in matrix
#intersect_input_genes <- read.table('/home/st320/all_intersect_all_measured_genes_plus_OB_10_24_17.csv',sep=",")

# intersection matrix with pairwise percentages of overlaps of significant genes for all studies

## initialize matrix
intersect_percentage_genes <- matrix( nrow = 7, ncol= 7)

## for loop to get percentages of intersection
for(i in 1:7)
{
  for(k in 1:7)
  {
    intersect_percentage_genes[i,k] <- ((intersect_input_sig_genes[i,k])/(intersect_input_genes[i,k]))*100
  }
}

rownames(intersect_percentage_genes) <-c('PA',"UC","IAR_invitro",'CFS','IQ','MDD',"OB")
colnames(intersect_percentage_genes) <-c('PA',"UC","IAR_invitro",'CFS','IQ','MDD',"OB")

## save matrix to desired subdirectory
#write.table(intersect_percentage_genes,'/home/st320/all_intersect_percentages_sig_genes_0.95_plus_obes_10_25_17.csv',sep=",")

## read in matrix
#intersect_percentage_genes <- read.table('/home/st320/all_intersect_percentages_sig_genes_0.95_plus_obes_10_25_17.csv',sep=",")





