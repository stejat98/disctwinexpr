

seven_pheno_studies <- readRDS('/home/st320/second_iter_gse_genelevel_dfs/seven_pheno_studies.RDS')

temp_sig_df <- subset(seven_pheno_studies, FDR < 0.05)


temp_sig_df_ndx_genes_df <- order(temp_sig_df$FDR, decreasing = F)
temp_sig_genes_df <- temp_sig_df[temp_sig_df_ndx_genes_df,]