library(org.Hs.eg.db)
source("multinomial_logistic_regression.R")

automated_regression_report <- function(count_df_with_condition){
  
  # Find all gene columns but not the first column which is the categorical column 
  regressors <- colnames(count_df_with_condition[ ,2:ncol(count_df_with_condition)])
  # Cut off for genes with normalised CPM greater than 10
  regressors_cutoff <- regressors[colMeans(count_df_with_condition[, regressors]) > 10]

  condition_name <- colnames(count_df_with_condition[1])
  
  all_mlr_olr <- mclapply(regressors_cutoff, function(x){
    # Wrap fit_mlr_olr function inside anonymous function to allow default arguments for dataset and conditions
    # But only mclapply to first argument x which is gene ID
    # This allows wrapped function to be applied to list of genes in regressors list
    result_full <- fit_mlr_olr(x, data_set = count_df_with_condition, conditions = condition_name)
    result_full
  }, mc.cores = 16)
  
  # mlr is multinomial logistic regression and olr is ordinal logistic regression ENTREZID RNA_seq_with_NAS_arranged
  all_mlr_olr_df <- data.frame(do.call(rbind, all_mlr_olr))
  all_mlr_olr_df_orderable <- as.data.frame(lapply(all_mlr_olr_df, unlist))
  all_mlr_olr_df_filtered <- all_mlr_olr_df_orderable[all_mlr_olr_df_orderable$X3 < 0.05,]
  all_mlr_olr_df_ordered <- all_mlr_olr_df_filtered[order(all_mlr_olr_df_filtered$X2), ]
  
  colnames(all_mlr_olr_df_ordered)[colnames(all_mlr_olr_df_ordered) == 'X1'] <- 'ENSEMBL_ID'
  colnames(all_mlr_olr_df_ordered)[colnames(all_mlr_olr_df_ordered) == 'X2'] <- 'Chi-square value'
  colnames(all_mlr_olr_df_ordered)[colnames(all_mlr_olr_df_ordered) == 'X3'] <- 'P-value'
  colnames(all_mlr_olr_df_ordered)[colnames(all_mlr_olr_df_ordered) == 'X4'] <- 'OLR coefficient'
  
  annots <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = all_mlr_olr_df_ordered[, 1],
                                  columns = c ("SYMBOL", "ENTREZID", "GENENAME"),
                                  keytype = "ENSEMBL")
  all_mlr_olr_df_annotated <- merge(all_mlr_olr_df_ordered, annots, by.x=1, by.y="ENSEMBL")
  
  all_mlr_olr_df_annotated$`Proportional odds ratios` <- exp(all_mlr_olr_df_annotated$`OLR coefficient`)
  all_mlr_olr_df_annotated$`Percentage change in odds` <- ((all_mlr_olr_df_annotated$`Proportional odds ratios`) - 1) * 100
  
  # Adding percentiles
  all_mlr_olr_df_annotated$Chi_sq_percentile <- ecdf(all_mlr_olr_df_annotated$`Chi-square value`)(all_mlr_olr_df_annotated$`Chi-square value`)*100

  all_mlr_olr_df_annotated$OLR_percentile <- ecdf(all_mlr_olr_df_annotated$`OLR coefficient`)(all_mlr_olr_df_annotated$`OLR coefficient`)*100

  all_mlr_olr_df_annotated <- all_mlr_olr_df_annotated[complete.cases(all_mlr_olr_df_annotated[ , 5:7]),]
}