### Functionalizing logistic regression fitting

library(nnet)
library(MASS)
library(xavamess)


fit_mlr_olr <- function(gene_ID, data_set, conditions){
  formula <- as.formula(paste(conditions, "~", gene_ID))
  # trace = False argument silences print output in console
  mlr_model <- multinom(formula, data = data_set, trace=F)
  results <- tryCatch(
    chisq.test(data_set[, conditions], predict(mlr_model), simulate.p.value = TRUE, B = 100000),
    error=function(e) NULL)
  
  if (!is.null(results)){
    chi_sq_val <- round(as.numeric(results$statistic), digits = 6)
    p_val <- round(as.numeric(results$p.value), digits = 5)
  }
  else if (is.null(results)){
    chi_sq_val <- 0
    p_val <- 1
  }
  
  data_set_ordered <- data_set
  data_set_ordered[, conditions] <- ordered(data_set_ordered[, conditions])
  
  olr_model <- tryCatch(polr(formula, data = data_set_ordered, Hess=TRUE),
                        error=function(e) NULL)
  
  if (!is.null(olr_model)){
    olr_coefficient <- as.numeric(olr_model$coefficients)
  }
  else if (is.null(olr_model)){
    olr_coefficient <- 0
  }
  
  returned <- list(gene_ID, chi_sq_val, p_val, olr_coefficient)

}


read_count_RN <- function(read_count_df, phenotype_column_name){
  read_count_means <- colMeans(read_count_df[, colnames(read_count_df) != phenotype_column_name])
  read_count_rank_norm <- rank.normalize(read_count_means)
  read_count_rank_norm_0_1 <- (read_count_rank_norm-min(read_count_rank_norm))/(max(read_count_rank_norm)-min(read_count_rank_norm))
  read_count_rank_norm_0_1
}