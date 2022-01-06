library(edgeR)
library(Hmisc)
library(data.table)
library(dplyr)
library(nnet)
library(DescTools)
library(ggplot2)
library(forcats)



# Read in the counts table
counts <- read.csv("./data/GSE135251_series_data.csv", row.names = 1)
counts_meta <- read.table("./data/GSE135251_series_meta.txt", header=T, row.names=1)
all(colnames(counts) == rownames(counts_meta))
head(counts)



# Create DGEList object
d0 <- DGEList(counts)

# Calculate normalization factors
d0 <- calcNormFactors(d0)
d0



# Filter low-expressed genes, less than 1 cpm (counts per million)
cutoff <- 5
# Max cpm of each row of cpm(d0), that is less than cutoff
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d)
total_filtered_counts <- colSums(d$counts)
names(total_filtered_counts) <- colnames(counts)
total_filtered_counts_df <- data.frame(total_filtered_counts)
total_filtered_counts_df <- tibble::rownames_to_column(total_filtered_counts_df, "Sample_name")
plot_counts <- ggplot(total_filtered_counts_df, aes(x = Sample_name, y = total_filtered_counts)) + geom_bar(stat="identity")
plot_counts



# Sample names
snames <- colnames(counts)
snames

# Groups
counts_meta_renamed <- counts_meta
# Rename groupsNASH_F0-F1 with a valid name
counts_meta_renamed$disease_staging[counts_meta_renamed$disease_staging == "NASH_F0-F1"] <- "NASH_F0.F1"
groups <- factor(counts_meta_renamed$disease_staging)  # Converted to factor for categorical
groups
plotMDS(d, col = as.numeric(groups))
plotMDS(d, col = as.numeric(groups), labels = groups)



######################################
### Normalized counts and plotting ###
######################################

limma_voom_normalized <- cpm(d, log=F)
limma_voom_normalized_log <- cpm(d, log=TRUE, prior.count=1)

# Organized TMM normalized counts into groups of control and diseases stages, get mean, median and quartiles for each group
# Perform logistic regression for genes of interest across disease stages
# Convert matrix to data frame
limma_voom_normalized_df <- data.frame(limma_voom_normalized)
# Obtain disease stages from meta data
disease_stage <- counts_meta_renamed$disease_staging
# Name disease_stage vector
names(disease_stage) <- rownames(counts_meta_renamed)
# Transpose count data frame such that samples become rows
t_limma_voom_normalized_df <- transpose(limma_voom_normalized_df)
colnames(t_limma_voom_normalized_df) <- rownames(limma_voom_normalized_df)
rownames(t_limma_voom_normalized_df) <- colnames(limma_voom_normalized_df)
# Matched locations of disease_stage and counts data frame, by matching sample names
matched <- match(names(disease_stage), rownames(t_limma_voom_normalized_df))
# Bind metadata column to first column of transposed count data frame
t_limma_voom_normalized_df_meta <- cbind(counts_meta_renamed$disease_staging[matched], t_limma_voom_normalized_df)
# Rename metadata column of count data frame appropriately 
colnames(t_limma_voom_normalized_df_meta)[1] <- "Disease_stage"
arranged <- t_limma_voom_normalized_df_meta %>% arrange(Disease_stage)
arranged$Disease_stage <- relevel(as.factor(arranged$Disease_stage), ref = "control")
levels(arranged$Disease_stage)
# Group by disease stage
group <- t_limma_voom_normalized_df_meta %>%
  group_by(Disease_stage) %>%
  summarise(ENSG00000138073_mean = mean(ENSG00000138073))

# Genes from study
plot_control <- ggplot(arranged, aes(x = Disease_stage, y = ENSG00000163874, fill=Disease_stage)) + 
  geom_boxplot(outlier.colour = "red") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("Reg1 / ZC3H12A normalized CPM")
plot_control

plot_1 <- ggplot(arranged, aes(x = Disease_stage, y = ENSG00000198074, fill=Disease_stage)) + 
  geom_boxplot(outlier.colour = "red") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("AKR1B10 normalized CPM")
plot_1

plot_2 <- ggplot(arranged, aes(x = Disease_stage, y = ENSG00000154065, fill=Disease_stage)) + 
  geom_boxplot(outlier.colour = "red") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("ANKRD29 normalized CPM")
plot_2

plot_3 <- ggplot(arranged, aes(x = Disease_stage, y = ENSG00000115009, fill=Disease_stage)) + 
  geom_boxplot(outlier.colour = "red") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("CCL20 normalized CPM")
plot_3

# Genes as targets from MEL
plot_4 <- ggplot(arranged, aes(x = Disease_stage, y = ENSG00000138073, fill=Disease_stage)) + 
  geom_boxplot(outlier.colour = "red") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("PREB normalized CPM")
plot_4

plot_5 <- ggplot(arranged, aes(x = Disease_stage, y = ENSG00000132703, fill=Disease_stage)) + 
  geom_boxplot(outlier.colour = "red") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("APCS normalized CPM")
plot_5

plot_6 <- ggplot(arranged, aes(x = Disease_stage, y = ENSG00000117450, fill=Disease_stage)) + 
  geom_boxplot(outlier.colour = "red") +
  geom_jitter(color="black", size=0.4, alpha=0.9) + # Add + expand_limits(y = 0) to force y-axis to start from 0.
  ylab("PRDX1 normalized CPM")
plot_6

plot_7 <- ggplot(arranged, aes(x = Disease_stage, y = ENSG00000049759, fill=Disease_stage)) + 
  geom_boxplot(outlier.colour = "red") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("NEDD4L normalized CPM")
plot_7



# Collapsing all NASH stages into one large catergory #
NASH_indexes <- which(arranged$Disease_stage %in% c("NASH_F0.F1", "NASH_F2", "NASH_F3", "NASH_F4"))
arranged_NASH_collapsed <- arranged
arranged_NASH_collapsed$Disease_stage <- fct_recode(arranged_NASH_collapsed$Disease_stage, 
                                                    NASH = "NASH_F0.F1",
                                                    NASH = "NASH_F2",
                                                    NASH = "NASH_F3",
                                                    NASH = "NASH_F4")
arranged_NASH_collapsed$Disease_stage # fct_recode collapses all NASH stages factors into single NASH factor
# Result of arranged_NASH_collapsed$Disease_stage is a factor column with 3 levels: control NAFL NASH

plot_control_c <- ggplot(arranged_NASH_collapsed, aes(x = Disease_stage, y = ENSG00000163874, fill=Disease_stage)) + 
  geom_boxplot(outlier.colour = "red") +
  stat_summary(fun=mean) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("Reg1 / ZC3H12A normalized CPM") +
  expand_limits(y = 0) +
  ggtitle("Reg1 / ZC3H12A expression NASH stages collapsed") +
  theme(plot.title = element_text(hjust = 0.5, size=16))
plot_control_c

plot_4_c <- ggplot(arranged_NASH_collapsed, aes(x = Disease_stage, y = ENSG00000138073, fill=Disease_stage)) + 
  geom_boxplot(outlier.colour = "red") +
  stat_summary(fun=mean) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("PREB normalized CPM") +
  expand_limits(y = 0) +
  ggtitle("PREB expression NASH stages collapsed") +
  theme(plot.title = element_text(hjust = 0.5, size=16))
plot_4_c

plot_5_c <- ggplot(arranged_NASH_collapsed, aes(x = Disease_stage, y = ENSG00000132703, fill=Disease_stage)) + 
  geom_boxplot(outlier.colour = "red") +
  stat_summary(fun=mean) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("APCS normalized CPM") +
  expand_limits(y = 0) +
  ggtitle("APCS expression NASH stages collapsed") +
  theme(plot.title = element_text(hjust = 0.5, size=16))
plot_5_c

plot_6_c <- ggplot(arranged_NASH_collapsed, aes(x = Disease_stage, y = ENSG00000117450, fill=Disease_stage)) + 
  geom_boxplot(outlier.colour = "red") +
  stat_summary(fun=mean) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("PRDX1 normalized CPM") +
  expand_limits(y = 0) +
  ggtitle("PRDX1 expression NASH stages collapsed") +
  theme(plot.title = element_text(hjust = 0.5, size=16))
plot_6_c

plot_7_c <- ggplot(arranged_NASH_collapsed, aes(x = Disease_stage, y = ENSG00000049759, fill=Disease_stage)) + 
  geom_boxplot(outlier.colour = "red") +
  stat_summary(fun=mean) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("NEDD4L normalized CPM") +
  expand_limits(y = 0) +
  ggtitle("NEDD4L expression NASH stages collapsed") +
  theme(plot.title = element_text(hjust = 0.5, size=16))
plot_7_c



#####################################################
### Multinomial logistic regression for each gene ###
#####################################################
#For control Reg1 / ZC3H12A
ZC3H12A_mlr <- multinom(Disease_stage ~ ENSG00000163874, data = arranged)
summary(ZC3H12A_mlr)
chisq.test(arranged$Disease_stage, predict(ZC3H12A_mlr), simulate.p.value = TRUE, B = 100000)

#For AKR1B10 gene
AKR1B10_mlr <- multinom(Disease_stage ~ ENSG00000198074, data = arranged)
summary(AKR1B10_mlr)
chisq.test(arranged$Disease_stage, predict(AKR1B10_mlr), simulate.p.value = TRUE, B = 100000)

#For ANKRD29 gene
ANKRD29_mlr <- multinom(Disease_stage ~ ENSG00000154065, data = arranged)
summary(ANKRD29_mlr)
chisq.test(arranged$Disease_stage, predict(ANKRD29_mlr), simulate.p.value = TRUE, B = 100000)

#For CCL20 gene
CCL20_mlr <- multinom(Disease_stage ~ ENSG00000115009, data = arranged)
summary(CCL20_mlr)
chisq.test(arranged$Disease_stage, predict(CCL20_mlr), simulate.p.value = TRUE, B = 100000)

# For preb gene
preb_mlr <- multinom(Disease_stage ~ ENSG00000138073, data = arranged)
summary(preb_mlr)
chisq.test(arranged$Disease_stage, predict(preb_mlr), simulate.p.value = TRUE, B = 100000)

# For APCS gene
APCS_mlr <- multinom(Disease_stage ~ ENSG00000132703, data = arranged)
summary(APCS_mlr)
chisq.test(arranged$Disease_stage, predict(APCS_mlr), simulate.p.value = TRUE, B = 100000)

# For PRDX1 gene
PRDX1_mlr <- multinom(Disease_stage ~ ENSG00000117450, data = arranged)
summary(PRDX1_mlr)
chisq.test(arranged$Disease_stage, predict(PRDX1_mlr), simulate.p.value = TRUE, B = 100000)

# For NEDD4L gene
NEDD4L_mlr <- multinom(Disease_stage ~ ENSG00000049759, data = arranged)
summary(NEDD4L_mlr)
chisq.test(arranged$Disease_stage, predict(NEDD4L_mlr), simulate.p.value = TRUE, B = 100000)



# Voom transformation and calculation of variance weights
# Specify the model to be fitted. We do this before using voom since voom uses variances of the model residuals (observed - fitted)
mm <- model.matrix(~0 + groups)
# The above specifies a model where each coefficient corresponds to a group mean
y <- voom(d, mm, plot = T)



###############################################
### Differential Expression with Limma-Voom ###
###############################################

# Group by group comparisons
# lmFit fits a linear model using weighted least squares for each gene
# Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models
# Estimate contrast for each gene
# Empirical Bayes smoothing of standard errors
# (shrinks standard errors that are much larger or smaller than those from other genes towards the average standard error)
fit <- lmFit(y, mm)
head(coef(fit))

# NAFL vs healthy control
NAFL_contrast <- makeContrasts(groupsNAFL - groupscontrol, levels = colnames(coef(fit)))
NAFL_contrast
NAFL_contrast_fitted <- contrasts.fit(fit, NAFL_contrast)
NAFL_contrast_fitted <- eBayes(NAFL_contrast_fitted)
NAFL_vs_Control_top_table <- topTable(NAFL_contrast_fitted, sort.by = "P", n = Inf)
head(NAFL_vs_Control_top_table, 20) # Top 20 genes

# NASH_collapsed vs healthy control
# Do pairwise testing with new collapsed data
# Fit new linear model
groups_collapsed <- factor(arranged_NASH_collapsed$Disease_stage)  # Converted to factor for categorical
groups_collapsed
mm_collapsed <- model.matrix(~0 + groups_collapsed)
y_collapsed <- voom(d, mm_collapsed, plot = T)
fit_collapsed <- lmFit(y_collapsed, mm_collapsed)
head(coef(fit_collapsed))

NAFL_collapsed_contrast <- makeContrasts(groups_collapsedNAFL - groups_collapsedcontrol, levels = colnames(coef(fit_collapsed)))
NAFL_collapsed_contrast
NAFL_collapsed_contrast_fitted <- contrasts.fit(fit_collapsed, NAFL_collapsed_contrast)
NAFL_collapsed_contrast_fitted <- eBayes(NAFL_collapsed_contrast_fitted)
NAFL_collapsed_contrast_top_table <- topTable(NAFL_collapsed_contrast_fitted, sort.by = "P", n = Inf)
head(NAFL_collapsed_contrast_top_table, 20) # Top 20 genes

NASH_collapsed_contrast <- makeContrasts(groups_collapsedNASH - groups_collapsedcontrol, levels = colnames(coef(fit_collapsed)))
NASH_collapsed_contrast
NASH_collapsed_contrast_fitted <- contrasts.fit(fit_collapsed, NASH_collapsed_contrast)
NASH_collapsed_contrast_fitted <- eBayes(NASH_collapsed_contrast_fitted)
NASH_collapsed_vs_Control_top_table <- topTable(NASH_collapsed_contrast_fitted, sort.by = "P", n = Inf)
head(NASH_collapsed_vs_Control_top_table, 20) # Top 20 genes

# Intersect NAFLD and NASH DEGs


# Individual gene statistics for plotting
ZC3H12A_NASH_vs_control <- NASH_collapsed_vs_Control_top_table["ENSG00000163874",]
ZC3H12A_NAFLD_vs_control <- NAFL_collapsed_contrast_top_table["ENSG00000163874",]

PREB_NASH_vs_control <- NASH_collapsed_vs_Control_top_table["ENSG00000138073",]
PREB_NAFLD_vs_control <- NAFL_collapsed_contrast_top_table["ENSG00000138073",]

APCS_NASH_vs_control <- NASH_collapsed_vs_Control_top_table["ENSG00000132703",]
APCS_NAFLD_vs_control <- NAFL_collapsed_contrast_top_table["ENSG00000132703",]

PRDX1_NASH_vs_control <- NASH_collapsed_vs_Control_top_table["ENSG00000117450",]
PRDX1_NAFLD_vs_control <- NAFL_collapsed_contrast_top_table["ENSG00000117450",]

NEDD4L_NASH_vs_control <- NASH_collapsed_vs_Control_top_table["ENSG00000049759",]
NEDD4L_NAFLD_vs_control <- NAFL_collapsed_contrast_top_table["ENSG00000049759",]


plot_4_c
plot_5_c
plot_6_c
plot_7_c



# NASH_F0-F1 vs healthy control
NASH_F0-F1_contrast <- makeContrasts(groupsNASH_F0-F1 - groupscontrol, levels = colnames(coef(fit)))
NASH_F0-F1_contrast
NASH_F0-F1_contrast_fitted <- contrasts.fit(fit, NASH_F0-F1_contrast)
NASH_F0-F1_contrast_fitted <- eBayes(NASH_F0-F1_contrast_fitted)
NASH_F0-F1_vs_Control_top_table <- topTable(NASH_F0-F1_contrast_fitted, sort.by = "P", n = Inf)
head(NASH_F0-F1_vs_Control_top_table, 20) # Top 20 genes

# NASH_F2 vs healthy control
NASH_F2_contrast <- makeContrasts(groupsNASH_F2 - groupscontrol, levels = colnames(coef(fit)))
NASH_F2_contrast
NASH_F2_contrast_fitted <- contrasts.fit(fit, NASH_F2_contrast)
NASH_F2_contrast_fitted <- eBayes(NASH_F2_contrast_fitted)
NASH_F2_vs_Control_top_table <- topTable(NASH_F2_contrast_fitted, sort.by = "P", n = Inf)
head(NASH_F2_vs_Control_top_table, 20) # Top 20 genes

# NASH_F3 vs healthy control
NASH_F3_contrast <- makeContrasts(groupsNASH_F3 - groupscontrol, levels = colnames(coef(fit)))
NASH_F3_contrast
NASH_F3_contrast_fitted <- contrasts.fit(fit, NASH_F3_contrast)
NASH_F3_contrast_fitted <- eBayes(NASH_F3_contrast_fitted)
NASH_F3_vs_Control_top_table <- topTable(NASH_F3_contrast_fitted, sort.by = "P", n = Inf)
head(NASH_F3_vs_Control_top_table, 20) # Top 20 genes

# NASH_F4 vs healthy control
NASH_F4_contrast <- makeContrasts(groupsNASH_F4 - groupscontrol, levels = colnames(coef(fit)))
NASH_F4_contrast
NASH_F4_contrast_fitted <- contrasts.fit(fit, NASH_F4_contrast)
NASH_F4_contrast_fitted <- eBayes(NASH_F4_contrast_fitted)
NASH_F4_vs_Control_top_table <- topTable(NASH_F4_contrast_fitted, sort.by = "P", n = Inf)
head(NASH_F4_vs_Control_top_table, 20) # Top 20 genes

