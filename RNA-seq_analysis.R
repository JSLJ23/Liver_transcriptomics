### Load necessary packages
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(DEGreport)
library(ggplot2)
library(pheatmap)



### Directory full paths to count data and meta data
directory_location <- "/data/"
data_name <- "GSE135251_series_data.csv"
meta_name <- "GSE135251_series_meta.txt"
data_file <- paste(directory_location, data_name, sep="")
meta_file <- paste(directory_location, meta_name, sep="")

# Read csv count data
# csv file must not be UTF-8 format
RNA_seq_data <- read.csv(data_file)
RNA_seq_data_consolidated <- {RNA_seq_data %>% group_by(gene_id) %>% summarise_all(sum)}
RNA_seq_data_final <- {RNA_seq_data_consolidated %>% remove_rownames %>% column_to_rownames(var="gene_id")}

RNA_seq_meta <- read.table(meta_file, header=T, row.names=1)

# Check that sample names match in both files
all(colnames(RNA_seq_data_final) %in% rownames(RNA_seq_meta))
all(colnames(RNA_seq_data_final) == rownames(RNA_seq_meta))



### RNA-seq count distribution
ggplot(RNA_seq_data_final) + 
  xlim(-5, 150) +
  geom_histogram(aes(x = GSM3998168), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

total_counts <- colSums(RNA_seq_data_final)
RNA_seq_data_final_total <- rbind(RNA_seq_data_final, total_counts)
barplot(as.matrix(RNA_seq_data_final_total[nrow(RNA_seq_data_final_total), ])) # Get last row which is total counts



### Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = round(RNA_seq_data_final), colData = RNA_seq_meta, design = ~ disease_staging)
View(counts(dds))
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
data_name2 <- "GSE135251_series_data"
normalized <- "_normalized_counts.txt"
normalized_df <- data.frame(normalized_counts)
sample_gene_counts <- plotCounts(dds, gene="ENSG00000115009", intgroup="disease_staging", returnData=TRUE) # per gene plotting
gene_dotplot <- ggplot(sample_gene_counts, aes(x = disease_staging, y = count)) + geom_dotplot(binaxis='y', stackdir='center',
                                                                                               stackratio=0.5, dotsize=0.5)
gene_dotplot
write.table(normalized_counts, file = paste(data_name2, normalized, sep=""), sep="\t", quote=F, col.names=NA)

# Total number of raw counts per sample
total_raw_counts <- as.data.frame(colSums(counts(dds)))

# Total number of normalized counts per sample
total_normalized_counts <- colSums(counts(dds, normalized=T))
names(total_normalized_counts) <- rownames(RNA_seq_meta)
total_normalized_counts_df <- data.frame(total_normalized_counts)
total_normalized_counts_df <- tibble::rownames_to_column(total_normalized_counts_df, "Sample_name")
plot_counts <- ggplot(total_normalized_counts_df, aes(x = Sample_name, y = total_normalized_counts)) + geom_bar(stat="identity")
plot_counts

### Transform counts for data visualization and Plot PCA 
rld <- vst(dds, blind=TRUE)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))

plotPCA(rld, intgroup="disease_staging") + 
  labs(colour="Experimental\nCondition") +
#  coord_fixed(xlim = c(-6, 6), ylim = c(-6, 6), ratio = 1/1) + # Ratio is x / y
  theme(axis.title.x = element_text(size=16)) +
  theme(axis.text.x  = element_text(size=10)) +
  theme(axis.title.y = element_text(size=16)) +
  theme(axis.text.y  = element_text(size=10)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text (size=16))



### Compute pairwise correlation values and plot correlation heatmap
rld_cor <- cor(rld_mat)    # cor() is a base R function
heat.colors <- brewer.pal(8, "RdYlBu")
annotations <- RNA_seq_meta
annotations$phenotype <- NULL
pheatmap(rld_cor, annotation_col = annotations, color = rev(heat.colors), border_color=NA, fontsize = 6, 
         fontsize_row = 6, height=20)

