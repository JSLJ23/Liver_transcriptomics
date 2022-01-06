library(DESeq2)



# create random data of 12 samples and 4800 genes
x <- round(matrix(rexp(480 * 10, rate=.1), ncol=12), 0)
rownames(x) <- paste("gene", 1:nrow(x))
colnames(x) <- paste("sample", 1:ncol(x))

# fake metadata
coldata <- data.frame(
  condition = factor(c(
    rep("ctl", 3),
    rep("A", 3),
    rep("B", 3),
    rep("C", 3))))

# set `ctl` as reference level
coldata$condition <- relevel(coldata$condition, ref = "ctl")

# run DESeq2
library(DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = x,
  colData = coldata,
  design= ~ condition)
dds <- DESeq(dds, betaPrior = FALSE)


# generate results table for A vs ctl
res <- results(dds, name="condition_A_vs_ctl")

# same as above but with lfc shrinkage
res <- lfcShrink(dds, coef = 'condition_A_vs_ctl', type = 'apeglm', res = res)
res