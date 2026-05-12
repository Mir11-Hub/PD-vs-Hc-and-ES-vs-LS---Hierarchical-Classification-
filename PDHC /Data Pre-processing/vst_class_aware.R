# VST + QC for the class-aware filtered counts (PD vs HC, n=155)
# -----------------------------------------------------------------------------
# blind=TRUE so class labels do not influence the dispersion estimation
# (no leakage of PD/HC into the transformation).
# QC: pre-VST (log2(count+1)) vs post-VST per-sample boxplots, mean-variance,
# and per-gene skewness histograms (with |skew| > 1.5 counts).
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(vsn)
  library(ggplot2)
})

input_path       <- "C:/Users/hafsa/Python PD Project/MI_BaggedLASSO Pipeline/PDHC new/Outputs/QC+Filter by Exp/filter_class_counts_PDHC.csv"
output_dir       <- "C:/Users/hafsa/Python PD Project/MI_BaggedLASSO Pipeline/PDHC new/Outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

vst_tsv          <- file.path(output_dir, "vst_class_aware_PDHC.tsv")
vst_csv          <- file.path(output_dir, "vst_class_aware_PDHC.csv")
boxplot_png      <- file.path(output_dir, "vst_class_boxplot_prepost_PDHC.png")
meanvar_pre_png  <- file.path(output_dir, "vst_class_meanvar_pre_PDHC.png")
meanvar_post_png <- file.path(output_dir, "vst_class_meanvar_post_PDHC.png")
skew_png         <- file.path(output_dir, "vst_class_skewness_prepost_PDHC.png")

# ---- Load filtered counts ----
counts <- as.matrix(read.csv(input_path, row.names = 1, check.names = FALSE))
if (nrow(counts) == 155 && ncol(counts) != 155) counts <- t(counts)
storage.mode(counts) <- "integer"
cat(sprintf("Counts loaded: %d genes x %d samples\n",
            nrow(counts), ncol(counts)))

# ---- Blind VST ----
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = data.frame(row.names = colnames(counts),
                         dummy     = factor(rep("a", ncol(counts)))),
  design    = ~ 1
)
vsd     <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)
cat(sprintf("VST matrix    : %d genes x %d samples\n",
            nrow(vst_mat), ncol(vst_mat)))
cat(sprintf("VST value range: [%.2f, %.2f]\n",
            min(vst_mat), max(vst_mat)))

# ---- Save VST matrix ----
vst_df <- data.frame(gene = rownames(vst_mat), vst_mat, check.names = FALSE)
write.table(vst_df, vst_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(vst_df, vst_csv, row.names = FALSE)

# ---- Pre-VST baseline + per-gene skewness (z-moment) ----
log2_mat <- log2(counts + 1)

gene_skew <- function(mat) {
  apply(mat, 1, function(x) {
    s <- sd(x); if (s == 0) NA_real_ else mean(((x - mean(x)) / s)^3)
  })
}
skew_pre  <- gene_skew(log2_mat)
skew_post <- gene_skew(vst_mat)

n_skew_pre  <- sum(abs(skew_pre)  > 1.5, na.rm = TRUE)
n_skew_post <- sum(abs(skew_post) > 1.5, na.rm = TRUE)
cat(sprintf("\nGenes with |skewness| > 1.5  pre-VST  (log2): %d / %d (%.1f%%)\n",
            n_skew_pre, length(skew_pre), 100 * n_skew_pre / length(skew_pre)))
cat(sprintf("Genes with |skewness| > 1.5  post-VST       : %d / %d (%.1f%%)\n",
            n_skew_post, length(skew_post), 100 * n_skew_post / length(skew_post)))

# ---- Per-sample boxplots: pre vs post ----
png(boxplot_png, width = 2400, height = 1800, res = 200)
par(mfrow = c(2, 1), mar = c(6, 4, 3, 1))
boxplot(log2_mat, las = 2, outline = FALSE,
        col = "#9DC3E6", border = "#1F4E79",
        cex.axis = 0.5, ylab = "log2(count + 1)",
        main = "Pre-VST per-sample distribution (class-aware filter, PD vs HC)")
boxplot(vst_mat, las = 2, outline = FALSE,
        col = "#4C8DAE", border = "#1F4E79",
        cex.axis = 0.5, ylab = "VST expression",
        main = "Post-VST per-sample distribution (class-aware filter, PD vs HC)")
dev.off()

# ---- Mean-variance ----
png(meanvar_pre_png, width = 1500, height = 1200, res = 200)
mv_pre <- meanSdPlot(log2_mat, ranks = FALSE, plot = FALSE)
print(mv_pre$gg + ggtitle("Mean-variance pre-VST (log2(count+1), class-aware PD vs HC)"))
dev.off()

png(meanvar_post_png, width = 1500, height = 1200, res = 200)
mv_post <- meanSdPlot(vst_mat, ranks = FALSE, plot = FALSE)
print(mv_post$gg + ggtitle("Mean-variance post-VST (class-aware PD vs HC)"))
dev.off()

# ---- Skewness pre/post side by side ----
skew_xlim <- range(c(skew_pre, skew_post), na.rm = TRUE)
png(skew_png, width = 2200, height = 1000, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 4, 1))
hist(skew_pre, breaks = 80, xlim = skew_xlim,
     col = "#9DC3E6", border = "white",
     xlab = "Per-gene skewness of log2(count+1)", ylab = "Number of genes",
     main = sprintf("Pre-VST  |skew| > 1.5 : %d / %d (%.1f%%)",
                    n_skew_pre, length(skew_pre),
                    100 * n_skew_pre / length(skew_pre)))
abline(v = c(-1.5, 0, 1.5), lty = c(3, 2, 3),
       col = c("red", "red", "red"), lwd = c(1, 2, 1))
hist(skew_post, breaks = 80, xlim = skew_xlim,
     col = "#4C8DAE", border = "white",
     xlab = "Per-gene skewness of VST values", ylab = "Number of genes",
     main = sprintf("Post-VST |skew| > 1.5 : %d / %d (%.1f%%)",
                    n_skew_post, length(skew_post),
                    100 * n_skew_post / length(skew_post)))
abline(v = c(-1.5, 0, 1.5), lty = c(3, 2, 3),
       col = c("red", "red", "red"), lwd = c(1, 2, 1))
dev.off()

cat(sprintf("\nVST matrix (CSV)      : %s\n", vst_csv))
cat(sprintf("Boxplot (pre+post)    : %s\n", boxplot_png))
cat(sprintf("Mean-variance pre/post: %s, %s\n", meanvar_pre_png, meanvar_post_png))
cat(sprintf("Skewness (pre+post)   : %s\n", skew_png))
