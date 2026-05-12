# QC: per-gene mean histogram (log10 scale) — PD vs HC (n=155)
# -----------------------------------------------------------------------------
# Subsets the main raw-count matrix to the 155 PD+HC samples listed in
# Meta_data_PDHC.csv (the full cohort: 102 PD = ES+LS, 53 HC), then plots a
# histogram of log10(rowMeans(counts)+1).
# -----------------------------------------------------------------------------

input_path    <- "C:/Users/hafsa/Python PD Project/MI_BaggedLASSO Pipeline/2. Dataset_Rationale_Refined_Features.csv"
metadata_path <- "C:/Users/hafsa/Python PD Project/MI_BaggedLASSO Pipeline/PDvsHC/Meta_data_PDHC.csv"
output_dir    <- "C:/Users/hafsa/Python PD Project/MI_BaggedLASSO Pipeline/PDHC new/Outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_png    <- file.path(output_dir, "qc_gene_mean_histogram_PDHC.png")

# ---- Load raw counts ----
counts <- as.matrix(read.csv(input_path, row.names = 1, check.names = FALSE))
if (nrow(counts) == 155 && ncol(counts) == 28844) counts <- t(counts)
stopifnot(nrow(counts) == 28844, ncol(counts) == 155)
storage.mode(counts) <- "integer"

# ---- Subset to PD+HC samples (full cohort, 155) ----
meta <- read.csv(metadata_path, stringsAsFactors = FALSE)
keep_samples <- intersect(colnames(counts), meta$sample_id)
stopifnot(length(keep_samples) == nrow(meta))
counts <- counts[, keep_samples]
meta   <- meta[match(keep_samples, meta$sample_id), ]
cat(sprintf("PD vs HC subset: %d genes x %d samples\n",
            nrow(counts), ncol(counts)))
print(table(meta$condition))

# ---- Plot ----
log10_mean <- log10(rowMeans(counts) + 1)

png(output_png, width = 1800, height = 1000, res = 200)
hist(log10_mean, breaks = 100,
     col = "#4C8DAE", border = "white",
     xlab = "log10(mean count + 1)",
     ylab = "Number of genes",
     main = sprintf("Per-gene mean counts (log10 scale) — PD vs HC, n=%d",
                    ncol(counts)))
abline(v = log10(2),  lty = 2, col = "red",    lwd = 1.5)
abline(v = log10(11), lty = 2, col = "orange", lwd = 1.5)
legend("topright", lty = 2, col = c("red", "orange"),
       legend = c("mean = 1", "mean = 10"), bty = "n")
dev.off()

cat(sprintf("\nQC plot: %s\n", output_png))
