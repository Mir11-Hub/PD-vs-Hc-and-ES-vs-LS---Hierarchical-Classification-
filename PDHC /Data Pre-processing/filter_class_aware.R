# Class-aware gene filter (PD vs HC, n=155)
# -----------------------------------------------------------------------------
# Wraps edgeR::filterByExpr with `group = condition` (PD vs HC). With the
# minority class (HC) at 53 samples and the majority class (PD) at 102,
# filterByExpr requires detection in ~70% of the smallest class plus a
# total-count floor.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(edgeR))

input_path        <- "C:/Users/hafsa/Python PD Project/MI_BaggedLASSO Pipeline/2. Dataset_Rationale_Refined_Features.csv"
metadata_path     <- "C:/Users/hafsa/Python PD Project/MI_BaggedLASSO Pipeline/PDvsHC/Meta_data_PDHC.csv"
output_dir        <- "C:/Users/hafsa/Python PD Project/MI_BaggedLASSO Pipeline/PDHC new/Outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_keep_tsv   <- file.path(output_dir, "filter_class_keep_PDHC.tsv")
output_keep_csv   <- file.path(output_dir, "filter_class_keep_PDHC.csv")
output_counts_tsv <- file.path(output_dir, "filter_class_counts_PDHC.tsv")
output_counts_csv <- file.path(output_dir, "filter_class_counts_PDHC.csv")
output_png        <- file.path(output_dir, "filter_class_zero_fraction_PDHC.png")

# ---- Load counts ----
counts <- as.matrix(read.csv(input_path, row.names = 1, check.names = FALSE))
if (nrow(counts) == 155 && ncol(counts) == 28844) counts <- t(counts)
stopifnot(nrow(counts) == 28844, ncol(counts) == 155)
storage.mode(counts) <- "integer"

# ---- Align metadata (full 155-sample cohort) ----
meta <- read.csv(metadata_path, stringsAsFactors = FALSE)
keep_samples <- intersect(colnames(counts), meta$sample_id)
stopifnot(length(keep_samples) == nrow(meta))
counts <- counts[, keep_samples]
meta   <- meta[match(keep_samples, meta$sample_id), ]
condition <- factor(meta$condition, levels = c("HC", "PD"))
cat("Sample counts by condition:\n"); print(table(condition))

# ---- filterByExpr ----
keep <- filterByExpr(counts, group = condition)
counts_filt <- counts[keep, ]
cat(sprintf("\nGenes kept: %d / %d (%.1f%%)\n",
            sum(keep), length(keep), 100 * mean(keep)))
cat(sprintf("Filtered matrix dim: %d genes x %d samples\n",
            nrow(counts_filt), ncol(counts_filt)))

# ---- QC: zero fraction before vs after ----
zf_before <- rowMeans(counts == 0)
zf_after  <- rowMeans(counts_filt == 0)

png(output_png, width = 1800, height = 1000, res = 200)
par(mfrow = c(1, 2))
hist(zf_before, breaks = seq(0, 1, by = 0.02),
     col = "#B85450", border = "white",
     main = "Before filter", xlab = "Zero fraction", ylab = "Genes")
hist(zf_after, breaks = seq(0, 1, by = 0.02),
     col = "#4C8DAE", border = "white",
     main = sprintf("After class-aware filter (n=%d)", nrow(counts_filt)),
     xlab = "Zero fraction", ylab = "Genes")
dev.off()

# ---- Save keep flags + filtered counts ----
keep_df <- data.frame(
  gene        = rownames(counts),
  keep        = keep,
  total_count = rowSums(counts),
  zero_frac   = zf_before
)
write.table(keep_df, output_keep_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

counts_filt_df <- data.frame(gene = rownames(counts_filt), counts_filt,
                             check.names = FALSE)
write.table(counts_filt_df, output_counts_tsv, sep = "\t",
            quote = FALSE, row.names = FALSE)

# ---- TSV -> CSV faithful re-encoding ----
write.csv(read.table(output_keep_tsv, header = TRUE, sep = "\t",
                     check.names = FALSE, stringsAsFactors = FALSE),
          output_keep_csv, row.names = FALSE)
write.csv(read.table(output_counts_tsv, header = TRUE, sep = "\t",
                     check.names = FALSE, stringsAsFactors = FALSE),
          output_counts_csv, row.names = FALSE)

cat(sprintf("\nKeep table     : %s\n", output_keep_csv))
cat(sprintf("Filtered counts: %s\n", output_counts_csv))
cat(sprintf("QC plot        : %s\n", output_png))
