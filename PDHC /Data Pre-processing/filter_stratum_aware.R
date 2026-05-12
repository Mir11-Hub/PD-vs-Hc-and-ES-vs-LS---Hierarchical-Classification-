# Stratum-aware gene filter (condition x cell_type x sex) — PD vs HC, n=155
# -----------------------------------------------------------------------------
# Per gene, per class, per (cell_type x sex) sub-stratum, compute fraction of
# samples with count >= min_count (CPM scaled). Gene passes if its best class
# has >= min_prop detection in EVERY sub-stratum, plus a total-count floor.
# More conservative than the class-aware filter; preserves class-discriminative
# genes while requiring confounder-robust detection (sex + cell_type).
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(edgeR))

input_path        <- "C:/Users/hafsa/Python PD Project/MI_BaggedLASSO Pipeline/2. Dataset_Rationale_Refined_Features.csv"
metadata_path     <- "C:/Users/hafsa/Python PD Project/MI_BaggedLASSO Pipeline/PDvsHC/Meta_data_PDHC.csv"
output_dir        <- "C:/Users/hafsa/Python PD Project/MI_BaggedLASSO Pipeline/PDHC new/Outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_keep       <- file.path(output_dir, "filter_stratum_keep_PDHC.tsv")
output_keep_csv   <- file.path(output_dir, "filter_stratum_keep_PDHC.csv")
output_counts_tsv <- file.path(output_dir, "filter_stratum_counts_PDHC.tsv")
output_counts_csv <- file.path(output_dir, "filter_stratum_counts_PDHC.csv")
output_png        <- file.path(output_dir, "filter_stratum_zero_fraction_PDHC.png")

# Tunables
min_count       <- 10
min_total_count <- 15
min_prop        <- 0.7
min_stratum_n   <- 3

# ---- Load + align to 155-sample cohort ----
counts <- as.matrix(read.csv(input_path, row.names = 1, check.names = FALSE))
if (nrow(counts) == 155 && ncol(counts) == 28844) counts <- t(counts)
stopifnot(nrow(counts) == 28844, ncol(counts) == 155)
storage.mode(counts) <- "integer"

meta <- read.csv(metadata_path, stringsAsFactors = FALSE)
keep_samples <- intersect(colnames(counts), meta$sample_id)
stopifnot(length(keep_samples) == nrow(meta))
counts <- counts[, keep_samples]
meta   <- meta[match(keep_samples, meta$sample_id), ]

condition <- factor(meta$condition, levels = c("HC", "PD"))
sex       <- factor(meta$sex)
cell_type <- factor(meta$cell_type)

cat("Sample counts by condition x cell_type x sex:\n")
print(addmargins(table(condition, cell_type, sex)))

# ---- CPM threshold equivalent to min_count at median lib size ----
dge        <- DGEList(counts = counts)
median_lib <- median(dge$samples$lib.size)
cpm_thresh <- min_count / (median_lib / 1e6)
cat(sprintf("\nMedian library size: %s reads\nCPM threshold (= count>=%d at median lib): %.3f\n",
            format(median_lib, big.mark = ","), min_count, cpm_thresh))

detected <- cpm(dge) >= cpm_thresh   # gene x sample logical

# ---- Per-class, per-(cell_type x sex) detection rate ----
class_min_rate <- matrix(0,
                         nrow = nrow(counts), ncol = nlevels(condition),
                         dimnames = list(rownames(counts), levels(condition)))

for (cls in levels(condition)) {
  cls_idx <- which(condition == cls)
  sub_stratum <- droplevels(interaction(cell_type[cls_idx], sex[cls_idx], drop = TRUE))

  rates <- vapply(levels(sub_stratum), function(s) {
    s_cols <- cls_idx[sub_stratum == s]
    if (length(s_cols) < min_stratum_n) {
      return(rep(NA_real_, nrow(counts)))
    }
    rowMeans(detected[, s_cols, drop = FALSE])
  }, numeric(nrow(counts)))

  class_min_rate[, cls] <- apply(rates, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA_real_ else min(x)
  })
}

best_class_rate <- apply(class_min_rate, 1, max, na.rm = TRUE)
total_count     <- rowSums(counts)
keep <- best_class_rate >= min_prop & total_count >= min_total_count

cat(sprintf("\nGenes kept: %d / %d (%.1f%%)\n",
            sum(keep), length(keep), 100 * mean(keep)))
cat(sprintf("Dropped by stratum rule  : %d\n",  sum(best_class_rate < min_prop)))
cat(sprintf("Dropped by total-count   : %d\n",  sum(total_count < min_total_count)))
cat("\nWhich class drove retention (argmax) for kept genes:\n")
print(table(apply(class_min_rate[keep, , drop = FALSE], 1, function(x) {
  levels(condition)[which.max(replace(x, is.na(x), -Inf))]
})))

# ---- QC plot ----
zf_before <- rowMeans(counts == 0)
zf_after  <- zf_before[keep]

png(output_png, width = 1800, height = 1000, res = 200)
par(mfrow = c(1, 2))
hist(zf_before, breaks = seq(0, 1, by = 0.02),
     col = "#B85450", border = "white",
     main = "Before filter", xlab = "Zero fraction", ylab = "Genes")
hist(zf_after, breaks = seq(0, 1, by = 0.02),
     col = "#4C8DAE", border = "white",
     main = sprintf("After stratum-aware filter (n=%d)", sum(keep)),
     xlab = "Zero fraction", ylab = "Genes")
dev.off()

# ---- Save keep table + filtered counts ----
out <- data.frame(
  gene             = rownames(counts),
  keep             = keep,
  total_count      = total_count,
  best_class_rate  = best_class_rate,
  class_min_rate
)
write.table(out, output_keep, sep = "\t", quote = FALSE, row.names = FALSE)

counts_filt    <- counts[keep, ]
counts_filt_df <- data.frame(gene = rownames(counts_filt), counts_filt,
                             check.names = FALSE)
write.table(counts_filt_df, output_counts_tsv, sep = "\t",
            quote = FALSE, row.names = FALSE)

write.csv(read.table(output_keep, header = TRUE, sep = "\t",
                     check.names = FALSE, stringsAsFactors = FALSE),
          output_keep_csv, row.names = FALSE)
write.csv(read.table(output_counts_tsv, header = TRUE, sep = "\t",
                     check.names = FALSE, stringsAsFactors = FALSE),
          output_counts_csv, row.names = FALSE)

cat(sprintf("\nKeep table     : %s\n", output_keep_csv))
cat(sprintf("Filtered counts: %s\n", output_counts_csv))
cat(sprintf("QC plot        : %s\n", output_png))
