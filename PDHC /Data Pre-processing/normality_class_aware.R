# Normality assessment for the class-aware filtered data — PD vs HC (n=155)
# -----------------------------------------------------------------------------
# Tests per-gene normality on:
#   (1) raw class-aware filtered counts  (filter_class_counts_PDHC.csv)  — expected NOT normal
#   (2) VST-transformed data             (vst_class_aware_PDHC.csv)      — expected close to normal
# so you can quantify how much VST helped.
#
# Per gene the script computes:
#   - Shapiro-Wilk p-value          (H0: gene is normally distributed across the 155 samples)
#   - skewness                      (0 = symmetric)
#   - excess kurtosis               (0 = normal-tailed; >0 = heavy tails)
#
# Aggregate views:
#   a) Histogram of Shapiro p-values + fraction of genes with p > 0.05
#   b) Histograms of per-gene skewness and excess kurtosis
#   c) 3x3 Q-Q-plot grid for genes spanning the expression spectrum
#   d) 3x3 histogram-with-normal-overlay grid for the same genes
# Plus a summary TSV/CSV.
# -----------------------------------------------------------------------------

# ---- 1. Inputs / outputs ----------------------------------------------------
raw_path <- "C:/Users/hafsa/Python PD Project/MI_BaggedLASSO Pipeline/PDHC new/Outputs/QC+Filter by Exp/filter_class_counts_PDHC.csv"
vst_path <- "C:/Users/hafsa/Python PD Project/MI_BaggedLASSO Pipeline/PDHC new/Outputs/VST/vst_class_aware_PDHC.csv"

output_dir <- "C:/Users/hafsa/Python PD Project/MI_BaggedLASSO Pipeline/PDHC new/Outputs/Normality Analysis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pval_png    <- file.path(output_dir, "normality_class_pvalues_PDHC.png")
skewkurt_png<- file.path(output_dir, "normality_class_skewkurt_PDHC.png")
qq_png      <- file.path(output_dir, "normality_class_qq_grid_PDHC.png")
hist_png    <- file.path(output_dir, "normality_class_hist_grid_PDHC.png")
summary_tsv <- file.path(output_dir, "normality_class_summary_PDHC.tsv")
summary_csv <- file.path(output_dir, "normality_class_summary_PDHC.csv")

# ---- 2. Helpers -------------------------------------------------------------
load_matrix <- function(path) {
  m <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  if (nrow(m) == 155 && ncol(m) != 155) m <- t(m)   # defensive re-orient
  m
}

# Per-gene metrics across the 155 samples (rows = genes, cols = samples)
per_gene_normality <- function(mat) {
  shapiro_p <- apply(mat, 1, function(x) {
    if (length(unique(x)) < 3) return(NA_real_)        # constant / near-constant
    tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
  })
  mu  <- rowMeans(mat)
  sdv <- apply(mat, 1, sd)
  z   <- (mat - mu) / sdv
  skew <- rowMeans(z^3)
  kurt <- rowMeans(z^4) - 3                            # excess kurtosis
  data.frame(mean_expr = mu, sd_expr = sdv,
             shapiro_p = shapiro_p, skew = skew, kurt = kurt,
             row.names = rownames(mat))
}

# Pick 9 representative genes spanning expression quantiles
pick_repr_genes <- function(mat, n = 9) {
  mu <- rowMeans(mat)
  q  <- quantile(mu, probs = seq(0.1, 0.9, length.out = n), na.rm = TRUE)
  vapply(q, function(target) {
    rownames(mat)[which.min(abs(mu - target))]
  }, character(1))
}

# ---- 3. Run on raw and (if present) VST -------------------------------------
cat("Loading raw class-aware counts ...\n")
raw_mat <- load_matrix(raw_path)
storage.mode(raw_mat) <- "double"  # Shapiro needs numeric, not int
cat(sprintf("  dim: %d genes x %d samples\n", nrow(raw_mat), ncol(raw_mat)))

cat("Computing per-gene normality on RAW ...\n")
raw_stats <- per_gene_normality(raw_mat)

have_vst <- file.exists(vst_path)
if (have_vst) {
  cat("Loading VST class-aware matrix ...\n")
  vst_mat <- load_matrix(vst_path)
  cat(sprintf("  dim: %d genes x %d samples\n", nrow(vst_mat), ncol(vst_mat)))
  cat("Computing per-gene normality on VST ...\n")
  vst_stats <- per_gene_normality(vst_mat)
} else {
  cat(sprintf("VST file not found at %s — running RAW-only.\n", vst_path))
  vst_stats <- NULL
}

# ---- 4. Plot a) Shapiro p-value histograms ---------------------------------
png(pval_png, width = 1800, height = 900, res = 200)
par(mfrow = c(1, if (have_vst) 2 else 1), mar = c(4, 4, 3, 1))

frac_norm_raw <- mean(raw_stats$shapiro_p > 0.05, na.rm = TRUE)
hist(raw_stats$shapiro_p, breaks = 50, col = "#B85450", border = "white",
     xlab = "Shapiro-Wilk p-value", ylab = "Genes",
     main = sprintf("RAW counts (PD vs HC) — %.1f%% of genes p > 0.05",
                    100 * frac_norm_raw))
abline(v = 0.05, lty = 2, col = "blue")

if (have_vst) {
  frac_norm_vst <- mean(vst_stats$shapiro_p > 0.05, na.rm = TRUE)
  hist(vst_stats$shapiro_p, breaks = 50, col = "#4C8DAE", border = "white",
       xlab = "Shapiro-Wilk p-value", ylab = "Genes",
       main = sprintf("VST values (PD vs HC) — %.1f%% of genes p > 0.05",
                      100 * frac_norm_vst))
  abline(v = 0.05, lty = 2, col = "blue")
}
dev.off()

# ---- 5. Plot b) Skewness + kurtosis histograms -----------------------------
png(skewkurt_png, width = 1800, height = if (have_vst) 1500 else 800, res = 200)
par(mfrow = c(if (have_vst) 2 else 1, 2), mar = c(4, 4, 3, 1))

hist(raw_stats$skew, breaks = 80, col = "#B85450", border = "white",
     xlab = "Skewness", main = "RAW (PD vs HC) — per-gene skewness")
abline(v = 0, lty = 2, col = "blue")
hist(raw_stats$kurt, breaks = 80, col = "#B85450", border = "white",
     xlab = "Excess kurtosis", main = "RAW (PD vs HC) — per-gene excess kurtosis")
abline(v = 0, lty = 2, col = "blue")

if (have_vst) {
  hist(vst_stats$skew, breaks = 80, col = "#4C8DAE", border = "white",
       xlab = "Skewness", main = "VST (PD vs HC) — per-gene skewness")
  abline(v = 0, lty = 2, col = "blue")
  hist(vst_stats$kurt, breaks = 80, col = "#4C8DAE", border = "white",
       xlab = "Excess kurtosis", main = "VST (PD vs HC) — per-gene excess kurtosis")
  abline(v = 0, lty = 2, col = "blue")
}
dev.off()

# ---- 6. Plot c+d) Q-Q + histogram grid for 9 representative genes ----------
diag_mat   <- if (have_vst) vst_mat else raw_mat
diag_label <- if (have_vst) "VST" else "RAW"
repr_genes <- pick_repr_genes(diag_mat)

png(qq_png, width = 1800, height = 1800, res = 200)
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))
for (g in repr_genes) {
  x <- diag_mat[g, ]
  qqnorm(x, main = sprintf("%s (%s)\nmean=%.2f, p=%.2g",
                           g, diag_label, mean(x),
                           tryCatch(shapiro.test(x)$p.value,
                                    error = function(e) NA_real_)),
         pch = 19, col = "#1F4E79", cex = 0.6)
  qqline(x, col = "red", lwd = 1.5)
}
dev.off()

png(hist_png, width = 1800, height = 1800, res = 200)
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))
for (g in repr_genes) {
  x <- diag_mat[g, ]
  hist(x, breaks = 20, freq = FALSE,
       col = "#4C8DAE", border = "white",
       xlab = sprintf("%s (%s)", g, diag_label),
       main = sprintf("mean=%.2f sd=%.2f", mean(x), sd(x)))
  curve(dnorm(x, mean(diag_mat[g, ]), sd(diag_mat[g, ])),
        add = TRUE, col = "red", lwd = 2)
}
dev.off()

# ---- 7. Summary table ------------------------------------------------------
mk_summary_row <- function(label, stats) {
  data.frame(
    matrix              = label,
    n_genes             = nrow(stats),
    mean_shapiro_p      = mean(stats$shapiro_p, na.rm = TRUE),
    median_shapiro_p    = median(stats$shapiro_p, na.rm = TRUE),
    pct_genes_p_gt_0.05 = 100 * mean(stats$shapiro_p > 0.05, na.rm = TRUE),
    median_abs_skew     = median(abs(stats$skew), na.rm = TRUE),
    median_excess_kurt  = median(stats$kurt, na.rm = TRUE)
  )
}
summary_df <- mk_summary_row("raw_filtered", raw_stats)
if (have_vst) summary_df <- rbind(summary_df, mk_summary_row("vst", vst_stats))

write.table(summary_df, summary_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(summary_df, summary_csv, row.names = FALSE)

cat("\n=========== Summary (PD vs HC, class-aware) ===========\n")
print(summary_df, row.names = FALSE)

verdict <- function(stats, label) {
  pct <- 100 * mean(stats$shapiro_p > 0.05, na.rm = TRUE)
  cat(sprintf("\n[%s]  %.1f%% of genes are not rejected as normal at alpha=0.05.\n",
              label, pct))
  if (pct >= 70) {
    cat("  -> Approximately normal. Methods that assume Gaussianity are reasonable.\n")
  } else if (pct >= 30) {
    cat("  -> Partially normal. Use rank-based or robust methods, or rely on the\n")
    cat("     central limit theorem if your downstream test averages many genes.\n")
  } else {
    cat("  -> Largely non-normal. Stick to non-parametric / count-based methods,\n")
    cat("     or apply a stronger transform (e.g. rank-INT) before parametric ML.\n")
  }
}
verdict(raw_stats, "RAW")
if (have_vst) verdict(vst_stats, "VST")

cat(sprintf("\nP-value plot   : %s\n", pval_png))
cat(sprintf("Skew/kurt plot : %s\n", skewkurt_png))
cat(sprintf("Q-Q grid       : %s\n", qq_png))
cat(sprintf("Histogram grid : %s\n", hist_png))
cat(sprintf("Summary (TSV)  : %s\n", summary_tsv))
cat(sprintf("Summary (CSV)  : %s\n", summary_csv))
