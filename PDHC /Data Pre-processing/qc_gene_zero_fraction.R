# QC: per-gene zero-fraction histogram — PD vs HC (n=155)
# -----------------------------------------------------------------------------
# Per gene, fraction of samples with count == 0; histogram of those fractions.
# -----------------------------------------------------------------------------

input_path    <- "C:/Users/hafsa/Python PD Project/MI_BaggedLASSO Pipeline/2. Dataset_Rationale_Refined_Features.csv"
metadata_path <- "C:/Users/hafsa/Python PD Project/MI_BaggedLASSO Pipeline/PDvsHC/Meta_data_PDHC.csv"
output_dir    <- "C:/Users/hafsa/Python PD Project/MI_BaggedLASSO Pipeline/PDHC new/Outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_png    <- file.path(output_dir, "qc_gene_zero_fraction_PDHC.png")
output_tsv    <- file.path(output_dir, "qc_gene_zero_fraction_PDHC.tsv")

counts <- as.matrix(read.csv(input_path, row.names = 1, check.names = FALSE))
if (nrow(counts) == 155 && ncol(counts) == 28844) counts <- t(counts)
stopifnot(nrow(counts) == 28844, ncol(counts) == 155)
storage.mode(counts) <- "integer"

meta <- read.csv(metadata_path, stringsAsFactors = FALSE)
keep_samples <- intersect(colnames(counts), meta$sample_id)
stopifnot(length(keep_samples) == nrow(meta))
counts <- counts[, keep_samples]

n_samples <- ncol(counts)
zero_fraction <- rowMeans(counts == 0)
cat(sprintf("PD vs HC subset: %d genes x %d samples\n",
            nrow(counts), n_samples))
cat(sprintf("Genes with zero fraction == 1.0 : %d\n",
            sum(zero_fraction == 1.0)))
cat(sprintf("Genes with zero fraction == 0.0 : %d\n",
            sum(zero_fraction == 0.0)))
cat(sprintf("Genes with zero fraction >= 0.50: %d\n",
            sum(zero_fraction >= 0.5)))

png(output_png, width = 1800, height = 1000, res = 200)
hist(zero_fraction, breaks = seq(0, 1, by = 0.02),
     col = "#B85450", border = "white",
     xlab = "Fraction of samples with count == 0",
     ylab = "Number of genes",
     main = sprintf("Per-gene zero fraction — PD vs HC (%d samples)",
                    n_samples))
abline(v = 0.5, lty = 2, col = "black")
text(0.51, par("usr")[4] * 0.9, "50% zeros", pos = 4, cex = 0.85)
dev.off()

write.table(data.frame(gene = rownames(counts), zero_fraction = zero_fraction),
            output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\nQC plot : %s\nQC table: %s\n", output_png, output_tsv))
