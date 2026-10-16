source("setup.R")
source("config.R")

src <- list()
pdf(file.path(out_dir, "FigS3_genes_per_cell.pdf"), height = 12, width = 16)
par(mfrow = c(4, 4))
for (i in seq_len(nrow(sample_info))) {
    id <- sample_info$sample_id[i]
    f <- list.files(GEO_DIR, paste0("^", sample_info$GSM[i], "_", id, "_cell_feature_matrix.h5$"), full.names = TRUE)
    counts <- Read10X_h5(f)[[1]]   # 'Gene Expression', as in the per-sample reports
    counts <- counts[rowSums(counts) >= sample_info$gene_filter_min_total_count[i], ]
    if (sample_info$batch[i] == 1) counts <- counts[, colSums(counts) > 0]
    n_genes <- Matrix::colSums(counts > 0)
    h <- hist(n_genes, n = 100, col = "dodgerblue", xlab = "Genes per cell", main = id)
    src[[id]] <- data.frame(sample = id, bin_start = head(h$breaks, -1), bin_end = h$breaks[-1], n_cells = h$counts)
}
dev.off()
save_source(do.call(rbind, src), "FigS3_genes_per_cell_histograms")

sessionInfo()
