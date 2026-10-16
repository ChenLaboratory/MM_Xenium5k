## 02_Pseudobulk_DE.R

source("setup.R")
source("config.R")

## ---- Build pseudo-bulk samples ---------------------------------------

cells <- readRDS(file.path(OUTPUT_DIR, "per_sample_cells_spe.rds"))
ct <- sapply(samples, function(id) {
    spe <- cells[[id]]
    Matrix::rowSums(counts(spe)[, spe$cell_type == "Plasma Cell"])
})
group <- factor(sample_info$analysis_group[match(samples, sample_info$sample_id)],
    levels = c("CTRL", "MGUS", "MM"))
y <- DGEList(counts = ct, group = group)

## ---- Gene annotation + filtering -------------------------------------
anno <- read.delim(file.path(annotation_dir, "240110_Homo_sapiens.gene_info.gz"))[, -1]
m <- match(alias2SymbolTable(rownames(y)), anno$Symbol)
y$genes <- setNames(anno[m, c("Symbol", "GeneID", "chromosome", "type_of_gene")],
    c("Symbol", "GeneID", "Chr", "Type"))
y <- y[filterByExpr(y, min.count = 30, min.total.count = 50), , keep.lib.sizes = FALSE]
y <- normLibSizes(y, method = "TMM")
plotMDS(y, col = col.p2[group])

## ---- Differential expression -----------------------------------------
design <- model.matrix(~ group)
colnames(design) <- c("Int", "MGUS", "MM")
fit <- glmQLFit(y, design, robust = TRUE)
tests <- list("MM_vs_MGUS" = glmQLFTest(fit, contrast = makeContrasts(MM - MGUS, levels = design)),
              "MGUS_vs_Ctrl" = glmQLFTest(fit, coef = 2))
for (k in names(tests)) {
    print(topTags(tests[[k]], n = 10))
    print(summary(decideTests(tests[[k]])))
    write.csv(topTags(tests[[k]], n = Inf),
        file.path(OUTPUT_DIR, "tables", paste0("pseudobulk_DEG-PC-", k, ".csv")), row.names = FALSE)
    write.csv(topKEGG(kegga(tests[[k]], geneid = "GeneID", species = "Hs"), n = Inf),
        file.path(OUTPUT_DIR, "tables", paste0("pseudobulk_KEGG-PC-", k, ".csv")), row.names = FALSE)
}
