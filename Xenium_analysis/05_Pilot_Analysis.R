## 05_Pilot_Analysis.R

source("setup.R")
source("config.R")
library(SingleR)

pilot <- list(
    "RM-1_377gene_panel" = list(prefix = geo_prefix("RM-1"), min_transcripts = 10,
        min_genes = 5, resolution = 0.2, reference = "human"),
    "mouse_BM_379gene_panel" = list(prefix = MOUSE_PREFIX, min_transcripts = 15,
        min_genes = 10, resolution = 0.2, reference = "mouse"))

for (name in names(pilot)) {
    p <- pilot[[name]]

    ## ---- Read in and QC --------------------------------------------------
    counts <- Read10X_h5(paste0(p$prefix, "_cell_feature_matrix.h5"))[[1]]
    cell_info <- read.csv(paste0(p$prefix, "_cells.csv.gz"))
    rownames(cell_info) <- cell_info$cell_id
    cell_info <- cell_info[, -1]
    so <- CreateSeuratObject(counts = counts, project = name, assay = "RNA",
        meta.data = cell_info, min.cells = 0, min.features = 0)
    so <- so[, so$nCount_RNA >= p$min_transcripts & so$nFeature_RNA >= p$min_genes]

    ## ---- Clustering (all panel genes are used for PCA) -------------------
    dimUsed <- 20
    so <- NormalizeData(so, scale.factor = 100, verbose = FALSE)
    so <- ScaleData(so, features = rownames(so), verbose = FALSE)
    so <- RunPCA(so, features = rownames(so), dims = 1:dimUsed, verbose = FALSE)
    so <- RunUMAP(so, dims = 1:dimUsed, seed.use = 2024, verbose = FALSE)
    so[["spatial"]] <- CreateDimReducObject(embeddings = cbind(Sp_1 = so$x_centroid, Sp_2 = so$y_centroid),
        key = "Sp_", assay = "RNA")
    so <- FindNeighbors(so, dims = 1:dimUsed, verbose = FALSE)
    so <- FindClusters(so, resolution = p$resolution, verbose = FALSE)

    ## ---- SingleR cell type annotation ------------------------------------
    if (p$reference == "human") {
        sce_ref <- readRDS(file.path(REF_DIR, "sce_ref_Jeremy.rds"))
        labels <- sce_ref$cluster_Jeremy
        label_column <- "SingleR_Jeremy"
        cols <- cols_ct
    } else {
        load(file.path(REF_DIR, "NicheData10x.rda"))
        sce_ref <- as.SingleCellExperiment(UpdateSeuratObject(NicheData10x))
        labels <- sce_ref$ident
        label_column <- "SingleR_labels_2"
        cols <- NicheDataColors
    }
    pred <- SingleR(test = as.SingleCellExperiment(so), ref = sce_ref, labels = labels,
        de.method = "wilcox", num.threads = 4)
    so[[label_column]] <- pred$labels
    print(table(pred$labels))
    so <- DietSeurat(so, layers = c("counts", "data"), dimreducs = c("pca", "umap", "spatial"))
    saveRDS(so, file.path(OUTPUT_DIR, "pilot", paste0(name, "_seurat.rds")))

    print(DimPlot(so, reduction = "umap", group.by = label_column, raster = FALSE, cols = cols) +
        DimPlot(so, reduction = "spatial", group.by = label_column, raster = FALSE, cols = cols,
            pt.size = 0.03) + coord_fixed())
}
