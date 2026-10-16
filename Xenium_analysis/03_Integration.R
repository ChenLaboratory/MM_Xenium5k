## 03_Integration.R

source("setup.R")
source("config.R")

samples_int <- c("MM-4", "MGUS-2", "MM-1", "MM-2", "MM-3", "Ctrl-1", "MGUS-1", "SM-1",
    "MM-10", "SM-5", "MM-9", "Ctrl-4", "Ctrl-3", "MM-8", "MM-7", "Ctrl-2", "SM-4", "MM-6",
    "MM-5", "SM-3", "SM-2")
annotation <- read.delim(file.path(annotation_dir, "per_sample_cluster_annotation.tsv"))
cells <- readRDS(file.path(OUTPUT_DIR, "per_sample_cells_spe.rds"))

integrate_cells <- function(name, cell_class, resolution_1, remove_1, nfeatures_2,
                            resolution_2, remove_2, label_column, labels, label_levels) {
    ## ---- Combine cells of one class from all samples ---------------------
    counts <- coord <- list()
    for (id in samples_int) {
        spe <- cells[[id]]
        clst <- annotation$Cluster[annotation$Sample == id & annotation$Celltype == cell_class]
        spe <- spe[, spe$seurat_clusters %in% clst]
        counts[[id]] <- counts(spe)
        colnames(counts[[id]]) <- paste(id, colnames(spe), sep = "_")
        coord[[id]] <- spatialCoords(spe)
    }
    so <- CreateSeuratObject(counts = do.call(cbind, counts), assay = "RNA", min.cells = 0, min.features = 0)
    so$orig.ident <- factor(rep(samples_int, times = sapply(counts, ncol)), levels = samples)
    so$x_centroid <- do.call(rbind, coord)[, 1]
    so$y_centroid <- do.call(rbind, coord)[, 2]
    so$group <- factor(sample_info$analysis_group[match(so$orig.ident, sample_info$sample_id)],
        levels = c("CTRL", "MGUS", "MM"))
    so <- so[rowSums(so[["RNA"]]$counts) >= 1000, ]

    ## ---- First clustering ------------------------------------------------
    dimUsed <- 30
    so <- NormalizeData(so, scale.factor = 100, verbose = FALSE)
    so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
    so <- ScaleData(so, features = VariableFeatures(so), verbose = FALSE)
    so <- RunPCA(so, dims = 1:dimUsed, features = VariableFeatures(so), verbose = FALSE)
    so <- RunUMAP(so, dims = 1:dimUsed, seed.use = 2024, verbose = FALSE)
    so[["spatial"]] <- CreateDimReducObject(embeddings = cbind(Sp_1 = so$x_centroid, Sp_2 = so$y_centroid),
        key = "Sp_", assay = "RNA")
    so <- FindNeighbors(so, dims = 1:dimUsed, verbose = FALSE)
    so <- FindClusters(so, resolution = resolution_1, verbose = FALSE)
    print(DimPlot(so, reduction = "umap", raster = FALSE, cols = col.p2))

    ## ---- Removal of unwanted clusters and re-clustering ------------------
    so <- so[, !so$seurat_clusters %in% remove_1]
    so <- FindVariableFeatures(so, nfeatures = nfeatures_2, verbose = FALSE)
    so <- ScaleData(so, features = VariableFeatures(so), verbose = FALSE)
    so <- RunPCA(so, dims = 1:dimUsed, features = VariableFeatures(so), verbose = FALSE)
    so <- RunUMAP(so, dims = 1:dimUsed, seed.use = 2024, verbose = FALSE)
    so <- FindNeighbors(so, dims = 1:dimUsed, verbose = FALSE)
    so <- FindClusters(so, resolution = resolution_2, verbose = FALSE)
    if (length(remove_2)) {
        so <- so[, !so$seurat_clusters %in% remove_2]
        so$seurat_clusters <- factor(as.integer(droplevels(so$seurat_clusters)) - 1)
        Idents(so) <- "seurat_clusters"
    }
    so[[label_column]] <- factor(labels[as.integer(so$seurat_clusters)], levels = label_levels)
    print(table(so$seurat_clusters, so$orig.ident))
    print(DimPlot(so, reduction = "umap", group.by = label_column, raster = FALSE,
        cols = setNames(col.p2[seq_along(labels)], labels)))

    ## ---- Cluster markers -------------------------------------------------
    markers <- FindAllMarkers(so, logfc.threshold = 1, min.pct = 0.1, only.pos = TRUE)
    markers <- markers[markers$p_val_adj < 0.05, ]
    markers[[label_column]] <- labels[as.integer(markers$cluster)]
    write.csv(markers, file.path(OUTPUT_DIR, "tables", paste0(name, "_cluster_markers.csv")),
        row.names = FALSE)
    write.csv(as.data.frame.matrix(table(so[[label_column, drop = TRUE]], so$orig.ident)),
        file.path(OUTPUT_DIR, "tables", paste0(name, "_cluster_counts_per_sample.csv")))
    so <- DietSeurat(so, layers = c("counts", "data"), dimreducs = c("pca", "umap", "spatial"))
    saveRDS(so, file.path(OUTPUT_DIR, "integration", paste0(name, "_cells_integrated_seurat.rds")))
    so
}

## ---- Plasma cells (Fig 3, Fig 4B) ------------------------------------
so_pc <- integrate_cells(name = "plasma", cell_class = "plasma",
    resolution_1 = 1.0, remove_1 = c(7, 10, 14, 15, 17, 19, 21, 24),
    nfeatures_2 = 5000, resolution_2 = 0.5, remove_2 = integer(0),
    label_column = "PC_cluster_paper",
    labels = paste0("PC", c(9, 2, 1, 3:8, 0, 10:16)), label_levels = paste0("PC", 0:16))

## ---- Stromal cells (Fig 5) -------------------------------------------
so_st <- integrate_cells(name = "stromal", cell_class = "non",
    resolution_1 = 0.5, remove_1 = c(2, 3, 6, 8, 9, 10, 11, 13, 14, 16, 17, 18, 19, 21),
    nfeatures_2 = 1500, resolution_2 = 0.4, remove_2 = 6,
    label_column = "stromal_cell_type",
    labels = c("MSC", "SEC2", "SEC1", "Osteoblast1", "VSMC", "Adipocyte", "AEC",
               "Osteoblast2", "Cartilage", "Osteo-MSC", "Osteoclast"),
    label_levels = c("SEC1", "SEC2", "AEC", "VSMC", "MSC", "Osteo-MSC", "Osteoblast1",
                     "Osteoblast2", "Osteoclast", "Cartilage", "Adipocyte"))
