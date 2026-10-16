## 01_Per_Sample_Analysis.R

source("setup.R")
source("config.R")
library(SingleR)

## Reference: human bone marrow atlas of Bandyopadhyay et al. (Cell 2024),
## with related cell types merged into 14 classes.
sce_ref <- readRDS(file.path(REF_DIR, "sce_ref_Jeremy.rds"))
GOI <- read.csv(file.path(annotation_dir, "cell_type_marker_genes.csv"))
dimUsed <- 20
cols_density <- c(scider:::col.spec, "grey80")
cells_spe <- grid_spe <- grid_seurat <- list()

for (id in samples) {
    p <- sample_info[sample_info$sample_id == id, ]

## ---- Read in and QC --------------------------------------------------
counts <- Read10X_h5(paste0(geo_prefix(id), "_cell_feature_matrix.h5"))[[1]]
cell_info <- read.csv(paste0(geo_prefix(id), "_cells.csv.gz"))
rownames(cell_info) <- cell_info$cell_id
cell_info <- cell_info[, -1]
so <- CreateSeuratObject(counts = counts, project = id, assay = "RNA",
    meta.data = cell_info, min.cells = 0, min.features = 0)
so <- so[rowSums(so[["RNA"]]$counts) >= p$gene_filter_min_total_count, ]

par(mfrow = c(1, 3))
hist(so$nCount_RNA, n = 100, col = "dodgerblue", xlab = "", main = "# transcripts per cell")
hist(so$nFeature_RNA, n = 100, col = "dodgerblue", xlab = "Genes per cell", main = "# genes per cell (Supp Fig 3)")
plot(so$nCount_RNA, so$nFeature_RNA, xlab = "# transcripts", ylab = "# genes", pch = 16, cex = 0.1, log = "xy")
abline(h = p$cell_qc_min_genes, v = p$cell_qc_min_transcripts, lwd = 2, lty = 2, col = "red")

keep.cell <- so$nCount_RNA >= p$cell_qc_min_transcripts & so$nFeature_RNA >= p$cell_qc_min_genes
if (p$cell_qc_single_nucleus_only) keep.cell <- keep.cell & so$nucleus_count == 1
print(table(keep.cell))
so <- so[, keep.cell]

## ---- Clustering ------------------------------------------------------
so <- NormalizeData(so, scale.factor = 100, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = 500, mean.cutoff = c(0.1, 100), verbose = FALSE)
so <- ScaleData(so, features = VariableFeatures(so), verbose = FALSE)
so <- RunPCA(so, dims = 1:dimUsed, features = VariableFeatures(so), verbose = FALSE)
so <- RunUMAP(so, dims = 1:dimUsed, seed.use = p$umap_seed, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:dimUsed, verbose = FALSE)
so <- FindClusters(so, resolution = p$per_sample_cluster_resolution, verbose = FALSE)

genes <- GOI$Gene[GOI$Gene %in% rownames(so)]
print(DimPlot(so, reduction = "umap", raster = FALSE, cols = col.p2) +
    DotPlot(so, features = genes, dot.scale = 6) + coord_flip())

## ---- SingleR cell type annotation ------------------------------------
pred <- SingleR(test = as.SingleCellExperiment(so), ref = sce_ref,
    labels = sce_ref$cluster_Jeremy, de.method = "wilcox", num.threads = 4)
so$SingleR_Jeremy <- pred$labels
print(table(so$SingleR_Jeremy))

print(DimPlot(so, reduction = "umap", group.by = "SingleR_Jeremy", raster = FALSE, cols = cols_ct))
xy <- data.frame(x = so$x_centroid, y = so$y_centroid, cell_type = so$SingleR_Jeremy)
print(ggplot(xy, aes(x, y, colour = cell_type)) + geom_point(size = 0.03) +
    scale_colour_manual(values = cols_ct) + coord_fixed() + theme_void())

    ## ---- scider: cell type densities -------------------------------------
    ## Cell-level SpatialExperiment with the raw counts of all genes; kernel
    ## density of each cell type on an 80-um hexagonal grid.
    spe <- SpatialExperiment(assays = list(counts = counts[, colnames(so)]), sample_id = id,
        colData = so[[]], spatialCoordsNames = c("x_centroid", "y_centroid"))
    spe$cell_type <- so$SingleR_Jeremy
    spe$umap_1 <- Embeddings(so, "umap")[, 1]
    spe$umap_2 <- Embeddings(so, "umap")[, 2]
    spe <- gridDensity(spe, grid.length.x = 80, bandwidth = 80, id = "cell_type")

    ## Hexagonal bin of each cell (scider grid node "<column>-<row>"), assigned
    ## as in scider::gridSPE()
    gi <- metadata(spe)$grid_info
    gd <- metadata(spe)$grid_density
    idx <- hexDensity::xy2hcell(x = spatialCoords(spe)[, 1], y = spatialCoords(spe)[, 2],
        xbins = gi$xbins, xbnds = gi$xlim, ybnds = gi$ylim, shape = gi$shape)
    spe$hex_id <- gd$node[idx]

    ct_names <- sort(unique(spe$cell_type))
    plist <- lapply(ct_names, function(ct) plotDensity(spe, coi = ct, probs = 0.3, cols = cols_density))
    gridExtra::grid.arrange(grobs = plist, ncol = 7)

    ## ---- scider: hexagon neighbourhoods ----------------------------------
    ## Cells are aggregated into 80-um hexagonal bins; bins with few transcripts
    ## are removed and the bins are clustered into neighbourhoods.
    spe_gr <- gridSPE(spe, id = "cell_type", cell.count = TRUE)
    spe_gr$LibSize <- colSums(counts(spe_gr))
    spe_gr$nGenes <- colSums(counts(spe_gr) != 0)
    spe_gr <- spe_gr[, spe_gr$LibSize >= p$grid_qc_min_transcripts & spe_gr$nGenes >= p$grid_qc_min_genes]
    spe_gr$hex_id <- gd$node[match(paste(spatialCoords(spe_gr)[, 1], spatialCoords(spe_gr)[, 2]),
        paste(gd$x_grid, gd$y_grid))]

    so_gr <- CreateSeuratObject(counts(spe_gr), min.cells = 0, min.features = 0)
    so_gr$hex_id <- spe_gr$hex_id
    so_gr <- so_gr[rowSums(so_gr[["RNA"]]$counts) >= p$grid_gene_filter_min_total_count, ]
    so_gr <- NormalizeData(so_gr, scale.factor = 1000, verbose = FALSE)
    so_gr <- FindVariableFeatures(so_gr, nfeatures = 500, verbose = FALSE)
    so_gr <- ScaleData(so_gr, features = rownames(so_gr), verbose = FALSE)
    so_gr <- RunPCA(so_gr, dims = 1:dimUsed, verbose = FALSE)
    so_gr <- RunUMAP(so_gr, dims = 1:dimUsed, seed.use = 2024, verbose = FALSE)
    so_gr <- FindNeighbors(so_gr, dims = 1:dimUsed, verbose = FALSE)
    so_gr <- FindClusters(so_gr, resolution = p$neighbourhood_resolution, verbose = FALSE)
    print(table(Idents(so_gr)))
    spe_gr$niche <- so_gr$seurat_clusters

    print(DimPlot(so_gr, reduction = "umap", raster = FALSE, cols = col.p2, pt.size = 1) +
        plotGrid(spe_gr, group.by = "niche", cols = col.pLight, pol.border = TRUE))

    cells_spe[[id]] <- spe
    grid_spe[[id]] <- spe_gr
    grid_seurat[[id]] <- so_gr
}

