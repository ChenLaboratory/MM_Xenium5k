source("setup.R")
source("config.R")

so <- readRDS(file.path(data_dir, "integration", "plasma_cells_integrated_seurat.rds"))
Idents(so) <- "PC_cluster_paper"
so$sample <- factor(so$orig.ident, levels = samples)

### Figure 3A
p <- DimPlot(so, reduction = "umap", raster = FALSE, cols = pc_colours)
ggsave(file.path(out_dir, "Fig3A_UMAP_PC_clusters.pdf"), p, height = 7.5, width = 9)

### Figure 3B: cells of each disease group highlighted in dark blue
umap <- Embeddings(so, "umap")
pdf(file.path(out_dir, "Fig3B_UMAP_by_group.pdf"), height = 6, width = 18)
par(mfrow = c(1, 3))
for (g in c("CTRL", "MGUS", "MM")) {
    ord <- order(so$group == g)
    plot(umap[ord, ], pch = 16, cex = 0.4, main = g,
         col = ifelse(so$group == g, "darkblue", "grey")[ord])
}
dev.off()

### Figure 3C and Supplemental Table 3: plasma cell clusters per sample
tab <- table(so$PC_cluster_paper, so$sample)
save_source(as.data.frame.matrix(tab), "TableS3_plasma_cluster_counts_per_sample", row.names = TRUE)
prop <- as.data.frame(prop.table(tab, 2) * 100)
colnames(prop) <- c("cluster", "sample", "percent")
prop$sample <- factor(prop$sample, levels = rev(samples))
save_source(prop, "Fig3C_plasma_cluster_proportions")
p <- ggplot(prop, aes(x = percent, y = sample, fill = cluster)) +
    geom_col(position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = pc_colours) + xlab("Cell proportion (%)") + ylab(NULL) + theme_classic()
ggsave(file.path(out_dir, "Fig3C_plasma_cluster_proportions.pdf"), p, height = 7, width = 7)

### Figure 3D: top 10 markers of each cluster (rows in Seurat cluster order, as in the paper)
markers <- read.csv(file.path(data_dir, "tables", "plasma_cluster_markers.csv"))
top <- do.call(rbind, lapply(split(markers, markers$cluster), head, n = 10))
genes <- unique(top$gene)
p <- DotPlot(so, features = genes, dot.scale = 8) + coord_flip()
save_source(p$data, "Fig3D_dotplot_data")
ggsave(file.path(out_dir, "Fig3D_dotplot_markers.pdf"), p, height = 30, width = 8)

### Figure 3E and Supplemental Figure 6: UMAP and spatial plot of each sample
pdf(file.path(out_dir, "Fig3E_FigS6_UMAP_spatial_per_sample.pdf"), height = 6, width = 12)
for (id in samples) {
    sel <- so$sample == id
    df <- data.frame(umap, cluster = so$PC_cluster_paper, sel = sel)
    p1 <- ggplot(df, aes(umap_1, umap_2)) + geom_point(colour = "grey92", size = 0.05) +
        geom_point(data = df[sel, ], aes(colour = cluster), size = 0.05) +
        scale_colour_manual(values = pc_colours) + ggtitle(id) + theme_classic() +
        theme(legend.position = "none")
    df <- data.frame(x = so$x_centroid, y = so$y_centroid, cluster = so$PC_cluster_paper)[sel, ]
    p2 <- ggplot(df, aes(x, y, colour = cluster)) + geom_point(size = 0.3, shape = 16) +
        scale_colour_manual(values = pc_colours) + coord_fixed() + ggtitle(id) + theme_spatial
    gridExtra::grid.arrange(p1, p2, ncol = 2)
}
dev.off()

### Figure 4B: FISH-associated genes in plasma cells, samples ordered by cytogenetics
fish_order <- c("Ctrl-1", "Ctrl-2", "Ctrl-3", "Ctrl-4", "MM-1", "MM-9", "MM-2", "MGUS-1",
    "SM-1", "MM-3", "MM-4", "MM-6", "MM-7", "MM-10", "SM-2", "SM-3", "SM-5", "MM-8",
    "MM-5", "SM-4", "MGUS-2")
fish_genes <- c("CD19", "CD79A", "CD27", "CCND1", "CCND2", "FGFR3", "MAF", "ITGB7", "MYC", "MAFB", "FRZB")
so$sample_fish <- factor(so$orig.ident, levels = fish_order)
p <- DotPlot(so, features = fish_genes, group.by = "sample_fish", dot.scale = 8) + coord_flip() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
save_source(p$data, "Fig4B_dotplot_data")
ggsave(file.path(out_dir, "Fig4B_dotplot_FISH_genes.pdf"), p, height = 5, width = 10)

### Supplemental Figure 8: MKI67 in plasma cells
mki67_plot <- function(id) {
    FeaturePlot(subset(so, sample == id), reduction = "spatial", features = "MKI67",
                raster = FALSE, pt.size = 0.8, order = TRUE) + coord_fixed() + ggtitle(id)
}
pdf(file.path(out_dir, "FigS8_MKI67_spatial.pdf"), height = 7, width = 6)
for (id in samples) print(mki67_plot(id))
dev.off()
pb <- rbind(
    cbind(comparison = "MGUS/SM vs Ctrl", subset(read.csv(file.path(data_dir, "tables", "pseudobulk_DEG-PC-MGUS_vs_Ctrl.csv")), Symbol == "MKI67")),
    cbind(comparison = "MM vs MGUS/SM", subset(read.csv(file.path(data_dir, "tables", "pseudobulk_DEG-PC-MM_vs_MGUS.csv")), Symbol == "MKI67")))
save_source(pb, "FigS8_MKI67_pseudobulk_tests")

### Supplemental Figure 9: spatially distinct plasma cell clusters and MKI67
sel_clusters <- list("MM-3" = c("PC2", "PC5", "PC14"), "MM-4" = c("PC3", "PC8"),
    "MM-6" = c("PC9", "PC11"), "MM-7" = c("PC4", "PC12", "PC13"), "MM-8" = c("PC0", "PC15"))
pdf(file.path(out_dir, "FigS9_spatial_PC_clusters_MKI67.pdf"), height = 8, width = 16)
for (id in names(sel_clusters)) {
    df <- data.frame(x = so$x_centroid, y = so$y_centroid, cluster = so$PC_cluster_paper)[so$sample == id, ]
    df <- df[df$cluster %in% sel_clusters[[id]], ]
    plist <- list(ggplot(df, aes(x, y, colour = cluster)) + geom_point(size = 0.4, shape = 16) +
        scale_colour_manual(values = pc_colours) + coord_fixed() + ggtitle(paste(id, "composite")) + theme_spatial)
    for (k in sel_clusters[[id]])
        plist[[k]] <- ggplot(df[df$cluster == k, ], aes(x, y)) +
            geom_point(size = 0.4, shape = 16, colour = pc_colours[k]) + coord_fixed() + ggtitle(k) + theme_spatial
    plist$MKI67 <- mki67_plot(id)
    gridExtra::grid.arrange(grobs = plist, nrow = 1)
}
dev.off()
fm <- list.files(file.path(data_dir, "tables", "plasma_cluster_FindMarkers"), full.names = TRUE)
fm <- do.call(rbind, lapply(fm, function(f) {
    x <- read.csv(f, row.names = 1)
    cbind(comparison = sub(".csv", "", basename(f)), gene = "MKI67", x["MKI67", ])
}))
save_source(fm, "FigS9_MKI67_FindMarkers")

