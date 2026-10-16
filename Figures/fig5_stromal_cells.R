source("setup.R")
source("config.R")

so <- readRDS(file.path(data_dir, "integration", "stromal_cells_integrated_seurat.rds"))
Idents(so) <- "stromal_cell_type"
so$sample <- factor(so$orig.ident, levels = samples)
so$group <- factor(sample_info$analysis_group[match(so$orig.ident, sample_info$sample_id)],
                   levels = c("CTRL", "MGUS", "MM"))

### Figure 5A
p <- DimPlot(so, reduction = "umap", raster = FALSE, cols = stroma_colours)
ggsave(file.path(out_dir, "Fig5A_UMAP_stromal_clusters.pdf"), p, height = 7.5, width = 9)

### Figure 5B: top 10 markers of each cluster (rows in Seurat cluster order, as in the paper)
markers <- read.csv(file.path(data_dir, "tables", "stromal_cluster_markers.csv"))
top <- do.call(rbind, lapply(split(markers, markers$cluster), head, n = 10))
genes <- unique(top$gene)
p <- DotPlot(so, features = genes, dot.scale = 8) + coord_flip() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
save_source(p$data, "Fig5B_dotplot_data")
ggsave(file.path(out_dir, "Fig5B_dotplot_markers.pdf"), p, height = 18, width = 8)

### Figure 5C: cells of each disease group highlighted in dark blue
umap <- Embeddings(so, "umap")
pdf(file.path(out_dir, "Fig5C_UMAP_by_group.pdf"), height = 6, width = 18)
par(mfrow = c(1, 3))
for (g in c("CTRL", "MGUS", "MM")) {
    ord <- order(so$group == g)
    plot(umap[ord, ], pch = 16, cex = 0.4, main = g,
         col = ifelse(so$group == g, "darkblue", "grey")[ord])
}
dev.off()

### Figure 5D: stromal clusters per sample
tab <- table(so$stromal_cell_type, so$sample)
save_source(as.data.frame.matrix(tab), "Fig5D_stromal_cluster_counts_per_sample", row.names = TRUE)
prop <- as.data.frame(prop.table(tab, 2) * 100)
colnames(prop) <- c("cluster", "sample", "percent")
prop$sample <- factor(prop$sample, levels = rev(samples))
save_source(prop, "Fig5D_stromal_cluster_proportions")
p <- ggplot(prop, aes(x = percent, y = sample, fill = cluster)) +
    geom_col(position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = stroma_colours) + xlab("Cell proportion (%)") + ylab(NULL) + theme_classic()
ggsave(file.path(out_dir, "Fig5D_stromal_cluster_proportions.pdf"), p, height = 7, width = 7)

### Figure 5E and Supplemental Figure 10: UMAP and spatial plot of each sample
pdf(file.path(out_dir, "Fig5E_FigS10_UMAP_spatial_per_sample.pdf"), height = 6, width = 12)
for (id in samples) {
    sel <- so$sample == id
    df <- data.frame(umap, cluster = so$stromal_cell_type, sel = sel)
    p1 <- ggplot(df, aes(umap_1, umap_2)) + geom_point(colour = "grey92", size = 0.05) +
        geom_point(data = df[sel, ], aes(colour = cluster), size = 0.1) +
        scale_colour_manual(values = stroma_colours) + ggtitle(id) + theme_classic() +
        theme(legend.position = "none")
    df <- data.frame(x = so$x_centroid, y = so$y_centroid, cluster = so$stromal_cell_type)[sel, ]
    p2 <- ggplot(df, aes(x, y, colour = cluster)) + geom_point(size = 0.6, shape = 16) +
        scale_colour_manual(values = stroma_colours) + coord_fixed() + ggtitle(id) + theme_spatial
    gridExtra::grid.arrange(p1, p2, ncol = 2)
}
dev.off()

