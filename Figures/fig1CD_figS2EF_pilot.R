source("setup.R")
source("config.R")

### Figure 1C-D: SingleR labels from Bandyopadhyay et al. (2024)
so <- readRDS(file.path(data_dir, "pilot", "RM-1_377gene_panel_seurat.rds"))
p1 <- DimPlot(so, reduction = "umap", group.by = "SingleR_Jeremy", raster = FALSE, cols = cols_ct) +
    ggtitle("RM-1")
p2 <- DimPlot(so, reduction = "spatial", group.by = "SingleR_Jeremy", raster = FALSE,
              cols = cols_ct, pt.size = 0.03) + coord_fixed() + ggtitle("RM-1")
ggsave(file.path(out_dir, "Fig1C_UMAP_RM-1.pdf"), p1, height = 7, width = 8.5)
ggsave(file.path(out_dir, "Fig1D_spatial_RM-1.pdf"), p2, height = 6, width = 10)
save_source(as.data.frame(table(cell_type = so$SingleR_Jeremy)), "Fig1CD_RM-1_cell_type_counts")

### Supplemental Figure 2E-F: SingleR labels from Baccin et al. (2020).
so <- readRDS(file.path(data_dir, "pilot", "mouse_BM_379gene_panel_seurat.rds"))
p1 <- DimPlot(so, reduction = "umap", group.by = "SingleR_labels_2", raster = FALSE) + ggtitle("Mouse BM")
p2 <- DimPlot(so, reduction = "spatial", group.by = "SingleR_labels_2", raster = FALSE, pt.size = 0.03) +
    coord_fixed() + ggtitle("Mouse BM")
ggsave(file.path(out_dir, "FigS2E_UMAP_mouse.pdf"), p1, height = 7.5, width = 11)
ggsave(file.path(out_dir, "FigS2F_spatial_mouse.pdf"), p2, height = 8, width = 10)
save_source(as.data.frame(table(cell_type = so$SingleR_labels_2)), "FigS2EF_mouse_cell_type_counts")

