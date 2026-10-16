suppressPackageStartupMessages({
    library(Seurat)
    library(SpatialExperiment)
    library(scider)
    library(ggplot2)
})

annotation_dir <- "../annotation"  # path to annotation files

## ---- Colour palettes -------------------------------------------------
col.pMedium <- c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D", "#AD8BC9",
    "#A8786E", "#ED97CA", "#A2A2A2", "#CDCC5D", "#6DCCDA")
col.pDark <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")
col.pLight <- c("#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
    "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5")
col.p2 <- rep(c(col.pDark, col.pLight, col.pMedium), 2)
col.group <- c(CTRL = "darkblue", MGUS = "gold3", MM = "salmon")

## Cell type colours (SingleR labels transferred from Bandyopadhyay et al. 2024)
cols_ct <- readRDS(file.path(annotation_dir, "cell_type_colours.rds"))

## ---- Cluster labels --------------------------------------------------
pc_levels <- paste0("PC", 0:16)
pc_colours <- setNames(col.p2[c(10, 3, 2, 4:9, 1, 11:17)], pc_levels)

## Manuscript names of the integrated stromal cell clusters
stroma_levels <- c("SEC1", "SEC2", "AEC", "VSMC", "MSC", "Osteo-MSC",
    "Osteoblast1", "Osteoblast2", "Osteoclast", "Cartilage", "Adipocyte")
stroma_colours <- setNames(col.p2[c(3, 2, 7, 5, 1, 10, 4, 8, 11, 9, 6)], stroma_levels)

theme_spatial <- theme_void() + theme(plot.title = element_text(hjust = 0.5))
