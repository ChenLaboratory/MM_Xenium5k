## 04_Plasma_Cluster_Markers.R

source("setup.R")
source("config.R")

so <- readRDS(file.path(OUTPUT_DIR, "integration", "plasma_cells_integrated_seurat.rds"))
Idents(so) <- "PC_cluster_paper"
comparisons <- list("MM-3" = c("PC2", "PC5", "PC14"), "MM-4" = c("PC3", "PC8"),
    "MM-6" = c("PC9", "PC11"), "MM-7" = c("PC4", "PC12", "PC13"), "MM-8" = c("PC0", "PC15"))
dir.create(file.path(OUTPUT_DIR, "tables", "plasma_cluster_FindMarkers"), showWarnings = FALSE)
mki67 <- list()
for (id in names(comparisons)) {
    so_id <- subset(so, orig.ident == id)
    for (pair in combn(comparisons[[id]], 2, simplify = FALSE)) {
        name <- paste0(id, "_", pair[1], "_vs_", pair[2])
        markers <- FindMarkers(so_id, ident.1 = pair[1], ident.2 = pair[2])
        write.csv(markers, file.path(OUTPUT_DIR, "tables", "plasma_cluster_FindMarkers", paste0(name, ".csv")))
        mki67[[name]] <- markers["MKI67", c("avg_log2FC", "p_val_adj")]
    }
}
print(do.call(rbind, mki67))
