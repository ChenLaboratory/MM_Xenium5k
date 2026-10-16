source("setup.R")
source("config.R")

col.manual <- list(
    "MM-4" = c(1, 3, 5, 7), "MGUS-2" = c(1, 2), "MM-1" = c(3, 2, 5), "MM-2" = c(3, 1),
    "MM-3" = c(3, 4, 5), "Ctrl-1" = c(1, 6, 2), "MGUS-1" = c(1, 2, 6), "SM-1" = c(1, 2),
    "MM-10" = c(3, 6, 4, 2), "SM-5" = c(6, 1, 2), "MM-9" = c(3, 6, 1, 2), "Ctrl-4" = c(1, 2),
    "Ctrl-3" = c(1, 6, 2), "MM-8" = c(1, 2, 3), "MM-7" = c(1, 3, 5), "Ctrl-2" = c(1, 4, 2),
    "SM-4" = c(1, 2), "MM-6" = c(3, 5, 1, 6), "MM-5" = c(1, 3), "SM-3" = c(1), "SM-2" = c(1, 4))

ct_names <- names(cols_ct)
ct_cols <- paste0("cell_counts_", gsub("[^a-z0-9]+", "_", tolower(ct_names)))
reord <- c(setdiff(ct_names, "Plasma Cell"), "Plasma Cell")

comp_plot <- function(ct_count, id, title) {
    perc <- 100 * ct_count / rowSums(ct_count)
    dat <- data.frame(Neighbourhood = factor(rep(rownames(perc), ncol(perc)), levels = rownames(perc)),
                      Cell_type = factor(rep(colnames(perc), each = nrow(perc)), levels = reord),
                      Percent = as.vector(perc))
    p <- ggplot(dat, aes(fill = Cell_type, x = Percent, y = Neighbourhood)) +
        geom_col(position = position_stack(reverse = TRUE), colour = "black", linewidth = 0.1) +
        scale_fill_manual(values = cols_ct[reord]) + ggtitle(paste(id, title)) +
        xlab("Cell proportion (%)") + theme_classic()
    list(plot = p, data = cbind(sample = id, dat))
}

src_with <- src_without <- list()
n_cells_grid <- list()
grid_so <- read_per_sample("grid_seurat")
grid_spe <- read_per_sample("grid_spe")
pdf(file.path(out_dir, "Fig6_FigS11_FigS12_neighbourhoods.pdf"), height = 12, width = 13)
for (id in samples) {
    so_gr <- grid_so[[id]]
    spe_gr <- grid_spe[[id]]
    cols <- col.manual[[id]]
    p1 <- DimPlot(so_gr, reduction = "umap", raster = FALSE, cols = col.p2[cols], pt.size = 1) + ggtitle(id)
    p2 <- plotGrid(spe_gr, group.by = "niche", cols = col.pLight[cols], pol.border = TRUE)

    ct_count <- rowsum(as.matrix(colData(spe_gr)[, ct_cols]), as.character(spe_gr$niche))
    colnames(ct_count) <- ct_names
    w <- comp_plot(ct_count, id, "- with plasma cells (Fig 6B, Supp Fig 12A)")
    wo <- comp_plot(ct_count[, colnames(ct_count) != "Plasma Cell", drop = FALSE], id,
                    "- without plasma cells (Supp Fig 12B)")
    n_cells_grid[[id]] <- spe_gr$cell_counts_overall
    src_with[[id]] <- w$data
    src_without[[id]] <- wo$data
    gridExtra::grid.arrange(grobs = list(p1, p2, w$plot, wo$plot),
        layout_matrix = rbind(c(1, 2), c(3, 3), c(4, 4)), heights = c(3, 1, 1))
}
dev.off()
rm(grid_so, grid_spe)
save_source(do.call(rbind, src_with), "Fig6B_FigS12A_neighbourhood_composition_with_PC")
save_source(do.call(rbind, src_without), "FigS12B_neighbourhood_composition_without_PC")

### Figure 6A (plasma cells, megakaryocytes, macrophages, MSCs) and
### Supplemental Figure 13 (other cell types): kernel density of each cell type
cols_density <- c(scider:::col.spec, "grey80")
cells <- read_per_sample("cells_spe")
pdf(file.path(out_dir, "Fig6A_FigS13_cell_type_density.pdf"), height = 12, width = 22)
for (id in samples) {
    spe <- cells[[id]]
    plist <- lapply(ct_names, function(ct)
        plotDensity(spe, coi = ct, probs = 0.3, cols = cols_density) + ggtitle(ct))
    gridExtra::grid.arrange(grobs = plist, ncol = 7, top = id)
}
dev.off()

sessionInfo()
