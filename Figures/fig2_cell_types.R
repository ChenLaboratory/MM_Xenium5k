source("setup.R")
source("config.R")

cells <- read_per_sample("cells_spe")[samples]

### Figure 2A: UMAP of each sample coloured by cell type
umap_plot <- function(spe, id) {
    df <- data.frame(UMAP_1 = spe$umap_1, UMAP_2 = spe$umap_2, cell_type = spe$cell_type)
    ggplot(df, aes(UMAP_1, UMAP_2, colour = cell_type)) +
        geom_point(size = 0.05) + scale_colour_manual(values = cols_ct) +
        ggtitle(id) + theme_classic() + theme(legend.position = "none")
}
pdf(file.path(out_dir, "Fig2A_UMAP_cell_types.pdf"), height = 16, width = 20)
gridExtra::grid.arrange(grobs = Map(umap_plot, cells, samples), ncol = 5)
dev.off()

### Figure 2B: cell type proportions
ct <- sapply(cells, function(spe) table(factor(spe$cell_type, levels = names(cols_ct))))
save_source(ct, "Fig2B_cell_type_counts", row.names = TRUE)
prop <- data.frame(sample = factor(rep(samples, each = nrow(ct)), levels = rev(samples)),
                   cell_type = factor(rep(rownames(ct), length(samples)), levels = names(cols_ct)),
                   percent = as.vector(100 * t(t(ct) / colSums(ct))))
save_source(prop, "Fig2B_cell_type_proportions")
p <- ggplot(prop, aes(x = percent, y = sample, fill = cell_type)) +
    geom_col(position = position_stack(reverse = TRUE), colour = "black", linewidth = 0.1) +
    scale_fill_manual(values = cols_ct) + xlab("Cell proportion (%)") + ylab(NULL) + theme_classic()
ggsave(file.path(out_dir, "Fig2B_cell_type_proportions.pdf"), p, height = 7, width = 8)

### Quasi-binomial GLMs (Results, paragraph on Figure 2B)
group <- factor(sample_info$analysis_group, levels = c("MGUS", "CTRL", "MM"))
cp <- t(t(ct) / colSums(ct))
res <- sapply(c("Plasma Cell", "Red Cell", "Megakaryocyte"), function(x) {
    fit <- glm(cp[x, ] ~ group, family = quasibinomial(), weights = colSums(ct))
    coef(summary(fit))[c("groupCTRL", "groupMM"), 4]
})
rownames(res) <- c("P_CTRL_vs_MGUS/SM", "P_MM_vs_MGUS/SM")
avg <- sapply(split(seq_along(group), group), function(j) 100 * rowMeans(cp[, j, drop = FALSE]))
sink(file.path(out_dir, "Fig2B_cell_type_proportion_tests.txt"))
print(signif(res, 3))
print(round(avg, 2))
sink()
save_source(signif(res, 3), "Fig2B_quasibinomial_P_values", row.names = TRUE)

### Figure 2C (SM-1, MM-3) and Supplemental Figure 5 (all samples): spatial maps
spatial_plot <- function(spe, id, size = 0.03) {
    df <- data.frame(spatialCoords(spe), cell_type = spe$cell_type)
    ggplot(df, aes(x_centroid, y_centroid, colour = cell_type)) +
        geom_point(size = size, shape = 16) + scale_colour_manual(values = cols_ct) +
        coord_fixed() + ggtitle(id) + theme_spatial + theme(legend.position = "none")
}
pdf(file.path(out_dir, "Fig2C_spatial_cell_types_SM-1_MM-3.pdf"), height = 8, width = 12)
gridExtra::grid.arrange(spatial_plot(cells[["SM-1"]], "SM-1"), spatial_plot(cells[["MM-3"]], "MM-3"), ncol = 2)
dev.off()
pdf(file.path(out_dir, "FigS5_spatial_cell_types.pdf"), height = 10, width = 14)
for (id in samples) print(spatial_plot(cells[[id]], id) + theme(legend.position = "right") +
    guides(colour = guide_legend(override.aes = list(size = 3))))
dev.off()

