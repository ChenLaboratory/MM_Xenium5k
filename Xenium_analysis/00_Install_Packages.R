## 00_Install_Packages.R

# CRAN
install.packages(c(
    "Seurat",
    "hdf5r",
    "ggplot2",
    "gridExtra",
    "Matrix"
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "SpatialExperiment",
    "scider",
    "SingleR",
    "edgeR",
    "limma"
))
