## config.R

if (!exists("annotation_dir")) source("setup.R")

## ---- Directories -----------------------------------------------------
GEO_DIR      <- "GSE299193_RAW"
REF_DIR      <- "Reference"
MOUSE_PREFIX <- "Mouse/Region_4"
OUTPUT_DIR   <- "output"

for (d in c("integration", "pilot", "tables")) {
    dir.create(file.path(OUTPUT_DIR, d), showWarnings = FALSE, recursive = TRUE)
}

## ---- Sample sheet ----------------------------------------------------
sample_info <- read.delim(file.path(annotation_dir, "sample_info.tsv"))
samples <- sample_info$sample_id[sample_info$panel == "Xenium Prime 5K"]

geo_prefix <- function(id) {
    gsm <- sample_info$GSM[sample_info$sample_id == id]
    file.path(GEO_DIR, paste0(gsm, "_", sub("^RM-1$", "RM", id)))
}
