## config.R

if (!exists("annotation_dir")) source("setup.R")

## ---- Sample sheet ----------------------------------------------------
sample_info <- read.delim(file.path(annotation_dir, "sample_info.tsv"))
sample_info <- sample_info[sample_info$panel == "Xenium Prime 5K", ]
samples <- c(paste0("Ctrl-", 1:4), paste0("MGUS-", 1:2), paste0("SM-", 1:5), paste0("MM-", 1:10))
sample_info <- sample_info[match(samples, sample_info$sample_id), ]

## ---- Helpers ---------------------------------------------------------
## Per-sample objects: named list of the 21 samples, one file per object type
## ("cells_spe", "grid_spe" or "grid_seurat")
read_per_sample <- function(type) readRDS(file.path(data_dir, paste0("per_sample_", type, ".rds")))

## Write the numbers plotted in a panel to out_dir/source_data/<name>.csv
save_source <- function(x, name, row.names = FALSE)
    write.csv(x, file.path(out_dir, "source_data", paste0(name, ".csv")), row.names = row.names)
