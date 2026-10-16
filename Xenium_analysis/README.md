# Xenium Analysis Pipeline

Analysis of the 21 Xenium Prime 5K bone marrow trephines (Ctrl-1–4, MGUS-1–2, SM-1–5, MM-1–10) and the two pilot data sets (RM-1, mouse bone). The outputs have the same layout as the Zenodo record, so the scripts in [`../Figures`](../Figures) can be run on them.

## Workflow

Scripts, run in order:

| Script | Description | Paper |
|---|---|---|
| `00_Install_Packages.R` | R package installation | |
| `01_Per_Sample_Analysis.R` | Per-sample QC, clustering, SingleR cell type annotation, scider cell type densities and hexagon neighbourhoods | Fig 2, 6; Supp Figs 3, 5, 11–13 |
| `02_Pseudobulk_DE.R` | Pseudo-bulk edgeR DE of plasma cells between groups, KEGG analysis | Supp Fig 8 |
| `03_Integration.R` | Combined clustering of the plasma cells and of the stromal cells | Fig 3, 4B, 5; Supp Figs 6, 10; Supp Table 3 |
| `04_Plasma_Cluster_Markers.R` | Markers between plasma cell clusters occupying distinct regions of the same sample | Supp Fig 9 |
| `05_Pilot_Analysis.R` | RM-1 and mouse bone pilot data sets | Fig 1C–D; Supp Fig 2E–F |

- `setup.R` loads the packages and defines the annotation directory and colour palettes;
- `config.R` defines paths, the sample sheet and helper functions.


## Running the pipeline


```r
setwd("path/to/this/repository/Xenium_analysis")

source("00_Install_Packages.R")
source("01_Per_Sample_Analysis.R")
source("02_Pseudobulk_DE.R")
source("03_Integration.R")
source("04_Plasma_Cluster_Markers.R")
source("05_Pilot_Analysis.R")
```
