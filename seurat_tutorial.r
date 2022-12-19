library(dplyr)
library(Seurat)
library(patchwork)

# Set working directory
setwd('/gpfs/home/acs9950/singlecell/seurat_tutorial/')

# Load the dataset
sc_data <- Read10X(data.dir = "raw/")