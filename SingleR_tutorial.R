dyn.load("/opt/ohpc/pub/apps/glpk/5.0/lib/libglpk.so.40")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib/libproj.so.19")
dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")
dyn.load("/opt/ohpc/pub/apps/libpng/1.6.37/lib/libpng16.so.16")
dyn.load("/usr/lib64/atlas/libsatlas.so.3")

#myPaths = .libPaths()
#myPaths[1] = "/home/u12/darrenc/R/library_R_v4.2_test"
#.libPaths(myPaths)
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(celldex)
library(SingleR)
library(Seurat)
#source("~/ccsgab_projects/cmm523_code/singleR_patch_a.R")

# We'll use a PBMC dataset from the R package scRNAseq
sce <- scRNAseq::KotliarovPBMCData(mode = c('rna'))
seu <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))
rm(sce)
seu <- NormalizeData(object = seu)
# Additional Seurat preprocessing steps - 
# This is not necessary for this tutorial but I will use it for visualisation later
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:10) 

#saveRDS(seu,"/xdisk/darrenc/darrenc/singler_seurat_obj.rds")
seu = readRDS("/xdisk/darrenc/darrenc/singler_seurat_obj.rds")

raw_counts <- LayerData(seu, assay = "RNA", layer = 'counts') #
raw_counts[c('VIM', 'BCL2', 'TP53', 'CD4'),1:5]
norm_counts <- LayerData(seu, assay = "RNA", layer = 'data') #
norm_counts[c('VIM', 'BCL2', 'TP53', 'CD4'),1:5]

# 2. Get reference dataset
ref <- celldex::HumanPrimaryCellAtlasData()
unique(ref$label.main)
unique(ref$label.fine)

# Subset to include only relevant cell types (CAREFUL!)
ref <- ref[,grepl('DC|B_cell|Neutrophils|T_cells|Monocyte|Erythroblast|
                 Macrophage|NK_cell|Platelets|Myelocyte', ref$label.main)]
unique(ref$label.main)

# # 3. Run SingleR
# 
# ct_ann <- SingleR(sc_data = norm_counts, # we could also use sce or raw_counts
#                   ref_data = ref@assays@data@listData$logcounts, 
#                   types = ref$label.main)
# 
# str(ct_ann)

# 3. Run SingleR
ct_ann <- SingleR(test = norm_counts, # we could also use sce or raw_counts
                  ref = ref, 
                  labels = ref$label.main,
                  de.method = 'wilcox')

plotScoreHeatmap(ct_ann)
plotDeltaDistribution(ct_ann, ncol = 4, dots.on.top = FALSE)

# Add to seurat object
rownames(ct_ann)[1:5] # make sure you have cell IDs
seu <- AddMetaData(seu, ct_ann$pruned.labels, col.name = 'SingleR_HCA')
# Visualise them on the UMAP
seu <- SetIdent(seu, value = "SingleR_HCA")
DimPlot(seu, label = T , repel = T, label.size = 3) + NoLegend()

pdf("~/ccsgab_projects/singler_umap.pdf")
DimPlot(seu) + NoLegend()
dev.off()
