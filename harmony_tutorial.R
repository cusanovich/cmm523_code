dyn.load("/opt/ohpc/pub/apps/glpk/5.0/lib/libglpk.so.40")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib/libproj.so.19")
dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")
dyn.load("/opt/ohpc/pub/apps/libpng/1.6.37/lib/libpng16.so.16")
dyn.load("/usr/lib64/atlas/libsatlas.so.3")
library(Seurat)
library(SeuratDisk)
library(harmony)
library(dplyr)
library(cowplot)
library(curl)

#system("mkdir /xdisk/darrenc/cmm_523/darrenc/harmony_tutorial")
curl_download("https://zenodo.org/records/8164711/files/pbmc_stim.RData?download=1","/xdisk/darrenc/cmm_523/darrenc/harmony_tutorial/pbmc_stim.RData")

## Source required data
data("pbmc_stim")
pbmc <- CreateSeuratObject(counts = cbind(pbmc.stim, pbmc.ctrl), project = "PBMC", min.cells = 5)

## Separate conditions

pbmc@meta.data$stim <- c(rep("STIM", ncol(pbmc.stim)), rep("CTRL", ncol(pbmc.ctrl)))

pbmc <- pbmc %>%
  NormalizeData(verbose = FALSE)

VariableFeatures(pbmc) <- split(row.names(pbmc@meta.data), pbmc@meta.data$stim) %>% lapply(function(cells_use) {
  pbmc[,cells_use] %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    VariableFeatures()
}) %>% unlist %>% unique

pbmc <- pbmc %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(features = VariableFeatures(pbmc), npcs = 20, verbose = FALSE)

pbmc <- pbmc %>% 
  RunHarmony("stim", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)

harmony.embeddings <- Embeddings(pbmc, reduction = "harmony")

p1 <- DimPlot(object = pbmc, reduction = "harmony", pt.size = .1, group.by = "stim")
p2 <- VlnPlot(object = pbmc, features = "harmony_1", group.by = "stim",  pt.size = .1)
plot_grid(p1,p2)

DimHeatmap(object = pbmc, reduction = "harmony", cells = 500, dims = 1:3)

pbmc <- pbmc %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.5) 

pbmc <- pbmc %>%
  RunTSNE(reduction = "harmony")


p1 <- DimPlot(pbmc, reduction = "tsne", group.by = "stim", pt.size = .1)
p2 <- DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = .1)
plot_grid(p1, p2)

FeaturePlot(object = pbmc, features= c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", "CCL2", "PPBP"), 
            min.cutoff = "q9", cols = c("lightgrey", "blue"), pt.size = 0.5)

pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony",  dims = 1:20)

p1 <- DimPlot(pbmc, reduction = "umap", group.by = "stim", pt.size = .1)
p2 <- DimPlot(pbmc, reduction = "umap", label = TRUE,  pt.size = .1)
plot_grid(p1, p2)
