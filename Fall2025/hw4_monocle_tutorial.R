dyn.load("/opt/ohpc/pub/apps/glpk/5.0/lib/libglpk.so.40")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib64/libproj.so.19")
dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/libs/gnu13/hdf5/1.14.0/lib/libhdf5_hl.so.310")
library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)


# #InstallData("hcabm40k")
# data("hcabm40k")
# hcabm40k = UpdateSeuratObject(object = hcabm40k)
# hcabm40k <- SplitObject(hcabm40k, split.by = "orig.ident")
# for (i in seq_along(hcabm40k)) {
#   hcabm40k[[i]] <- NormalizeData(hcabm40k[[i]]) %>% FindVariableFeatures()
# }
# features <- SelectIntegrationFeatures(hcabm40k)
# for (i in seq_along(along.with = hcabm40k)) {
#   hcabm40k[[i]] <- ScaleData(hcabm40k[[i]], features = features) %>% RunPCA(features = features)
# }
# anchors <- FindIntegrationAnchors(hcabm40k, reference = c(1, 2), reduction = "rpca", dims = 1:30)
# integrated <- IntegrateData(anchors, dims = 1:30)
# integrated <- ScaleData(integrated)
# integrated <- RunPCA(integrated)
# integrated <- RunUMAP(integrated, dims = 1:30, reduction.name = "UMAP")
# integrated <- FindNeighbors(integrated, dims = 1:30)
# integrated <- FindClusters(integrated)

con <- url("https://seurat.nygenome.org/monocle3/hcabm40k_integrated.Rds")
integrated <- readRDS(file = con)
close(con = con)
integrated = UpdateSeuratObject(integrated)
DimPlot(integrated, group.by = c("orig.ident", "ident"))

pdf("~/ccsgab_projects/cmm523_code/Fall2025/monocle_umap.pdf")
DimPlot(integrated, group.by = "orig.ident") + NoLegend()
dev.off()

cds <- as.cell_data_set(integrated)
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

integrated.sub <- subset(as.Seurat(cds,assay=NULL), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

max.avp <- which.max(unlist(FetchData(integrated.sub, "AVP")))
max.avp <- colnames(integrated.sub)[max.avp]
cds <- order_cells(cds, root_cells = max.avp)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)

# Set the assay back as 'integrated'
integrated.sub <- as.Seurat(cds, assay = NULL)
FeaturePlot(integrated.sub, "monocle3_pseudotime")
