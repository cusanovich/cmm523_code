dyn.load("/opt/ohpc/pub/apps/glpk/5.0/lib/libglpk.so.40")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib64/libproj.so.19")
dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/libs/gnu13/hdf5/1.14.0/lib/libhdf5_hl.so.310")

#options(future.globals.maxSize = 1000 * 1024^2)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)

library(curl)
system("mkdir /xdisk/darrenc/cmm_523/darrenc/wnn_tutorial/")
curl_download("https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5",
              "/xdisk/darrenc/cmm_523/darrenc/wnn_tutorial/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
curl_download("https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
              "/xdisk/darrenc/cmm_523/darrenc/wnn_tutorial/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")
curl_download("https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi",
              "/xdisk/darrenc/cmm_523/darrenc/wnn_tutorial/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi")

# the 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5("/xdisk/darrenc/cmm_523/darrenc/wnn_tutorial/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object
pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/xdisk/darrenc/cmm_523/darrenc/wnn_tutorial/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
pbmc[["ATAC"]] <- chrom_assay

VlnPlot(pbmc, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

# RNA analysis
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

# # perform sub-clustering on cluster 6 to find additional structure
# pbmc <- FindSubCluster(pbmc, cluster = 6, graph.name = "wsnn", algorithm = 3)
# Idents(pbmc) <- "sub.cluster"
# 
# # add annotations
# pbmc <- RenameIdents(pbmc, '19' = 'pDC','20' = 'HSPC','15' = 'cDC')
# pbmc <- RenameIdents(pbmc, '0' = 'CD14 Mono', '9' ='CD14 Mono', '5' = 'CD16 Mono')
# pbmc <- RenameIdents(pbmc, '10' = 'Naive B', '11' = 'Intermediate B', '17' = 'Memory B', '21' = 'Plasma')
# pbmc <- RenameIdents(pbmc, '7' = 'NK')
# pbmc <- RenameIdents(pbmc, '4' = 'CD4 TCM', '13'= "CD4 TEM", '3' = "CD4 TCM", '16' ="Treg", '1' ="CD4 Naive", '14' = "CD4 Naive")
# pbmc <- RenameIdents(pbmc, '2' = 'CD8 Naive', '8'= "CD8 Naive", '12' = 'CD8 TEM_1', '6_0' = 'CD8 TEM_2', '6_1' ='CD8 TEM_2', '6_4' ='CD8 TEM_2')
# pbmc <- RenameIdents(pbmc, '18' = 'MAIT')
# pbmc <- RenameIdents(pbmc, '6_2' ='gdT', '6_3' = 'gdT')
# pbmc$celltype <- Idents(pbmc)

##############Annotations are wrong on Seurat Tutorial###############
library(tidyverse)
library(celldex)
library(SingleR)
library(pheatmap)

ref <- celldex::HumanPrimaryCellAtlasData()
unique(ref$label.main)
unique(ref$label.fine)

ref <- ref[,grepl('DC|B_cell|T_cells|Monocyte|NK_cell|^HSC_CD34',
                  ref$label.main)]
unique(ref$label.main)

norm_counts <- LayerData(pbmc, assay = "SCT", layer = 'data')

ct_ann <- SingleR(test = norm_counts, # we could also use sce or raw_counts
                  ref = ref, 
                  labels = ref$label.fine,
                  de.method = 'wilcox')

plotScoreHeatmap(ct_ann)
plotDeltaDistribution(ct_ann, ncol = 4, dots.on.top = FALSE)

tab <- table(pbmc@meta.data$seurat_clusters, ct_ann$pruned.labels)
tab_mat = as.data.frame(tab)
pheatmap(log10(tab+1), color=viridis::viridis(100), silent=TRUE)
apply(tab,1,function(x) which(x == max(x)))
tab_mat = as.data.frame(tab)
wide_tab <- reshape(tab_mat, idvar = "Var2", timevar = "Var1", direction = "wide")

pbmc <- AddMetaData(pbmc, ct_ann$pruned.labels, col.name = 'SingleR_HCA')
#DimPlot(pbmc,group.by = "SingleR_HCA",label=T,reduction = "wnn.umap")
#DimPlot(pbmc,split.by = "seurat_clusters",label=T,reduction = "wnn.umap")
#DimPlot(pbmc,group.by = "seurat_clusters",label=T,reduction = "wnn.umap", label.size = 6)

Idents(pbmc) <- pbmc@meta.data$seurat_clusters
pbmc <- RenameIdents(pbmc, '18' = 'pDC','21' = 'HSPC','14' = 'cDC')
pbmc <- RenameIdents(pbmc, '0' = 'CD14 Mono', '8' ='CD14 Mono', '6' = 'CD16 Mono')
pbmc <- RenameIdents(pbmc, '9' = 'Naive B', '10' = 'Intermediate B', '16' = 'Memory B', '22' = 'Plasma')
pbmc <- RenameIdents(pbmc, '7' = 'NK')
pbmc <- RenameIdents(pbmc, '4' = 'CD4 TCM', '12'= "CD4 TEM", '5' = "CD4 TCM", '19' ="Treg", '1' ="CD4 Naive", '13' = "CD4 Naive")
pbmc <- RenameIdents(pbmc, '2' = 'CD8 Naive', '15'= "CD8 Naive", '11' = 'CD8 TEM_1', '3' = 'CD8 TEM_2', '20'= "CD8 Naive")
pbmc <- RenameIdents(pbmc, '17' = 'MAIT')
pbmc$celltype <- Idents(pbmc)
#########################################################

p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

## to make the visualization easier, subset T cell clusters
celltype.names <- levels(pbmc)
tcell.names <- grep("CD4|CD8|Treg", celltype.names,value = TRUE)
tcells <- subset(pbmc, idents = tcell.names)
CoveragePlot(tcells, region = 'CD8A', features = 'CD8A', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE)

#pdf("~/ccsgab_projects/wnn_browser_track.pdf")
#CoveragePlot(tcells, region = 'CD8A', features = 'CD8A', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE)
#dev.off()

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
DefaultAssay(pbmc) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(pbmc), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
pbmc <- SetAssayData(pbmc, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

# Note that this step can take 30-60 minutes 
pbmc <- RunChromVAR(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

#returns MA0466.2
motif.name <- ConvertMotifID(pbmc, name = 'CEBPB')
gene_plot <- FeaturePlot(pbmc, features = "sct_CEBPB", reduction = 'wnn.umap')
motif_plot <- FeaturePlot(pbmc, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
gene_plot | motif_plot

markers_rna <- presto:::wilcoxauc.Seurat(X = pbmc, group_by = 'celltype', assay = 'data', seurat_assay = 'SCT')
markers_motifs <- presto:::wilcoxauc.Seurat(X = pbmc, group_by = 'celltype', assay = 'data', seurat_assay = 'chromvar')
motif.names <- markers_motifs$feature
colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
markers_rna$gene <- markers_rna$RNA.feature
markers_motifs$gene <- ConvertMotifID(pbmc, id = motif.names)

# a simple function to implement the procedure above
topTFs <- function(celltype, padj.cutoff = 1e-2) {
  ctmarkers_rna <- dplyr::filter(
    markers_rna, RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0) %>% 
    arrange(-RNA.auc)
  ctmarkers_motif <- dplyr::filter(
    markers_motifs, motif.group == celltype, motif.padj < padj.cutoff, motif.logFC > 0) %>% 
    arrange(-motif.auc)
  top_tfs <- inner_join(
    x = ctmarkers_rna[, c(2, 11, 6, 7)], 
    y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
  )
  top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
  top_tfs <- arrange(top_tfs, -avg_auc)
  return(top_tfs)
}

# identify top markers in NK and visualize
head(topTFs("NK"), 3)

motif.name <- ConvertMotifID(pbmc, name = 'TBX21')
gene_plot <- FeaturePlot(pbmc, features = "sct_TBX21", reduction = 'wnn.umap')
motif_plot <- FeaturePlot(pbmc, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
gene_plot | motif_plot

# identify top markers in pDC and visualize
head(topTFs("pDC"), 3)

motif.name <- ConvertMotifID(pbmc, name = 'TCF4')
gene_plot <- FeaturePlot(pbmc, features = "sct_TCF4", reduction = 'wnn.umap')
motif_plot <- FeaturePlot(pbmc, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
gene_plot | motif_plot

# identify top markers in HSPC and visualize
head(topTFs("CD16 Mono"),3)

motif.name <- ConvertMotifID(pbmc, name = 'SPI1')
gene_plot <- FeaturePlot(pbmc, features = "sct_SPI1", reduction = 'wnn.umap')
motif_plot <- FeaturePlot(pbmc, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
gene_plot | motif_plot

# identify top markers in other cell types
head(topTFs("Naive B"), 3)

head(topTFs("HSPC"), 3)

head(topTFs("Plasma"), 3)
