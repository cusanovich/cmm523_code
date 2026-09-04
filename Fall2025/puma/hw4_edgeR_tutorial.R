dyn.load("/opt/ohpc/pub/apps/glpk/5.0/lib/libglpk.so.40")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib64/libproj.so.19")
dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/libs/gnu13/hdf5/1.14.0/lib/libhdf5_hl.so.310")
library(edgeR)
library(Seurat)
source("~/ccsgab_projects/cmm523_code/new_edgeR_patch.R")

download.file("https://bioinf.wehi.edu.au/edgeR/UserGuideData/SeuratObj.rds", "~/ccsgab_projects/SeuratObj.rds")
so <- readRDS("~/ccsgab_projects/SeuratObj.rds")
head(so@meta.data)

p1 <- Seurat::DimPlot(so, reduction="tsne", cols=2:8)
p2 <- Seurat::DimPlot(so, reduction="tsne", group.by="group")
p1 | p2

y <- Seurat2PB(so, sample="group", cluster="seurat_clusters")
dim(y)
head(y$samples, n=10L)

summary(y$samples$lib.size)
keep.samples<-y$samples$lib.size>5e4
table(keep.samples)
y<-y[, keep.samples]

keep.genes<-filterByExpr(y, group=y$samples$cluster)
table(keep.genes)
y<-y[keep.genes,, keep=FALSE]

y<-normLibSizes(y) 
head(y$samples,n=10L)

summary(y$samples$norm.factors)

cluster<-as.factor(y$samples$cluster)
plotMDS(y, pch=16,col=c(2:8)[cluster], main="MDS")
legend("topleft", legend=paste0("cluster",levels(cluster)),
       pch=16, col=2:8, cex=0.8)

donor<-factor(y$samples$sample)
design<-model.matrix(~cluster+ donor)
colnames(design)<-gsub("donor","",colnames(design))
colnames(design)[1]<-"Int"
head(design)

dim(design)

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

ncls <- nlevels(cluster)
contr <- rbind( matrix(1/(1-ncls), ncls, ncls),
                matrix(0, ncol(design)-ncls, ncls) )
diag(contr) <- 1
contr[1,] <- 0
rownames(contr) <- colnames(design)
colnames(contr) <- paste0("cluster", levels(cluster))
contr

qlf<-list()
for(i in 1:ncls){
  qlf[[i]]<-glmQLFTest(fit, contrast=contr[,i])
  qlf[[i]]$comparison <-paste0("cluster",levels(cluster)[i],"_vs_others")
}

topTags(qlf[[1]], n=10L)

dt<-lapply(lapply(qlf,decideTests),summary)
dt.all<-do.call("cbind", dt)
dt.all

top<-20
topMarkers<-list()
for(i in 1:ncls){
  ord<-order(qlf[[i]]$table$PValue, decreasing=FALSE)
  up<-qlf[[i]]$table$logFC[ord] > 0
  topMarkers[[i]]<-rownames(y)[ord[up][1:top]]
}
topMarkers<-unique(unlist(topMarkers))
topMarkers

lcpm<-cpm(y, log=TRUE)
annot<-data.frame(cluster=paste0("cluster",cluster))
rownames(annot)<-colnames(y)
ann_colors<-list(cluster=2:8)
names(ann_colors$cluster)<-paste0("cluster",levels(cluster))
pheatmap::pheatmap(lcpm[topMarkers,],breaks=seq(-2,2,length.out=101),
                   color=colorRampPalette(c("blue","white","red"))(100), scale="row",
                   cluster_cols=TRUE, border_color="NA",fontsize_row=5,
                   treeheight_row=70, treeheight_col=70,cutree_cols=7,
                   clustering_method="ward.D2", show_colnames=FALSE,
                   annotation_col=annot, annotation_colors=ann_colors)
ann_colors

pheatmap::pheatmap(lcpm[topMarkers,],breaks=seq(-2,2,length.out=101),
                                         color=colorRampPalette(c("blue","white","red"))(100), scale="row",
                                         cluster_cols=TRUE, border_color="NA",fontsize_row=5,
                                         treeheight_row=70, treeheight_col=70,cutree_cols=7,
                                         clustering_method="ward.D2", show_colnames=FALSE,
                                         annotation_col=annot, annotation_colors=ann_colors,
                                         show_rownames = F,legend=F,filename="~/ccsgab_projects/edger_heatmap.pdf")



