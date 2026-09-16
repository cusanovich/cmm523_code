#' Convert a Seurat object to a monocle3 cell_data_set
#'
#' Does what SeuratWrappers::as.cell_data_set() does, without the dependency.
#' SeuratWrappers pulls in a large tree of spatial-transcriptomics packages
#' (Banksy, SpatialExperiment, leidenAlg, sccore) that this course never uses.
#'
#' A cell_data_set is not a mysterious object. It is a counts matrix, a table
#' of cell metadata, a table of gene metadata, and optionally some dimensional
#' reductions. Everything below is moving those four things across.
#'
#' @param obj         A Seurat object.
#' @param assay       Which assay to take counts from. Default "RNA".
#' @param reductions  Named vector mapping Seurat reduction names to the names
#'                    monocle3 expects, e.g. c(pca = "PCA", umap = "UMAP").
#'                    Pass NULL to skip and let monocle3 compute its own.
#'
#' @return A monocle3 cell_data_set.
seurat_to_cds <- function(obj,
                          assay = "RNA",
                          reductions = c(pca = "PCA", umap = "UMAP")) {

  stopifnot(requireNamespace("monocle3", quietly = TRUE))

  # --- 1. counts ------------------------------------------------------------
  # monocle3 wants raw counts, not normalized data -- it does its own
  # normalization internally. In Seurat 5 this is a "layer"; older code and
  # tutorials call it a "slot".
  counts <- SeuratObject::LayerData(obj, assay = assay, layer = "counts")

  # --- 2. cell metadata -----------------------------------------------------
  cell_meta <- obj@meta.data

  # Carry the active identities across as well. Seurat stores the "current"
  # labels separately from the metadata table, and it is easy to lose them.
  cell_meta$seurat_ident <- as.character(Seurat::Idents(obj))

  # --- 3. gene metadata -----------------------------------------------------
  # monocle3 REQUIRES a column literally named gene_short_name. Several of its
  # plotting functions look for it by name and fail with an unhelpful error if
  # it is absent.
  gene_meta <- data.frame(
    gene_short_name = rownames(counts),
    row.names = rownames(counts),
    stringsAsFactors = FALSE
  )

  cds <- monocle3::new_cell_data_set(
    expression_data = counts,
    cell_metadata   = cell_meta,
    gene_metadata   = gene_meta
  )

  # --- 4. size factors ------------------------------------------------------
  # Normally set by preprocess_cds(). If you are bringing your own reductions
  # you will skip that step, so set them here or downstream functions complain.
  cds <- monocle3::estimate_size_factors(cds)

  # --- 5. dimensional reductions -------------------------------------------
  # Optional. Transferring Seurat's PCA and UMAP means the trajectory is drawn
  # on the embedding the student already knows, rather than a new one monocle3
  # computed. That makes the result easier to interpret and easier to compare.
  if (!is.null(reductions)) {
    for (seurat_name in names(reductions)) {
      monocle_name <- reductions[[seurat_name]]
      if (seurat_name %in% SeuratObject::Reductions(obj)) {
        SingleCellExperiment::reducedDims(cds)[[monocle_name]] <-
          SeuratObject::Embeddings(obj, reduction = seurat_name)
      } else {
        warning("Reduction '", seurat_name, "' not found in the Seurat object; skipping.")
      }
    }
  }

  cds
}


# ============================================================================
# NOTES AND LIMITATIONS
#
# What this does NOT copy across:
#
#   Clusters and partitions. monocle3 keeps these in an internal slot
#   (cds@clusters) with a structure that is not part of its documented API, so
#   writing into it directly is fragile. Run monocle3::cluster_cells(cds)
#   instead -- it takes a few seconds and gives monocle3 the partition
#   information that learn_graph() needs.
#
#   Normalized or scaled data. monocle3 normalizes internally from counts.
#   Handing it pre-normalized values would normalize them twice.
#
#   Multiple assays. Only the assay you name is transferred. If you need both
#   RNA and ATAC in one cell_data_set, this function is not enough.
#
# Typical use:
#
#   cds <- seurat_to_cds(pbmc)
#   cds <- monocle3::cluster_cells(cds)
#   cds <- monocle3::learn_graph(cds)
#   cds <- monocle3::order_cells(cds)
#
# If any of that fails in a way that points at the object rather than your
# data, SeuratWrappers::as.cell_data_set() is the fallback -- it handles some
# edge cases this does not.
# ============================================================================
