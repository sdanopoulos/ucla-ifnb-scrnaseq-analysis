#
# ==============================================================================
# Title: IFNb Explant Single-Cell RNA-Seq Analysis (Seurat 5.3.0)
# Author: Mariani Lab
# Date: May 2026
# R Version: 4.4.3 | Seurat Version: 5.3.0
# Description: This script performs quality control, normalization, RPCA
# integration, lineage subclustering, differential-expression analysis, and
# CellChat cell-cell communication analysis.
## Samples:
## Sample_H29934-IFN
## Sample_H29934-CTL
## Sample_H29925-IFN
## Sample_H29925-CTL
## Sample_H29921-IFN
## Sample_H29921-CTL
## Sample_H29912-IFN
## Sample_H29912-CTL
# ==============================================================================

# 1. ENVIRONMENT & DEPENDENCIES
rm(list = ls())

library(Seurat)
library(CellChat)
library(tidyverse)
library(here)
library(patchwork)
library(dplyr)
library(ggplot2)

set.seed(1234)

# Replace this placeholder with the folder where all output files should be saved.
output_dir <- "path to output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. LOAD DATA 
## Sample_H29912-CTL
H29912.CTL.data <- Read10X(data.dir = "path to data/Sample_H29912-CTL/filtered_feature_bc_matrix/")
H29912.CTL <- CreateSeuratObject(counts = H29912.CTL.data, project = "H29912_CTL" ,min.cells = 3)
H29912.CTL@meta.data$sampleID <- "H29912.CTL"
H29912.CTL@meta.data$sample <- "H29912"
H29912.CTL@meta.data$group <- "CTL"

## Sample_H29912-IFN
H29912.IFN.data <- Read10X(data.dir = "path to data/Sample_H29912-IFN/filtered_feature_bc_matrix")
H29912.IFN <- CreateSeuratObject(counts = H29912.IFN.data, project = "H29912_IFN" ,min.cells = 3)
H29912.IFN@meta.data$sampleID <- "H29912.IFN"
H29912.IFN@meta.data$sample <- "H29912"
H29912.IFN@meta.data$group <- "IFN"

## Sample_H29921-CTL
H29921.CTL.data <- Read10X(data.dir = "path to data/Sample_H29921-CTL/filtered_feature_bc_matrix/")
H29921.CTL <- CreateSeuratObject(counts = H29921.CTL.data, project = "H29921_CTL" ,min.cells = 3)
H29921.CTL@meta.data$sampleID <- "H29921.CTL"
H29921.CTL@meta.data$sample <- "H29921"
H29921.CTL@meta.data$group <- "CTL"

## Sample_H29921-IFN
H29921.IFN.data <- Read10X(data.dir = "path to data/Sample_H29921-IFN/filtered_feature_bc_matrix/")
H29921.IFN <- CreateSeuratObject(counts = H29921.IFN.data, project = "H29921_IFN" ,min.cells = 3)
H29921.IFN@meta.data$sampleID <- "H29921.IFN"
H29921.IFN@meta.data$sample <- "H29921"
H29921.IFN@meta.data$group <- "IFN"

## Sample_H29925-CTL
H29925.CTL.data <- Read10X(data.dir = "path to data/Sample_H29925-CTL/filtered_feature_bc_matrix/")
H29925.CTL <- CreateSeuratObject(counts = H29925.CTL.data, project = "H29925_CTL" ,min.cells = 3)
H29925.CTL@meta.data$sampleID <- "H29925.CTL"
H29925.CTL@meta.data$sample <- "H29925"
H29925.CTL@meta.data$group <- "CTL"

## Sample_H29925-IFN
H29925.IFN.data <- Read10X(data.dir = "path to data/Sample_H29925-IFN/filtered_feature_bc_matrix/")
H29925.IFN <- CreateSeuratObject(counts = H29925.IFN.data, project = "H29925_IFN" ,min.cells = 3)
H29925.IFN@meta.data$sampleID <- "H29925.IFN"
H29925.IFN@meta.data$sample <- "H29925"
H29925.IFN@meta.data$group <- "IFN"

## Sample_H29934-CTL
H29934.CTL.data <- Read10X(data.dir = "path to data/Sample_H29934-CTL/filtered_feature_bc_matrix/")
H29934.CTL <- CreateSeuratObject(counts = H29934.CTL.data, project = "H29934_CTL" ,min.cells = 3)
H29934.CTL@meta.data$sampleID <- "H29934.CTL"
H29934.CTL@meta.data$sample <- "H29934"
H29934.CTL@meta.data$group <- "CTL"

## Sample_H29934-IFN
H29934.IFN.data <- Read10X(data.dir = "path to data/Sample_H29934-IFN/filtered_feature_bc_matrix/")
H29934.IFN <- CreateSeuratObject(counts = H29934.IFN.data, project = "H29934_IFN" ,min.cells = 3)
H29934.IFN@meta.data$sampleID <- "H29934.IFN"
H29934.IFN@meta.data$sample <- "H29934"
H29934.IFN@meta.data$group <- "IFN"

## combine 8 samples
obj <- merge(H29912.CTL, y = c(H29912.IFN, H29921.CTL, H29921.IFN, H29925.CTL, H29925.IFN, H29934.CTL, H29934.IFN), 
                 add.cell.ids = c("H29912_CTL","H29912_IFN", "H29921_CTL", "H29921_IFN", "H29925_CTL", "H29925_IFN", "H29934_CTL", "H29934_IFN"), project = "IFNB")

# Record the number of cells retained at each QC stage. 
qc_counts_raw <- obj@meta.data %>%
  count(sampleID, name = "n_cells") %>%
  mutate(stage = "Before QC")

# The merged object is the only object needed from this point onward.
rm(
  H29912.CTL.data, H29912.IFN.data, H29921.CTL.data, H29921.IFN.data,
  H29925.CTL.data, H29925.IFN.data, H29934.CTL.data, H29934.IFN.data,
  H29912.CTL, H29912.IFN, H29921.CTL, H29921.IFN,
  H29925.CTL, H29925.IFN, H29934.CTL, H29934.IFN
)
invisible(gc())

# 3. QUALITY CONTROL (QC)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Apply manuscript-specific filters across all batches
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)

qc_counts_after_thresholds <- obj@meta.data %>%
  count(sampleID, name = "n_cells") %>%
  mutate(stage = "After feature/mitochondrial QC")

# Removing doublets
meta_scdblfinder <- read.table(
  file = "path to data/meta_scdblfinder.result.sct.04-2025.txt",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

if (!"scDblFinder.class" %in% colnames(meta_scdblfinder)) {
  stop("The doublet metadata file does not contain 'scDblFinder.class'.")
}

missing_doublet_cells <- setdiff(colnames(obj), rownames(meta_scdblfinder))
if (length(missing_doublet_cells) > 0) {
  stop(
    "Doublet metadata is missing ", length(missing_doublet_cells),
    " retained cells. Confirm that its row names match the merged Seurat cell names."
  )
}

meta_scdblfinder <- meta_scdblfinder[
  colnames(obj),
  "scDblFinder.class",
  drop = FALSE
]
obj <- AddMetaData(object = obj, metadata = meta_scdblfinder)
obj <- subset(obj, subset=scDblFinder.class == 'singlet')

qc_counts_after_doublet_removal <- obj@meta.data %>%
  count(sampleID, name = "n_cells") %>%
  mutate(stage = "After doublet removal")

write_csv(
  bind_rows(
    qc_counts_raw,
    qc_counts_after_thresholds,
    qc_counts_after_doublet_removal
  ) %>% select(stage, sampleID, n_cells),
  file.path(output_dir, "QC_cell_counts_by_sample_and_stage.csv")
)

# QC plots
p_nfeature <- VlnPlot(
  obj,
  features = "nFeature_RNA",
  group.by = "sampleID",
  pt.size = 0
) +
  geom_hline(yintercept = c(200, 8000), linetype = "dashed", linewidth = 0.45) +
  coord_cartesian(ylim = c(0, 8500)) +
  ggtitle("nFeature_RNA")

p_percent_mt <- VlnPlot(
  obj,
  features = "percent.mt",
  group.by = "sampleID",
  pt.size = 0
) +
  geom_hline(yintercept = 10, linetype = "dashed", linewidth = 0.45) +
  coord_cartesian(ylim = c(0, 11)) +
  ggtitle("percent.mt")

p_ncount <- VlnPlot(
  obj,
  features = "nCount_RNA",
  group.by = "sampleID",
  pt.size = 0
) +
  ggtitle("nCount_RNA")

p_qc_after_filtering <- p_nfeature + p_percent_mt + p_ncount

ggsave(
  filename = file.path(
    output_dir,
    "QC_violin_plots_after_filtering.png"
  ),
  plot = p_qc_after_filtering,
  width = 14,
  height = 4.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

# 4. NORMALIZATION
## SCT transform
obj <- SCTransform(obj, verbose = TRUE)
obj <- RunPCA(obj, npcs = 50, verbose =T)
ElbowPlot(obj, ndims = 50)

# 5. INTEGRATION
## RPCA integration
integrated_seurat <- IntegrateLayers(
  object = obj,
  method = RPCAIntegration,
  normalization.method = "SCT",
  assay = "SCT",
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = TRUE)

## identify clusters
integrated_seurat <- FindNeighbors(integrated_seurat, reduction = "integrated.rpca", dims = 1:25)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.2, cluster.name = "rpca_clusters.02")
integrated_seurat <- RunUMAP(integrated_seurat, reduction = "integrated.rpca", dims = 1:25, 
                             reduction.name = "umap.rpca",
                             n.neighbors = 35L, min.dist = 0.2)

# 6. GENERATE VISUALIZATIONS
# Plot A: UMAP colored by Cluster
DimPlot(integrated_seurat, reduction = "umap.rpca", label = TRUE, repel = TRUE)

# Plot B: UMAP colored by condition to inspect condition mixing after integration
DimPlot(integrated_seurat, reduction = "umap.rpca", group.by = "group", alpha = 0.5)

# 7. ANALYSIS: IDENTIFYING CELL-TYPE MARKERS 
integrated_seurat[["RNA"]] <- JoinLayers(integrated_seurat[["RNA"]])
integrated_seurat <- PrepSCTFindMarkers(integrated_seurat)
all_cluster_markers <- FindAllMarkers(object = integrated_seurat, only.pos = TRUE)
write.csv(
  all_cluster_markers,
  file.path(output_dir, "ifnb_all_cluster_markers.csv"),
  row.names = FALSE
)

# 8. Azimuth annotation
library(Azimuth)
obj.lung <- RunAzimuth(integrated_seurat, reference = "lungref", query.modality = "SCT", umap.name = "ref.umap.lung", assay = "SCT")
obj.fetus <- RunAzimuth(obj.lung, reference = "fetusref", query.modality = "SCT", umap.name = "ref.umap.fetus", assay = "SCT")
integrated_seurat <- obj.fetus

saveRDS(
  integrated_seurat,
  file = file.path(output_dir, "ifnb.explants.big.lineage.2026.rds")
)

# 9. Subclustering: Mesenchymal lineage
## subset mesenchymal cells from cluster
obj <- integrated_seurat
Idents(object = obj) <- "rpca_clusters.02"
mesenchymal <- subset(x = obj, idents = c("0", "2", "4", "6", "8", "9", "10"))

## split layers
obj <- mesenchymal
DefaultAssay(object = obj) <- "RNA"
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sampleID)

## normalization
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 3000)
obj <- ScaleData(obj, features = VariableFeatures(obj))
obj <- RunPCA(obj, npcs = 50)
ElbowPlot(obj, ndims = 50)
obj <- FindNeighbors(obj, dims = 1:20, reduction = "pca")
obj <- FindClusters(obj, resolution = 0.2, cluster.name = "unintegrated_clusters")
obj <- RunUMAP(obj, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")

## RPCA integration
obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = TRUE) 

## identify clusters
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:20)
obj <- FindClusters(obj, resolution = 0.20, cluster.name = "rpca_clusters.20")
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:20,
               n.neighbors = 30L, min.dist = 0.1,
               reduction.name = "umap.rpca",return.model = TRUE)
DimPlot(obj, reduction = "umap.rpca", group.by = "rpca_clusters.20", label = TRUE)

## Find markers
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
Idents(obj) <- "rpca_clusters.20"
obj.markers <- FindAllMarkers(obj, only.pos = TRUE)
write.csv(
  obj.markers,
  file.path(output_dir, "mesenchymal_subcluster_markers.csv"),
  row.names = FALSE
)

## 10. ANALYSIS: CONDITION-SPECIFIC DE
# Purpose: Find genes changed by IFNb treatment *within* a specific cell type
# Example: What changes in Cluster 2 when you treat the cells with IFNb
##
obj$cluster.group <- paste(obj$rpca_clusters.20, obj$group, sep = "_")
Idents(obj) <- "cluster.group"

## cluster 2
cluster2_de_results <- FindMarkers(
  object = obj,
  ident.1 = "2_IFN",
  ident.2 = "2_CTL",
  test.use = "MAST"
)

# Convert row names to an explicit 'gene' column before exporting
cluster2_de_results <- cluster2_de_results %>% rownames_to_column(var = "gene")
write.csv(
  cluster2_de_results,
  file.path(output_dir, "mesenchymal_cluster2_IFN_vs_CTL_MAST.csv"),
  row.names = FALSE
)

# 11. Plot cell-cycle proliferation genes
Idents(obj) <- "group"

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

obj <- CellCycleScoring(object = obj, s.features = s.genes, g2m.features = g2m.genes)
DotPlot(obj, features = g2m.genes, split.by = "group", dot.scale=8, cols="RdBu")
DotPlot(obj, features = s.genes, split.by = "group", dot.scale=8, cols="RdBu")

mesenchymal <- obj
saveRDS(
  mesenchymal,
  file = file.path(output_dir, "ifnb.explants.mesenchymal.lineage.2026.rds")
)

# 12. Subclustering: Epithelial lineage
## subset epithelial cells from cluster
obj <- integrated_seurat
Idents(object = obj) <- "rpca_clusters.02"
# These parent clusters reproduce the original marker/Azimuth-based lineage
# assignment. Re-check this list if the full-object clustering parameters change.
epithelial <- subset(x = obj, idents = c("1", "12", "13"))
obj <- epithelial
DefaultAssay(object = obj) <- "RNA"
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sampleID)

## normalization
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 3000)
obj <- ScaleData(obj, features = VariableFeatures(obj))
obj <- RunPCA(obj, npcs = 50)
ElbowPlot(obj, ndims = 50)
obj <- FindNeighbors(obj, dims = 1:12, reduction = "pca")
obj <- FindClusters(obj, resolution = 0.2, cluster.name = "unintegrated_clusters")
obj <- RunUMAP(obj, dims = 1:12, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(
  obj,
  reduction = "umap.unintegrated",
  group.by = "unintegrated_clusters",
  label = TRUE
)

## RPCA integration
obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = TRUE) 

## identify cluster
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:12)
obj <- FindClusters(obj, resolution = 0.20, cluster.name = "rpca_clusters.20")
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:12, reduction.name = "umap.rpca",
               n.neighbors = 35L, min.dist = 0.2, return.model = TRUE)
DimPlot(obj, reduction = "umap.rpca", group.by = "rpca_clusters.20", label = TRUE)

## Find markers
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
Idents(obj) <- "rpca_clusters.20"
obj.markers <- FindAllMarkers(obj, only.pos = TRUE)
write.csv(
  obj.markers,
  file.path(output_dir, "epithelial_subcluster_markers.csv"),
  row.names = FALSE
)

## 13. ANALYSIS: CONDITION-SPECIFIC DE
# Purpose: Find genes changed by IFNb treatment *within* a specific cell type
# Example: What changes in Cluster 2 when you treat the cells with IFNb
##
obj$cluster.group <- paste(obj$rpca_clusters.20, obj$group, sep = "_")

Idents(obj) <- "cluster.group"

## cluster 2
cluster2_de_results <- FindMarkers(object = obj, 
                                   ident.1 = "2_IFN", 
                                   ident.2 = "2_CTL",
                                   test.use = "MAST")

# Convert row names to an explicit 'gene' column before exporting
cluster2_de_results <- cluster2_de_results %>% rownames_to_column(var = "gene")
write.csv(
  cluster2_de_results,
  file.path(output_dir, "epithelial_cluster2_IFN_vs_CTL_MAST.csv"),
  row.names = FALSE
)

# 14. Plot gene of interest
Idents(obj) <- "group"
genes <- c("ABCA3", "ADGRF5", "NAPSA", "SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SLC34A2")
DotPlot(object = obj, features = genes, dot.scale=8, cols="RdBu") +coord_flip()


# 15. plot cell-cycle proliferation genes
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

obj <- CellCycleScoring(object = obj, s.features = s.genes, g2m.features = g2m.genes)
DotPlot(obj, features = g2m.genes, split.by = "group", dot.scale=8, cols="RdBu")
DotPlot(obj, features = s.genes, split.by = "group", dot.scale=8, cols="RdBu")

epithelial <- obj
saveRDS(
  epithelial,
  file = file.path(output_dir, "ifnb.explants.epithelial.lineage.2026.rds")
)


# 16. CELLCHAT CELL-CELL COMMUNICATION ANALYSIS

# Assign EC labels to the epithelial lineage subclusters
my_data <- epithelial
Idents(my_data) <- "rpca_clusters.20"
epithelial_cluster_map <- setNames(
  sprintf("EC%02d", 0:7),
  as.character(0:7)
)
my_data <- RenameIdents(my_data, epithelial_cluster_map)
my_data$cluster_name <- Idents(my_data)
obj.1 <- my_data

# Assign MC labels to the mesenchymal lineage subclusters
my_data <- mesenchymal
Idents(my_data) <- "rpca_clusters.20"
mesenchymal_cluster_map <- setNames(
  sprintf("MC%02d", 0:9),
  as.character(0:9)
)
my_data <- RenameIdents(my_data, mesenchymal_cluster_map)
my_data$cluster_name <- Idents(my_data)
obj.2 <- my_data

epi.cells <- rownames(obj.1@meta.data)[obj.1$group == "IFN"]
mes.cells <- rownames(obj.2@meta.data)[obj.2$group == "IFN"]

epi.data <- LayerData(obj.1, assay = "RNA", layer = "data")[, epi.cells, drop = FALSE]
mes.data <- LayerData(obj.2, assay = "RNA", layer = "data")[, mes.cells, drop = FALSE]

common.genes <- intersect(rownames(epi.data), rownames(mes.data))
data.input <- cbind(
  epi.data[common.genes, , drop = FALSE],
  mes.data[common.genes, , drop = FALSE]
)

meta <- rbind(
  obj.1@meta.data[epi.cells, , drop = FALSE],
  obj.2@meta.data[mes.cells, , drop = FALSE]
)
meta <- meta[colnames(data.input), , drop = FALSE]
meta$samples <- meta$sample

cellchat <- createCellChat(
  object = data.input,
  meta = meta,
  group.by = "cluster_name"
)

# Select the Human database
cellchat@DB <- CellChatDB.human

# Subset and identify overexpressed features
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication probability
cellchat <- computeCommunProb(cellchat, type = "triMean")

# Filter out weak interactions (minimum 10 cells per group)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Compute pathway-level signaling and network centrality
cellchat <- computeCommunProbPathway(cellchat)

## Calculate the aggregated cell-cell communication network 
cellchat <- aggregateNet(cellchat)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
netAnalysis_signalingRole_scatter(cellchat)
cellchat.ifn <- cellchat
saveRDS(cellchat.ifn, file.path(output_dir, "cellchat_IFN.rds"))

# Create the CTL CellChat object from the normalized RNA matrix
epi.cells <- rownames(obj.1@meta.data)[obj.1$group == "CTL"]
mes.cells <- rownames(obj.2@meta.data)[obj.2$group == "CTL"]

epi.data <- LayerData(obj.1, assay = "RNA", layer = "data")[, epi.cells, drop = FALSE]
mes.data <- LayerData(obj.2, assay = "RNA", layer = "data")[, mes.cells, drop = FALSE]

common.genes <- intersect(rownames(epi.data), rownames(mes.data))
data.input <- cbind(
  epi.data[common.genes, , drop = FALSE],
  mes.data[common.genes, , drop = FALSE]
)

meta <- rbind(
  obj.1@meta.data[epi.cells, , drop = FALSE],
  obj.2@meta.data[mes.cells, , drop = FALSE]
)
meta <- meta[colnames(data.input), , drop = FALSE]
meta$samples <- meta$sample

cellchat <- createCellChat(
  object = data.input,
  meta = meta,
  group.by = "cluster_name"
)

# Select the Human database
cellchat@DB <- CellChatDB.human

# Subset and identify overexpressed features
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication probability
cellchat <- computeCommunProb(cellchat, type = "triMean")

# Filter out weak interactions (minimum 10 cells per group)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Compute pathway-level signaling and network centrality
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network 
cellchat <- aggregateNet(cellchat)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
netAnalysis_signalingRole_scatter(cellchat)
cellchat.ctl <- cellchat
saveRDS(cellchat.ctl, file.path(output_dir, "cellchat_CTL.rds"))

# Merge the condition-specific CellChat results for comparison
object.list <- list(CTRL = cellchat.ctl, IFN = cellchat.ifn)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(cellchat, file.path(output_dir, "cellchat_CTRL_IFN_merged.rds"))

# Compare the total number of inferred interactions
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(2,1),
                           title.name = "Total Number of Interactions",
                           size.text = 10, 
                           color.use = c("#B2182B","#2166AC" ))

# Compare total inferred interaction strength
gg2 <- compareInteractions(cellchat, show.legend = F, 
                           group = c(2,1), measure = "weight", 
                           color.use = c("#B2182B","#2166AC" ),
                           title.name = "Total Interaction Strength (Weight)")
gg1 + gg2

# Information Flow Comparison (Pathway Level)
rankNet(cellchat, slot.name = "netP", mode = "comparison", 
        measure = "weight", sources.use = NULL, targets.use = NULL, 
        stacked = T, do.stat = FALSE,
        #color.use = c("#00BFC4", "#F8766D"),
        color.use = c("#2166AC", "#B2182B"),
        font.size = 8,
        do.flip = FALSE)

## Pathways of interest: WNT, BMP, and FGF
for (pathway in c("WNT", "BMP", "FGF")) {
  pathway_present <- vapply(
    object.list,
    function(x) pathway %in% x@netP$pathways,
    logical(1)
  )
  
  if (!all(pathway_present)) {
    warning(
      pathway,
      " signaling was not detected in every condition and was not plotted."
    )
    next
  }
  
  weight.max <- getMaxWeight(
    object.list,
    slot.name = "netP",
    attribute = pathway
  )
  
  png(
    filename = file.path(
      output_dir,
      paste0("CellChat_", pathway, "_CTRL_IFN_circle.png")
    ),
    width = 3600,
    height = 1800,
    res = 300
  )
  
  par(mfrow = c(1, 2), xpd = TRUE)
  
  for (i in seq_along(object.list)) {
    netVisual_aggregate(
      object.list[[i]],
      signaling = pathway,
      layout = "circle",
      edge.weight.max = weight.max[1],
      edge.width.max = 10,
      signaling.name = paste(pathway, names(object.list)[i])
    )
  }
  
  dev.off()
}

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)
