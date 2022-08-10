source("source_satac.R")

# load seurat object with denoised gene activity counts
combined <- readRDS("results/combined_brca_denoised_rna.rds")
DefaultAssay(combined) <- "RNA_dca"
Idents(combined) <- "peaks_snn_res.0.5"
combined <- RunPCA(combined)

# link peaks to genes
markers <- FindAllMarkers(combined, logfc.threshold = 0.1, only.pos = T, assay = "RNA") %>%
  filter(p_val_adj < 0.05)
DefaultAssay(combined) <- "peaks"
combined <- RegionStats(combined, genome = BSgenome.Hsapiens.UCSC.hg38)

combined <- LinkPeaks(
  object = combined,
  peak.assay = "peaks",
  expression.assay = "RNA_dca",
  genes.use = unique(markers$gene),
  min.distance = 2000
)
linked_peaks <- combined@assays[["peaks"]]@links@elementMetadata %>% as.data.frame()

## visium - create seurat object with stutility 
infoTable <- data.frame(
  samples = "Visium_BCSA4/human_feature_bc_matrix.h5",
  spotfiles = "Visium_BCSA4/spatial/tissue_positions_list.csv",
  imgs = "Visium_BCSA4/spatial/tissue_hires_image.png",
  json = "Visium_BCSA4/spatial/scalefactors_json.json"
)

visium <- InputFromTable(infotable = infoTable, 
                     min.gene.count = 0, 
                     min.gene.spots = 0,
                     min.spot.count = 0,
                     platform =  "Visium")
visium <- LoadImages(visium)
visium <- visium %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:10) %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.4)

#ON terminal:
## run dca on visium count matrix using default paramaters

# load denoised visium data
mtx <- read.table("visium_denoised/mean.tsv")
colnames(mtx) <- gsub("\\.", "-", colnames(mtx))
visium <- subset(visium, cells = colnames(mtx))
visium[["RNA_dca"]] <- CreateAssayObject(counts = as.matrix(mtx))
rm(mtx)

# plot representative genes
DefaultAssay(combined) <- "RNA_dca"
p1 <- ST.FeaturePlot(combined, "ERBB2", cols = magenta_scale, max.cutoff = "q99", indices = 1, min.cutoff = "q1")
p2 <- ST.FeaturePlot(combined, "C1QB", cols = magenta_scale, max.cutoff = "q99", indices = 1, min.cutoff = "q1")

DefaultAssay(combined) <- "dca"
p3 <- ST.FeaturePlot(combined, "chr17-39723461-39723962", cols = magenta_scale, max.cutoff = "q99", indices = 1, min.cutoff = "q1")
p4 <- ST.FeaturePlot(combined, "chr1-22652746-22653247", cols = magenta_scale, max.cutoff = "q99", indices = 1, min.cutoff = "q1")

DefaultAssay(visium) <- "RNA_dca"
p5 <- ST.FeaturePlot(visium, "ERBB2", cols = cyan_scale, max.cutoff = "q95", min.cutoff = "q1")
p6 <- ST.FeaturePlot(visium, "C1QB", cols = cyan_scale, max.cutoff = "q99", min.cutoff = "q1")

p1|p3|p5
p2|p4|p6

saveRDS(visium, "visium_denoised_rna.rds")


# compare sATAC with scRNA-seq
## scRNA-seq atlas of BRCA obtained from https://www.nature.com/articles/s41588-021-00911-1
meta_sc <- read.csv("data/Wu_etal_2021_BRCA_scRNASeq/metadata.csv", row.names = 1)
mtx <- Read10X("data/Wu_etal_2021_BRCA_scRNASeq/", gene.column = 1)

# preprocess single cell data
sc_obj <- CreateSeuratObject(mtx, meta.data = meta_sc) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

# label transfer using the original authors' classification of cell types
transfer.anchors <- FindTransferAnchors(reference = sc_obj, 
                                        query = combined, 
                                        features = VariableFeatures(sc_obj), 
                                        reference.assay = "RNA", 
                                        query.assay = "RNA_dca", 
                                        reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = sc_obj$celltype_major, 
                                     weight.reduction = combined[["pca"]],
                                     dims = 1:15)

combined <- AddMetaData(combined, 
                        metadata = celltype.predictions)

# plot prediction scores as an indication of the type of cells found in the sATAC clusters 
prediction_scores <- paste0("prediction.score.", levels(as.factor(obj$celltype_major)))
DotPlot(combined, features = prediction_scores, cols = c("gray93", "#AE017E")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# repeat analysis for the visium section
transfer.anchors <- FindTransferAnchors(reference = sc_obj, 
                                        query = visium, 
                                        features = VariableFeatures(sc_obj), 
                                        reference.assay = "RNA", 
                                        query.assay = "RNA", 
                                        reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = sc_obj$celltype_major, 
                                     weight.reduction = visium[["pca"]],
                                     dims = 1:15)

visium <- AddMetaData(visium, metadata = celltype.predictions)

DotPlot(visium, features =prediction_scores, cols = c("gray93", "#AE017E")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# save object with the single cell predictions
saveRDS(visium, "visium_brca_scpredicted.rds")
saveRDS(combined, "atac_brca_scpredicted.rds")
