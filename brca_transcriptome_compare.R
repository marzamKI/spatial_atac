source("source_satac.R")

meta_sc <- read.csv("data/Wu_etal_2021_BRCA_scRNASeq/metadata.csv", row.names = 1)
mtx <- Read10X("data/Wu_etal_2021_BRCA_scRNASeq/", gene.column = 1)
sc_obj <- CreateSeuratObject(mtx, meta.data = meta_sc) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

combined <- readRDS("results/combined_brca_denoised_rna.rds")
DefaultAssay(combined) <- "RNA_dca"
Idents(combined) <- "peaks_snn_res.0.5"
combined <- RunPCA(combined)

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

prediction_scores <- paste0("prediction.score.", levels(as.factor(obj$celltype_major)))
DotPlot(combined, features = prediction_scores, cols = c("gray93", "#AE017E")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

saveRDS(combined, "atac_brca_scpredicted.rds")
