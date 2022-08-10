source("source_satac.R")

# compare sATAC with visium

#process all visium visiumctions
infoTable <- data.frame(
  samples = c("V10B01-135_D1_manual_alignment_v100/mouse_filtered_feature_bc_matrix.h5",
              "V10S29-086_D1_spaceranger/outs/mouse_feature_bc_matrix.h5",
              "220421_Visium_e12/220421_e12_B1/outs/filtered_feature_bc_matrix.h5"),
  spotfiles = c("V10B01-135_D1_manual_alignment_v100/spatial/tissue_positions_list.csv",
                "V10S29-086_D1_spaceranger/outs/spatial/tissue_positions_list.csv",
                "220421_Visium_e12/220421_e12_B1/outs/spatial/tissue_positions_list.csv"),
  imgs = c("V10B01-135_D1_manual_alignment_v100/spatial/tissue_hires_image.png",
           "V10S29-086_D1_spaceranger/outs/spatial/tissue_hires_image.png",
           "220421_Visium_e12/220421_e12_B1/outs/spatial/tissue_hires_image.png"),
  json = c("V10B01-135_D1_manual_alignment_v100/spatial/scalefactors_json.json",
           "V10S29-086_D1_spaceranger/outs/spatial/scalefactors_json.json",
           "220421_Visium_e12/220421_e12_B1/outs/spatial/scalefactors_json.json")
)

visium <- InputFromTable(infotable = infoTable, 
                     min.gene.count = 0, 
                     min.gene.spots = 0,
                     min.spot.count = 0,
                     platform =  "Visium")

visium <- LoadImages(visium, time.resolve = FALvisium, verbovisium = TRUE)
visium$visiumction <- visium@tools$Staffli@meta.data$sample

# Perform standard analysis of each modality independently 
# RNA analysis
visium <- NormalizeData(visium) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(visium)

visium <- RunHarmony(visium, 
                 group.by.vars = "section", 
                 reduction = "pca", dims.use = 1:15, 
                 assay.use = "RNA", 
                 project.dim = F,
                 verbovisium = T) %>%
  RunUMAP(reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony") %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.3)
visium$seurat_clusters_harmony <- visium$seurat_clusters

saveRDS(visium, "visium_allages.rds")
write.csv(visium@assays$RNA@counts, "visium_counts.csv")

### denoise visium with DCA and default parameters ###
mtx <- read.table("visium_denoised/mean.tsv")
colnames(mtx) <- gsub("\\.", "-", colnames(mtx))
visium[["RNA_dca"]] <- CreateAssayObject(counts = as.matrix(mtx))
rm(mtx)
saveRDS(visium, "visium_denoised_rna.rds")

DefaultAssay(visium) <- "RNA"
DefaultAssay(combined) <- "RNA"

markers <- FindAllMarkers(combined, logfc.threshold = 0.2, only.pos = T)
table(markers$cluster)

# calculate gene signatures for each cluster - 
# first from sATAC and check scores on Visium, then from Visium and check scores on sATAC clusters
# DEA followed by seurat's module score 
genes_signatures <- list()
for(i in 0:10){
  genes <- markers %>%
    filter(p_val_adj < 0.05) %>%
    filter(cluster %in% i) %>%
    select(gene)
  cluster <- as.character(i)
  genes_signatures[[cluster]] <- genes$gene
}

visium <- AddModuleScore(visium, genes_signatures)
combined <- AddModuleScore(combined, genes_signatures)

pdf("signaturescore_atactovisium.pdf", height = 8, width = 7)
for(i in 1:11){
  signature <- paste0("Cluster",i)
  p1 <- ST.FeaturePlot(visium, signature, 
                       cols = magenta_scale, 
                       min.cutoff = "q1", max.cutoff = "q99", 
                       pt.size = 0.9)
  p2 <- ST.FeaturePlot(combined, signature, 
                       cols = magenta_scale, 
                       indices = c(1,3,5), 
                       min.cutoff = "q1", max.cutoff = "q99", 
                       pt.size = 1)
  print(p1 | p2)
}
dev.off()

#gene signature from visium
Idents(visium) <- "seurat_clusters_harmony"
markers <- FindAllMarkers(visium, logfc.threshold = 0.2, only.pos = T)
table(markers$cluster)

#calculate signatures
genes_signatures <- list()
for(i in 0:11){
  genes <- markers %>%
    filter(p_val_adj < 0.05) %>%
    filter(cluster %in% i) %>%
    select(gene)
  cluster <- as.character(i)
  genes_signatures[[cluster]] <- genes$gene
}

visium <- AddModuleScore(visium, genes_signatures, name = "Visium_cluster")
combined <- AddModuleScore(combined, genes_signatures, name = "Visium_cluster")

pdf("signaturescore_visiumtoatac.pdf", height = 8, width = 7)
for(i in 1:12){
  signature <- paste0("Visium_cluster",i)
  p1 <- ST.FeaturePlot(visium, signature, 
                       cols = cyan_scale, 
                       min.cutoff = "q1", max.cutoff = "q99", 
                       pt.size = 0.9)
  p2 <- ST.FeaturePlot(combined, signature, 
                       cols = cyan_scale, 
                       indices = c(1,3,5), 
                       min.cutoff = "q1", max.cutoff = "q99", 
                       pt.size = 1)
  print(p1 | p2)
}
dev.off()

# plot module score using dotplot
levels(visium) <- c("3", "1", "2","4","8", "9", "6","0","5","10","11","7")
signatures <- paste0("Visium_cluster",c("4","3","5","10", "9", "12","2","7","6","1","11", "8"))
levels(combined) <- rev(c("5", "7", "2","9","3", "4", "8","0","1","10","6"))
DotPlot(combined, features = signatures, cols = c(cyan_scale[1], cyan_scale[8]))
signatures <- paste0("Cluster",c("6","8","3","1", "11", "2","4","5","9","10","7"))
DotPlot(visium, features = signatures, cols = c(magenta_scale[1], magenta_scale[8]))


# perform clustering of visium data based on DCA and compare with original clustering
DefaultAssay(visium) <- "RNA_dca"
visium <- NormalizeData(visium) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(visium)

visium <- RunHarmony(visium, 
                 group.by.vars = "section", 
                 reduction = "pca", dims.use = 1:15, 
                 assay.use = "RNA_dca", 
                 project.dim = F,
                 verbovisium = T) %>%
  RunUMAP(reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony") %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.3)
visium$seurat_clusters_harmony_dca <- visium$seurat_clusters

Idents(visium) <- "RNA_snn_res.0.3"
levels(visium) <- c("3", "2", "1","4","9","11", "6","8","0", "5","10","7")

Idents(visium) <- "seurat_clusters_harmony_dca"
levels(visium) <- c("4", "2", "3","7","5", "11","1", "6","0","12","8","10","9")
visium$seurat_clusters_harmony_dca <- Idents(visium)

p1 <- ST.FeaturePlot(visium, "ident", ncol = 3, cols = bright_visium)
p2 <- ST.FeaturePlot(visium, "seurat_clusters_harmony_dca", ncol = 3, cols = visium_dca)
p1+p2

Idents(visium) <- "seurat_clusters_harmony_dca"
prop <- prop.table(table(visium$seurat_clusters_harmony_dca, visium@active.ident), 1)*100
heatmap.2(prop[order(nrow(prop):1),], 
          Colv = NA, Rowv = NA, scale="none", 
          xlab="peaks cluster", ylab="dca clusters", 
          col = hm_colors, RowSideColors=rev(visium_dca), ColSideColors = bright_visium)
