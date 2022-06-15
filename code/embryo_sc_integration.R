#
source("source_satac.R")

#scRNAseq from http://mousebrain.org/development/downloads.html
#sc preprocessing
lfile <- connect(filename = "~/Downloads/dev_all.loom", mode = "r+", skip.validate = TRUE)
names(lfile$col.attrs)

cells <- lfile$col.attrs$CellID[]
regions <- lfile$col.attrs$Tissue[] %in% "ForebrainDorsal"
age <- lfile$col.attrs$Age[] %in% c("e15.0", "e15.5")
cell_types <- lfile$col.attrs$Class[] %in% c("Neuroblast", "Neuron", "Radial glia")
unique_cells <- !cells %in% cells[duplicated(cells)]
data.neurons <- lfile[["matrix"]][unique_cells & regions & cell_types & age, ]

# get metadata
attrs <- names(lfile$col.attrs)
attrs <- attrs[c(1, 3:11, 26:37)]

meta <- list()
for(i in attrs){
  meta[[i]] <- lfile$col.attrs[[i]][][unique_cells & regions & cell_types & age]
}
meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta_df$CellID

# get genes
genes <- lfile$row.attrs$Gene[]
dim(data.neurons)
cells <- lfile$col.attrs[["CellID"]][unique_cells & regions & cell_types & age]
rownames(data.neurons) <- cells
colnames(data.neurons) <- genes
data.neurons <- t(data.neurons)

obj_dev <- CreateSeuratObject(data.neurons, meta.data = meta_df)

# filter other cells
classes <- c("Cortical or hippocampal glutamatergic", "Forebrain",
             "Forebrain glutamatergic", "Neuronal intermediate progenitor",
             "Dorsal forebrain")
Idents(obj_dev) <- "Subclass"
obj_dev <- subset(obj_dev, idents = classes, invert = F)

lfile$close_all()
rm(lfile)
rm(data.neurons)

obj_dev <- obj_dev %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(obj_dev, features = rownames(obj_dev)) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.1) %>%
  RunUMAP(dims = 1:10)

new.cluster.ids <- c("Radial glia", "Neuroblasts", "Neurons", "Cycling neuroblasts")

names(new.cluster.ids) <- levels(obj_dev)
obj_dev <- RenameIdents(obj_dev, new.cluster.ids)
obj_dev$new_clusters <- Idents(obj_dev)
saveRDS(obj_dev, "data/dev_e15.rds")

set.seed(123)
obj_dev <- obj_dev[, sample(colnames(obj_dev), size =1500, replace=F)]

# subset sATAC data
combined <- readRDS("combined_denoised_rna_nmf.rds") 

# load ctx spot coordinates and subset data
samples <- c("220327_A1", "220327_A2")
table <- list()
for(i in seq_along(samples)){
  sample <- samples[i]
  table[[sample]] <- read.table(paste0("meta/", sample, "_ctx.tsv"), sep = "\t", header = T)
  table[[sample]]$barcode_1 <- paste0(table[[sample]]$barcode, "-1_", i)
}
infoTable <- do.call("rbind", table) %>% as.data.frame()

cortex <- SubsetSTData(combined, expression = barcode %in% infoTable$barcode_1)
saveRDS(cortex, "cortex_denoised_rna.rds")
rm(combined)

DefaultAssay(cortex) <- "RNA_dca"
cortex <- cortex %>% 
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

cortex$mode <- "sATAC"
obj_dev$mode <- "scRNA"
obj_list <- list(obj_dev, cortex)

obj_list <- lapply(X = obj_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = obj_list)
anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features)
combined_e15 <- IntegrateData(anchorset = anchors)
DefaultAssay(combined_e15) <- "integrated"

combined_e15 <- combined_e15 %>%
  ScaleData() %>% 
  RunPCA() %>%
  RunUMAP(dims = 1:7) %>%
  FindNeighbors(dims = 1:7) %>%
  FindClusters(resolution = 0.1)
DimPlot(combined_e15, reduction = "umap", split.by = "mode", ncol = 2)

#monocle
rgs <- rownames(combined_e15@meta.data)[combined_e15$integrated_snn_res.0.1 == "3"]
cds <- as.cell_data_set(combined_e15) %>%
  cluster_cells(reduction_method = "UMAP") %>%
  learn_graph(use_partition = T) %>%
  order_cells(reduction_method = "UMAP", root_cells = rgs)

# plot trajectories colored by pseudotime
plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)

cortex <- AddMetaData(
  object = cortex,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime_cortex_e15")

FeatureOverlay(cortex, features = "Pseudotime_cortex_e15", sampleids = 1:2, ncol = 2, value.scale = "all") +
  plot_layout(guides = "collect")


combined_e15 <- AddMetaData(
  object = combined_e15,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime_e15")
FeaturePlot(combined_e15, "Pseudotime_e15", cols = viridis(10, option = "magma"), split.by = "mode")

saveRDS(cortex, "e15_ctx_subset.rds")
saveRDS(combined_e15, "coembed_e15_atac_sc.rds")


