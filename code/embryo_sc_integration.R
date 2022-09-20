<<<<<<< HEAD
# integration with snATAC-seq and scRNA-seq
source("source_satac.R")

# scATAC-seq (generated in the lab)
# load ENCODE peaks
gr_12 <- read.table(
  file = paste0(path, "E12_peaks_ALL.bed"),
  col.names = c("chr", "start", "end")
) %>% makeGRangesFromDataFrame()

gr_13 <- read.table(
  file = paste0(path, "ALL_E13_PEAKS.bed"),
  col.names = c("chr", "start", "end")
) %>% makeGRangesFromDataFrame()

gr_15 <- read.table(
  file = paste0(path, "E15_ALL.bed"),
  col.names = c("chr", "start", "end")
) %>% makeGRangesFromDataFrame()

gr <- GenomicRanges::reduce(c(gr_12, gr_13, gr_15)) 
rm(gr_12, gr_13, gr_15)

# load fragments and metadata files
setwd("~/Documents/sequencing/snatac2208/")
fragments <- list.files(pattern = "\\.tsv.gz$")
singlecell <- list.files(pattern = "\\.csv$")

object <- list()
# create objects with common peak set
for(i in seq_along(fragments)){
  frag_name = strsplit(fragments[i], "_fragments.tsv.gz") %>% unlist()
  
  # load metadata
  object$md[[frag_name]] <- read.table(
    file = singlecell[i],
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] %>%
    filter(is__cell_barcode == 1)
  
  # create fragment objects
  object$frag[[frag_name]] <- CreateFragmentObject(
    path = fragments[i],
    cells = rownames(object$md[[frag_name]])
  )
  
  # make count matrix
  object$counts[[frag_name]] <- FeatureMatrix(
    fragments = object$frag[[frag_name]],
    features = gr,
    cells = rownames(object$md[[frag_name]])
  )
  
  chrom_assay <- CreateChromatinAssay(
    counts = object$counts[[frag_name]],
    sep = c(":", "-"),
    genome = 'mm10',
    fragments = object$frag[[frag_name]]
  )
  
  object$obj[[frag_name]] <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = object$md[[frag_name]])
  
  object$obj[[frag_name]]$sample <- frag_name
  
}

# merge in unique object
snatac <- merge(object$obj$e12,
                c(object$obj$e13, object$obj$e15))

snatac$logunique <- log10(snatac$passed_filters) # log #unique fragments
#normalize, dim red, clustering
snatac <- snatac %>% 
  RunTFIDF() %>%
  FindTopFeatures(min.cutoff = 'q0') %>%
  RunSVD() %>%
  RunUMAP(reduction = 'lsi', dims = 1:10) %>%
  FindNeighbors(reduction = 'lsi', dims = 1:10) %>%
  FindClusters(algorithm = 3, resolution = 0.7)

# integrate with sATAC
satac <- readRDS("combined_lsi_q0.rds")

# integrate sATAC and scRNA datasets with seurat's CCA strategy
satac$mode <- "sATAC"
snatac$mode <- "snATAC"

DefaultAssay(satac) <- "peaks"
DefaultAssay(snatac) <- "peaks"

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(satac, snatac),
  anchor.features = rownames(satac),
  reduction = "rlsi",
  dims = 1:10
)

# merge
combined <- merge(satac, snatac) %>%
  FindTopFeatures(min.cutoff = 10)%>%
  RunTFIDF() %>%
  RunSVD() %>%
  RunUMAP(reduction = "lsi", dims = 1:10)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:10
)

# create a new UMAP using the integrated embeddings
integrated <- integrated %>%
  RunHarmony(group.by.vars = "sample", 
             reduction = "integrated_lsi", dims.use = 1:10, 
             assay.use = "peaks", 
             project.dim = F,
             verbose = T) %>%
  RunUMAP(reduction = "harmony", dims = 1:10, reduction.name = "umap.harmony") 

# cluster integrated data on umap
DefaultAssay(integrated) <- "peaks"
integrated <- integrated %>% 
  FindNeighbors(reduction = "umap.harmony", dims = 1:2) %>%
  FindClusters(resolution = 0.1, algorithm = 1)

integrated$integrated_harmony <- Idents(integrated)
integrated$cluster_mode <- paste0(integrated$integrated_harmony, 
                                  "_",
                                  integrated$mode)

# log normalize peaks data for correlation plots
Idents(integrated) <- "cluster_mode"
integrated[["peaks_lognorm"]] <- integrated[["peaks"]]

DefaultAssay(integrated) <- "peaks_lognorm"
integrated <- NormalizeData(integrated)

# calculate gene activity in the integrated object
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(integrated) <- annotations
gene.activities <- GeneActivity(integrated)
integrated[['RNA']] <- CreateAssayObject(counts = gene.activities)
integrated <- NormalizeData(
  object = integrated,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(integrated$nCount_RNA)
)

# make average expression object for peaks and RNA assays to compare signal across modalities
avg_cl_mode <- AverageExpression(integrated, return.seurat = T)

df_rna <- avg_cl_mode@assays$RNA@data %>% 
  as.data.frame()
cor_rna <- cor(df_rna, method = "spearman")

df_peaks <- avg_cl_mode@assays$peaks@data %>% 
  as.data.frame()
cor_peaks <- cor(df_peaks, method = "spearman")

# make scatterplots
cl_order <- c("5", "9", "13", "2", #cns
              "0", "11", "3", "10", #pns
              "1", "4", "19", "8", "16", "12", # face
              "14", "17", # olfactory
              "6", # muscle
              "7", "15", "18") # liver

plots <- list()
for(i in cl_order){
  satac_cl <- paste0(i, "_sATAC")
  snatac_cl <- paste0(i, "_snATAC")
  name <- as.character(i)
  plots[[name]] <-
    ggplot(df_peaks, aes_(as.name(satac_cl), as.name(snatac_cl))) +
    scattermore::geom_scattermore(pointsize = 2, alpha = 0.4) + xlab(satac_cl) + ylab(snatac_cl) + theme_bw() +
    annotate("text", y = 2, x = 2, 
             label = round(cor(df_peaks[satac_cl], df_peaks[snatac_cl],method = "spearman"), digits = 4)) + theme_bw()
  
}

ggpubr::ggarrange(plotlist = plots, ncol = 5, nrow = 4)


# do heatmap with top differentially accessible genes
DefaultAssay(integrated) <- "RNA"
markers_top <- FindAllMarkers(integrated, only.pos = T, logfc.threshold = 0.2, min.pct = 0.05) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = T) %>%
  arrange(factor(cluster, levels = cl_order))

cl_order_subset <- c("5", "9", "13", "2", "0", "3", "10", "1", "4", "8",
                     "14", "6",  "7",  "15", "18")
my_levels <- paste0(rep(cl_order_subset, each = 2), c("_sATAC", "_snATAC"))

Idents(integrated) <- "mode"
integrated_subset <- subset(integrated, idents = "sATAC", invert = F)
Idents(integrated_subset) <- "integrated_harmony"
avg_subset <- AverageExpression(integrated_subset, return.seurat = T, assays = "RNA")
levels(avg_subset) <- cl_order
avg_subset <- ScaleData(avg_subset, features = rownames(avg_subset))
p1 <- DoHeatmap(avg_subset, markers_top$gene, draw.lines = F, slot = "scale.data") +
  scale_fill_viridis(option = "magma")

Idents(integrated) <- "mode"
integrated_subset <- subset(integrated, idents = "sATAC", invert = T)
Idents(integrated_subset) <- "integrated_harmony"
avg_subset <- AverageExpression(integrated_subset, return.seurat = T, assays = "RNA")
levels(avg_subset) <- cl_order
avg_subset <- ScaleData(avg_subset, features = rownames(avg_subset))
p2 <- DoHeatmap(avg_subset, markers_top$gene, draw.lines = F, slot = "scale.data") +
  scale_fill_viridis(option = "magma")
p1 | p2



# ENCODE snATAC-seq dataset of forebrain (Preissl et al.)
setwd("~/Documents/sequencing/preissl")
fragments <- list.files(pattern = "\\.tsv.gz$")

object <- list()
# create common peak set
for(i in seq_along(fragments)){
  frag_name = strsplit(fragments[i], "_fragments.tsv.gz") %>% unlist()
  
  # create fragment objects
  object$frag[[frag_name]] <- CreateFragmentObject(
    path = fragments[i]
  )
  
  # make count matrix
  object$counts[[frag_name]] <- FeatureMatrix(
    fragments = object$frag[[frag_name]],
    features = gr
  )
  
  # figure out number of unique fragments - 
  # should be able to count how many times each spot barcode occurs in the fragment file
  frag <- read.table(fragments[i])
  passed_filters <- table(frag[,4]) %>% as.data.frame()
  colnames(passed_filters) <- c("barcode", "passed_filter")
  total <- aggregate(frag$V5, by=list(barcode=frag$V4), FUN=sum)
  colnames(total) <- c("barcode", "total")
  md <- merge(passed_filters, total, by = "barcode")
  rownames(md) <- md$barcode
  
  object$md[[frag_name]] <- md
}

fragments_rep <- c("e15.5_rep1", "e15.5_rep2",
                   "e12.5_rep1", "e12.5_rep2",
                   "e13.5_rep1", "e13.5_rep2")
fragments_df <- data.frame(x = fragments,
                           y = fragments_rep)

for(i in seq_along(fragments)){
  frag_name = strsplit(fragments[i], "_fragments.tsv.gz") %>% unlist()
  rownames(object$md[[frag_name]]) <- object$md[[frag_name]]$barcode
  
  chrom_assay <- CreateChromatinAssay(
    counts = object$counts[[frag_name]],
    sep = c(":", "-"),
    genome = 'mm10',
    fragments = object$frag[[frag_name]]
  )
  
  object$obj[[frag_name]] <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = object$md[[frag_name]])
  
  object$obj[[frag_name]]$sample <- frag_name
  object$obj[[frag_name]] <- RenameCells(object$obj[[frag_name]], fragments_df[i,2])
}

#merge objects and remove columns that are not listed as cells in xgi file
snatac <- merge(object$obj$ENCSR374QZJ_1_15.5, 
                y = object$obj[!names(object$obj) %in% "ENCSR374QZJ_1_15.5"] )

e12_cells <- read.table("GSE100033_Preissl/GSM2668118_e12.5.nchrM.merge.sel_cell.xgi.txt.gz")
e13_cells <- read.table("GSE100033_Preissl/GSM2668119_e13.5.nchrM.merge.sel_cell.xgi.txt.gz")
e15_cells <- read.table("GSE100033_Preissl/GSM2668121_e15.5.nchrM.merge.sel_cell.xgi.txt.gz")

snatac <- subset(snatac, cells = c(e12_cells$V1, e13_cells$V1, e15_cells$V1))
snatac$logunique <- log10(snatac$passed_filter)


#subset satac dataset to only contain forebrain spots (manually selected)
meta <- list()
cells <- NULL
fb_spots <- list.files("meta/")
for(i in seq_along(fb_spots)){
  meta[[i]] <- read.table(paste0("meta/",fb_spots[i]))
  meta[[i]]$barcode <- paste0(meta[[i]]$barcode, "-1_", i)
  cells <- c(cells, meta[[i]]$barcode)
}

satac <- SubsetSTData(satac, spots = cells)

# integrate sATAC and scRNA datasets with seurat's CCA strategy
satac$mode <- "sATAC"
snatac$mode <- "snATAC"

DefaultAssay(satac) <- "peaks"
DefaultAssay(snatac) <- "peaks"

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(satac, snatac),
  anchor.features = rownames(satac),
  reduction = "rlsi",
  dims = 1:15
)

# merge
combined <- merge(satac, snatac)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:15
)

# create a new UMAP using the integrated embeddings
integrated <- integrated %>%
  RunHarmony(group.by.vars = "sample", 
             reduction = "integrated_lsi", dims.use = 1:15, 
             assay.use = "peaks", 
             project.dim = F,
             verbose = T) %>%
  RunUMAP(reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony") %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.1)
integrated$integrated_harmony <- Idents(integrated)

# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = snatac,
  query = satac,
  reference.reduction = "lsi", 
  reduction = "lsiproject",
  dims = 1:15
)

snatac$query_clusters <- Idents(snatac)
predictions <- TransferData(
  anchorset = transfer.anchors, 
  refdata = snatac$query_clusters,
  weight.reduction = "lsiproject",
  dims = 1:15
)

satac <- AddMetaData(
  object = satac,
  metadata = predictions
)


=======
source("source_satac.R")

# compare sATAC with scRNAseq - 
# select dorsal forebrain neural progenitors and neurons from the transcriptomics data 
# and integrate with our sATAC spots 
>>>>>>> b4fc63364f096ce1ab6698b8d47bfc19f222f4ac
#scRNAseq from http://mousebrain.org/development/downloads.html

#sc preprocessing
lfile <- connect(filename = "~/Downloads/dev_all.loom", mode = "r+", skip.validate = TRUE)
names(lfile$col.attrs)

cells <- lfile$col.attrs$CellID[]
regions <- lfile$col.attrs$Tissue[] %in% "ForebrainDorsal" # select dorsal forebrain
age <- lfile$col.attrs$Age[] %in% c("e15.0", "e15.5") # look at E15
cell_types <- lfile$col.attrs$Class[] %in% c("Neuroblast", "Neuron", "Radial glia") # only include neurogenesis-related cells
unique_cells <- !cells %in% cells[duplicated(cells)]
data.neurons <- lfile[["matrix"]][unique_cells & regions & cell_types & age, ]

# get metadata
attrs <- names(lfile$col.attrs)
attrs <- attrs[c(1, 3:11, 26:37)] # keep relevant metadata (e.g., cluster IDs)

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

# filter cells not implicated in cortical neurogenesis
classes <- c("Cortical or hippocampal glutamatergic", "Forebrain",
             "Forebrain glutamatergic", "Neuronal intermediate progenitor",
             "Dorsal forebrain")
Idents(obj_dev) <- "Subclass"
obj_dev <- subset(obj_dev, idents = classes, invert = F) # keep above IDs

lfile$close_all()
rm(lfile)
rm(data.neurons)

# preprocess single cell data
obj_dev <- obj_dev %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(obj_dev, features = rownames(obj_dev)) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.1) %>%
  RunUMAP(dims = 1:10)

# rename clusters
new.cluster.ids <- c("Radial glia", "Neuroblasts", "Neurons", "Cycling neuroblasts")
names(new.cluster.ids) <- levels(obj_dev)
obj_dev <- RenameIdents(obj_dev, new.cluster.ids)
obj_dev$new_clusters <- Idents(obj_dev)

# save subset object
saveRDS(obj_dev, "data/dev_e15.rds")

# randomly sample cells from the object to match the size of our dataset better
set.seed(123)
obj_dev <- obj_dev[, sample(colnames(obj_dev), size =1500, replace=F)]

# subset sATAC data to only include E15 data from the dorsal forebrain
combined <- readRDS("combined_denoised_rna_nmf.rds") 

# manually select spots overlaying cortex
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

# preprocess cortex data
DefaultAssay(cortex) <- "RNA_dca"
cortex <- cortex %>% 
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

# integrate sATAC and scRNA datasets with seurat's CCA strategy
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

# dim red and clustering of the merged dataset
combined_e15 <- combined_e15 %>%
  ScaleData() %>% 
  RunPCA() %>%
  RunUMAP(dims = 1:7) %>%
  FindNeighbors(dims = 1:7) %>%
  FindClusters(resolution = 0.1)
DimPlot(combined_e15, reduction = "umap", split.by = "mode", ncol = 2)

# run pseudotemporal analysis using monocle3
rgs <- rownames(combined_e15@meta.data)[combined_e15$integrated_snn_res.0.1 == "3"] # set radial glia cluster as root
cds <- as.cell_data_set(combined_e15) %>%
  cluster_cells(reduction_method = "UMAP") %>%
  learn_graph(use_partition = T) %>%
  order_cells(reduction_method = "UMAP", root_cells = rgs)

# plot trajectories colored by pseudotime
plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)

# add monocle's pseudotime ordering to sATAC object and plot on spatial
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

# save subset objects with pseudotime ordering metadata
saveRDS(cortex, "e15_ctx_subset.rds")
saveRDS(combined_e15, "coembed_e15_atac_sc.rds")


