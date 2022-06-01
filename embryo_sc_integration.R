#
source("source_satac.R")

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
DefaultAssay(obj_dev) <- "RNA"

Idents(cortex) <- "sample"
cortex$sample_section <- cortex$sample
table(cortex$sample)
#cortex <- SubsetSTData(cortex, expression = sample %in% c("1", "2", "7", "8", "9", "10"))
#table(cortex$sample)

obj_dev <- FindVariableFeatures(obj_dev)
cortex <- FindVariableFeatures(cortex, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cortex)
cortex <- ScaleData(cortex, vars.to.regress = "passed_filters")
cortex <- RunPCA(cortex, features = VariableFeatures(object = cortex))

cortex$Age <- "e15_ATAC"

combined_pseudot <- merge(obj_dev, cortex)
obj_list <- SplitObject(combined_pseudot, split.by = "Age")
rm(combined_pseudot)
rm(obj_dev)



