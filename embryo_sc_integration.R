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

