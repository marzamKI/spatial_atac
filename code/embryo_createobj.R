source("source_satac.R")

# set paths to raw data and metadata
dirs <- list.dirs("data", recursive = F, full.names = F)
table <- list()

# create a table containing paths to raw and meta data (i.e., output from cellranger + spatial info)
for(i in dirs){
  table[[i]] <- c(samples = paste0("data/", i, "/outs/raw_peak_bc_matrix.h5"),
                  singlecell = paste0("data/", i, "/outs/singlecell.csv"),
                  fragments = paste0("data/", i, "/outs/fragments.tsv.gz"),
                  tissue_paths = paste0("meta/", i, "_tissue.tsv"),
                  spotfiles = paste0("meta/", i, "_tissue.csv"))
}

infoTable <- do.call("rbind", table) %>% as.data.frame()

object <- list(md = list(),
               frag = list(),
               counts = list())

# build new matrices with new peak set
# peaks obtained from https://www.science.org/doi/10.1126/science.aav1898
# load peaks and create genomicranges object
gr_12 <- read.table(
  file = "E12_peaks_ALL.bed",
  col.names = c("chr", "start", "end")
) %>% makeGRangesFromDataFrame()

gr_13 <- read.table(
  file = "ALL_E13_PEAKS.bed",
  col.names = c("chr", "start", "end")
) %>% makeGRangesFromDataFrame()

gr_15 <- read.table(
  file = "E15_ALL.bed",
  col.names = c("chr", "start", "end")
) %>% makeGRangesFromDataFrame()

gr <- GenomicRanges::reduce(c(gr_12, gr_13, gr_15)) 

# create fragment objects and count matrices for each section separately
for(i in seq_along(dirs)){
  folder = strsplit(infoTable$samples[i], "raw_peak_bc_matrix.h5") %>% unlist()
  # load metadata
  object$md[[i]] <- read.table(
    file = infoTable$singlecell[i],
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row (NO_BARCODE)
  
  # create fragment objects
  object$frag[[i]] <- CreateFragmentObject(
    path = infoTable$fragments[i],
    cells = rownames(object$md[[i]])
  )
  
  if (file.exists(paste0(folder, "encode_peak_bc_matrix_gr12_13_15.h5"))) {
    message("File exists.")
    object$counts[[i]] <- Read10X_h5(paste0(folder, "encode_peak_bc_matrix_gr12_13_15.h5"))
    
  } else {
    
    # make count matrix
    object$counts[[i]] <- FeatureMatrix(
      fragments = object$frag[[i]],
      features = gr,
      cells = rownames(object$md[[i]])
    )
    
    write10xCounts(path = paste0(folder, "encode_peak_bc_matrix_gr12_13_15.h5"),
                   x = object$counts[[i]],
                   type = "HDF5")
  }
  
}

# save unmerged objects and metadata
saveRDS(object, "data/embryo_unmerged.rds")

# make tissue position files
for(i in seq_along(dirs)){
  spotfile <- read.table(infoTable$tissue_paths[i], sep = "\t", header = T)
  
  spotfile_tissue <- cbind(
    spotfile$row, #array_row
    spotfile$col, #array_col
    spotfile$x, #pxl_row_in_fullres
    spotfile$y #pxl_col_in_fullres
  ) %>% as.data.frame()
  tissue <- rep(1, nrow(spotfile_tissue))
  spotfile_tissue <- cbind(tissue, spotfile_tissue)
  rownames(spotfile_tissue) <- paste0(spotfile$barcode, "-1")
  
  # save new tissue position file with the correct columns for creating STutility object
  write.csv(spotfile_tissue, 
            paste0(strsplit(infoTable$tissue_paths[i], "tsv"), "csv"))
  
}

# load spotfiles
tissue_md <- list()
object$mtx <- list()

for(i in seq_along(dirs)){
  spotfile <- read.csv(infoTable$spotfiles[i],
                       col.names = c("barcode", "tissue", "y", "x", "pixel_y", "pixel_x")
  )
  
  mtx <- object$counts[[i]]
  mtx <- mtx[,colnames(mtx) %in% spotfile$barcode]
  # filter count matrix to retain only spots overlaying tissue
  
  object$mtx[[i]] <- mtx
  
  #save spot info in a list
  tissue_md[[i]] <- spotfile
}

signac_object <- list()

# create signac objects for each section and run normalization and dimensionality reduction
for(i in seq_along(dirs)){
  assay <- CreateChromatinAssay(object$mtx[[i]], 
                                fragments = object$frag[[i]])
  
  signac_object[[i]] <- CreateSeuratObject(assay, 
                                           assay = "peaks", 
                                           meta.data=object$md[[i]])
  signac_object[[i]]$section <- rownames(infoTable)[i]
  signac_object[[i]]$sample <- i
  # compute LSI
  signac_object[[i]] <- FindTopFeatures(signac_object[[i]], min.cutoff = 10)
  signac_object[[i]] <- RunTFIDF(signac_object[[i]])
  signac_object[[i]] <- RunSVD(signac_object[[i]])
  signac_object[[i]] <- RunUMAP(signac_object[[i]], reduction = 'lsi', dims = 1:7) %>%
    FindNeighbors(reduction = 'lsi', dims = 1:7) %>%
    FindClusters(algorithm = 3, resolution = 0.5)
  signac_object[[i]]$individual_clusters <- signac_object[[i]]$peaks_snn_res.0.5
}

# merge objects
combined <- merge(
  x = signac_object[[1]],
  y = c(signac_object[2:length(dirs)])
)

# add spatial data to meta
for(i in 1:length(dirs)){
  tissue_md[[i]]$barcode <- paste0(tissue_md[[i]]$barcode, "_", i)
  rownames(tissue_md[[i]]) <- tissue_md[[i]]$barcode 
}

tissue_md_combined <- do.call("rbind", tissue_md) 
combined <- AddMetaData(combined, tissue_md_combined)

# add image
# image is manually cropped so that it only shows the capture area
table <- list()
for(i in dirs){
  table[[i]] <- c(samples = paste0("data/", i, "/outs/encode_peak_bc_matrix_gr12_13_15.h5"),
                  spotfiles = paste0("meta/", i, "_tissue.csv"),
                  imgs = paste0("images/", i, "_cropped.jpg"))
}

infoTable <- do.call("rbind", table) %>% as.data.frame()

for (i in 1:length(dirs)) {
  # Read *_tissue.csv file
  xy.raw <- setNames(read.csv(file = infoTable$spotfiles[i]), 
                     nm = c("barcode", "tissue", "y", "x", "pixel_y", "pixel_x"))
  xy <- xy.raw[, c("x", "y")]
  
  img_path <- infoTable$imgs[i]
  img <- readJPEG(img_path)
  
  sf <- c(ncol(img)/128, nrow(img)/78)
  xy$x <- xy$x*sf[1]
  xy$y <- xy$y*sf[2]
  
  # Create a new spot selection table with proper image pixel coordinates which match the cropped images
  spotfile <- data.frame(xy.raw$barcode, xy.raw$tissue, xy.raw$y, xy.raw$x, round(xy$y), round(xy$x))
  write.table(spotfile, file = paste0(strsplit(infoTable$spotfiles[i], ".csv"), "_positions_list.csv"), 
              sep = ",", quote = F, row.names = F, col.names = F)
}

infoTable$spotfiles <- paste0(strsplit(infoTable$spotfiles, ".csv"), "_positions_list.csv")
se <- InputFromTable(infoTable, scaleVisium = 1)
se <- LoadImages(se, time.resolve = F)
#combined <- readRDS("data/combined_denoised_rna.rds")
combined@tools[["Staffli"]] <- se@tools$Staffli

#normalize, dim red, clustering
combined <- RunTFIDF(combined) %>%
  FindTopFeatures(min.cutoff = 'q0') %>%
  RunSVD() %>% 
  RunHarmony(group.by.vars = "section", 
             reduction = "lsi", dims.use = 1:7, 
             assay.use = "peaks", 
             project.dim = F,
             verbose = T) %>%
  RunUMAP(reduction = "harmony", dims = 1:7, reduction.name = "umap.harmony") %>%
  FindNeighbors(reduction = "harmony", dims = 1:7) %>%
  FindClusters(resolution = 0.7)
combined$seurat_clusters_harmony <- combined$seurat_clusters

# create matrix for DCA
mtx <- as_matrix(combined@assays$peaks@data)

#save matrix for denoising
write.csv(mtx, "combined_q0_peaks.csv")

#run dca on terminal
# dca for peaks
######## 
# dca combined_q0_peaks.csv dca_peaks_q0 \
#  --threads 3 \
#  --nosizefactors --nonorminput --nologinput --nocheckcounts \
#  --saveweights
#
#########

# load denoised matrices and save objects - separately for peaks and gene activity due to large size
# peaks
denoised_counts <- read.table("dca_peaks_q0/mean.tsv", row.names = 1) %>% 
  as.matrix()
colnames(denoised_counts) <- sub("\\.", "-", colnames(denoised_counts))

# add denoised matrix to object
# filter to retain same spots and features
combined[["dca"]] <- combined[["peaks"]]
combined@assays$dca@data <- denoised_counts
combined@assays$dca@counts <- combined@assays$dca@counts[rownames(combined@assays$dca@counts) %in%
                                                           rownames(combined@assays$dca@data),]
combined@assays$dca@var.features <- combined@assays$dca@var.features[combined@assays$dca@var.features %in%
                                                                       rownames(combined@assays$dca@data)]
rm(denoised_counts)
saveRDS(combined, "results/combined_denoised_peaks.rds")

# calculate gene activity
DefaultAssay(combined) <- "peaks"
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(combined) <- annotations
gene.activities <- GeneActivity(combined)
genes <- rownames(gene.activities)
genes_filter <- c(grep("Pcdh", genes), grep("PCDH", genes), 
                  grep("UGT", genes), grep("Ugt", genes))
gene.activities_filtered <- gene.activities[-genes_filter,]

# add the gene activity matrix to the Seurat object as a new assay and normalize it
#remove PCDH and UGT genes
combined[['RNA']] <- CreateAssayObject(counts = gene.activities_filtered) %>%
  NormalizeData(assay = 'RNA', normalization.method = 'LogNormalize')

#save raw counts for DCA
write.csv(combined@assays$RNA@counts, "gene_activity_counts.csv")

saveRDS(combined, "combined_lsi_q0.rds")

# dca for gene activity
######## 
# dca gene_activity_counts.csv dca_gene_activity \
#  --threads 3 --saveweights
#
#########

# gene activity
combined[["dca"]] <- NULL
mtx <- read.table("dca_gene_activity/mean.tsv")
colnames(mtx) <- gsub("\\.", "-", colnames(mtx))
DefaultAssay(combined) <- "RNA"
combined <- subset(combined, cells = colnames(mtx))

combined[["RNA_dca"]] <- CreateAssayObject(counts = as.matrix(mtx))
rm(mtx)
saveRDS(combined, "results/combined_denoised_rna.rds")
