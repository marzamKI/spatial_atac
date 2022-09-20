source("source_satac.R")

# set paths to raw data and metadata
dirs <- list.dirs("data/brca", recursive = F, full.names = F)
table <- list()

# create a table containing paths to raw and meta data (i.e., output from cellranger + spatial info)
for(i in dirs){
  table[[i]] <- c(samples = paste0("data/brca/", i, "/raw_peak_bc_matrix.h5"),
                  singlecell = paste0("data/brca/", i, "/singlecell.csv"),
                  fragments = paste0("data/brca/", i, "/fragments.tsv.gz"),
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
gr <- read.csv(
  file = "BRCA_peakCalls.csv",
  col.names = c("chr", "start", "end"),
  sep = ";"
) %>% makeGRangesFromDataFrame()

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
  
  if (file.exists(paste0(folder, "brca_peak_bc_matrix.h5"))) {
    message("File exists.")
    object$counts[[i]] <- Read10X_h5(paste0(folder, "brca_peak_bc_matrix.h5"))
    
  } else {
    
    # make count matrix
    object$counts[[i]] <- FeatureMatrix(
      fragments = object$frag[[i]],
      features = gr,
      cells = rownames(object$md[[i]])
    )
    
    write10xCounts(path = paste0(folder, "brca_peak_bc_matrix.h5"),
                   x = object$counts[[i]],
                   type = "HDF5") # save new count matrix
  }
  
}

# save unmerged objects and metadata
saveRDS(object, "data/brca/brca_unmerged.rds")

# make tissue position files
for(i in seq_along(infoTable$tissue_paths)){
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
            paste0(strsplit(infoTable$tissue_paths[i], "tsv"), "csv"), 
            sep = ",")
  
}

# load spotfiles
tissue_md <- list()
object$mtx <- list()

for(i in seq_along(infoTable$spotfiles)){
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
for(i in seq_along(infoTable$spotfiles)){
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
}

# merge objects
combined <- merge(
  x = signac_object[[1]],
  y = c(signac_object[2:3])
)

table(combined$section)

# add spatial data to meta
for(i in seq_along(infoTable$tissue_paths)){
  tissue_md[[i]]$barcode <- paste0(tissue_md[[i]]$barcode, "_", i)
  rownames(tissue_md[[i]]) <- tissue_md[[i]]$barcode 
}

tissue_md_combined <- do.call("rbind", tissue_md) 
combined <- AddMetaData(combined, tissue_md_combined)

# add image
# image is manually cropped so that it only shows the capture area
table <- list()
for(i in dirs){
  table[[i]] <- c(samples = paste0("data/brca", i, "/brca_peak_bc_matrix.h5"),
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
combined@tools[["Staffli"]] <- se@tools$Staffli # create STutility object for spatial plots

#preprocess raw data
combined <- combined %>% 
  RunTFIDF() %>%
  FindTopFeatures(min.cutoff = 'q0') %>%
  RunSVD() %>%
  RunUMAP(reduction = 'lsi', dims = 2:7) %>%
  FindNeighbors(reduction = 'lsi', dims = 2:7) %>%
  FindClusters(algorithm = 3, resolution = 0.5)

#save matrix for denoising
write.csv(combined@assays$peaks@data,
          "data/combined_brca_q0_peak_bc_matrix.csv")

#run DCA on terminal
# dca for peaks
######## 
# dca combined_brca_q0_peak_bc_matrix.csv dca_peaks_q0_brca \
#  --threads 3 \
#  --nosizefactors --nonorminput --nologinput --nocheckcounts \
#  --saveweights
#
#########

# load denoised matrices and save objects - separately for peaks and gene activity due to large size
# peaks
denoised_counts <- read.table("dca_peaks_q0_brca/mean.tsv", row.names = 1) %>% 
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

## compare DCA clustering
# clustering between DCA and peaks
# no need for harmony integration as the sections come from the same sample
DefaultAssay(combined) <- "dca"
combined <- RunTFIDF(combined) %>%
  FindTopFeatures(min.cutoff = 'q0') %>%
  RunSVD() %>%
  RunUMAP(dims = 2:9) %>%
  FindNeighbors(dims = 2:9) %>%
  FindClusters(resolution = 0.5)
combined$dca_snn_res.0.5 <- combined$seurat_clusters

levels(combined$peaks_snn_res.0.5) <- c("2", "4", "1", "3", "0")
levels(combined$dca_snn_res.0.5) <- c("2", "3", "0", "4", "1")

# Check the proportion of spots assigned to each cluster when using denoised vs original data
prop <- prop.table(table(combined$peaks_snn_res.0.5, combined$dca_snn_res.0.5), 1)*100
heatmap(prop[order(nrow(prop):1),], 
        Colv = NA, Rowv = NA, scale="none", 
        xlab="peaks cluster", ylab="dca clusters", 
        col = hm_colors, RowSideColors=rev(colors_okate), ColSideColors = colors_okate)

saveRDS(combined, "results/combined_brca_denoised_peaks.rds")

# calculate gene activity
DefaultAssay(combined) <- "peaks"
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(combined) <- annotations
gene.activities <- GeneActivity(combined)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
#remove PCDH and UGT genes
gene.activities <- gene.activities[-grep("PCDH", rownames(gene.activities)),]
gene.activities <- gene.activities[-grep("UGT", rownames(gene.activities)),]
combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
write.csv(combined@assays$RNA@counts, "data/gene_activity_counts_brca.csv")

#run DCA on terminal
# dca for gene activity
######## 
# dca gene_activity_counts_brca.csv dca_gene_activity_brca \
#  --threads 3 --saveweights
#
#########

# gene activity
combined[["dca"]] <- NULL
mtx <- read.table("dca_gene_activity_brca/mean.tsv")
colnames(mtx) <- gsub("\\.", "-", colnames(mtx))
DefaultAssay(combined) <- "RNA"
combined <- subset(combined, cells = colnames(mtx))

# add denoised matrix to object
combined[["RNA_dca"]] <- CreateAssayObject(counts = as.matrix(mtx))
rm(mtx)
saveRDS(combined, "results/combined_brca_denoised_rna.rds")





