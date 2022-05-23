source("source_satac.R")

# set paths
dirs <- list.dirs("data", recursive = F, full.names = F)
table <- list()

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

# make new matrices with encode peaks
# load peaks
# load peaks
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

# create common peak set
for(i in seq_along(dirs)){
  folder = strsplit(infoTable$samples[i], "raw_peak_bc_matrix.h5") %>% unlist()
  # load metadata
  object$md[[i]] <- read.table(
    file = infoTable$singlecell[i],
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
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

saveRDS(object, "results/samples_unmerged.rds")
#object <- readRDS("results/samples_unmerged.rds")


# make tissue files
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
  
  object$mtx[[i]] <- mtx
  
  #save spot info in a list
  tissue_md[[i]] <- spotfile
}

signac_object <- list()

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
write.csv(mtx, "combined_q0_peaks.csv")

#run dca on terminal
######## 
# dca combined_q0_peaks.csv dca_peaks_q0 \
#   --threads 3 \
#   --nosizefactors --nonorminput --nologinput --nocheckcounts \
#   --saveweights
#
#########
