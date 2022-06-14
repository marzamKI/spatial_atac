source("source_satac.R")

# set paths
dirs <- list.dirs("data/brca", recursive = F, full.names = F)
table <- list()

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


# make new matrices with encode peaks
# load peaks
gr <- read.csv(
  file = "BRCA_peakCalls.csv",
  col.names = c("chr", "start", "end"),
  sep = ";"
) %>% makeGRangesFromDataFrame()

# create common peak set
for(i in 1:3){
  
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
  
  # make count matrix
  object$counts[[i]] <- FeatureMatrix(
    fragments = object$frag[[i]],
    features = gr,
    cells = rownames(object$md[[i]])
  )
  
  # save sparse matrix
  folder = strsplit(infoTable$samples[i], "raw_peak_bc_matrix.h5") %>% unlist()
  
  if (file.exists(paste0(folder, "brca_peak_bc_matrix.h5"))) {
    message("File exists.")
  } else {
    write10xCounts(path = paste0(folder, "brca_peak_bc_matrix.h5"),
                   x = object$counts[[i]],
                   type = "HDF5")
  }
  
}

saveRDS(object, "data/brca_unmerged.rds")

