#QC metrics 

#ArchR
dirs <- list.dirs("data", recursive = F, full.names = F)
table <- list()
arrow <- list()
for(i in dirs){
  table[[i]] <- c(barcodes = paste0("meta/", i, "_tissue.tsv")$barcode,
                  fragments = paste0("data/", i, "/outs/fragments.tsv.gz"))
  
  arrow[[i]] <- createArrowFiles(
    inputFiles = table[[i]]$fragments,
    sampleNames = i,
    validBarcodes = table[[i]]$barcodes,
    geneAnnotation = getGeneAnnotation(),
    genomeAnnotation = getGenomeAnnotation(),
    minTSS = 1,
    minFrags = 100,
    maxFrags = 1e+05,
  )
  
  archr_obj[[i]] <- ArchRProject(
    ArrowFiles = arrow[[i]], 
    outputDirectory = "ArchR",
    copyArrows = FALSE 
  )
  
}

plots_tss <- list()
plots_frag <- list()
for(i in dirs){
  plots_tss[[i]] <- plotTSSEnrichment(ArchRProj = archr_obj[[i]]) + ylim(0,8)
  plots_frag[[i]] <- plotFragmentSizes(archr_obj[[i]]) + ylim(0,1)
}

ggarrange(plotlist = plots_tss)
ggarrange(plotlist = plots_frag)


# QC plots using signac
embryo_den$pct_mito <- embryo_den$mitochondrial / embryo_den$total 
embryo_den$pct_TSS <- embryo_den$TSS_fragments / embryo_den$passed_filters 

VlnPlot(embryo_den, 'pct_TSS', group.by = 'section_id', cols = section_cols) + ylim(0, 0.5) 
VlnPlot(embryo_den, 'pct_mito', group.by = 'section_id', cols = section_cols) + ylim(0, 0.5)
embryo_den$unique_log <- log(embryo_den$passed_filters, base = 10)
embryo_den$tss_log <- log(embryo_den$TSS_fragments, base = 10)
VlnPlot(embryo_den,c('unique_log'), pt.size = 0.01, group.by = 'cluster_ids') & scale_fill_manual(values = graf_bright_re1) & ylim(0,6)  
VlnPlot(embryo_den,c('tss_log'), pt.size = 0.01, group.by = 'cluster_ids') & scale_fill_manual(values = graf_bright_re1) & ylim(0,6) 


# Combine with 10X snATACseq data, E18 brain
counts <- Read10X_h5(filename = "data/atac_v1_E18_brain_flash_5k_filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "data/atac_v1_E18_brain_flash_5k_singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = 'data/atac_v1_E18_brain_flash_5k_fragments.tsv.gz',
  min.cells = 1,
  min.features = 1
)

E18_10X <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

E18_10X$pct_mito <- E18_10X$mitochondrial / E18_10X$total 
E18_10X$pct_TSS <- E18_10X$TSS_fragments / E18_10X$passed_filters 

VlnPlot(E18_10X, 'pct_TSS') + ylim(0, 0.5)
VlnPlot(E18_10X, 'pct_mito') + ylim(0, 0.5)

E18_10X$unique_log <- log(E18_10X$passed_filters, base = 10)
E18_10X$tss_log <- log(E18_10X$TSS_fragments, base = 10)
VlnPlot(E18_10X, c("unique_log")) + ylim(0,6) 
VlnPlot(E18_10X, c("tss_log")) + ylim(0,6) 

#distribution plots across sections
df <- table(combined$section, combined$peaks_snn_res.0.7) %>% as.data.frame()
prop_df <- prop.table(table(combined$section, combined$peaks_snn_res.0.7), 1) %>% as.data.frame()
levels(prop_df$Var1) <- c("220403_C2", "220403_D2", "220403_A2","220403_B2","220327_A1", "220327_A2") %>% rev()
levels(prop_df$Var2) <- rev(my_levels)

pdf("results/barplot_clacrosssections.pdf", height = 4, width = 8)
ggplot(prop_df, aes(fill=Var2, y=Freq, x=Var1)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=cols) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank())
dev.off()


# peak calling and genome binning
DefaultAssay(combined) <- "peaks"

# call peaks on the merged dataset using MACS2
peaks_all <- CallPeaks(
  object = combined,
  macs2.path = "~/opt/anaconda3/envs/atac_env/bin/macs2"
)

# check overlap between peaks called with MACS2 and ENCODE peaks
overlap <- findOverlaps(combined@assays$peaks@ranges, peaks_all)
overlaps <- findOverlaps(peaks_all, combined@assays$peaks@ranges, select = "arbitrary") 
table(is.na(overlaps))

# create feature matrices with MACS2 peaks
object <- list()
for(i in seq_along(new.path)){
  object$frag[[i]] <- CreateFragmentObject(new.path[i])
  
  object$counts_macs[[i]] <- FeatureMatrix(
    object$frag[[i]],
    peaks_all,
    verbose = TRUE
  )
  
  object$assay_macs[[i]] <- CreateChromatinAssay(object$counts_macs[[i]], 
                                                 fragments = object$frag[[i]])
  
  object$signac_macs[[i]] <- CreateSeuratObject(object$assay_macs[[i]], 
                                                assay = "peaks")
  
  object$signac_macs[[i]]$sample <- i
  
}

# merge objects and process data same way as before
combined_macs <- merge(object$signac_macs[[1]], c(object$signac_macs[[2]],
                                                  object$signac_macs[[3]],
                                                  object$signac_macs[[4]],
                                                  object$signac_macs[[5]],
                                                  object$signac_macs[[6]]))

combined_macs <- subset(combined_macs, cells = colnames(combined))
combined_macs@tools <- combined@tools
ST.FeaturePlot(combined_macs, "nCount_peaks")

combined_macs <- combined_macs %>%
  RunTFIDF() %>%
  FindTopFeatures(min.cutoff = 'q0') %>%
  RunSVD() %>% 
  RunHarmony(group.by.vars = "sample", 
             reduction = "lsi", dims.use = 1:7, 
             assay.use = "peaks", 
             project.dim = F,
             verbose = T) %>%
  RunUMAP(reduction = "harmony", dims = 1:7, reduction.name = "umap.harmony") %>%
  FindNeighbors(reduction = "harmony", dims = 1:7) %>%
  FindClusters(resolution = 0.7)
combined_macs$seurat_clusters_harmony <- combined_macs$seurat_clusters
levels(combined) <- c("5", "7", "2","9","3", "4", "8","0","1","10","6")
cols <-  c("#4477AA","#58AAD2","#58BEC8", "#2F9558","#669C39","#CCBB44", "#A8A8A8","#E05B77","#B73D77", "#B06992","#E08762") #embryo clusters res 0.7
cols_macs <- c("#4477AA","#58AAD2","#669C39", "#2F9558","#CCBB44", "#A8A8A8","#E05B77","#B73D77","#B27699", "#B06992","#E08762")

Idents(combined_macs) <- "seurat_clusters_harmony"
levels(combined_macs) <- c("4", "1", "2","3","6", "7", "9","0","5","10","8")

# compare clusters obtained with MACS2 and ENCODE peaks
p1 <- ST.FeaturePlot(combined_macs, "ident", ncol = 2, pt.size = 0.7, cols = cols_macs) +
  ggtitle("MACS2 peaks")
p2 <- ST.FeaturePlot(combined, "ident", ncol = 2, pt.size = 0.7, cols = cols) +
  ggtitle("ENCODE peaks")
ggpubr::ggarrange(p1, p2)

combined <- AddMetaData(combined, combined_macs@active.ident, col.name = "macs_cl")
prop <- prop.table(table(combined@active.ident, combined$macs_cl), 1)*100
heatmap.2(prop[order(nrow(prop):1),], 
          Colv = NA, Rowv = NA, scale="none", 
          xlab="peaks cluster", ylab="dca clusters", 
          col = hm_colors, RowSideColors=rev(cols), ColSideColors = cols_macs)



## compare metrics with other spatial ATAC method (Deng et al.)
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171943

# make new matrices with encode peaks
# load peaks
path_frag <- "~/Documents/sequencing/sATAC/Deng/GSE171943_RAW/embryo/fragments/"
fragments <- list.files(path_frag, pattern = "\\.tsv.gz$")

object <- list()
# create common peak set
for(i in seq_along(fragments)){
  frag_name = strsplit(fragments[i], ".fragments.tsv.gz") %>% unlist()
  
  # create fragment objects
  object$frag[[frag_name]] <- CreateFragmentObject(
    path = paste0(path_frag, fragments[i])
  )
  
  # make count matrix
  object$counts[[frag_name]] <- FeatureMatrix(
    fragments = object$frag[[frag_name]],
    features = gr
  )
}

# figure out number of unique fragments - 
# do all versions - even those without the spatial file (some spots won't be on tissue tho)
object_nospatial <- list()
for(i in seq_along(fragments)){
  frag_name = strsplit(fragments[i], ".fragments.tsv.gz") %>% unlist()
  
  frag <- read.table(paste0(path_frag, fragments[i]))
  passed_filters <- table(frag[,4]) %>% as.data.frame()
  colnames(passed_filters) <- c("barcode", "passed_filter")
  total <- aggregate(frag$V5, by=list(barcode=frag$V4), FUN=sum)
  colnames(total) <- c("barcode", "total")
  md <- merge(passed_filters, total, by = "barcode")
  rownames(md) <- md$barcode
  
  object_nospatial$md[[frag_name]] <- md
  
  chrom_assay <- CreateChromatinAssay(
    counts = object$counts[[frag_name]],
    sep = c(":", "-"),
    genome = 'mm10',
    fragments = object$frag[[frag_name]]
  )
  
  object_nospatial$obj[[frag_name]] <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks_",
    meta.data = object_nospatial$md[[frag_name]][colnames(chrom_assay),])
  
  object_nospatial$obj[[frag_name]]$sample <- frag_name

  sig_obj <- object_nospatial$obj[[frag_name]] %>%
    RenameCells(new.names = paste0(colnames(sig_obj), "_", frag_name))
  object_nospatial$obj[[frag_name]] <- sig_obj
}

satac <- merge(object_nospatial$obj$GSM5238385_ME11_50um, 
               y = object_nospatial$obj[!names(object_nospatial$obj) %in% "GSM5238385_ME11_50um"] )
satac$logunique <- log10(satac$passed_filter)
satac$passed_filters <- satac$passed_filter

Idents(satac) <- "V2"
satac <- subset(satac, ident = "1")
satac_all <- merge(combined, satac)
satac_all$logunique <- log10(satac_all$passed_filters)
VlnPlot(satac_all, "logunique", group.by = "sample", pt.size = 0) + NoLegend()





