combined <- readRDS("combined_q0_denoised_peaks.rds")
DefaultAssay(combined) <- "dca"
combined <- FindTopFeatures(combined, min.cutoff = "q75")
combined <- ScaleData(combined, features = VariableFeatures(combined))
spatgenes_dca <- CorSpatialGenes(combined, assay = "dca")

combined <- RunNMF_fix(combined, nfactors = 8, rescale = T, assay = "dca", slot = "scale.data", n.cores = 100) # Specificy nfactors to choose the number of factors, default=20.
saveRDS(combined, "combined_nmf_peaks.rds")

cscale <- c("lightgray", "mistyrose", "red", "darkred", "black")

for(i in 1:8){
  print(
    ST.DimPlot(combined, 
               dims = i,
               ncol = 2, # Sets the number of columns at dimensions level
               grid.ncol = 1, # Sets the number of columns at sample level
               reduction = "NMF", 
               pt.size = 1, 
               center.zero = F, 
               cols = cscale, 
               show.sb = FALSE)
  )
}
dev.off()

plots <- list()
for(i in 1:8){
  plots[[i]] <- 
    ST.DimPlot(combined, 
               dims = i,
               ncol = 2, # Sets the number of columns at dimensions level
               grid.ncol = 1, # Sets the number of columns at sample level
               reduction = "NMF", 
               pt.size = 0.3, 
               center.zero = F, 
               cols = cscale, 
               show.sb = FALSE)
  
}
ggarrange(plotlist = plots, ncol = 4, nrow = 2)

## nmf compare
my_levels <- c("5", "7", "2","9","3", "4", "8","0","1","10","6")
levels(combined$peaks_snn_res.0.7) <- my_levels
fun_color_range <- colorRampPalette(c("gray93", "#AE017E"))  # Create color generating function
my_colors <- fun_color_range(100)                         
combined$peaks_id <- combined$peaks_snn_res.0.7

combined <- RunHarmony(combined, 
                       group.by.vars = "section", 
                       reduction = "NMF", dims.use = 1:8, 
                       assay.use = "dca", 
                       project.dim = F,
                       verbose = T) %>%
  RunUMAP(reduction = "harmony", dims = 1:8, reduction.name = "umap.harmony") %>%
  FindNeighbors(reduction = "harmony", dims = 1:8) %>%
  FindClusters(resolution = 0.7)
combined$seurat_clusters_harmony <- combined$seurat_clusters
p2 <- ST.FeaturePlot(object = combined, features = "peaks_id", pt.size = 1, ncol = 6) + ggtitle("peaks clusters")
p1 <- ST.FeaturePlot(object = combined, features = "seurat_clusters_harmony", pt.size = 1, ncol = 6) + ggtitle(ncol(combined@reductions$NMF))
print(p1+p2)

my_levels <- c("5", "7", "2","9","3", "4", "0","1","10","6", "8")
Idents(combined) <- "peaks_snn_res.0.7"
levels(combined) <- my_levels
combined$peaks_snn_res.0.7 <- Idents(combined)

my_levels <- c("3", "0", "6", "10", "4", "2","9",   "1",  "11", "5", "12","7","8")
Idents(combined) <- "seurat_clusters_harmony"
levels(combined) <- my_levels
combined$seurat_clusters_harmony <- Idents(combined)

cols_nmf <-c("#4477AA","#55A1CC", "#66CCEE",
             "#43AA90", "#228833", "#76A13B",
             "#CCBB44",  "#BBBBBB",
             "#EE6677", "#CB4C77", "#AA3377", "#B27699","#DD905D")

cols <-c("#4477AA","#55A1CC", "#66CCEE",
         "#228833", "#76A13B", 
         "#CCBB44",  "#BBBBBB",
         "#CB4C77", "#AA3377", "#B27699","#DD905D") 
prop <- prop.table(table(combined$peaks_snn_res.0.7, combined$seurat_clusters_harmony), 1)*100
heatmap(prop[order(nrow(prop):1),], 
        Colv = NA, Rowv = NA, scale="none", 
        xlab="peaks cluster", ylab="nmf clusters", 
        col = my_colors, ColSideColors = cols_nmf,
        RowSideColors = cols)

## weights peaks 
head(combined@reductions$NMF@feature.loadings)
FactorGeneLoadingPlot(combined, factor = 4)
top_20 <- list()
for(i in c(2,4,5,7)){
  ftr <- paste0("factor_", i)
  nmf <- combined@reductions$NMF@feature.loadings[, ftr]
  gene <- names(nmf)
  df <- data.frame(gene, val = nmf, stringsAsFactors = F)
  df <- df[order(df$val, decreasing = T), ]
  df <- df[1:20, ]
  df$gene <- factor(df$gene, levels = df$gene)
  df$factor <- ftr
  top_20[[i]] <- df
  
}

DefaultAssay(combined) <- "dca"
p1 <- ST.FeaturePlot(combined, features =as.character(top_20[[4]]$gene[5]), ncol = 3, indices = c(1,3,5), pt.size = 0.7, cols = magenta_scale)
p2 <- ST.FeaturePlot(combined, features =as.character(top_20[[5]]$gene[1]), ncol = 3, indices = c(1,3,5), pt.size = 0.7, cols = magenta_scale)
p3 <- ST.FeaturePlot(combined, features =as.character(top_20[[7]]$gene[10]), ncol = 3, indices = c(1,3,5), pt.size = 0.7, cols = magenta_scale)
p4 <- ST.FeaturePlot(combined, features =as.character(top_20[[2]]$gene[1]), ncol = 3, indices = c(1,3,5), pt.size = 0.7, cols = magenta_scale)
p1+p2+p3+p4

