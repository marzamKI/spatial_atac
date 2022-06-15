source("source_satac.R")

#correlation peaks
combined <- readRDS("combined_lsi_q0.rds")
Idents(combined) <- "peaks_snn_res.0.7"
combined <- RunTFIDF(combined) %>%
  FindTopFeatures(min.cutoff = 'q75') %>%
  RunSVD()

# peaks in promoter regions vs distal elements
DefaultAssay(combined) <- "peaks"
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'

tss.positions <- GetTSSPositions(ranges = annotations)
tss.positions <- Extend(
  x = tss.positions,
  upstream = 1000,
  downstream = 100,
  from.midpoint = TRUE
)

closest_tss <- ClosestFeature(combined, rownames(combined), annotation = tss.positions)
tss_peaks <- closest_tss %>% filter(distance == 0)
distal_peaks <- closest_tss %>% filter(distance > 0)

col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256) %>% rev()

avg_tss <- AverageExpression(combined, assays = "peaks", slot = "data", features = intersect(tss_peaks$query_region, VariableFeatures(combined))) #12886 peaks
cor_mat <- cor(avg_tss$peaks, method = "spearman")

pdf("promoter_peaks_corplot_q75_12886peaks.pdf", height = 7, width = 7)
heatmap.2(cor_mat, dendrogram="none", Rowv=FALSE, Colv=FALSE,
          col = col, scale="none", key=TRUE, density.info="none",
          trace="none", symkey=T,symbreaks=T,
          ColSideColors = cols)
dev.off()

avg_distal <- AverageExpression(combined, assays = "peaks", slot = "data", features = intersect(distal_peaks$query_region, VariableFeatures(combined))) #49768 features
cor_mat <- cor(avg_distal$peaks, method = "spearman")

pdf("distal_peaks_corplot_q75_49768peaks.pdf", height = 7, width = 7)
heatmap.2(cor_mat, dendrogram="none", Rowv=FALSE, Colv=FALSE,
          col = col, scale="none", key=TRUE, density.info="none",
          trace="none", symkey=T,symbreaks=T,
          ColSideColors = cols)
dev.off()



# nmf
DefaultAssay(combined) <- "peaks"
combined <- FindTopFeatures(combined, min.cutoff = "q75")
combined <- ScaleData(combined, features = VariableFeatures(combined))
spatgenes_dca <- CorSpatialGenes(combined, assay = "peaks")

combined <- RunNMF_fix(combined, nfactors = 8, rescale = T, assay = "peaks", slot = "scale.data", n.cores = 100) # Specificy nfactors to choose the number of factors, default=20.
saveRDS(combined, "combined_nmf8_peaks.rds")

for(i in 1:8){
  print(
    ST.DimPlot(combined, 
               dims = i,
               ncol = 2, # Sets the number of columns at dimensions level
               grid.ncol = 1, # Sets the number of columns at sample level
               reduction = "NMF", 
               center.zero = F, 
               cols = magenta_scale, 
               show.sb = FALSE)
  )
}

combined <- RunHarmony(combined, 
                       group.by.vars = "section", 
                       reduction = "NMF", dims.use = 1:8, 
                       assay.use = "peaks", 
                       project.dim = F,
                       verbose = T) %>%
  RunUMAP(reduction = "harmony", dims = 1:8, reduction.name = "umap.harmony") %>%
  FindNeighbors(reduction = "harmony", dims = 1:8) %>%
  FindClusters(resolution = 0.7)
combined$seurat_nmf_harmony <- combined$seurat_clusters

Idents(combined) <- "seurat_nmf_harmony"
levels(combined) <- c("3", "2", "9", "5", "6", "1","0","8","7", "10", "4","11")
combined$seurat_nmf_harmony <- Idents(combined)

prop <- prop.table(table(combined$peaks_id, combined$seurat_nmf_harmony), 1)*100
heatmap(prop[order(nrow(prop):1),], 
        Colv = NA, Rowv = NA, scale="none", 
        xlab="peaks cluster", ylab="nmf clusters", 
        col = hm_colors, ColSideColors = cols_nmf,
        RowSideColors = cols)

p1 <- ST.FeaturePlot(combined, "ident", indices = c("1", "3", "5"), ncol = 3, cols = cols)
p2 <- ST.FeaturePlot(combined, "seurat_nmf_harmony", indices = c("1", "3", "5"), ncol = 3, cols = cols_nmf)
p1+p2

## weights peaks 
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



## Gene ontology
DefaultAssay(combined) <- "RNA"
markers <- FindAllMarkers(combined, only.pos = T, logfc.threshold = 0.1, min.pct = 0.05)

ST.FeaturePlot(combined, "ident", indices = 6, split.labels = T, ncol = 3)
my_levels <- c("5", "7", "2","9","3", "4", "8","0","1","10","6")
go <- list()
for(i in my_levels){
  fct_genes <- markers %>%
    filter(cluster == i) %>%
    select(gene)
  
  gostres <- gost(query = fct_genes$gene,
                  organism = "mmusculus")
  gostres$result$log_p <- -log10(gostres$result$p_value)
  gostres$result$cluster <- i
  
  go[[i]] <- gostres$result
}

go_df <- do.call(rbind, go)

table(go_df[go_df$source == "GO:BP", colnames(go_df) == "cluster"])

top_go <- go_df %>% 
  filter(source == "GO:BP") %>%
  filter(term_size < 2927) %>% #mean term size - to avoid very general terms
  group_by(cluster) %>%
  top_n(5, log_p) %>%
  arrange(factor(cluster, levels = my_levels))
top_go$fct_term_name <- paste0(top_go$cluster, "_", top_go$term_name)
top_go$fct_term_name <- factor(top_go$fct_term_name, levels = top_go$fct_term_name)

# Lollipop - horizontal version
ggplot(top_go, aes(x=fct_term_name, y=log_p)) +
  geom_segment( aes(x=fct_term_name, xend=fct_term_name, y=0, yend=log_p, color = cluster)) +
  geom_point( size=2, alpha=1, aes(color = cluster)) +
  scale_x_discrete(limits=rev) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  NoLegend()

# motif analysis
# MOTIF ANALYSIS PER CLUSTER
# Curated motif lists and archetypes were obtained from https://www.vierstra.org/resources/motif_clustering and RDS file downloaded from  https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/Vierstra_Individual_Motifs.rds  or https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/Vierstra_Archetype_Motifs_v2.1.rds
Vierstra <- readRDS("/Users/enric.llorens/Downloads/Vierstra_Individual_Motifs.rds")
# add motif information
DefaultAssay(combined) <- 'peaks'
combined <- AddMotifs(
  object = combined,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = Vierstra
)

# Identify peaks linked to differentially accessible genes based on co-accessibility and example with cluster 5:
combined <- LinkPeaks(
  combined,
  peak.assay = "peaks",
  expression.assay = "RNA", min.cells = 3, score_cutoff = 0.05,
  genes.use = rownames(da_genes_clust_01), min.distance = 1000
)

c5_genes <- read.csv2('/.../cl5_da_genes.csv', 
                      header = FALSE )$V1
linked_c5 <- GetLinkedPeaks(combined, c5_genes, min.abs.score = 0.1, assay = 'peaks')

open.peaks <- AccessiblePeaks(combined, idents = c('5'))
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[linked_c5, ],
  n = 13000
)
Vierstra.motifs.linked.5 <- FindMotifs(
  object = combined,
  features = linked_c5,
  background = peaks.matched
)

# Generate rank plot for enriched motifs
df <- data.frame(Vierstra.motifs.linked.5$motif.name, Vierstra.motifs.linked.5$pvalue)
df$Vierstra.motifs.linked.5.pvalue<- -log10(df$Vierstra.motifs.linked.5.pvalue)
df <- df[order(df$Vierstra.motifs.linked.5.pvalue, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggplot(df, aes(rank, Vierstra.motifs.linked.5.pvalue)) + geom_point(size = 1) 

# MOTIF ANALYSIS ON E15 CORTEX SPOTS (RG referes to SOX2+ and N refers to SOX2-)
cortex <- readRDS('/.../e15_ctx_subset.rds')
DefaultAssay(cortex) <- 'peaks'
cortex <- AddMotifs(
  object = cortex,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = Vierstra
)

# Obtain the top 500 peaks with largest fold change differences between SOX2+ and SOX2- spots
RG_fc <- FoldChange(cortex, ident.1 = 'RGCTX', ident.2 = NULL)
N_fc <- FoldChange(cortex, ident.1 = 'NCTX', ident.2 = NULL)

top500.RG.peaks <- RG_fc %>% top_n(500, avg_log2FC)
top500.RG.peaks <- rownames(top500.RG.peaks)
top500.N.peaks <- N_fc %>% top_n(500, avg_log2FC)
top500.N.peaks <- rownames(top500.N.peaks)

# match the overall GC content in the peak set
open.peaks <- AccessiblePeaks(cortex, idents = c("RGCTX", "NCTX"))
meta.feature <- GetAssayData(cortex, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top500.RG.peaks, ],
  n = 13000
)
enriched.motifs.RG.500 <- FindMotifs(
  object = cortex,
  features = top500.RG.peaks,
  background = peaks.matched
)

open.peaks <- AccessiblePeaks(cortex, idents = c("RGCTX", "NCTX"))
meta.feature <- GetAssayData(cortex, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top500.N.peaks, ],
  n = 13000
)
enriched.motifs.N.500 <- FindMotifs(
  object = cortex,
  features = top500.N.peaks,
  background = peaks.matched
)

df <- data.frame(enriched.motifs.RG.500$motif.name, enriched.motifs.RG.500$pvalue)
df$enriched.motifs.RG.500.pvalue<- -log10(df$enriched.motifs.RG.500.pvalue)
df <- df[order(df$enriched.motifs.RG.500.pvalue, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggplot(df, aes(rank, enriched.motifs.RG.500.pvalue)) + geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(90)), ], aes(x = rank, y = enriched.motifs.RG.500.pvalue, label = enriched.motifs.RG.500.motif.name), 
    size = 1.5,
    nudge_x = 2,
    color = "black",
    max.overlaps = 35
  ) 

df <- data.frame(enriched.motifs.N.500$motif.name, enriched.motifs.N.500$pvalue)
df$enriched.motifs.N.500.pvalue<- -log10(df$enriched.motifs.N.500.pvalue)
df <- df[order(df$enriched.motifs.N.500.pvalue, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggplot(df, aes(rank, enriched.motifs.N.500.pvalue)) + geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(90)), ], aes(x = rank, y = enriched.motifs.N.500.pvalue, label = enriched.motifs.N.500.motif.name), 
    size = 1.5,
    nudge_x = 2,
    color = "black",
    max.overlaps = 35
  ) 

### clustering concordance
# compare clustering original vs dca
DefaultAssay(combined) <- "dca"
combined <- RunTFIDF(combined) %>%
  FindTopFeatures(min.cutoff = 'q0') %>%
  RunSVD() %>%
  RunHarmony(group.by.vars = "section", 
             reduction = "lsi", dims.use = 1:7, 
             assay.use = "dca", 
             project.dim = F,
             verbose = T) %>%
  RunUMAP(reduction = "harmony", dims = 1:7, reduction.name = "umap.harmony") %>%
  FindNeighbors(reduction = "harmony", dims = 1:7) %>%
  FindClusters(resolution = 0.7)
combined$dca_snn_res.0.7 <- combined$seurat_clusters

# compare denoised and peaks
for(i in 1:6){
  p1 <- ST.FeaturePlot(combined, "peaks_snn_res.0.7", ncol = 2, split.labels = T, indices = i)
  p2 <- ST.FeaturePlot(combined, "dca_snn_res.0.7", ncol = 2, split.labels = T, indices = i)
  print(p1|p2)
}

Idents(combined) <- "peaks_snn_res.0.7"
levels(combined) <- c("5", "7", "2","9","3", "4", "8","0","1","10","6")
combined$peaks_snn_res.0.7 <- Idents(combined)

Idents(combined) <- "dca_snn_res.0.7"
levels(combined) <- c("5", "9", "2","7","3", "4", "0","1","8","6")
combined$dca_snn_res.0.7 <- Idents(combined)

prop <- prop.table(table(combined$dca_snn_res.0.7, combined$peaks_snn_res.0.7), 1)*100
heatmap.2(prop[order(nrow(prop):1),], 
          Colv = NA, Rowv = NA, scale="none", 
          xlab="peaks cluster", ylab="dca clusters", 
          col = hm_colors, RowSideColors=rev(cols_dca), ColSideColors = cols)

#forebrain
DefaultAssay(combined) <- "peaks"
p1 <- ST.FeaturePlot(combined, "chr3-88205246-88207924", ncol = 2, pt.size = 0.7, cols = magenta_scale, max.cutoff = "q95")
DefaultAssay(combined) <- "dca"
p2 <- ST.FeaturePlot(combined, "chr3-88205246-88207924", ncol = 2, pt.size = 0.7, cols = magenta_scale, max.cutoff = "q95")
p1 | p2

#limb
DefaultAssay(combined) <- "peaks"
p1 <- ST.FeaturePlot(combined, "chr19-5071232-5073427", ncol = 2, pt.size = 0.7, cols = magenta_scale, max.cutoff = "q95")
DefaultAssay(combined) <- "dca"
p2 <- ST.FeaturePlot(combined, "chr19-5071232-5073427", ncol = 2, pt.size = 0.7, cols = magenta_scale, max.cutoff = "q95")
p1 | p2

# liver
DefaultAssay(combined) <- "peaks"
p1 <- ST.FeaturePlot(combined, "chr11-32283339-32285942", ncol = 2, pt.size = 0.7, cols = magenta_scale, max.cutoff = "q95")
DefaultAssay(combined) <- "dca"
p2 <- ST.FeaturePlot(combined, "chr11-32283339-32285942", ncol = 2, pt.size = 0.7, cols = magenta_scale, max.cutoff = "q95")
p1 | p2

