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

#####
# ADD MOTIF ANALYSIS!!
#####


