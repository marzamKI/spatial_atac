#QC metrics 

#ArchR
# Load fragments file from Cell Ranger
fragments_C1 <- '/.../fragments.tsv.gz'

# Obtain tissue or region of interest barcodes from Loupe Browser
aligned_bc_C1 <- read.csv2('/.../C1_barcodes.csv', 
                           header = FALSE )$V1

fragments_C2 <- '/.../fragments.tsv.gz'
aligned_bc_C2 <- read.csv2('/.../C2_barcodes.csv', 
                           header = FALSE )$V1
fragments_D2 <- '/.../fragments.tsv.gz'
aligned_bc_D2 <- read.csv2('/.../D2_barcodes.csv', 
                           header = FALSE )$V1

# Create arrow files per section
D2T <- createArrowFiles(
  inputFiles = fragments_D2,
  sampleNames = 'D2T',
  validBarcodes = aligned_bc_D2,
  geneAnnotation = getGeneAnnotation(),
  genomeAnnotation = getGenomeAnnotation(),
  minTSS = 1,
  minFrags = 100,
  maxFrags = 1e+05,
  
)

# Create ArchR project 

C1T <- ArchRProject(
  ArrowFiles = C1T, 
  outputDirectory = "/.../ArchR",
  copyArrows = FALSE 
)

C2T <- ArchRProject(
  ArrowFiles = C2T, 
  outputDirectory = "/.../ArchR",
  copyArrows = FALSE 
)

D2T <- ArchRProject(
  ArrowFiles = D2T, 
  outputDirectory = "/.../ArchR",
  copyArrows = FALSE 
)

plotTSSEnrichment(ArchRProj = C1T) + ylim(0,8)
plotTSSEnrichment(ArchRProj = C2T) + ylim(0,8)
plotTSSEnrichment(ArchRProj = D2T) + ylim(0,8)
plotFragmentSizes(C1T) + ylim(0,1)
plotFragmentSizes(C2T) + ylim(0,1)
plotFragmentSizes(D2T) + ylim(0,1)


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
counts <- Read10X_h5(filename = "/.../atac_v1_E18_brain_flash_5k_filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "/.../atac_v1_E18_brain_flash_5k_singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = '/.../atac_v1_E18_brain_flash_5k_fragments.tsv.gz',
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





