### load packages

# omics
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)
library(Signac)
library(Seurat)
library(DropletUtils)
library(harmony)
library(gprofiler2)
library(loomR)

#spatial
library(STutility)
library(parallel)
library(ica)
library(spdep)
library(jpeg)

#plots
library(ggplot2)
library(patchwork)
library(scales)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(gplots)

#other
library(dplyr)
library(purrr)
library(Matrix)

###color palettes
magenta_scale <- c("gray93", "#FFF7F3","#FDE0DD","#FCC5C0","#FA9FB5","#F768A1","#DD3497","#AE017E")
cyan_scale <- c("gray93","#F7FCF0","#E0F3DB","#CCEBC5","#7BCCC4","#4EB3D3","#2B8CBE","#0868AC","#084081","#084081","#084081")
fun_color_range <- colorRampPalette(c("gray93", "#AE017E"))  # Create color generating function
hm_colors <- fun_color_range(100)                         

#embryo
bright_visium <- c("#4477AA",  "#9ECAE1","#669C39", "#58BEC8", "#ccddaa", "#CCBB44", "#2F9558", "#B73D77","#E05B77", "#B06992", "#A8A8A8", "#E08762")
visium_dca <- c("#4477AA","#58BEC8","#669C39","#9ECAE1","#ccddaa","#CCBB44","#2F9558","#cc6677","#E05B77","#B73D77","#B06992","#A8A8A8","#E08762")
color_sections <- c("#E69F00", "#D55E00", "#aaaa00", "#009E73","#56B4E9", "#0072B2")
cols <-  c("#4477AA","#58AAD2","#58BEC8", "#2F9558","#669C39","#CCBB44", "#A8A8A8","#E05B77","#B73D77", "#B06992","#E08762") #embryo clusters res 0.7
cols_dca <-  c("#4477AA","#58AAD2","#58BEC8", "#2F9558","#669C39","#CCBB44", "#E05B77","#B73D77", "#B06992","#E08762")
cols_nmf <-c("#4477AA", "#66CCEE", "#228833", "#76A13B", "#CCBB44", "#BBBBBB","#EE6677", "#CB4C77", "#B27699","#AA3377", "#E69F00", "#D55E00")

#brca

#functions
source("as_matrix.R")
source("RunNMF_fix.R")


