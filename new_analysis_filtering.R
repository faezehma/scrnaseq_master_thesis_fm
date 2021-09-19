## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
setwd("/Volumes/SD Card/Thesis/")
source("/Volumes/SD Card/Thesis/my_functions.R")

#pbmc_1k_v3 #neuron_1k_v3 #heart_1k_v3 #CISE13 #neuron_5k_v3
## ------------------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(DropletUtils)
library(dplyr)
library(scater)
library(scran)
library(cowplot)
#expression matrix for both filtered data and raw data comes from 10x genomics. 
#fname_filtered <- "/Volumes/SD Card/Thesis/pbmc_10k_v3filtered_feature_bc_matrix/" 
#fname_filtered <- "/Volumes/SD Card/Thesis_2/output_5k_neuron_v3/filtered_feature_bc_matrix/"
fname_filtered <- "/Volumes/SD Card/Thesis_2/CISE13_counts/filtered_feature_bc_matrix/"
#fname_raw <- "/Volumes/SD Card/Thesis_2/output_5k_neuron_v3/raw_feature_bc_matrix/"
#fname_raw <- "/Volumes/SD Card/Thesis/pbmc_10k_v3_raw_feature_bc_matrix/"
fname_raw <- "/Volumes/SD Card/Thesis_2/CISE13_counts/raw_feature_bc_matrix/"
expression_matrix_filtered<-Read10X(data.dir =  fname_filtered)
expression_matrix_raw<-Read10X(data.dir =  fname_raw)

dim(expression_matrix_filtered$`Gene Expression`)
#33538  1222 #31053  1301 #31053  1011 #31053  5810 #31053  6997
dim(expression_matrix_raw$`Gene Expression`)
#33538 6794880 #31053 6794880 #31053 6794880 #31053 737280 #31053 6794880
expression_matrix_filtered<- expression_matrix_filtered$`Gene Expression`
expression_matrix_raw<-expression_matrix_raw$`Gene Expression`
## ------------------------------------------------------------------------
#create singlecellexperiment object:
expression_matrix_filtered<-as(expression_matrix_filtered, "dgCMatrix")
s_filtered_cellranger3<-SingleCellExperiment(list(counts=expression_matrix_filtered))
dim(s_filtered_cellranger3)
#33538  1222 #31053  1301 #31053  1011 #31053  5810 #31053  6997

## ------------------------------------------------------------------------
#create singlecellexperiment object:
expression_matrix_raw<-as(expression_matrix_raw, "dgCMatrix")
s_raw_emptydrops<-SingleCellExperiment(list(counts=expression_matrix_raw))
dim(s_raw_emptydrops)
e.out_s_raw_emptydrops <- emptyDrops(counts(s_raw_emptydrops))
is.cell <- e.out_s_raw_emptydrops$FDR <= 0.01
sum(e.out_s_raw_emptydrops$FDR <= 0.01, na.rm=TRUE)
s_raw_emptydrops <- s_raw_emptydrops[,which(e.out_s_raw_emptydrops$FDR <= 0.01)]
dim(s_raw_emptydrops) # 33538  1206 #31053  1497 #31053   954 #31053  7918 #31053  6935
table(Sig=e.out_s_raw_emptydrops$FDR <= 0.01, Limited=e.out_s_raw_emptydrops$Limited)
## ------------------------------------------------------------------------
#create seurat object for raw
s_raw_200 <- CreateSeuratObject(counts = expression_matrix_raw, min.features = 200)
dim(s_raw_200) #33538  1202 #31053  1479 #941 #31053  7129 #31053  9125
## ------------------------------------------------------------------------

library(VennDiagram)
library(grid)
venn.plot <- venn.diagram(
  x = list(
    "CellRanger 3" = colnames(x = s_filtered_cellranger3),
    "EmptyDrops" = colnames(x = s_raw_emptydrops),
    "Seurat 200" = colnames(x = s_raw_200)
    
  ),
  euler.d = TRUE,
  #filename=NULL,
  filename = "cell_calling_pbmc1k_venn.tiff",
  fil= c("red","blue","green"),
  cex = 2,
  cat.cex = 0.7,
  reverse = TRUE
);

## ------------------------------------------------------------------------
library("UpSetR")
my_list <- list(
  "CellRanger 3" = colnames(x = s_filtered_cellranger3),
  "EmptyDrops" = colnames(x = s_raw_emptydrops),
  "Seurat 200" = colnames(x = s_raw_200)
  )
png("cell_calling_pbmc1k_upset.png", res = 100)
upset(fromList(my_list),sets.bar.color = c("red","blue","green"),
      main.bar.color = "black")
dev.off()
## ------------------------------------------------------------------------
library(R.utils)
xx.1 <- list(s_filtered_cellranger3 = colnames(x = s_filtered_cellranger3),
             s_raw_emptydrops = colnames(x = s_raw_emptydrops),
             s_raw_3_200 = colnames(x = s_raw_3_200))
             
Intersect <- function (x) {  
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function (x) {  
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function (x, y) {
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}

Setdiff.s_filtered_cellranger3 <- Setdiff(xx.1["s_filtered_cellranger3"], xx.1[c("s_raw_3_200","s_raw_emptydrops")])
Is.Setdiff.s_filtered_cellranger3 <- colnames(s_filtered_cellranger3) %in% Setdiff.s_filtered_cellranger3
summary(Is.Setdiff.s_filtered_cellranger3)

Setdiff.s_raw_emptydrops <- Setdiff(xx.1["s_raw_emptydrops"], xx.1[c("s_raw_3_200","s_filtered_cellranger3")])
Is.Setdiff.s_raw_emptydrops <- colnames(s_raw_emptydrops) %in% Setdiff.s_raw_emptydrops
summary(Is.Setdiff.s_raw_emptydrops)

Setdiff.s_raw_3_200_emptydrops <- Setdiff(xx.1[c("s_raw_3_200","s_raw_emptydrops")], xx.1["s_filtered_cellranger3"])
Is.Setdiff.s_raw_3_200_emptydrops<- colnames(s_raw_3_200) %in% Setdiff.s_raw_3_200_emptydrops
summary(Is.Setdiff.s_raw_3_200_emptydrops)

Setdiff.s_raw_3_200_filtered_cellranger3 <- Setdiff(xx.1[c("s_raw_3_200","s_filtered_cellranger3")], xx.1["s_raw_emptydrops"])
Is.Setdiff.s_raw_3_200_filtered_cellranger3<- colnames(s_filtered_cellranger3) %in% Setdiff.s_raw_3_200_filtered_cellranger3
summary(Is.Setdiff.s_raw_3_200_filtered_cellranger3)

Intersect.all <- Intersect(xx.1[c("s_raw_emptydrops","s_filtered_cellranger3","s_raw_3_200")])
Is.Intersect.all<- colnames(s_filtered_cellranger3) %in% Intersect.all
summary(Is.Intersect.all)

Union.all <- Union(xx.1[c("s_raw_emptydrops","s_filtered_cellranger3","s_raw_3_200")])
Is.Union.all<- colnames(expression_matrix_raw) %in% Union.all
summary(Is.Union.all)
## ------------------------------------------------------------------------
cells.all <- expression_matrix_raw[ ,Union.all]
dim(cells.all)
## ------------------------------------------------------------------------
cells.Intersect.all <- cells.all[,Intersect.all]
cells.Setdiff.s_filtered_cellranger3 <- cells.all[,Setdiff.s_filtered_cellranger3]
cells.Setdiff.s_raw_emptydrops <- cells.all[,Setdiff.s_raw_emptydrops]
cells.Setdiff.s_raw_3_200_emptydrops<- cells.all[,Setdiff.s_raw_3_200_emptydrops]
cells.Setdiff.s_raw_3_200_filtered_cellranger3 <- cells.all[,Setdiff.s_raw_3_200_filtered_cellranger3]
## ------------------------------------------------------------------------
region0 = cells.all
region1=cells.Intersect.all 
region2=cells.Setdiff.s_filtered_cellranger3 
region3=cells.Setdiff.s_raw_emptydrops 
region4=cells.Setdiff.s_raw_3_200_emptydrops
region5=cells.Setdiff.s_raw_3_200_filtered_cellranger3 

cells.all <- CreateSeuratObject(cells.all)
region1 <- CreateSeuratObject(region1)
region2 <- CreateSeuratObject(region2)
region3 <- CreateSeuratObject(region3)
region4 <- CreateSeuratObject(region4)
region5 <- CreateSeuratObject(region5)

region1$original <- region1$orig.ident
region1$orig.ident <- "Intersect all initial methods"
region2$original <- region2$orig.ident
region2$orig.ident <- "cellranger"
region3$original <- region3$orig.ident
region3$orig.ident <- "DropletUtils"
region4$original <- region4$orig.ident
region4$orig.ident <- "Seurat3-200 & DropletUtils"
region5$original <- region5$orig.ident
region5$orig.ident <- "Seurat3-200 & cellranger"

combined.object.pbmc <- merge(region1, y =c(region2, region3, region4, region5),add.cell.ids = c("intersect all initial methods","cellranger", "DropletUtils","Seurat3-200 & DropletUtils","Seurat3-200 & cellranger"),project = "initialfiltering")
unique <-unique(sapply(X = strsplit(colnames(combined.object.pbmc), split = "_"), FUN = "[", 1))
table(combined.object.pbmc$orig.ident)
## ------------------------------------------------------------------------
# saveRDS(combined.object.pbmc, file = "new_analysis/combined.object.pbmc.rds")
# combined.object.pbmc$log10GenesPerUMI <- log10(combined.object.pbmc$nFeature_RNA)/log10(combined.object.pbmc$nCount_RNA)
# combined.object.pbmc$mitoRatio <- PercentageFeatureSet(object = combined.object.pbmc, pattern = "^MT-")
# combined.object.pbmc$mitoRatio <- combined.object.pbmc@meta.data$mitoRatio/100  
# 
# metadata <- combined.object.pbmc@meta.data
# # Add cell IDs to metadata
# metadata$cells <- rownames(metadata)
# Rename columns
# metadata <- metadata %>%
#   dplyr::rename(seq_folder = orig.ident,
#                 nUMI = nCount_RNA,
#                 nGene = nFeature_RNA)
# metadata %>% 
#   ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
#   geom_density(alpha = 0.2) + 
#   scale_x_log10() + 
#   theme_classic() +
#   ylab("Cell density")
# 
# metadata %>% 
#   ggplot(aes(color=sample, x=nGene, fill= sample)) + 
#   geom_density(alpha = 0.2) + 
#   theme_classic() +
#   scale_x_log10()

#merge 3 filtering methods in seurat object
intersect_all <- region1
Idents(intersect_all) <- "All three"
intersect_all$orig.ident <- "All three"

cellranger <- as(counts(s_filtered_cellranger3), "dgCMatrix")
cellranger <- CreateSeuratObject(cellranger)
Idents(cellranger) <- "cellranger"
cellranger$orig.ident <- "cellranger"

emptydrops <- as(counts(s_raw_emptydrops), "dgCMatrix")
emptydrops <- CreateSeuratObject(emptydrops)
Idents(emptydrops) <- "emptydrops"
emptydrops$orig.ident <- "emptydrops"

Idents(s_raw_200) <- "seurat200"
s_raw_200$orig.ident <- "seurat200"
###########################################################
methods_withintersect <- merge(cellranger, y = c(emptydrops, s_raw_3_200, intersect_all), add.cell.ids = c("cellranger", "emptydrops", "seurat200", "All three"))
methods<- rep("cellranger",ncol(methods_withintersect))
methods[Idents(methods_withintersect) == "emptydrops"] <- "emptydrops"
methods[Idents(methods_withintersect) == "seurat200"] <- "seurat200"
methods[Idents(methods_withintersect) == "All three"] <- "All three"
merge.methods <- AddMetaData(methods_withintersect, methods, col.name = "methods")
table(Idents(methods_withintersect))
table(methods_withintersect$orig.ident)
saveRDS(methods_withintersect, file = "new_analysis/methods_withintersect.neuron5kv3.rds")
##########################################################
merge.methods <- merge(cellranger, y = c(emptydrops, s_raw_200), add.cell.ids = c("cellranger", "emptydrops", "seurat200"))

methods<- rep("cellranger",ncol(merge.methods))
methods[Idents(merge.methods) == "emptydrops"] <- "emptydrops"
methods[Idents(merge.methods) == "seurat200"] <- "seurat200"


merge.methods <- AddMetaData(merge.methods, methods, col.name = "methods")
table(Idents(merge.methods))
table(merge.methods$orig.ident)
#save merged methods:
saveRDS(merge.methods, file = "/Volumes/SD Card/Thesis/new_analysis/merge.methods.neuron_5k_v3_test_again.rds")
#add more variables to control QC:
# Add number of genes per UMI for each cell to metadata
merge.methods$log10GenesPerUMI <- log10(merge.methods$nFeature_RNA) / log10(merge.methods$nCount_RNA)
merge.methods$UMIperGenes <- merge.methods$nCount_RNA/merge.methods$nFeature_RNA
merge.methods$log10nUMI <- log10(merge.methods$nCount_RNA)
merge.methods$log10nGene <- log10(merge.methods$nFeature_RNA)
merge.methods$sqrtnUMI <- sqrt(merge.methods$nCount_RNA)

merge.methods$mitoRatio <- PercentageFeatureSet(object = merge.methods, pattern = "^MT-")
merge.methods$mitoRatio <- merge.methods@meta.data$mitoRatio / 100

merge.methods$pct_mt <- merge.methods@meta.data$mitoRatio * 100

metadata <- merge.methods@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns

metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

metadata$filtering <- NA
metadata$filtering[which(str_detect(metadata$cells, "^cellranger_"))] <- "cellranger"
metadata$filtering[which(str_detect(metadata$cells, "^emptydrops_"))] <- "emptydrops"
metadata$filtering[which(str_detect(metadata$cells, "^seurat200_"))] <- "seurat200"

# Add metadata back to Seurat object
merge.methods@meta.data <- metadata
head(merge.methods@meta.data)



pdf(file='new_analysis/plots_pbmc_1k_v3.pdf', width = 15, height = 10)
#density plot

metadata %>% 
  ggplot(aes(color=filtering, x=nUMI, fill= filtering)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density")+
  # geom_vline(xintercept = 3647.539, color = "red", linetype = "dashed")+
  # annotate("text", x = 3647.539, y = Inf,
  #          label = "2SD",
  #          colour = "red", angle = 90, hjust = 1, vjust = -0.5)+
  # geom_vline(xintercept = 3047.895, color = "blue", linetype = "dashed")+
  # annotate("text", x = 3047.895, y = Inf,
  #          label = "2SD",
  #          colour = "blue", angle = 90, hjust = 1, vjust = -0.5)+
  # geom_vline(xintercept = 4008.667, color = "#009a00", linetype = "dashed")+
  # annotate("text", x = 4008.667, y = Inf,
  #          label = "2SD",
  #          colour = "#009a00", angle = 90, hjust = 1, vjust = -0.5)



cellranger %>%
  ggplot(aes(x=nUMI), fill = "gray") +
  geom_histogram(alpha = 0.1, bins =30) +
  scale_x_log10() +
  theme_classic() +
  # geom_vline(xintercept = median(metadata$nUMI), colour="red", size=1)+
  # facet_wrap(~filtering) +
  ylab("frequency")+
  geom_vline(xintercept = exp(threshold), color="red")

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=filtering, x=nGene, fill= filtering)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10()+
  ylab("Cell density")+
  # geom_vline(xintercept = 2084.491, color = "red", linetype = "dashed")+
  # annotate("text", x = 2084.491, y = Inf,
  #          label = "2SD",
  #          colour = "red", angle = 90, hjust = 1, vjust = -0.5)+
  # geom_vline(xintercept = 1892.344, color = "blue", linetype = "dashed")+
  # annotate("text", x = 1892.344, y = Inf,
  #          label = "2SD",
  #          colour = "blue", angle = 90, hjust = 1, vjust = -0.5)+
  # geom_vline(xintercept = 2051.162, color = "#009a00", linetype = "dashed")+
  # annotate("text", x = 2051.162, y = Inf,
  #          label = "2SD",
  #          colour = "#009a00", angle = 90, hjust = 1, vjust = -0.5)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=filtering, y=log10(nGene), fill=filtering)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes","pbmc_1k_v3")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  ggtitle("pbmc_1k_v3")+
  facet_wrap(~filtering)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=filtering, x=mitoRatio, fill=filtering)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ggtitle("pbmc_1k_v3") 

metadata %>%
  ggplot(aes(x=UMIperGenes, color = filtering, fill=filtering)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  ggtitle("pbmc_1k_v3")

ggplot(metadata, aes(x = nUMI,y=nGene,colour=mitoRatio)) +
  geom_point(size=0.7) +
  scale_colour_gradientn(colours = c("darkblue","darkgreen","green")) +
  xlab("Count depth (total UMI count)") +
  ylab("Nr of Genes") +
  # geom_vline(xintercept=ul.cd,
  #            color = "red", size=1) +
  # geom_vline(xintercept=ll.cd,
  #            color = "red", size=1) +
  # geom_hline(yintercept=ul.ng,
  #            color = "red", size=1) +
  # geom_hline(yintercept=ll.ng,
  #            color = "red", size=1) +
  geom_rug(col=rgb(0,0,0.5,alpha=.1)) +
  facet_wrap(~filtering) +
  ggtitle("pbmc_1k_v3")+
  theme_bw()
dev.off()
###########
cellranger <- subset(metadata, methods =="cellranger")
ggplot(ranger, aes(x = nUMI,y=nGene,colour=mitoRatio)) +
  geom_point(size=0.5) +
  scale_color_gradientn(colours = c("darkblue","darkgreen","green","yellow")) +
  xlab("Count depth (total UMI count)") +
  ylab("Nr of Genes") +
  geom_abline(intercept = 0, slope =1, color = "green", size=0.5) + 
  geom_abline(intercept = 0, slope =1/2, color = "cyan", size=0.5) +

  geom_rug(col=rgb(0,0,0.5,alpha=.1)) +
  facet_wrap(~filtering) +
  ggtitle("CISE13")+
  theme_bw()

#calculation of MADS:

med    <- median(cellranger$log10nUMI)
MAD    <- mad(cellranger$log10nUMI, center = med, na.rm = TRUE)
n <- 1
lower  <- med - n * MAD
lower
higher <- med + n * MAD
list(n = n, lower = lower, higher = higher,
      n_low = n_low, n_high = n_high)

#Dimensionality reduction for background method
sce <- as.SingleCellExperiment(methods_withintersect)
dimred_factors <- c(
  #"count depth" = "nUMI",
  #"Total features" = "nGene",
  #"Mitochondrial genes" = "mitoRatio",
  "Selection method" = "ident")

### PCA
# plot_list <- lapply(names(dimred_factors), function(fct_name) {
#    plotPCA(sce, colour_by = dimred_factors[fct_name], shape_by = "ident"
#           )+
#     ggtitle(fct_name)+
#     geom_point(alpha = 0.0000001)
#  })
# plot_grid(plotlist = plot_list, ncol = 1)


DimPlot(methods_withintersect,
        reduction = "umap",
        split.by = "ident",
        label = TRUE,
        label.size = 6,
        plot.title = "UMAP")
  
##############################To show data in PCA plot:
sce <- as.SingleCellExperiment(methods_withintersect)

colnames(colData(sce))
#ranger <- subset(colData(sce), ident == "cellranger")
patterns <- c("^hg19_MT-", "^mm10_mt-","^MT-", "^mt-")
is.mito<- grepl(paste(patterns, collapse="|"), rownames(sce))
sce<- calculateQCMetrics(sce,feature_controls=list(Mt=is.mito)) 


sce <- runPCA(sce, use_coldata = TRUE, exprs_values = "log10_total_counts")
plotPCA(
  sce,
  colour_by = "ident"
)+
  geom_point(alpha = 0.00000000000000000000001)+
  scale_fill_manual(values = c("grey44","Red", "Green", "Blue"))


p<- plotReducedDim(sce, use_dimred="PCA_coldata", colour_by = "ident")+
  geom_rug(col=rgb(0,0,0.5,alpha=.1))+
  ggtitle("CISE13")+ 
  geom_point(alpha = 0.00000000000000000000001)+
  scale_fill_manual(values = c("grey44","Red", "Green", "Blue"))
print(p)

colnames(colData(sce))[colnames(colData(sce))=="log10_total_counts"] <- "logcounts"


set.seed(123456)
sce <- runUMAP(sce)

sce <- runTSNE(sce, perplexity=50, 
                       dimred="PCA", n_dimred=10)
reducedDim(sce, "UMAP")
plotReducedDim(sce, use_dimred = "UMAP", colour_by = "ident")

plotTSNE(sce, colour_by = "ident")
# plotUMAP(
#   sce,
#   colour_by = "ident"
# )

# # run PCA with 1000 top variable genes
# sce <- runPCA(sce, ntop = 1000, exprs_values = "logcounts", ncomponents = 20)
# 
# # PCA - with different coloring, first 4 components
# # first by sample
# plotPCA(sce,ncomponents=4,colour_by="ident")+ggtitle("CISE13")+ 
#   geom_point(alpha = 0.0000001)
#   
# plotPCA(sce,ncomponents=4,colour_by= "mitoRatio")



####apply MADs or manual cells filtering on each dataset: at  first apply ablines for QC plots:
##count depth vs number of genes plot:

#MADs threshod
n <- 1
med    <- median(log(metadata$nUMI))
MAD    <- mad(log(metadata$nUMI), center = med)
Q <- Qn(log(metadata$nUMI))
s <- scaleTau2(log(metadata$nUMI))
lower  <- med - n * MAD
lower <- med - n *Q
lower <- med - n *s
lower

med    <- median(log(metadata$nGene))
MAD    <- mad(log(metadata$nGene), center = med)
lowerg  <- med - n * MAD
lowerg

ll.cd = exp(lower)
ll.ng = exp(lowerg)
#doublet threshold should be apply after normalization and dimensionality reduction
# set.seed(1)
# doublet_scores <- doubletCells(sce)
# colData(sce)$DoubletScore <- doublet_scores
# doublet_out <- isOutlier(doublet_scores, nmads = 15, type = "higher")
# doublet_thresh <- attr(doublet_out, "thresholds")["higher"]

########################the other datasets are in visualization.R file:
#CISE dataset:
ll.cd = 10^2.763
ll.ng = 10^2.422
ul.cd = 10^3.729
ul.ng = 10^3.258

ll.cd2 = 10^2.597
ll.ng2 = 10^2.497
ul.cd2 = 10^3.756
ul.ng2 = 10^3.247

ll.cd3 = 10^2.563
ll.ng3 = 10^2.374
ul.cd3 = 10^3.763
ul.ng3 = 10^3.275

cellranger <- subset(metadata, methods =="cellranger")
emptydrops <- subset(metadata, methods =="emptydrops")
seurat200 <- subset(metadata, methods =="seurat200")
pdf(file='new_analysis/new_cutoff/SD_cutoff/all_joint_cutoff_CISE.pdf', width = 15, height = 10)
ggplot(cellranger, aes(x = nUMI,y=nGene,colour=mitoRatio)) +
  geom_point(size=1.5, alpha = 0.8) +
  scale_colour_gradientn(colours = c("darkblue","darkgreen","green", "darkorange2")) +
  xlab("Count depth (total UMI count)") +
  ylab("Nr of Genes") +
  geom_vline(xintercept=ul.cd,
             color = "red", size=0.6, linetype="dashed") +
  geom_vline(xintercept=ll.cd,
              color = "red", size=0.6, linetype="dashed") +
  geom_hline(yintercept=ul.ng,
             color = "red", size=0.6, linetype="dashed") +
  geom_hline(yintercept=ll.ng,
              color = "red", size=0.6,linetype="dashed") +
  
  geom_rug(col=rgb(0,0,0.5,alpha=.1)) +
  facet_wrap(~filtering) +
  ggtitle("CISE13")+
  theme_bw()
ggplot(emptydrops, aes(x = nUMI,y=nGene,colour=mitoRatio)) +
  geom_point(size=1.5, alpha = 0.8) +
  scale_colour_gradientn(colours = c("darkblue","darkgreen","green","darkorange2")) +
  xlab("Count depth (total UMI count)") +
  ylab("Nr of Genes") +
  geom_vline(xintercept=ul.cd2,
             color = "red", size=0.6, linetype="dashed") +
  geom_vline(xintercept=ll.cd2,
             color = "red", size=0.6, linetype="dashed") +
  geom_hline(yintercept=ul.ng2,
             color = "red", size=0.6, linetype="dashed") +
  geom_hline(yintercept=ll.ng2,
             color = "red", size=0.6,linetype="dashed") +
  
  geom_rug(col=rgb(0,0,0.5,alpha=.1)) +
  facet_wrap(~filtering) +
  ggtitle("CISE13")+
  theme_bw()

ggplot(seurat200, aes(x = nUMI,y=nGene,colour=mitoRatio)) +
  geom_point(size=1.5, alpha = 0.8) +
  scale_colour_gradientn(colours = c("darkblue","darkgreen","green","darkorange2")) +
  xlab("Count depth (total UMI count)") +
  ylab("Nr of Genes") +
  geom_vline(xintercept=ul.cd3,
             color = "red", size=0.6, linetype="dashed") +
  geom_vline(xintercept=ll.cd3,
             color = "red", size=0.6, linetype="dashed") +
  geom_hline(yintercept=ul.ng3,
             color = "red", size=0.6, linetype="dashed") +
  geom_hline(yintercept=ll.ng3,
             color = "red", size=0.6,linetype="dashed") +
  
  geom_rug(col=rgb(0,0,0.5,alpha=.1)) +
  facet_wrap(~filtering) +
  ggtitle("CISE13")+
  theme_bw()
dev.off()
#PBMC:
##for PBMC I prefer to take 3MADs for anual threshold based on this plot to check the high mito ratio in lower part of the plot
#4MADs is still some points with high mito ratio remained

####barcode rank plot: for each background methods output
#sce of raw count data:
cellranger <- subset(sce, ,ident=="cellranger")
emptydrops <- subset(sce, ,ident=="emptydrops")
seurat200 <- subset(sce, ,ident=="seurat200")
final = counts(emptydrops)

###############or on raw data find knee: Not necessary now
fname_raw <- "/Volumes/SD Card/Thesis_2/output_5k_neuron_v3/raw_feature_bc_matrix/"
expression_matrix_raw<-Read10X(data.dir =  fname_raw)
expression_matrix_raw<-as(expression_matrix_raw, "dgCMatrix")
sce_raw<-SingleCellExperiment(list(counts=expression_matrix_raw))
final = counts(sce_raw)
###############################################################################
#final <- counts(sce)
br.out <- barcodeRanks(final)
write.csv(br.out,'new_analysis/bar.out.cellranger.n5k.csv')
## Create the knee plot
plot(log10(br.out$rank), log10(br.out$total))
abline(h = log10(metadata(br.out)$knee))
#define knee for applying abline on hist plot
knee <- log10(metadata(br.out)$knee)
knee <- metadata(br.out)$knee
###
# # Making a plot.
# plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total count")
# o <- order(br.out$rank)
# lines(br.out$rank[o], br.out$fitted[o], col="red")
# 
# abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
# abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
# legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
#        legend=c("knee", "inflection"))

#histo plot for knee
#histogram
emptydrops <- subset(metadata, methods =="emptydrops")

hist(x = cellranger$log10nUMI,
     breaks = 100)
abline(v = med , col = "red")

abline(v = knee,
       col = "red", size=0.7)
abline(v = log10(f_low),
       col = "blue", size = 0.7)


####check different univariate outliers detection methods
library(sn)
library(univOutl)
mc<- mc(emptydrops$nUMI)
Q0 <- IQR(emptydrops$nUMI)
Q3 <- quantile(emptydrops$nUMI, 0.75) 
Q1 <- quantile(emptydrops$nUMI, 0.25)  
f_low = Q1-(1.5*exp(-4*mc)*Q0)

mad(cellranger$nUMI)
Sn(cellranger$nUMI)
Qn(cellranger$nUMI)
scaleTau2(cellranger$nUMI)
a1 <- LocScaleB(x=cellranger$nUMI,
                method = "IQR")
a2 <- LocScaleB(x=cellranger$nUMI,
                method = "dq")
a2 <- boxB(x= cellranger$nUMI, k=1.5, method='asymmetric')
a3 <- boxB(x=cellranger$nUMI, k=1.5, method='adjbox')

####tsne plot for 3 background methods:
library(scran)
sce <- as.SingleCellExperiment(methods_withintersect)
set.seed(1)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters=clusters)
sce <- normalize(sce)

library(scater)
sizeFactors(sce) <- librarySizeFactors(sce)
sce <- normalize(sce)

sce <- runPCA(sce)
sce <- runTSNE(sce, dimred="PCA")




#subset1 <- subset(methods_withintersect, idents = c("emptydrops", "All three"))
subset1 <- subset(sce, ,ident==c("emptydrops", "All three"))
subset2 <- subset(sce, ,ident==c("seurat200", "All three"))
plotTSNE(subset2, colour_by="ident")+ 
  scale_fill_manual(values = c("grey", "Green"))+
  xlab("t-SNE1")+ ylab("t-SNE2")+
  ggtitle("")


# Defining arrow coordinates.
coords <- reducedDim(sce, "TSNE")
plate.clust <- metadata(sce)$Platelets
platelets <- colMeans(coords[sce$Cluster==plate.clust,]) 
FUN <- function(coloration, ...) {
  plot(coords[,1], coords[,2], col=coloration, pch=16, xlab="t-SNE1", ylab="t-SNE2", cex.axis=1.2, cex.lab=1.4, ...)
  SHIFT <- c(0, -7)
  WIDTH <- c(0, 4)
  arrows(platelets[1] - SHIFT[1] - WIDTH[1], platelets[2] - SHIFT[2],
         platelets[1] - SHIFT[1], platelets[2] - SHIFT[2] - WIDTH[2],
         angle=20, length=0.1, lwd=2, col="black")
}

# Making a plot of detection status.
coloration <- rep("grey", nrow(coords))
coloration[sce$ident=="emptydrops"] <- "green"
coloration[sce$ident=="cellranger"] <- "red"
coloration[sce$ident=="seurat200 "] <- "blue"

#pdf("pics/by_detection.pdf")
FUN(coloration)+
  geom_point(alpha = 0.00000000000000000000001)
legend("bottomright", legend=c("Three", "CellRanger", "EmptyDrops", "Seurat200"), col=c("grey", "red", "green", "blue"), pch=16)
#dev.off()

###fit dataset to find outliers and based on that find optimal parameter to determine cutoff

cellranger <- subset(metadata, methods =="cellranger")
emptydrops <- subset(metadata, methods =="emptydrops")
seurat200 <- subset(metadata, methods =="seurat200")
my_vector <- cellranger$log10nUMI
my_vector2 <- cellranger$log10nGene
my_vector <- cellranger$nUMI
my_vector <- sqrt(my_vector)
my_vector2 <- cellranger$nUMI
my_vector2 <- sign(my_vector) * abs(my_vector)^(1/3)

library("fitdistrplus")
fitd <- fitdist(my_vector, "lnorm", method = "mle")
p <- plot(fitd)

#mixture gaussian distribution
library(mixtools)
mixmdl = normalmixEM(my_vector, k=2,maxit=10000,epsilon=1e-08, maxrestarts = 1000)
summary(mixmdl)
mixmdl$mu[2]

mixmdl2 = normalmixEM(my_vector2, k=2,maxit=10000,epsilon=1e-08, maxrestarts = 1000)
summary(mixmdl2)
mixmdl2$mu[2] 

# hist(my_vector, freq = FALSE, breaks = 100, main = "cellranger PBMC")
# #show the respective curves
# lines(my_vector,mixmdl$lambda[1]*dnorm(my_vector,mixmdl$mu[1],mixmdl$sigma[1]), col = "green")
# lines(my_vector,mixmdl$lambda[2]*dnorm(my_vector,mixmdl$mu[2],mixmdl$sigma[2]), col = "red")

#identify cut-off:
library("devtools")
devtools::install_github("choisy/cutoff")
library(cutoff)
my_vector_out <- em(my_vector, "normal","normal")
confint(my_vector_out, nb=100, level=.95)
hist(my_vector,100,F,xlab="log10nUMI",ylab="density",ylim=c(0,3),xlim = c(2.5,4.5) ,
          main=NULL,col="grey")
lines(my_vector_out,lwd=1.5,col="red")
(cut_off <- cutoff(my_vector_out))

polygon(c(cut_off[-1],rev(cut_off[-1])),c(0,0,.55,.55),
                col=rgb(0,0,1,.2),border=NA)
abline(v=cut_off[-1],lty=2,col="blue")
abline(v=cut_off[1],col="blue")

plot(hist(my_vector,breaks=101),col="grey",border="grey",freq=FALSE,
     xlab="log10nUMI",main="cellranger PBMC")
lines(density(my_vector),lty=2)

#####correct one for mixture guassian distribution:
plot.normal.components <- function(mixture,component.number,...) {
  curve(mixture$lambda[component.number] *
          dnorm(x,mean=mixture$mu[component.number],
                sd=mixture$sigma[component.number]), add=TRUE, ...)
}

plot(hist(my_vector,breaks=101),freq=FALSE,
     xlab="sqrt(nUMI)",main=" ")
lines(density(my_vector),lty=2, col ="red")

# abline(v= mixmdl$mu[2], col = "blue", lty = "dashed")
# text(x = mixmdl$mu[2],y = Inf, label= "mu", srt = 90)
sapply(1:2,plot.normal.components,mixture=mixmdl)
abline(v = mixmdl$mu[2], col = "blue")
abline(v = mixmdl$mu[2]-2.5*mixmdl$sigma[2], col = "red")
####
#making gaussian distribution in ggplot
sdnorm =
  function(x, mean=0, sd=1, lambda=1){lambda*dnorm(x, mean=mean, sd=sd)}

gg <- ggplot(cellranger, aes(x=sqrt(nUMI))) +
  geom_histogram(bins = 100, color = "black", fill = "gray",  
                 aes(y=..density.., fill=..count..)) +
  geom_vline(xintercept = mixmdl$mu[2], colour = "blue", linetype = "dashed") +
  annotate("text", x = mixmdl$mu[2], y = 0, label = "Mean", colour = "blue",
           angle = 90, hjust = -0.1, vjust = -0.5) +
  annotate("text", x = mixmdl$mu[2], y = Inf, label = round(mixmdl$mu[2], 3), colour = "blue",
           angle = 90, hjust = 1, vjust = -0.5)+
  stat_function(fun=sdnorm,args=list(mean=mixmdl$mu[1],sd=mixmdl$sigma[1], lambda=mixmdl$lambda[1]), color ="black")+
  stat_function(fun=sdnorm,args=list(mean=mixmdl$mu[2],sd=mixmdl$sigma[2], lambda=mixmdl$lambda[2]),color =  "black")


(lower_cutoff = mixmdl$mu[2]-3*mixmdl$sigma[2])
######make marginal plot -------------------------------------------------------------------------------------------------- 
ll.cd = 2.763
ul.cd = 3.729
ll.ng = 2.422
ul.ng = 3.258
g1 <- ggplot(cellranger, aes(x = log10nUMI,y=log10nGene,colour=mitoRatio)) +
  geom_point(size=0.7) +
  scale_colour_gradientn(colours = c("darkblue","darkgreen","green","yellow")) +
  xlab("log10nUMI") +
  ylab("log10nGene") +
  geom_vline(xintercept=ul.cd,
              color = "red", size= 0.6, linetype="dashed") +
  geom_vline(xintercept=ll.cd,
             color = "red", size=0.6, linetype="dashed") +
  
  geom_hline(yintercept=ul.ng,
             color = "red", size=0.6, linetype="dashed") +
  geom_hline(yintercept=ll.ng,
             color = "red", size=0.6,linetype="dashed") +
  geom_rug(col=rgb(0,0,0.5,alpha=.1)) +
  theme_bw()+
  ggtitle("cellranger CISE13")

# c1 <- sdnorm(cellranger$log10nUMI, mean=mixmdl$mu[1], sd=mixmdl$sigma[1], lambda=mixmdl$lambda[1])
# c2 <- sdnorm(cellranger$log10nUMI, mean=mixmdl$mu[2], sd=mixmdl$sigma[2], lambda=mixmdl$lambda[2])
# d <- c(c1,c2)
# df <- data.frame(d, id = as.factor(rep(c(1,2), each = 5810)))
# g2 <- ggplot(cellranger, aes(x=log10nUMI)) +
  # geom_histogram(bins = 100, color = "black", fill = "gray",  
  #                aes(y=..density.., fill=..count..)) +
  # geom_vline(xintercept = mixmdl$mu[2], colour = "blue", linetype = "dashed") +
  # annotate("text", x = mixmdl$mu[2], y = 0, label = "Mean", colour = "blue",
  #          angle = 90, hjust = -0.1, vjust = -0.5) +
  # annotate("text", x = mixmdl$mu[2], y = Inf, label = round(mixmdl$mu[2], 3), colour = "blue",
  #          angle = 90, hjust = 1, vjust = -0.5)+
xdens <- cowplot::axis_canvas(g1, axis = "x")+
  geom_density(data = cellranger, aes(x = log10nUMI))

  # stat_function(fun=sdnorm,args=list(mean=mixmdl$mu[1],sd=mixmdl$sigma[1], lambda=mixmdl$lambda[1]), color ="black")+
  # stat_function(fun=sdnorm,args=list(mean=mixmdl$mu[2],sd=mixmdl$sigma[2], lambda=mixmdl$lambda[2]),color =  "black")

# g3 <- ggplot(cellranger, aes(x=log10nGene)) +
  # geom_histogram(bins = 100, color = "black", fill = "gray",  
  #                aes(y=..density.., fill=..count..)) +
  # geom_vline(xintercept = mixmdl$mu[2], colour = "blue", linetype = "dashed") +
  # annotate("text", x = mixmdl$mu[2], y = 0, label = "Mean", colour = "blue",
  #          angle = 90, hjust = -0.1, vjust = -0.5) +
  # annotate("text", x = mixmdl$mu[2], y = Inf, label = round(mixmdl$mu[2], 3), colour = "blue",
  #          angle = 90, hjust = 1, vjust = -0.5)+
ydens <- cowplot::axis_canvas(g1, axis = "y", coord_flip = TRUE)+
  geom_density(data = cellranger, aes(x=log10nGene))+
  stat_function(fun=sdnorm,args=list(mean=mixmdl2$mu[1],sd=mixmdl2$sigma[1], lambda=mixmdl2$lambda[1]), color ="black")+
  stat_function(fun=sdnorm,args=list(mean=mixmdl2$mu[2],sd=mixmdl2$sigma[2], lambda=mixmdl2$lambda[2]),color =  "black")+
  coord_flip()

p1 <- cowplot::insert_xaxis_grob(g1, xdens, grid::unit(.2, "null"), position = "top")

p2 <- cowplot::insert_yaxis_grob(p1, ydens,
                                 grid::unit(.2, "null"), position = "right")
# plot_grid(
#   plot_grid(plotlist=list(g2, ggplot(), g1, g3), ncol=2, 
#             rel_widths=c(4,1), rel_heights=c(1,4), scale=1.1), rel_widths=c(5,1))
ggdraw(p2)

########################plotpca, plottsne and automated outlier detection ***before filtering***:
library(scran)
methods_withintersect <- read_rds("new_analysis/methods_withintersect.neuron5kv3.rds")

sce <- as.SingleCellExperiment(methods_withintersect)
set.seed(1)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters=clusters)
sce <- normalize(sce)

sce <- runPCA(sce)
sce <- runTSNE(sce, dimred="PCA")
sce <- runUMAP(sce)


subset0 <- subset(sce, ,ident==c("cellranger", "All three"))
subset1 <- subset(sce, ,ident==c("emptydrops", "All three"))
subset2 <- subset(sce, ,ident==c("seurat200", "All three"))

p000 <- plotUMAP(subset0, colour_by="ident")+ 
  scale_fill_manual(values = c("grey", "red"))+
  xlab("UMAP1")+ ylab("UMAP2")+
  ggtitle("")

p111 <- plotUMAP(subset1, colour_by="ident")+ 
  scale_fill_manual(values = c("grey", "green"))+
  xlab("UMAP1")+ ylab("UMAP2")+
  ggtitle("")

p222 <- plotUMAP(subset2, colour_by="ident")+ 
  scale_fill_manual(values = c("grey", "blue"))+
  xlab("UMAP1")+ ylab("UMAP2")+
  ggtitle("")

fig <- plot_grid(p000, p111, p222,
                 ncol = 2, labels = "AUTO")
ggsave(("new_analysis/plotUMAP_3backgroundmethods_CISE.pdf"), fig, width = 7,height=6, scale = 1.5)

plotPCA(sce, colour_by = "ident")+ 
  scale_fill_manual(values = c("grey", "red","Green", "blue"))

p00 <- plotPCA(subset0,colour_by = "ident")+ 
  scale_fill_manual(values = c("grey", "red"))

p11 <- plotPCA(subset1,colour_by = "ident")+ 
  scale_fill_manual(values = c("grey", "green"))

p22 <- plotPCA(subset2,colour_by = "ident")+ 
  scale_fill_manual(values = c("grey", "blue"))

fig <- plot_grid(p00, p11, p22,
          ncol = 2, labels = "AUTO")


ggsave(("new_analysis/plotPCA_3backgroundmethods_CISE.pdf"), fig, width = 7,height=6, scale = 1.5)


p0 <- plotTSNE(subset0, colour_by="ident")+ 
  scale_fill_manual(values = c("grey", "red"))+
  xlab("t-SNE1")+ ylab("t-SNE2")+
  ggtitle("")

p1 <- plotTSNE(subset1, colour_by="ident")+ 
  scale_fill_manual(values = c("grey", "green"))+
  xlab("t-SNE1")+ ylab("t-SNE2")+
  ggtitle("")

p2 <- plotTSNE(subset2, colour_by="ident")+ 
  scale_fill_manual(values = c("grey", "blue"))+
  xlab("t-SNE1")+ ylab("t-SNE2")+
  ggtitle("")

fig <- plot_grid(p0, p1, p2,
          ncol = 2, labels = "AUTO")

ggsave(("new_analysis/plottsne_3backgroundmethods_CISE.pdf"), fig, width = 7,height=6, scale = 1.5)

#colnames(colData(sce))
sce0 <- subset(sce, ,ident=="cellranger")
#ranger <- subset(colData(sce), ident == "cellranger")
patterns <- c("^hg19_MT-", "^mm10_mt-","^MT-", "^mt-")
is.mito<- grepl(paste(patterns, collapse="|"), rownames(sce))
sce<- calculateQCMetrics(sce,feature_controls=list(Mt=is.mito)) 
#for neuron 5k I have to use selected_variables and remove one of the defualt variables from coldata: total_features_by_counts_feature_control
sce <- runPCA(sce, use_coldata = TRUE, detect_outliers = TRUE, 
              selected_variables = c("pct_counts_in_top_100_features", "total_features_by_counts", "pct_counts_feature_control",
                                     "log10_total_counts_endogenous", "log10_total_counts_feature_control"))
p_out <- plotReducedDim(sce, use_dimred = "PCA_coldata", colour_by = "outlier")+
  geom_rug(col=rgb(0,0,0.5,alpha=.1))

p<- plotReducedDim(sce, use_dimred="PCA_coldata", colour_by = "ident")+ 
  scale_fill_manual(values = c("grey44","Red", "Green", "Blue"))+
  geom_rug(col=rgb(0,0,0.5,alpha=.1))+ 
  geom_point(alpha = 0.00000000000000000000001)
fig <- plot_grid(p_out, p,
          ncol = 2, labels = "AUTO")
ggsave(("new_analysis/plotPCAoutlier_3backgroundmethods_neuron5k.pdf"), fig, width = 9,height=6, scale = 1.5)

#################testing another transformation on Heart dataset:


cairo_pdf("new_analysis/transformation_heartdataset.sqrt.pdf")
plot(hist(my_vector,breaks=101),freq=FALSE,
     xlab="sqrt(nUMI)",main= "cellranger")
lines(density(my_vector),lty=2, col ="red")
dev.off()
cairo_pdf("new_analysis/transformation_heartdataset.sqrt1.pdf")
fitd <- fitdist(my_vector, "norm", method = "mle")
plot(fitd)
dev.off()
cairo_pdf("new_analysis/transformation_heartdataset.sqrt2.pdf")
mixmdl = normalmixEM(my_vector, k=2,maxit=10000,epsilon=1e-08, maxrestarts = 1000)
plot.normal.components <- function(mixture,component.number,...) {
  curve(mixture$lambda[component.number] *
          dnorm(x,mean=mixture$mu[component.number],
                sd=mixture$sigma[component.number]), add=TRUE, ...)
}

plot(hist(my_vector,breaks=101),freq=FALSE,
     xlab="sqrt(nUMI)",main=" ")
lines(density(my_vector),lty=2, col ="red")
sapply(1:2,plot.normal.components,mixture=mixmdl)
dev.off()
cairo_pdf("new_analysis/transformation_heartdataset.sqrt3.pdf")
plot(hist(my_vector2,breaks=101),freq=FALSE,
     xlab="cube root(nUMI)",main="cellranger")
lines(density(my_vector2),lty=2, col ="red")
dev.off()
cairo_pdf("new_analysis/transformation_heartdataset.sqrt4.pdf")
fitd <- fitdist(my_vector2, "norm", method = "mle")
plot(fitd)
dev.off()
cairo_pdf("new_analysis/transformation_heartdataset.sqrt5.pdf")
mixmdl = normalmixEM(my_vector2, k=2,maxit=10000,epsilon=1e-08, maxrestarts = 1000)
plot.normal.components <- function(mixture,component.number,...) {
  curve(mixture$lambda[component.number] *
          dnorm(x,mean=mixture$mu[component.number],
                sd=mixture$sigma[component.number]), add=TRUE, ...)
}

plot(hist(my_vector2,breaks=101),freq=FALSE,
     xlab="cube root(nUMI)",main="")+
lines(density(my_vector2),lty=2, col ="red")
sapply(1:2,plot.normal.components,mixture=mixmdl)
dev.off()