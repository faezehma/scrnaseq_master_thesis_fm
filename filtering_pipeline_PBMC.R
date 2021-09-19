## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
setwd("/Volumes/SD Card/Thesis")
source("/Volumes/SD Card/Thesis/my_functions.R")

#pbmc_1k_v3
## ------------------------------------------------------------------------
library(Seurat)
library(SingleCellExperiment)
library(DropletUtils)
#expression matrix for both filtered data and raw data comes from 10x genomics. 
fname_filtered <- "/Volumes/SD Card/Thesis/pbmc_1k_v3_filtered_feature_bc_matrix/filtered_feature_bc_matrix" #(cellrangerv3)
fname_raw <- "/Volumes/SD Card/Thesis/pbmc_1k_v3_raw_feature_bc_matrix/raw_feature_bc_matrix"
expression_matrix_filtered<-Read10X(data.dir =  fname_filtered)
expression_matrix_raw<-Read10X(data.dir =  fname_raw)

dim(expression_matrix_filtered)
#33538  1222
dim(expression_matrix_raw)
#33538 6794880

## ------------------------------------------------------------------------
#create singlecellexperiment object:
expression_matrix_filtered<-as(expression_matrix_filtered, "dgCMatrix")
s_filtered_cellranger3<-SingleCellExperiment(list(counts=expression_matrix_filtered))
dim(s_filtered_cellranger3)
#33538  1222

## ------------------------------------------------------------------------
#create singlecellexperiment object:
expression_matrix_raw<-as(expression_matrix_raw, "dgCMatrix")
s_raw_emptydrops<-SingleCellExperiment(list(counts=expression_matrix_raw))
dim(s_raw_emptydrops)
#33538 6794880
e.out_s_raw_emptydrops <- emptyDrops(counts(s_raw_emptydrops))
is.cell <- e.out_s_raw_emptydrops$FDR <= 0.01
sum(e.out_s_raw_emptydrops$FDR <= 0.01, na.rm=TRUE)
s_raw_emptydrops <- s_raw_emptydrops[,which(e.out_s_raw_emptydrops$FDR <= 0.01)]
dim(s_raw_emptydrops) # 33538  1206
table(Sig=e.out_s_raw_emptydrops$FDR <= 0.01, Limited=e.out_s_raw_emptydrops$Limited)



## ------------------------------------------------------------------------
#create seurat object for raw
#s_raw_3_200 <- CreateSeuratObject(counts = expression_matrix_raw, min.cells = 3, min.features = 200)
s_raw_200 <- CreateSeuratObject(counts = expression_matrix_raw, min.features = 200)
dim(s_raw_200) #15254  1202
## ------------------------------------------------------------------------
#create seurat object for raw 
# s_raw_2_100 <- CreateSeuratObject(counts = expression_matrix_raw, min.cells = 2, min.features = 100)
# dim(s_raw_2_100) #16298  1228
## ------------------------------------------------------------------------
library(VennDiagram)
library(grid)
venn.plot <- venn.diagram(
  x = list(
    "filtered_cellranger3" = colnames(x = s_filtered_cellranger3),
    "raw_emptydrops" = colnames(x = s_raw_emptydrops),
    "raw_3_200" = colnames(x = s_raw_3_200),
    "raw_2_100" = colnames(x = s_raw_2_100)
    
    ),
  euler.d = TRUE,
  #filename=NULL,
  filename = "output_PBMC_1k_v3/_venn_cells.tiff",
  fil= c("red","blue","green", "yellow"),
  cex = 2,
  cat.cex = 0.7,
  reverse = TRUE
);
# grid.newpage();
# grid.draw(venn.plot)

venn.plot <- venn.diagram(
  x = list(
    "filtered_cellranger3" = rownames(x = s_filtered_cellranger3),
    "raw_emptydrops" = rownames(x = s_raw_emptydrops),
    "raw_3_200" = rownames(x = s_raw_3_200),
    "raw_2_100" = rownames(x = s_raw_2_100)
    
    ),
  euler.d = TRUE,
  #filename=NULL,
  filename = "output_PBMC_1k_v3/_venn_genes.tiff",
  fil= c("red","blue","green", "yellow"),
  cex = 2,
  cat.cex = 0.7,
  reverse = TRUE
);

## ------------------------------------------------------------------------
library("UpSetR")
my_list <- list(
    "s_filtered_cellranger3" = colnames(x = s_filtered_cellranger3),
    "s_raw_emptydrops" = colnames(x = s_raw_emptydrops),
    "s_raw_3_200" = colnames(x = s_raw_3_200),
    "s_raw_2_100" = colnames(x = s_raw_2_100))
cairo_pdf("output_PBMC_1k_v3/upset_initialfiltering.pdf")
upset(fromList(my_list),sets.bar.color = c("red","blue","green", "yellow"),
      main.bar.color = "black")
dev.off()

## ------------------------------------------------------------------------
library(R.utils)
xx.1 <- list(s_filtered_cellranger3 = colnames(x = s_filtered_cellranger3),
    s_raw_emptydrops = colnames(x = s_raw_emptydrops),
    s_raw_3_200 = colnames(x = s_raw_3_200),
    s_raw_2_100 = colnames(x = s_raw_2_100))

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
#Intersect(xx.1)
#Setdiff(xx.1["s_raw_2_100"], xx.1[c("s_raw_3_200","s_raw_emptydrops", "s_filtered_cellranger3")])

## ------------------------------------------------------------------------
Setdiff.s_filtered_cellranger3 <- Setdiff(xx.1["s_filtered_cellranger3"], xx.1[c("s_raw_3_200","s_raw_emptydrops", "s_raw_2_100")])
Is.Setdiff.s_filtered_cellranger3 <- colnames(s_filtered_cellranger3) %in% Setdiff.s_filtered_cellranger3
summary(Is.Setdiff.s_filtered_cellranger3)

## ------------------------------------------------------------------------
Setdiff.s_raw_2_100 <- Setdiff(xx.1["s_raw_2_100"], xx.1[c("s_raw_3_200","s_raw_emptydrops", "s_filtered_cellranger3")])
Is.Setdiff.s_raw_2_100 <- colnames(s_raw_2_100) %in% Setdiff.s_raw_2_100
summary(Is.Setdiff.s_raw_2_100)

## ------------------------------------------------------------------------
Setdiff.s_raw_emptydrops <- Setdiff(xx.1["s_raw_emptydrops"], xx.1[c("s_raw_3_200","s_filtered_cellranger3", "s_raw_2_100")])
Is.Setdiff.s_raw_emptydrops <- colnames(s_raw_emptydrops) %in% Setdiff.s_raw_emptydrops
summary(Is.Setdiff.s_raw_emptydrops)
## ------------------------------------------------------------------------
Setdiff.s_raw_2_100_filtered_cellranger3 <- Setdiff(xx.1[c("s_raw_2_100", "s_filtered_cellranger3")], xx.1[c("s_raw_emptydrops","s_raw_3_200")])
Is.Setdiff.s_raw_2_100_filtered_cellranger3 <- colnames(s_raw_2_100) %in% Setdiff.s_raw_2_100_filtered_cellranger3
summary(Is.Setdiff.s_raw_2_100_filtered_cellranger3)

## ------------------------------------------------------------------------
Setdiff.s_raw_emptydrops_2_100 <- Setdiff(xx.1[c("s_raw_emptydrops", "s_raw_2_100")], xx.1[c("s_raw_3_200","s_filtered_cellranger3")])
Is.Setdiff.s_raw_emptydrops_2_100 <- colnames(s_raw_2_100) %in% Setdiff.s_raw_emptydrops_2_100
summary(Is.Setdiff.s_raw_emptydrops_2_100)

## ------------------------------------------------------------------------
Setdiff.s_raw_2_100_3_200_emptydrops <- Setdiff(xx.1[c("s_raw_2_100","s_raw_3_200","s_raw_emptydrops")], xx.1["s_filtered_cellranger3"])
Is.Setdiff.s_raw_2_100_3_200_emptydrops<- colnames(s_raw_2_100) %in% Setdiff.s_raw_2_100_3_200_emptydrops
summary(Is.Setdiff.s_raw_2_100_3_200_emptydrops)

## ------------------------------------------------------------------------
Setdiff.s_raw_2_100_3_200_filtered_cellranger3 <- Setdiff(xx.1[c("s_raw_2_100","s_raw_3_200","s_filtered_cellranger3")], xx.1["s_raw_emptydrops"])
Is.Setdiff.s_raw_2_100_3_200_filtered_cellranger3<- colnames(s_filtered_cellranger3) %in% Setdiff.s_raw_2_100_3_200_filtered_cellranger3
summary(Is.Setdiff.s_raw_2_100_3_200_filtered_cellranger3)

## ------------------------------------------------------------------------
#intersect cells of all four initial filtering methods
Intersect.all <- Intersect(xx.1[c("s_raw_2_100","s_raw_emptydrops","s_filtered_cellranger3","s_raw_3_200")])
Is.Intersect.all<- colnames(s_filtered_cellranger3) %in% Intersect.all
summary(Is.Intersect.all)

## ------------------------------------------------------------------------
Union.all <- Union(xx.1[c("s_raw_2_100","s_raw_emptydrops","s_filtered_cellranger3","s_raw_3_200")])
Is.Union.all<- colnames(expression_matrix_raw) %in% Union.all
summary(Is.Union.all)
## ------------------------------------------------------------------------
cells.all <- expression_matrix_raw [ ,Union.all]
dim(cells.all)
## ------------------------------------------------------------------------
cells.Intersect.all <- cells.all[,Intersect.all]
cells.Setdiff.s_raw_2_100_3_200_filtered_cellranger3 <- cells.all[,Setdiff.s_raw_2_100_3_200_filtered_cellranger3]
cells.Setdiff.s_raw_2_100_3_200_emptydrops <- cells.all[,Setdiff.s_raw_2_100_3_200_emptydrops]
cells.Setdiff.s_raw_emptydrops_2_100 <- cells.all[,Setdiff.s_raw_emptydrops_2_100]
cells.Setdiff.s_raw_2_100_filtered_cellranger3 <- cells.all[,Setdiff.s_raw_2_100_filtered_cellranger3]
cells.Setdiff.s_raw_emptydrops <- cells.all[,Setdiff.s_raw_emptydrops]
cells.Setdiff.s_raw_2_100 <- cells.all[,Setdiff.s_raw_2_100]
cells.Setdiff.s_filtered_cellranger3 <- cells.all[,Setdiff.s_filtered_cellranger3]

## ------------------------------------------------------------------------
region0 = cells.all
region1=cells.Intersect.all 
region2=cells.Setdiff.s_raw_2_100_3_200_filtered_cellranger3 
region3=cells.Setdiff.s_raw_2_100_3_200_emptydrops 
region4=cells.Setdiff.s_raw_emptydrops_2_100 
region5=cells.Setdiff.s_raw_2_100_filtered_cellranger3 
region6=cells.Setdiff.s_raw_emptydrops 
region7=cells.Setdiff.s_raw_2_100 
region8=cells.Setdiff.s_filtered_cellranger3 

## ------------------------------------------------------------------------
cells.all <- CreateSeuratObject(cells.all)
region1 <- CreateSeuratObject(region1)
region2 <- CreateSeuratObject(region2)
region3 <- CreateSeuratObject(region3)
region4 <- CreateSeuratObject(region4)
region5 <- CreateSeuratObject(region5)
region6 <- CreateSeuratObject(region6)
region7 <- CreateSeuratObject(region7)
region8 <- CreateSeuratObject(region8)

region1$original <- region1$orig.ident
region1$orig.ident <- "Intersect all initial methods"
region2$original <- region2$orig.ident
region2$orig.ident <- "Seurat2-100,3-200 & cellranger"
region3$original <- region3$orig.ident
region3$orig.ident <- "Seurat2-100,3-200 & DropletUtils"
region4$original <- region4$orig.ident
region4$orig.ident <- "Seurat2-100 & DropletUtils"
region5$original <- region5$orig.ident
region5$orig.ident <- "Seurat2-100 & cellranger"
region6$original <- region6$orig.ident
region6$orig.ident <- "DropletUtils"
region7$original <- region7$orig.ident
region7$orig.ident <- "Seurat2-100"
region8$original <- region8$orig.ident
region8$orig.ident <- "cellranger"
cells.all$original <- cells.all$orig.ident
cells.all$orig.ident <- "cells.all"
#it includes all cells without any initial filtering of genes:
# combined.object <- merge(cells.all, y =c(region1,region2, region3, region4, region5, region6, region7,region8), add.cell.ids = c("cells.all","r1","r2", "r3","r4","r5","r6","r7","r8"),project = "initialfiltering")
combined.object <- merge(region1, y =c(region2, region3, region4, region5, region6, region7,region8),add.cell.ids = c("intersect all initial methods","Seurat2-100,3-200 & cellranger", "Seurat2-100,3-200 & DropletUtils","Seurat2-100 & DropletUtils","Seurat2-100 & cellranger","DropletUtils","Seurat2-100","cellranger"),project = "initialfiltering")
unique <-unique(sapply(X = strsplit(colnames(combined.object), split = "_"), FUN = "[", 1))
table(combined.object$orig.ident)

## ------------------------------------------------------------------------
saveRDS(combined.object, file = "output_PBMC_1k_v3/combined.object.rds")

## ------------------------------------------------------------------------
combined.object <- readRDS(file = "combined.object.rds")
## ------------------------------------------------------------------------

threshold = median(log(region1@meta.data$nFeature_RNA)) - 2*mad(log(region1@meta.data$nFeature_RNA))
threshold1 = median(log(region1@meta.data$nFeature_RNA)) - 3*mad(log(region1@meta.data$nFeature_RNA))
threshold2 = median(log(region1@meta.data$nFeature_RNA)) - 4*mad(log(region1@meta.data$nFeature_RNA))

plot1 = VlnPlot(object = region1, features = ("nFeature_RNA"), group.by = "orig.ident", ncol = 2, pt.size = 0.1)
plot1 + geom_hline(yintercept = exp(threshold1), color = "blue")
plot2 =  VlnPlot(object = region1, features = ("nCount_RNA"), group.by = "orig.ident", ncol = 2, pt.size = 0.1)
plot2 + geom_hline(yintercept = c(exp(threshold1), exp(threshold), exp(threshold2)), color = "blue")

## ------------------------------------------------------------------------
# plotExpression(b,rownames(b)[1:6], ncol = 3, exprs_values = "counts",
#                 colour_by = "orig.id")

combined.object[["percent.mt"]] <- PercentageFeatureSet(object = combined.object, pattern = "^MT-")
## ------------------------------------------------------------------------
library(ggplot2)
threshold = median(log(combined.object@meta.data$nFeature_RNA)) - 2*mad(log(combined.object@meta.data$nFeature_RNA))
threshold1 = median(log(combined.object@meta.data$nFeature_RNA)) - 3*mad(log(combined.object@meta.data$nFeature_RNA))
threshold2 = median(log(combined.object@meta.data$nFeature_RNA)) - 4*mad(log(combined.object@meta.data$nFeature_RNA))
thresholda = median(log(combined.object@meta.data$nCount_RNA)) - 2*mad(log(combined.object@meta.data$nCount_RNA))
thresholdb = median(log(combined.object@meta.data$nCount_RNA)) - 3*mad(log(combined.object@meta.data$nCount_RNA))
thresholdc = median(log(combined.object@meta.data$nCount_RNA)) - 4*mad(log(combined.object@meta.data$nCount_RNA))
#nmads for percent.mt???
thresholdi = median(log(combined.object@meta.data$percent.mt)) + 2*mad(log(combined.object@meta.data$percent.mt))
thresholdii = median(log(combined.object@meta.data$percent.mt)) + 3*mad(log(combined.object@meta.data$percent.mt))
thresholdiii = median(log(combined.object@meta.data$percent.mt)) + 4*mad(log(combined.object@meta.data$percent.mt))

plot1 <- FeatureScatter(object = combined.object, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident", pt.size = 2)+
  geom_point(alpha = 0.1)+
  geom_hline(yintercept = c(exp(thresholdi), exp(thresholdii), exp(thresholdiii)), linetype="dashed", color = "black")+
  annotate("text", x = Inf, y = c(exp(thresholdi), exp(thresholdii), exp(thresholdiii)), label = c("2mads", "3mads","4mads"), colour = "black",
           angle = 0, hjust = 1, vjust = -0.5)

  plot1 + geom_vline(xintercept = c(exp(thresholda), exp(thresholdb), exp(thresholdc)), color = "black")

plot2 <- FeatureScatter(object = combined.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", pt.size = 2)+
  geom_point(alpha = 0.1)+
  geom_hline(yintercept = c(exp(threshold), exp(threshold1), exp(threshold2)), linetype="dashed", color = "black")+
  annotate("text", x = Inf, y = c(exp(threshold), exp(threshold1), exp(threshold2)), label = c("2mads", "3mads","4mads"), colour = "black",
           angle = 0, hjust = 1, vjust = -0.5)+
  theme(legend.position = "none")

  geom_vline(xintercept = c(exp(thresholda), exp(thresholdb), exp(thresholdc)), linetype="dashed", color = "black")+
  annotate("text", y = 6600, x = c(exp(thresholda), exp(thresholdb), exp(thresholdc)), label = c("2mads", "3mads","4mads"), colour = "black",
           angle = 90, hjust =  -0.1, vjust = -0.5)

geom_hline(yintercept = c(exp(threshold), exp(threshold1), exp(threshold2)), color = "black")
## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
QC_metrics <-VlnPlot(object = combined.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
QC_metrics

# p = VlnPlot(object = combined.object, features = ("nFeature_RNA"), ncol = 2, pt.size =0.1 )+
#   geom_point(alpha = 0.1, aes(color = combined.object@meta.data$orig.ident))
# p
# DotPlot(object = combined.object, features = "nFeature_RNA", split.by = 'orig.ident', cols = )

## ------------------------------------------------------------------------
table(combined.object@meta.data$orig.ident)
## ------------------------------------------------------------------------
test<- SingleCellExperiment(list(counts=combined.object[["RNA"]]@counts))
# test <- as(counts(test), "dgCMatrix")
# test <- CreateSeuratObject(test)
## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
#nmads = 2,3 and 4 
library(scater)
#function:
marioni_filtering <- function(x){
patterns <- c("^hg19-MT-", "^mm10-mt-","^MT-", "^mt-")
is.mito<- grepl(paste(patterns, collapse="|"), rownames(x))
x<- calculateQCMetrics(x,feature_controls=list(Mt=is.mito)) 
libsize.drop <- isOutlier(x$total_counts, nmads= 2, type="lower", log=TRUE)
feature.drop <- isOutlier(x$total_features_by_counts, nmads= 2, type="lower", log=TRUE)
high.mito <- isOutlier(x$pct_counts_Mt, nmads= 2, type="higher")
keep <- !(libsize.drop | feature.drop | high.mito)
outliers <- (libsize.drop | feature.drop | high.mito)
outhigh.mito <- (high.mito)
s <- x[,keep]
s@metadata$PassQC <- keep
s@metadata$outliers.nmads <- outliers
s@metadata$outhigh.mito <- outhigh.mito
return(s)}
## ------------------------------------------------------------------------
test_nmads4 <- marioni_filtering(test)
test_nmads4 #33538 1114
nmads4 <- as(counts(test_nmads4), "dgCMatrix")
nmads4 <- CreateSeuratObject(nmads4)
test_nmads3 <- marioni_filtering(test)
test_nmads3 #33538 1085
nmads3 <- as(counts(test_nmads3), "dgCMatrix")
nmads3 <- CreateSeuratObject(nmads3)
test_nmads2 <- marioni_filtering(test)
test_nmads2 #33538 1042
nmads2 <- as(counts(test_nmads2), "dgCMatrix")
nmads2 <- CreateSeuratObject(nmads2)
## ------------------------------------------------------------------------
#filter_by_MT <- colData(test)$pct_counts_Mt < 10
## ------------------------------------------------------------------------

#compare 3 different nmads in one umap plot with different colors
# s_test_nmads2_genes.markers$original.1 <- s_test_nmads2_genes.markers$orig.ident
# s_test_nmads2_genes.markers$orig.ident <- "nmads2"
# 
# s_test_nmads3_genes.markers$original.1 <- s_test_nmads3_genes.markers$orig.ident
# s_test_nmads3_genes.markers$orig.ident <- "nmads3"
# 
# s_test_nmads4_genes.markers$original.1 <- s_test_nmads4_genes.markers$orig.ident
# s_test_nmads4_genes.markers$orig.ident <- "nmads4"


 # nmads_merge.plot <- DimPlot(object = nmads_merge, reduction = "umap", group.by = "orig.ident")
## ------------------------------------------------------------------------
# x4 <- marioni_filtering_genes(test_nmads4)
# x4
# x3 <- marioni_filtering_genes(test_nmads3)
# x3
# x2 <- marioni_filtering_genes(test_nmads2)
# x2
# x4 <- as(counts(x4), "dgCMatrix")
# x4 <- CreateSeuratObject(x4)
# x3 <- as(counts(x3), "dgCMatrix")
# x3 <- CreateSeuratObject(x3)
# x2 <- as(counts(x2), "dgCMatrix")
# x2 <- CreateSeuratObject(x2)
# 
# x2@meta.data[, "nmads"] <- "nmads2"
# x3@meta.data[, "nmads"] <- "nmads3"
# x4@meta.data[, "nmads"] <- "nmads4"
# # 
# nmads_merge <- merge(x = x4,
#                      y = c(x3, x2),
#                      add.cell.ids = c("nmads4","nmads3", "nmads2"),
#                      project = "filtering")
# # 
# # 
# # tail(colnames(nmads_merge))
# # head(colnames(nmads_merge))
# head(nmads_merge@meta.data)
# nmads_merge[["percent.mt"]] <- PercentageFeatureSet(object = nmads_merge, pattern = "^MT-")
## ------------------------------------------------------------------------
#normalization
##lognormalize
# nmads.genes.log <- nmads_merge
# NormalizeData(object = nmads.genes.log, normalization.method = "LogNormalize", 
#               scale.factor = 10000)
# nmads.genes.log <- FindVariableFeatures(object = nmads.genes.log,selection.method = 'vst')
# nmads.genes.log <- ScaleData(object = nmads.genes.log, vars.to.regress = 'percent.mt')
# nmads.genes.log <- RunPCA(object = nmads.genes.log)
# nmads.genes.log<- FindNeighbors(object = nmads.genes.log, dims = 1:20)
# nmads.genes.log <- FindClusters(object = nmads.genes.log, resolution = c(0.4,0.5,0.6,0.8,1.2))
# nmads.genes.log <- RunUMAP(object = nmads.genes.log, dims = 1:20)
# 
# nmads.genes.log.plot <- DimPlot(object = nmads.genes.log, reduction = "umap", group.by = "nmads")
# nmads.genes.log.plot
# 
# plot.list <- list()
# for (i in unique(x = nmads.genes.log@meta.data$nmads)) {
#   plot.list[[i]] <- DimPlot(
#     object = nmads.genes.log, 
#     cells.highlight = WhichCells(object = nmads.genes.log, expression = nmads == i)
#   ) + NoLegend() + ggtitle(i)
# }
# myy_plot <- CombinePlots(plots = plot.list, ncol = 3)
## ------------------------------------------------------------------------
#function
marioni_filtering_genes <- function(x){
  ave.counts <- calcAverage(object = x, use_size_factors=FALSE)
  rowData(x)$ave.count <- ave.counts
  demo.keep <- ave.counts >0
  s <- x[demo.keep,]
  return(s)
}
## ------------------------------------------------------------------------
test_nmads3_genes <- marioni_filtering_genes(test_nmads3)
test_nmads3_genes

## ------------------------------------------------------------------------
test_nmads4_genes <- marioni_filtering_genes(test_nmads4)
test_nmads4_genes

## ------------------------------------------------------------------------
test_nmads2_genes <- marioni_filtering_genes(test_nmads2)
test_nmads2_genes

## ------------------------------------------------------------------------
test_nmads4_genes <- as(counts(test_nmads4_genes), "dgCMatrix")
test_nmads4_genes <- CreateSeuratObject(test_nmads4_genes)
test_nmads4_genes[["percent.mt"]] <- PercentageFeatureSet(object = test_nmads4_genes, pattern = "^MT-")
test_nmads4_genes

test_nmads3_genes <- as(counts(test_nmads3_genes), "dgCMatrix")
test_nmads3_genes <- CreateSeuratObject(test_nmads3_genes)
test_nmads3_genes[["percent.mt"]] <- PercentageFeatureSet(object = test_nmads3_genes, pattern = "^MT-")
test_nmads3_genes

test_nmads2_genes <- as(counts(test_nmads2_genes), "dgCMatrix")
test_nmads2_genes <- CreateSeuratObject(test_nmads2_genes)
test_nmads2_genes[["percent.mt"]] <- PercentageFeatureSet(object = test_nmads2_genes, pattern = "^MT-")
test_nmads2_genes
## ------------------------------------------------------------------------

dim(test_nmads3_genes)
saveRDS(test_nmads3_genes, file = "test_nmads3_genes.rds")

dim(test_nmads4_genes)
saveRDS(test_nmads4_genes, file = "test_nmads4_genes.rds")

dim(test_nmads2_genes)
saveRDS(test_nmads2_genes, file = "test_nmads2_genes.rds")
## ------------------------------------------------------------------------
library(sctransform)
test_nmads2_genes_sctransform <- test_nmads2_genes
test_nmads3_genes_sctransform <- test_nmads3_genes
test_nmads4_genes_sctransform <- test_nmads4_genes

## ------------------------------------------------------------------------
test_nmads2_genes_sctransform <- SCTransform(object = test_nmads2_genes_sctransform,vars.to.regress = 'percent.mt')
test_nmads2_genes_sctransform <- RunPCA(object = test_nmads2_genes_sctransform)
test_nmads2_genes_sctransform <- RunUMAP(object = test_nmads2_genes_sctransform, dims = 1:15)
test_nmads2_genes_sctransform <- RunTSNE(object = test_nmads2_genes_sctransform, dims = 1:15)
test_nmads2_genes_sctransform <- FindNeighbors(object = test_nmads2_genes_sctransform, dims = 1:15)
test_nmads2_genes_sctransform <- FindClusters(object = test_nmads2_genes_sctransform, resolution = c(0.4,0.5,0.6, 0.8,1.2))
umap_nmads2.sc <- DimPlot(object = test_nmads2_genes_sctransform, reduction = 'umap', pt.size = 1, label = TRUE)+ggtitle(paste0("nmads2\n",'SCTransform'))
find.all.markers.test_nmads2_genes_sctransform <- FindAllMarkers(object = test_nmads2_genes_sctransform, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
write.csv(find.all.markers.test_nmads2_genes_sctransform,'output_PBMC_1k_v3/all_test_nmads2_genes_sctransform.markers.csv')
## ------------------------------------------------------------------------
test_nmads3_genes_sctransform <- SCTransform(object = test_nmads3_genes_sctransform, vars.to.regress = 'percent.mt')
test_nmads3_genes_sctransform <- RunPCA(object = test_nmads3_genes_sctransform)
test_nmads3_genes_sctransform <- RunUMAP(object = test_nmads3_genes_sctransform, dims = 1:15)
test_nmads3_genes_sctransform <- RunTSNE(object = test_nmads3_genes_sctransform, dims = 1:15)
test_nmads3_genes_sctransform <- FindNeighbors(object = test_nmads3_genes_sctransform, dims = 1:15)
test_nmads3_genes_sctransform <- FindClusters(object = test_nmads3_genes_sctransform, resolution = c(0.4,0.5,0.6, 0.8,1.2))
umap_nmads3.sc <- DimPlot(object = test_nmads3_genes_sctransform, reduction = 'umap', pt.size = 1, label = TRUE)+ggtitle(paste0("nmads3\n",'SCTransform'))
find.all.markers.test_nmads3_genes_sctransform <- FindAllMarkers(object = test_nmads3_genes_sctransform, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
write.csv(find.all.markers.test_nmads3_genes_sctransform,'output_PBMC_1k_v3/all_test_nmads3_genes_sctransform.markers.csv')
## ------------------------------------------------------------------------
test_nmads4_genes_sctransform <- SCTransform(object = test_nmads4_genes_sctransform, vars.to.regress = 'percent.mt')
test_nmads4_genes_sctransform <- RunPCA(object = test_nmads4_genes_sctransform)
test_nmads4_genes_sctransform <- RunUMAP(object = test_nmads4_genes_sctransform, dims = 1:15)
test_nmads4_genes_sctransform <- RunTSNE(object = test_nmads4_genes_sctransform, dims = 1:15)
test_nmads4_genes_sctransform <- FindNeighbors(object = test_nmads4_genes_sctransform, dims = 1:15)
test_nmads4_genes_sctransform <- FindClusters(object = test_nmads4_genes_sctransform, resolution = c(0.4,0.5,0.6, 0.8,1.2))
umap_nmads4.sc <- DimPlot(object = test_nmads4_genes_sctransform, reduction = 'umap', pt.size = 1, label = TRUE)+ggtitle(paste0("nmads4\n",'SCTransform'))
find.all.markers.test_nmads4_genes_sctransform <- FindAllMarkers(object = test_nmads4_genes_sctransform, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
write.csv(find.all.markers.test_nmads4_genes_sctransform,'output_PBMC_1k_v3/all_test_nmads4_genes_sctransform.markers.csv')
## ------------------------------------------------------------------------
#for all test_nmads4_genes, test_nmads3_genes and test_nmads2_genes:
setwd("~/Documents/Thesis")
library(dplyr)
x <- test_nmads3_genes 
s_1 <-  NormalizeData(object = x, normalization.method = "LogNormalize", 
                    scale.factor = 10000)
#feature selection
#s_1 <- FindVariableFeatures(object = s_1,selection.method = 'vst', nfeatures = 2000)
s_1 <- FindVariableFeatures(object = s_1,selection.method = 'vst')
#all.genes <- rownames(x = s_1)
#length(s_1@assays$RNA@var.features)
s_1 <- ScaleData(object = s_1, vars.to.regress = 'percent.mt')
#runpca only on those variable genes:
s_1 <- RunPCA(object = s_1)
print(x = s_1[['pca']], dims = 1:5, nfeatures = 10)
#VizDimLoadings(object = s_1, dims = 1:2, reduction = 'pca')
cairo_pdf("output_PBMC_1k_v3/pca_test_nmads2_genes.pdf")
DimPlot(object = s_1, reduction = 'pca', pt.size = 2)
dev.off()
# s_1 <- JackStraw(object = s_1, num.replicate = 100)
# s_1 <- ScoreJackStraw(object = s_1, dims = 1:20)
# cairo_pdf("output_PBMC_1k_v3/JackStraw_test_nmads3_genes.pdf")
# JackStrawPlot(object = s_1, dims = 1:15)
# dev.off()
#ElbowPlot(object = s_1)
s_1 <- FindNeighbors(object = s_1, dims = 1:15)
#set a range of resolution for clustering:
s_1 <- FindClusters(object = s_1, resolution = c(0.4,0.5,0.6,0.8,1.2))
# Look at cluster IDs of the first 5 cells
head(x = Idents(object = s_1), 5)
levels(x = s_1)
sapply(grep("^RNA_snn_res",colnames(s_1@meta.data),value = TRUE),function(x) length(unique(s_1@meta.data[,x])))
#runumap:
#reticulate::py_install(packages = "umap-learn")
s_1 <- RunTSNE(object = s_1, dims = 1:15)
s_1 <- RunUMAP(object = s_1, dims = 1:15)
cairo_pdf("output_PBMC_1k_v3/umap_s_test_nmads2_genes.pdf")
DimPlot(object = s_1, reduction = 'umap', pt.size = 1)
dev.off()
#track cells from upstream analysis
s_1test <- s_1
Idents(s_1test) <- "orig.ident"
levels(s_1test)
cairo_pdf("output_PBMC_1k_v3/umap_s_test_nmads2_genes_track.pdf")
DimPlot(object = s_1test, reduction = 'umap',cols = c("grey", "red"), pt.size = 1) 
dev.off()
###``` Finding differentially expressed features (cluster biomarkers)```
# find markers for every cluster compared to all remaining cells, report only the positive ones #wilcox test is used
find.all.markers <- FindAllMarkers(object = s_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
write.csv(find.all.markers,'output_PBMC_1k_v3/all_s_test_nmads2_genes.markers.csv')
find.all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#find.all.markers$gene
top10.s_test_nmads4_genes.markers<- find.all.markers %>% group_by(cluster) %>% top_n(10, wt = avg_logFC)

s_test_nmads2_genes.markers <- s_1
## ------------------------------------------------------------------------
saveRDS(s_test_nmads2_genes.markers, file = "output_PBMC_1k_v3/s_test_nmads2_genes.markers.rds")
saveRDS(s_test_nmads3_genes.markers, file = "output_PBMC_1k_v3/s_test_nmads3_genes.markers.rds")
saveRDS(s_test_nmads4_genes.markers, file = "output_PBMC_1k_v3/s_test_nmads4_genes.markers.rds")

s_test_nmads2_genes.markers <- readRDS("s_test_nmads2_genes.markers.rds")
s_test_nmads3_genes.markers <- readRDS("s_test_nmads3_genes.markers.rds")
s_test_nmads4_genes.markers <- readRDS("s_test_nmads4_genes.markers.rds")

s_test_nmads2_genes.markers
s_test_nmads3_genes.markers
s_test_nmads4_genes.markers
## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
umap_nmads2.log <- DimPlot(object = s_test_nmads2_genes.markers, reduction = 'umap', pt.size = 1, label =TRUE)+ggtitle(paste0("nmads2\n",'LogNormalize'))
umap_nmads3.log <- DimPlot(object = s_test_nmads3_genes.markers, reduction = 'umap', pt.size = 1, label = TRUE)+ggtitle(paste0("nmads3\n",'LogNormalize'))
umap_nmads4.log <- DimPlot(object = s_test_nmads4_genes.markers, reduction = 'umap', pt.size = 1, label = TRUE)+ggtitle(paste0("nmads4\n",'LogNormalize'))
## ------------------------------------------------------------------------
library(gridExtra)
res.dir <- paste0("/Volumes/SD Card/Thesis/output_PBMC_1k_v3/")
dir.create(res.dir) 
ggsave(
  grid.arrange(umap_nmads2.log, umap_nmads2.sc,
               umap_nmads3.log, umap_nmads3.sc,
               umap_nmads4.log, umap_nmads4.sc,ncol=2),
  file=paste0(res.dir,"01_umapclustering_normd.pdf"), width = 15,height=20)

## ------------------------------------------------------------------------
s_2.sc <- test_nmads2_genes_sctransform
Idents(s_2.sc) <- "orig.ident"
levels(s_2.sc)
s_2.sc <- DimPlot(object = s_2.sc, reduction = 'umap',cols = c("grey", "red"), pt.size = 1)+ ggtitle(paste0("nmads2\n",'sctransform'))

s_3.sc <- test_nmads3_genes_sctransform
Idents(s_3.sc) <- "orig.ident"
levels(s_3.sc)
s_3.sc <- DimPlot(object = s_3.sc, reduction = 'umap',cols = c("grey", "red"), pt.size = 1)+ ggtitle(paste0("nmads3\n",'sctransform'))

s_4.sc <- test_nmads4_genes_sctransform
Idents(s_4.sc) <- "orig.ident"
levels(s_4.sc)
s_4.sc <- DimPlot(object = s_4.sc, reduction = 'umap',cols = c("grey", "red"), pt.size = 1)+ ggtitle(paste0("nmads4\n",'sctransform'))

s_2.log <- s_test_nmads2_genes.markers
Idents(s_2.log) <- "orig.ident"
levels(s_2.log)
s_2.log <- DimPlot(object = s_2.log, reduction = 'umap',cols = c("grey", "red"), pt.size = 1)+ ggtitle(paste0("nmads2\n",'LogNormalize'))

s_3.log <- s_test_nmads3_genes.markers
Idents(s_3.log) <- "orig.ident"
levels(s_3.log)
s_3.log <- DimPlot(object = s_3.log, reduction = 'umap',cols = c("grey", "red"), pt.size = 1)+ ggtitle(paste0("nmads3\n",'LogNormalize'))

s_4.log <- s_test_nmads4_genes.markers
Idents(s_4.log) <- "orig.ident"
levels(s_4.log)
s_4.log<- DimPlot(object = s_4.log, reduction = 'umap',cols = c("grey", "red"), pt.size = 1)+ ggtitle(paste0("nmads4\n",'LogNormalize'))

library(gridExtra)
res.dir <- paste0("/Volumes/SD Card/Thesis/output_PBMC_1k_v3/")
dir.create(res.dir) 
ggsave(
  grid.arrange(s_2.log, s_2.sc,
               s_3.log, s_3.sc,
               s_4.log, s_4.sc,ncol=2),
  file=paste0(res.dir,"01_umaptrack_normdiag.pdf"), width = 15,height=20)
## ------------------------------------------------------------------------
#COMPARE 2 NORMALIZATION METHODS
## ------------------------------------------------------------------------
s_test_nmads2_genes.markers$clusterID <- Idents(test_nmads2_genes_sctransform)
Idents(s_test_nmads2_genes.markers) <- 'clusterID'
plot1 <- DimPlot(object = test_nmads2_genes_sctransform, label = TRUE, cols = rainbow(13))+ ggtitle(paste0("nmads2\n",'sctransform'))  
plot2 <- DimPlot(object = s_test_nmads2_genes.markers, label = TRUE, cols = rainbow(13))
plot2 <- plot2 + ggtitle(paste0("nmads2\n",'Log-normalization')) 
CombinePlots(list(plot1,plot2))
ggsave("output_PBMC_1k_v3/s_test_nmads2_genes.markers.norm.pdf",width = 15,height=10, dpi = 320)
## ------------------------------------------------------------------------
s_test_nmads3_genes.markers$clusterID <- Idents(test_nmads3_genes_sctransform)
Idents(s_test_nmads3_genes.markers) <- 'clusterID'
plot1 <- DimPlot(object = test_nmads3_genes_sctransform, label = TRUE, cols = rainbow(12))+ ggtitle(paste0("nmads3\n",'sctransform'))  
plot2 <- DimPlot(object = s_test_nmads3_genes.markers, label = TRUE, cols = rainbow(12))
plot2 <- plot2 + ggtitle(paste0("nmads3\n",'Log-normalization')) 
CombinePlots(list(plot1,plot2))
ggsave("output_PBMC_1k_v3/s_test_nmads3_genes.markers.norm.pdf",width = 15,height=10, dpi = 320)
## ------------------------------------------------------------------------
s_test_nmads4_genes.markers$clusterID <- Idents(test_nmads4_genes_sctransform)
Idents(s_test_nmads4_genes.markers) <- 'clusterID'
plot1 <- DimPlot(object = test_nmads4_genes_sctransform, label = TRUE, cols = rainbow(12))+ ggtitle(paste0("nmads4\n",'sctransform'))  
plot2 <- DimPlot(object = s_test_nmads4_genes.markers, label = TRUE, cols = rainbow(12))
plot2 <- plot2 + ggtitle(paste0("nmads4\n",'Log-normalization')) 
CombinePlots(list(plot1,plot2))
ggsave("output_PBMC_1k_v3/s_test_nmads4_genes.markers.norm.pdf",width = 15,height=10, dpi = 320)
## ------------------------------------------------------------------------
pca.umi.ln2 <- diagplot(s_test_nmads2_genes.markers,'pca','umi')
pca.umi.sct2 <- diagplot(test_nmads2_genes_sctransform,'pca','umi')
CombinePlots(list(pca.umi.ln2 ,pca.umi.sct2))
ggsave("output_PBMC_1k_v3/s_test_nmads2_genes.markers.normpca.pdf",width = 15,height=10, dpi = 320)

umap.umi.ln2 <- diagplot(s_test_nmads2_genes.markers,'umap','umi')
umap.umi.sct2 <- diagplot(test_nmads2_genes_sctransform,'umap','umi')
CombinePlots(list(umap.umi.ln2 ,umap.umi.sct2))
ggsave("output_PBMC_1k_v3/s_test_nmads2_genes.markers.normumap.pdf",width = 15,height=10, dpi = 320)

umap.mito.ln2 <- diagplot_(s_test_nmads2_genes.markers,'umap','percent.mt')
umap.mito.sct2 <- diagplot_(test_nmads2_genes_sctransform,'umap','percent.mt')
CombinePlots(list(umap.mito.ln2 ,umap.mito.sct2))
ggsave("output_PBMC_1k_v3/s_test_nmads2_genes.markers.mitoumap.pdf",width = 15,height=10, dpi = 320)

pca.mito.ln2 <- diagplot_(s_test_nmads2_genes.markers,'pca','percent.mt')
pca.mito.sct2 <- diagplot_(test_nmads2_genes_sctransform,'pca','percent.mt')
CombinePlots(list(pca.mito.ln2 ,pca.mito.sct2))
ggsave("output_PBMC_1k_v3/s_test_nmads2_genes.markers.mitopca.pdf",width = 15,height=10, dpi = 320)

## ------------------------------------------------------------------------
pca.umi.ln3 <- diagplot(s_test_nmads3_genes.markers,'pca','umi')
pca.umi.sct3 <- diagplot(test_nmads3_genes_sctransform,'pca','umi')
CombinePlots(list(pca.umi.ln3 ,pca.umi.sct3))
ggsave("output_PBMC_1k_v3/s_test_nmads3_genes.markers.normpca.pdf",width = 15,height=10, dpi = 320)

umap.umi.ln3 <- diagplot(s_test_nmads3_genes.markers,'umap','umi')
umap.umi.sct3 <- diagplot(test_nmads3_genes_sctransform,'umap','umi')
CombinePlots(list(umap.umi.ln3 ,umap.umi.sct3))
ggsave("output_PBMC_1k_v3/s_test_nmads3_genes.markers.normumap.pdf",width = 15,height=10, dpi = 320)

umap.mito.ln3 <- diagplot_(s_test_nmads3_genes.markers,'umap','percent.mt')
umap.mito.sct3 <- diagplot_(test_nmads3_genes_sctransform,'umap','percent.mt')
CombinePlots(list(umap.mito.ln3 ,umap.mito.sct3))
ggsave("output_PBMC_1k_v3/s_test_nmads3_genes.markers.mitoumap.pdf",width = 15,height=10, dpi = 320)

pca.mito.ln3 <- diagplot_(s_test_nmads3_genes.markers,'pca','percent.mt')
pca.mito.sct3 <- diagplot_(test_nmads3_genes_sctransform,'pca','percent.mt')
CombinePlots(list(pca.mito.ln3 ,pca.mito.sct3))
ggsave("output_PBMC_1k_v3/s_test_nmads3_genes.markers.mitopca.pdf",width = 15,height=10, dpi = 320)
## ------------------------------------------------------------------------
pca.umi.ln4 <- diagplot(s_test_nmads4_genes.markers,'pca','umi')
pca.umi.sct4 <- diagplot(test_nmads4_genes_sctransform,'pca','umi')
CombinePlots(list(pca.umi.ln4 ,pca.umi.sct4))
ggsave("output_PBMC_1k_v3/s_test_nmads4_genes.markers.normpca.pdf",width = 15,height=10, dpi = 320)

umap.umi.ln4 <- diagplot(s_test_nmads4_genes.markers,'umap','umi')
umap.umi.sct4 <- diagplot(test_nmads4_genes_sctransform,'umap','umi')
CombinePlots(list(umap.umi.ln4 ,umap.umi.sct4))
ggsave("output_PBMC_1k_v3/s_test_nmads4_genes.markers.normumap.pdf",width = 15,height=10, dpi = 320)

umap.mito.ln4 <- diagplot_(s_test_nmads4_genes.markers,'umap','percent.mt')
umap.mito.sct4 <- diagplot_(test_nmads4_genes_sctransform,'umap','percent.mt')
CombinePlots(list(umap.mito.ln4 ,umap.mito.sct4))
ggsave("output_PBMC_1k_v3/s_test_nmads4_genes.markers.mitoumap.pdf",width = 15,height=10, dpi = 320)

pca.mito.ln4 <- diagplot_(s_test_nmads4_genes.markers,'pca','percent.mt')
pca.mito.sct4 <- diagplot_(test_nmads4_genes_sctransform,'pca','percent.mt')
CombinePlots(list(pca.mito.ln4 ,pca.mito.sct4))
ggsave("output_PBMC_1k_v3/s_test_nmads4_genes.markers.mitopca.pdf",width = 15,height=10, dpi = 320)
## ------------------------------------------------------------------------

library(gridExtra)
res.dir <- paste0("/Volumes/SD Card/Thesis/output_PBMC_1k_v3/")
dir.create(res.dir) 
ggsave(
  grid.arrange(pca.umi.ln2, pca.umi.sct2,
               pca.umi.ln3, pca.umi.sct3,
               pca.umi.ln4, pca.umi.sct4,ncol=2),
  file=paste0(res.dir,"01_pca_normdiag.pdf"), width = 15,height=20)

ggsave(
  grid.arrange(pca.mito.ln2, pca.mito.sct2,
               pca.mito.ln3, pca.mito.sct3,
               pca.mito.ln4, pca.mito.sct4,ncol=2),
  file=paste0(res.dir,"01_pca_mitodiag.pdf"), width = 15,height=20)

ggsave(
  grid.arrange(umap.umi.ln2, umap.umi.sct2,
               umap.umi.ln3, umap.umi.sct3,
               umap.umi.ln4, umap.umi.sct4,ncol=2),
  file=paste0(res.dir,"01_umap_normdiag.pdf"), width = 15,height=20)

ggsave(
  grid.arrange(umap.mito.ln2, umap.mito.sct2,
               umap.mito.ln3, umap.mito.sct3,
               umap.mito.ln4, umap.mito.sct4,ncol=2),
  file=paste0(res.dir,"01_umap_mitodiag.pdf"), width = 15,height=20)
## ------------------------------------------------------------------------
s_test_nmads4_genes.markers <- RenameCells(
    object = s_test_nmads4_genes.markers,
    new.names = sub("^r1_", "", colnames(s_test_nmads4_genes.markers))
)
s_test_nmads4_genes.markers <- RenameCells(
    object = s_test_nmads4_genes.markers,
    new.names = sub("^r3_", "", colnames(s_test_nmads4_genes.markers))
)
head(colnames(s_test_nmads4_genes.markers))
tail(colnames(s_test_nmads4_genes.markers))

## ------------------------------------------------------------------------
#devtools::install_github('dviraran/SingleR')
## ------------------------------------------------------------------------
s = seq(1,length(colnames(s_test_nmads4_genes.markers)),by=20000)
s.s = seq(1,length(rownames(s_test_nmads4_genes.markers)),by=20000)
for (i in s) {
  A = seq(i,min(i+20000-1,length(colnames(s_test_nmads4_genes.markers))))
  print(A)
}
for (i in s.s) {
  B = seq(i,min(i+20000-1,length(rownames(s_test_nmads4_genes.markers))))
  print(B)
}

## ------------------------------------------------------------------------
#singleR to help interpret cell clusters and annotate cells. it will be done for s_test_nmads4_genes.markers, s_test_nmads3_genes.markers and s_test_nmads2_genes.markers
library(SingleR)
singler4.log = CreateSinglerObject(s_test_nmads4_genes.markers[["RNA"]]@counts, annot = NULL,project.name= "initial_filtering",min.genes = 0,
  technology = "10X", species = "Human", citation = "",
  ref.list = list(), normalize.gene.length = F, variable.genes = "de",
  fine.tune = T, do.signatures = T, clusters = s_test_nmads4_genes.markers@active.ident, do.main.types = T, 
  reduce.file.size = T, numCores = SingleR.numCores)

## ------------------------------------------------------------------------
singler3.log = CreateSinglerObject(s_test_nmads3_genes.markers[["RNA"]]@counts, annot = NULL,project.name= "initial_filtering",min.genes = 0,
  technology = "10X", species = "Human", citation = "",
  ref.list = list(), normalize.gene.length = F, variable.genes = "de",
  fine.tune = T, do.signatures = T, clusters = s_test_nmads3_genes.markers@active.ident, do.main.types = T, 
  reduce.file.size = T, numCores = SingleR.numCores)

## ------------------------------------------------------------------------
singler2.log = CreateSinglerObject(s_test_nmads2_genes.markers[["RNA"]]@counts, annot = NULL,project.name= "initial_filtering",min.genes = 0,
  technology = "10X", species = "Human", citation = "",
  ref.list = list(), normalize.gene.length = F, variable.genes = "de",
  fine.tune = T, do.signatures = T, clusters = s_test_nmads2_genes.markers@active.ident, do.main.types = T, 
  reduce.file.size = T, numCores = SingleR.numCores)

## ------------------------------------------------------------------------
singler4.sc = CreateSinglerObject(test_nmads4_genes_sctransform[["RNA"]]@counts, annot = NULL,project.name= "initial_filtering",min.genes = 0,
                                   technology = "10X", species = "Human", citation = "",
                                   ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                                   fine.tune = T, do.signatures = T, clusters = test_nmads4_genes_sctransform@active.ident, do.main.types = T, 
                                   reduce.file.size = T, numCores = SingleR.numCores)
## ------------------------------------------------------------------------
singler3.sc = CreateSinglerObject(test_nmads3_genes_sctransform[["RNA"]]@counts, annot = NULL,project.name= "initial_filtering",min.genes = 0,
                                   technology = "10X", species = "Human", citation = "",
                                   ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                                   fine.tune = T, do.signatures = T, clusters = test_nmads3_genes_sctransform@active.ident, do.main.types = T, 
                                   reduce.file.size = T, numCores = SingleR.numCores)
## ------------------------------------------------------------------------
singler2.sc = CreateSinglerObject(test_nmads2_genes_sctransform[["RNA"]]@counts, annot = NULL,project.name= "initial_filtering",min.genes = 0,
                                   technology = "10X", species = "Human", citation = "",
                                   ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                                   fine.tune = T, do.signatures = T, clusters = test_nmads2_genes_sctransform@active.ident, do.main.types = T, 
                                   reduce.file.size = T, numCores = SingleR.numCores)
## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
#transfering singler cell labels to seurat object Idents to show cells names(annotated cells) in UMAP plot:
singler4.loglabels <- s_test_nmads4_genes.markers
head(s_test_nmads4_genes.markers$orig.ident)
Idents(singler4.loglabels)<- singler4.log$singler[[2]]$SingleR.single.main$labels
head(Idents(singler4.loglabels))
singler4.loglabels <- DimPlot(object = singler4.loglabels, reduction = 'umap', pt.size = 1)+ ggtitle(paste0("nmads4\n",'logNormalize'))
# ggsave("test.pdf",width = 15,height=10, dpi = 320)
## ------------------------------------------------------------------------
singler3.loglabels <- s_test_nmads3_genes.markers
head(s_test_nmads3_genes.markers$orig.ident)
Idents(singler3.loglabels)<- singler3.log$singler[[2]]$SingleR.single.main$labels
head(Idents(singler3.loglabels))
singler3.loglabels <- DimPlot(object = singler3.loglabels, reduction = 'umap', pt.size = 1)+ ggtitle(paste0("nmads3\n",'logNormalize'))
## ------------------------------------------------------------------------
singler2.loglabels <- s_test_nmads2_genes.markers
head(s_test_nmads2_genes.markers$orig.ident)
Idents(singler2.loglabels)<- singler2.log$singler[[2]]$SingleR.single.main$labels
head(Idents(singler2.loglabels))
singler2.loglabels <- DimPlot(object = singler2.loglabels, reduction = 'umap', pt.size = 1)+ ggtitle(paste0("nmads2\n",'logNormalize'))

## ------------------------------------------------------------------------
ggsave(
  grid.arrange(umap_nmads2.log, singler2.loglabels,
               umap_nmads3.log, singler3.loglabels,
               umap_nmads4.log, singler4.loglabels,ncol=2),
  file=paste0(res.dir,"01_umapsingler_lognorm.pdf"), width = 20,height=20)
## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
#transfering singler cell labels to seurat object Idents to show cells names(annotated cells) in UMAP plot: for sctransform normalization:
singler4.sclabels <- test_nmads4_genes_sctransform
head(test_nmads4_genes_sctransform$orig.ident)
Idents(singler4.sclabels)<- singler4.sc$singler[[2]]$SingleR.single.main$labels
head(Idents(singler4.sclabels))
singler4.sclabels <- DimPlot(object = singler4.sclabels, reduction = 'umap', pt.size = 1)+ ggtitle(paste0("nmads4\n",'SCTransform'))
## ------------------------------------------------------------------------
singler3.sclabels <- test_nmads3_genes_sctransform
head(test_nmads3_genes_sctransform$orig.ident)
Idents(singler3.sclabels)<- singler3.sc$singler[[2]]$SingleR.single.main$labels
head(Idents(singler3.sclabels))
singler3.sclabels <- DimPlot(object = singler3.sclabels, reduction = 'umap', pt.size = 1)+ ggtitle(paste0("nmads3\n",'SCTransform'))
## ------------------------------------------------------------------------
singler2.sclabels <- test_nmads2_genes_sctransform
head(test_nmads2_genes_sctransform$orig.ident)
Idents(singler2.sclabels)<- singler2.sc$singler[[2]]$SingleR.single.main$labels
head(Idents(singler2.sclabels))
singler2.sclabels <- DimPlot(object = singler2.sclabels, reduction = 'umap', pt.size = 1)+ ggtitle(paste0("nmads2\n",'SCTransform'))
## ------------------------------------------------------------------------
ggsave(
  grid.arrange(umap_nmads2.sc, singler2.sclabels,
               umap_nmads3.sc, singler3.sclabels,
               umap_nmads4.sc, singler4.sclabels,ncol=2),
  file=paste0(res.dir,"01_umapsingler_scnorm.pdf"), width = 20,height=20)
## ------------------------------------------------------------------------
ggsave(
  grid.arrange(singler2.loglabels, singler2.sclabels,
               singler3.loglabels, singler3.sclabels,
               singler4.loglabels, singler4.sclabels,ncol=2),
  file=paste0(res.dir,"01_umapsingler_lognorm_scnorm.pdf"), width = 20,height=20)

## ------------------------------------------------------------------------
# out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,
#                        singler$meta.data$xy,do.label = F,
#                        do.letters = T,labels = singler$meta.data$orig.ident, 
#                        dot.size = 3.3,alpha=0.5,label.size = 6)
# out$p
# ## ------------------------------------------------------------------------
# out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,
#                        singler$meta.data$xy,do.label = T,
#                        do.letters = F,labels=singler$singler[[2]]$SingleR.single.main$labels, 
#                        dot.size = 3.3,label.size = 5,alpha=0.5, font.size = 4)
# out$p
# ## ------------------------------------------------------------------------
# out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,
#                        singler$meta.data$xy,do.label = T,
#                        do.letters = F,labels=singler$seurat@active.ident, 
#                        dot.size = 1.3,label.size = 5,alpha=0.5, font.size = 6)
# out$p
## ------------------------------------------------------------------------

# kable(table(singler$singler[[2]]$SingleR.single.main$labels,singler$seurat@active.ident))
# 
# ## ------------------------------------------------------------------------
# SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single,top.n=25,
#                     clusters = singler$meta.data$orig.ident)
# ## ------------------------------------------------------------------------
# 
# SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single,top.n=25,
#                     clusters = singler$singler[[2]]$SingleR.single.main$labels)
## ------------------------------------------------------------------------
# out = SingleR.PlotTsne(singler$singler[[2]]$SingleR.single,
#                        singler$meta.data$xy,do.label=FALSE,
#                        do.letters =T,labels=singler$singler[[2]]$SingleR.single$labels, 
#                        dot.size = 5.3, font.size = 7)
# out$p

## ------------------------------------------------------------------------

# genes.use = c("GNLY", "IGLC2", "IGLC3", "FCGR3A", "IGKC", "CDKN1C", "GZMB", "ITM2C", "JCHAIN", "FCER1A")
# 
# df = data.frame(x=singler$meta.data$xy[,1],
#                 y=singler$meta.data$xy[,2],
#                 t(as.matrix(singler$seurat[["RNA"]]@data[genes.use,])))
# df = melt(df,id.vars = c('x','y'))
# ggplot(df,aes(x=x,y=y,color=value)) + 
#   geom_point(size=0.9)+scale_color_gradient(low="dark gray", high=" dark blue") + 
#   facet_wrap(~variable,ncol=4) +theme_classic()+xlab('')+ylab('')

## ------------------------------------------------------------------------

library(knitr)
singler2.log$seurat = s_test_nmads2_genes.markers
singler2.log$meta.data$orig.ident = s_test_nmads2_genes.markers@meta.data$orig.ident
singler2.log$meta.data$xy = s_test_nmads2_genes.markers@reductions$tsne@cell.embeddings #tsne coordinates
singler2.log$meta.data$clusters =s_test_nmads2_genes.markers@active.ident

singler2.sc$seurat = test_nmads2_genes_sctransform
singler2.sc$meta.data$orig.ident = test_nmads2_genes_sctransform@meta.data$orig.ident
singler2.sc$meta.data$xy = test_nmads2_genes_sctransform@reductions$tsne@cell.embeddings #tsne coordinates
singler2.sc$meta.data$clusters = test_nmads2_genes_sctransform@active.ident

my_list <- list(kable(table(singler2.log$meta.data$orig.ident,singler2.log$seurat@active.ident), "latex"),
kable(table(singler2.log$singler[[2]]$SingleR.single.main$labels,singler2.log$seurat@active.ident), "latex"))
write2pdf(my_list, paste0(res.dir, "singler2.log.pdf"), quiet = TRUE)

my_list<- list(kable(table(singler2.sc$meta.data$orig.ident,singler2.sc$seurat@active.ident), "latex"),
               kable(table(singler2.sc$singler[[2]]$SingleR.single.main$labels,singler2.sc$seurat@active.ident), "latex"))
write2pdf(my_list, paste0(res.dir, "singler2.sc.pdf"), quiet = TRUE)

singler3.log$seurat = s_test_nmads3_genes.markers
singler3.log$meta.data$orig.ident = s_test_nmads3_genes.markers@meta.data$orig.ident
singler3.log$meta.data$xy = s_test_nmads3_genes.markers@reductions$tsne@cell.embeddings #tsne coordinates
singler3.log$meta.data$clusters =s_test_nmads3_genes.markers@active.ident

singler3.sc$seurat = test_nmads3_genes_sctransform
singler3.sc$meta.data$orig.ident = test_nmads3_genes_sctransform@meta.data$orig.ident
singler3.sc$meta.data$xy = test_nmads3_genes_sctransform@reductions$tsne@cell.embeddings #tsne coordinates
singler3.sc$meta.data$clusters = test_nmads3_genes_sctransform@active.ident

my_list <- list(kable(table(singler3.log$meta.data$orig.ident,singler3.log$seurat@active.ident), "latex"),
                kable(table(singler3.log$singler[[2]]$SingleR.single.main$labels,singler3.log$seurat@active.ident), "latex"))
write2pdf(my_list, paste0(res.dir, "singler3.log.pdf"), quiet = TRUE)

my_list<- list(kable(table(singler3.sc$meta.data$orig.ident,singler3.sc$seurat@active.ident), "latex"),
               kable(table(singler3.sc$singler[[2]]$SingleR.single.main$labels,singler3.sc$seurat@active.ident), "latex"))
write2pdf(my_list, paste0(res.dir, "singler3.sc.pdf"), quiet = TRUE)

singler4.log$seurat = s_test_nmads4_genes.markers
singler4.log$meta.data$orig.ident = s_test_nmads4_genes.markers@meta.data$orig.ident
singler4.log$meta.data$xy = s_test_nmads4_genes.markers@reductions$tsne@cell.embeddings #tsne coordinates
singler4.log$meta.data$clusters =s_test_nmads4_genes.markers@active.ident

singler4.sc$seurat = test_nmads4_genes_sctransform
singler4.sc$meta.data$orig.ident = test_nmads4_genes_sctransform@meta.data$orig.ident
singler4.sc$meta.data$xy = test_nmads4_genes_sctransform@reductions$tsne@cell.embeddings #tsne coordinates
singler4.sc$meta.data$clusters = test_nmads4_genes_sctransform@active.ident

my_list <- list(kable(table(singler4.log$meta.data$orig.ident,singler4.log$seurat@active.ident), "latex"),
                kable(table(singler4.log$singler[[2]]$SingleR.single.main$labels,singler4.log$seurat@active.ident), "latex"))
write2pdf(my_list, paste0(res.dir, "singler4.log.pdf"), quiet = TRUE)

my_list<- list(kable(table(singler4.sc$meta.data$orig.ident,singler4.sc$seurat@active.ident), "latex"),
               kable(table(singler4.sc$singler[[2]]$SingleR.single.main$labels,singler4.sc$seurat@active.ident), "latex"))
write2pdf(my_list, paste0(res.dir, "singler4.sc.pdf"), quiet = TRUE)
## ------------------------------------------------------------------------
# singler_2labels_test <- s_test_nmads3_genes.markers
# head(s_test_nmads3_genes.markers$orig.ident)
# Idents(singler_2labels_test)<- singler_2$singler[[2]]$SingleR.single.main$labels
# head(Idents(singler_2labels_test))
# DimPlot(object = singler_2labels_test, reduction = 'umap', pt.size = 1, cols = rainbow(7), label = TRUE, label.size = 4)
# ggsave("test2.pdf",width = 15,height=10, dpi = 320)
## ------------------------------------------------------------------------
# out = SingleR.PlotTsne(singler_2$singler[[1]]$SingleR.single,
#                        singler_2$meta.data$xy,do.label = F,
#                        do.letters = T,labels = singler_2$meta.data$orig.ident, 
#                        dot.size = 4.3,alpha=0.5,label.size = 6)
# out$p
# ## ------------------------------------------------------------------------
# 
# out = SingleR.PlotTsne(singler_2$singler[[1]]$SingleR.single,
#                        singler_2$meta.data$xy,do.label = T,
#                        do.letters = F,labels=singler_2$singler[[2]]$SingleR.single.main$labels, 
#                        dot.size = 1.3,label.size = 5,alpha=0.5, font.size = 4)
# out$p
# ## ------------------------------------------------------------------------
# out = SingleR.PlotTsne(singler_2$singler[[1]]$SingleR.single,
#                        singler_2$meta.data$xy,do.label = T,
#                        do.letters = F,labels=singler_2$seurat@active.ident, 
#                        dot.size = 1.3,label.size = 5,alpha=0.5, font.size = 6)
# out$p
## ------------------------------------------------------------------------
# 
# SingleR.DrawHeatmap(singler_2$singler[[2]]$SingleR.single,top.n=25,
#                     clusters = singler_2$meta.data$orig.ident)
# ## ------------------------------------------------------------------------
# 
# SingleR.DrawHeatmap(singler_2$singler[[2]]$SingleR.single,top.n=25,
#                     clusters = singler_2$singler[[2]]$SingleR.single.main$labels)
# ## ------------------------------------------------------------------------
# out = SingleR.PlotTsne(singler_2$singler[[2]]$SingleR.single,
#                        singler_2$meta.data$xy,do.label=FALSE,
#                        do.letters =T,labels=singler_2$singler[[2]]$SingleR.single$labels, 
#                        dot.size = 5.3, font.size = 8)
# out$p
## ------------------------------------------------------------------------
# library(reshape2)
# library(data.table)
# genes.use = c("GNLY", "IGLC2", "IGLC3", "FCGR3A", "IGKC", "CDKN1C", "GZMB", "ITM2C", "JCHAIN", "FCER1A")
# 
# df = data.frame(x=singler_2$meta.data$xy[,1],
#                 y=singler_2$meta.data$xy[,2],
#                 t(as.matrix(singler_2$seurat[["RNA"]]@data[genes.use,])))
# df = melt(df,id.vars = c('x','y'))
# ggplot(df,aes(x=x,y=y,color=value)) + 
#   geom_point(size=0.9)+scale_color_gradient(low="dark gray", high=" dark blue") + 
#   facet_wrap(~variable,ncol=4) +theme_classic()+xlab('')+ylab('')
#   #theme(strip.background = element_blank())
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
# singler_3labels_test <- s_test_nmads2_genes.markers
# head(s_test_nmads3_genes.markers$orig.ident)
# Idents(singler_3labels_test)<- singler_3$singler[[2]]$SingleR.single.main$labels
# head(Idents(singler_3labels_test))
# DimPlot(object = singler_3labels_test, reduction = 'umap', pt.size = 1, cols = rainbow(7), label = TRUE, label.size = 4)
# ggsave("test3.pdf",width = 15,height=10, dpi = 320)
## ------------------------------------------------------------------------
# out = SingleR.PlotTsne(singler_3$singler[[1]]$SingleR.single,
#                        singler_3$meta.data$xy,do.label = F,
#                        do.letters = T,labels = singler_3$meta.data$orig.ident, 
#                        dot.size = 4.3,alpha=0.5,label.size = 6)
# out$p
# ## ------------------------------------------------------------------------
# out = SingleR.PlotTsne(singler_3$singler[[1]]$SingleR.single,
#                        singler_3$meta.data$xy,do.label = T,
#                        do.letters = F,labels=singler_3$singler[[2]]$SingleR.single.main$labels, 
#                        dot.size = 1.3,label.size = 5,alpha=0.5, font.size = 4)
# out$p
# ## ------------------------------------------------------------------------
# out = SingleR.PlotTsne(singler_3$singler[[1]]$SingleR.single,
#                        singler_3$meta.data$xy,do.label = T,
#                        do.letters = F,labels=singler_3$seurat@active.ident, 
#                        dot.size = 1.3,label.size = 5,alpha=0.5, font.size = 6)
# out$p
## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
# SingleR.DrawHeatmap(singler_3$singler[[2]]$SingleR.single,top.n=25,
#                     clusters = singler_3$meta.data$orig.ident)
# ## ------------------------------------------------------------------------
# 
# SingleR.DrawHeatmap(singler_3$singler[[2]]$SingleR.single,top.n=25,
#                     clusters = singler_3$singler[[2]]$SingleR.single.main$labels)
# ## ------------------------------------------------------------------------
# out = SingleR.PlotTsne(singler_3$singler[[2]]$SingleR.single,
#                        singler_3$meta.data$xy,do.label=FALSE,
#                        do.letters =T,labels=singler_3$singler[[2]]$SingleR.single$labels, 
#                        dot.size = 5.3, font.size = 8)
# out$p
# ## ------------------------------------------------------------------------
# library(reshape2)
# library(data.table)
# genes.use = c("GNLY", "IGLC2", "IGLC3", "FCGR3A", "IGKC", "CDKN1C", "GZMB", "ITM2C", "JCHAIN", "FCER1A")
# 
# df = data.frame(x=singler_3$meta.data$xy[,1],
#                 y=singler_3$meta.data$xy[,2],
#                 t(as.matrix(singler_3$seurat[["RNA"]]@data[genes.use,])))
# df = melt(df,id.vars = c('x','y'))
# ggplot(df,aes(x=x,y=y,color=value)) + 
#   geom_point(size=0.9)+scale_color_gradient(low="dark gray", high=" dark blue") + 
#   facet_wrap(~variable,ncol=4) +theme_classic()+xlab('')+ylab('')
# # +
# #   theme(strip.background = element_blank())

