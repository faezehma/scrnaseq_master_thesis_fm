library(Seurat)
library(SingleCellExperiment)

merge.methods <- readRDS("merge.methods.neuron_5k_v3.rds")
table(Idents(merge.methods))
merge.methods$log10GenesPerUMI <- log10(merge.methods$nFeature_RNA) / log10(merge.methods$nCount_RNA)
merge.methods$UMIperGenes <- merge.methods$nCount_RNA/merge.methods$nFeature_RNA
merge.methods$log10nUMI <- log10(merge.methods$nCount_RNA)
merge.methods$log10nGene <- log10(merge.methods$nFeature_RNA)


merge.methods$mitoRatio <- PercentageFeatureSet(object = merge.methods, pattern = "^mt-")
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

###################### Univariate outlier detection by using manual threshold: mito threshold,lower and upper threshold for nUMI & nGene
#for percent.mito, the MADs cut-off can be used and for nUMI and nGene, SDs cut-off based on mixture gaussian distribution:
#for percent.mito less stringent: 4MADs and for nUMI and nGene: 2.5 SD for all datasets and for all background removal methods.

###PBMC dataset:
selected_cellranger <- WhichCells(merge.methods, expression = pct_mt < 21.309 & methods == "cellranger")
selected_emptydrops <- WhichCells(merge.methods, expression = pct_mt < 20.526 & methods == "emptydrops")
selected_seurat200 <- WhichCells(merge.methods, expression = pct_mt < 20.645 & methods == "seurat200")
merge.methods.filtmito <- subset(merge.methods, cells = c(selected_cellranger,selected_emptydrops, selected_seurat200))
table(Idents(merge.methods.filtmito))

high.selected_cellranger <- WhichCells(merge.methods.filtmito, expression = nGene > 10^3.692 & methods == "cellranger")
high.selected_emptydrops <- WhichCells(merge.methods.filtmito, expression = nGene > 10^3.691 & methods == "emptydrops")
high.selected_seurat200 <- WhichCells(merge.methods.filtmito, expression = nGene > 10^3.693 & methods == "seurat200")

high.selected_cellranger_ <- WhichCells(merge.methods.filtmito, expression = nUMI > 10^4.321 & methods == "cellranger")
high.selected_emptydrops_ <- WhichCells(merge.methods.filtmito, expression = nUMI > 10^4.327  & methods == "emptydrops")
high.selected_seurat200_ <- WhichCells(merge.methods.filtmito, expression = nUMI > 10^4.32  & methods == "seurat200")
data.filt <- subset(merge.methods.filtmito, cells=setdiff(WhichCells(merge.methods.filtmito),c(high.selected_cellranger,high.selected_emptydrops,
                                                                                                    high.selected_seurat200,high.selected_cellranger_,high.selected_emptydrops_
                                                                                                    ,high.selected_seurat200_)))
table(Idents(data.filt))


low.selected_cellranger <- WhichCells(merge.methods.filtmito, expression = nGene <10^2.938 & methods == "cellranger")
low.selected_emptydrops <- WhichCells(merge.methods.filtmito, expression = nGene < 10^2.947& methods == "emptydrops")
low.selected_seurat200 <- WhichCells(merge.methods.filtmito, expression = nGene < 10^2.939& methods == "seurat200")

low.selected_cellranger_ <- WhichCells(merge.methods.filtmito, expression = nUMI < 10^3.417& methods == "cellranger")
low.selected_emptydrops_ <- WhichCells(merge.methods.filtmito, expression = nUMI < 10^3.402 & methods == "emptydrops")
low.selected_seurat200_ <- WhichCells(merge.methods.filtmito, expression = nUMI < 10^3.413 & methods == "seurat200")

data.filt.pbmc <- subset(data.filt, cells=setdiff(WhichCells(data.filt),c(low.selected_cellranger,low.selected_emptydrops,
                                                                                               low.selected_seurat200,
                                                                     low.selected_cellranger_,high.selected_emptydrops_
                                                                                               ,low.selected_seurat200_)))
table(Idents(data.filt.pbmc))
saveRDS(data.filt.pbmc, file = "new_analysis/data.filt.pbmc.rds")

#neuron_1k dataset:
selected_cellranger <- WhichCells(merge.methods, expression = pct_mt < 20.436 & methods == "cellranger")
selected_emptydrops <- WhichCells(merge.methods, expression = pct_mt <  21.733 & methods == "emptydrops")
selected_seurat200 <- WhichCells(merge.methods, expression = pct_mt <  22.344 & methods == "seurat200")
merge.methods.filtmito <- subset(merge.methods, cells = c(selected_cellranger,selected_emptydrops, selected_seurat200))
table(Idents(merge.methods.filtmito))

high.selected_cellranger <- WhichCells(merge.methods.filtmito, expression = nGene > 10^3.906 & methods == "cellranger")
high.selected_emptydrops <- WhichCells(merge.methods.filtmito, expression = nGene > 10^3.904 & methods == "emptydrops")
high.selected_seurat200 <- WhichCells(merge.methods.filtmito, expression = nGene > 10^3.909 & methods == "seurat200")

high.selected_cellranger_ <- WhichCells(merge.methods.filtmito, expression = nUMI > 10^4.697 & methods == "cellranger")
high.selected_emptydrops_ <- WhichCells(merge.methods.filtmito, expression = nUMI > 10^4.687  & methods == "emptydrops")
high.selected_seurat200_ <- WhichCells(merge.methods.filtmito, expression = nUMI > 10^4.705  & methods == "seurat200")
data.filt <- subset(merge.methods.filtmito, cells=setdiff(WhichCells(merge.methods.filtmito),c(high.selected_cellranger,high.selected_emptydrops,
                                                                                               high.selected_seurat200,high.selected_cellranger_,high.selected_emptydrops_
                                                                                               ,high.selected_seurat200_)))
table(Idents(data.filt))


low.selected_cellranger <- WhichCells(merge.methods.filtmito, expression = nGene < 10^3.21 & methods == "cellranger")
low.selected_emptydrops <- WhichCells(merge.methods.filtmito, expression = nGene < 10^3.22 & methods == "emptydrops")
low.selected_seurat200 <- WhichCells(merge.methods.filtmito, expression = nGene < 10^3.205 & methods == "seurat200")

low.selected_cellranger_ <- WhichCells(merge.methods.filtmito, expression = nUMI < 10^3.38 & methods == "cellranger")
low.selected_emptydrops_ <- WhichCells(merge.methods.filtmito, expression = nUMI <  10^3.421 & methods == "emptydrops")
low.selected_seurat200_ <- WhichCells(merge.methods.filtmito, expression = nUMI < 10^3.367 & methods == "seurat200")

data.filt.neuron1k <- subset(data.filt, cells=setdiff(WhichCells(data.filt),c(low.selected_cellranger,low.selected_emptydrops,
                                                                          low.selected_seurat200,
                                                                          low.selected_cellranger_,high.selected_emptydrops_
                                                                          ,low.selected_seurat200_)))
table(Idents(data.filt.neuron1k))
saveRDS(data.filt.neuron1k, file = "new_analysis/data.filt.neuron1k.rds")

#neuron5k dataset:
selected_cellranger <- WhichCells(merge.methods, expression = pct_mt < 13.99 & methods == "cellranger")
selected_emptydrops <- WhichCells(merge.methods, expression = pct_mt <  13.497 & methods == "emptydrops")
selected_seurat200 <- WhichCells(merge.methods, expression = pct_mt <  17.728 & methods == "seurat200")
merge.methods.filtmito <- subset(merge.methods, cells = c(selected_cellranger,selected_emptydrops, selected_seurat200))
table(Idents(merge.methods.filtmito))

high.selected_cellranger <- WhichCells(merge.methods.filtmito, expression = nGene > 10^3.836 & methods == "cellranger")
high.selected_emptydrops <- WhichCells(merge.methods.filtmito, expression = nGene >  10^3.838 & methods == "emptydrops")
high.selected_seurat200 <- WhichCells(merge.methods.filtmito, expression = nGene >10^3.858 & methods == "seurat200")

high.selected_cellranger_ <- WhichCells(merge.methods.filtmito, expression = nUMI >  10^4.592& methods == "cellranger")
high.selected_emptydrops_ <- WhichCells(merge.methods.filtmito, expression = nUMI >  10^4.547 & methods == "emptydrops")
high.selected_seurat200_ <- WhichCells(merge.methods.filtmito, expression = nUMI > 10^4.615  & methods == "seurat200")
data.filt <- subset(merge.methods.filtmito, cells=setdiff(WhichCells(merge.methods.filtmito),c(high.selected_cellranger,high.selected_emptydrops,
                                                                                               high.selected_seurat200,high.selected_cellranger_,high.selected_emptydrops_
                                                                                               ,high.selected_seurat200_)))
table(Idents(data.filt))


low.selected_cellranger <- WhichCells(merge.methods.filtmito, expression = nGene < 10^3.262 & methods == "cellranger")
low.selected_emptydrops <- WhichCells(merge.methods.filtmito, expression = nGene < 10^3.254 & methods == "emptydrops")
low.selected_seurat200 <- WhichCells(merge.methods.filtmito, expression = nGene < 10^3.212 & methods == "seurat200")

low.selected_cellranger_ <- WhichCells(merge.methods.filtmito, expression = nUMI < 10^3.408 & methods == "cellranger")
low.selected_emptydrops_ <- WhichCells(merge.methods.filtmito, expression = nUMI <  10^3.498 & methods == "emptydrops")
low.selected_seurat200_ <- WhichCells(merge.methods.filtmito, expression = nUMI < 10^3.359 & methods == "seurat200")

data.filt.neuron5k <- subset(data.filt, cells=setdiff(WhichCells(data.filt),c(low.selected_cellranger,low.selected_emptydrops,
                                                                              low.selected_seurat200,
                                                                              low.selected_cellranger_,high.selected_emptydrops_
                                                                              ,low.selected_seurat200_)))
table(Idents(data.filt.neuron5k))
saveRDS(data.filt.neuron5k, file = "new_analysis/data.filt.neuron5k.rds")

#CISE dataset:
selected_cellranger <- WhichCells(merge.methods, expression = pct_mt < 8.425 & methods == "cellranger")
selected_emptydrops <- WhichCells(merge.methods, expression = pct_mt < 10.614  & methods == "emptydrops")
selected_seurat200 <- WhichCells(merge.methods, expression = pct_mt <  9.801 & methods == "seurat200")
merge.methods.filtmito <- subset(merge.methods, cells = c(selected_cellranger,selected_emptydrops, selected_seurat200))
table(Idents(merge.methods.filtmito))

high.selected_cellranger <- WhichCells(merge.methods.filtmito, expression = nGene > 10^3.258 & methods == "cellranger")
high.selected_emptydrops <- WhichCells(merge.methods.filtmito, expression = nGene > 10^3.247  & methods == "emptydrops")
high.selected_seurat200 <- WhichCells(merge.methods.filtmito, expression = nGene > 10^3.275 & methods == "seurat200")

high.selected_cellranger_ <- WhichCells(merge.methods.filtmito, expression = nUMI > 10^3.729 & methods == "cellranger")
high.selected_emptydrops_ <- WhichCells(merge.methods.filtmito, expression = nUMI > 10^3.756 & methods == "emptydrops")
high.selected_seurat200_ <- WhichCells(merge.methods.filtmito, expression = nUMI > 10^3.763 & methods == "seurat200")
data.filt <- subset(merge.methods.filtmito, cells=setdiff(WhichCells(merge.methods.filtmito),c(high.selected_cellranger,high.selected_emptydrops,
                                                                                               high.selected_seurat200,high.selected_cellranger_,high.selected_emptydrops_
                                                                                               ,high.selected_seurat200_)))
table(Idents(data.filt))


low.selected_cellranger <- WhichCells(merge.methods.filtmito, expression = nGene < 10^2.422 & methods == "cellranger")
low.selected_emptydrops <- WhichCells(merge.methods.filtmito, expression = nGene <  10^2.497 & methods == "emptydrops")
low.selected_seurat200 <- WhichCells(merge.methods.filtmito, expression = nGene < 10^2.374 & methods == "seurat200")

low.selected_cellranger_ <- WhichCells(merge.methods.filtmito, expression = nUMI < 10^2.763 & methods == "cellranger")
low.selected_emptydrops_ <- WhichCells(merge.methods.filtmito, expression = nUMI < 10^2.597 & methods == "emptydrops")
low.selected_seurat200_ <- WhichCells(merge.methods.filtmito, expression = nUMI < 10^2.563 & methods == "seurat200")

data.filt.CISE <- subset(data.filt, cells=setdiff(WhichCells(data.filt),c(low.selected_cellranger,low.selected_emptydrops,
                                                                              low.selected_seurat200,
                                                                              low.selected_cellranger_,high.selected_emptydrops_
                                                                              ,low.selected_seurat200_)))
table(Idents(data.filt.CISE))
saveRDS(data.filt.CISE, file = "new_analysis/data.filt.CISE.rds")

###################### Multivariate outlier detection bu using runPCA from scater package to detect outliers and filter those outliers
#apply on 3 background removel method separately and then combined all togethor to build one combined object: seurat or sce 

# sce <- as.SingleCellExperiment(merge.methods)
# s0 <- subset(sce, ,ident=="cellranger")
# s1 <- subset(sce, ,ident=="emptydrops")
# s2 <- subset(sce, ,ident=="seurat200")
# my_s <- list(s0, s1, s2)
# result_s <- list()
# for (i in my_s){
#   patterns <- c("^MT-", "^mt-")
#   is.mito<- grepl(paste(patterns, collapse="|"), rownames(i))
#   out<- calculateQCMetrics(i,feature_controls=list(Mt=is.mito)) 
#   out <- runPCA(out, use_coldata = TRUE, detect_outliers = TRUE)
#                 
#   out <- scater::filter(out, outlier==FALSE)
#   result_s[[length(result_s)+1]] = out
# }
# 
# result_s2 <- list()
# for (i in result_s){
#   out <-as(counts(i), "dgCMatrix")
#   out <- CreateSeuratObject(out)
#   result_s2[[length(result_s2)+1]] = out
# }
# 
# combined <-  merge(result_s2[[1]], y = c(result_s2[[2]], result_s2[[3]]))
# table(Idents(combined))
# 
# #sce_final <- as.SingleCellExperiment(combined)
# #summary(sce_final$ident)
# saveRDS(combined, file = "new_analysis/data.filtPCA.CISE.rds")

###################### Multivariate outlier detection bu using runPCA from scater package to detect outliers and filter those outliers
#we don't subset the main object, only apply on main object that included all background removal methods:
methods_withintersect <- readRDS("new_analysis/methods_withintersect.neuron5kv3.rds")
sce <- as.SingleCellExperiment(methods_withintersect)
patterns <- c("^MT-", "^mt-")
is.mito<- grepl(paste(patterns, collapse="|"), rownames(sce))
sce <- calculateQCMetrics(sce,feature_controls=list(Mt=is.mito))
out <- runPCA(sce, use_coldata = TRUE, detect_outliers = TRUE, 
              selected_variables = c("pct_counts_in_top_100_features", "total_features_by_counts", "pct_counts_feature_control",
                                     "log10_total_counts_endogenous", "log10_total_counts_feature_control"))
out <- scater::filter(out, outlier==FALSE)
out_seurat <-as(counts(out), "dgCMatrix")
out_seurat <- CreateSeuratObject(out_seurat)
table(Idents(out_seurat))
saveRDS(out_seurat, file = "new_analysis/data_intersect.filtPCA.neuron5k.rds")

#visualization after outlier detection by pca
sce <- out
#sce <- as.SingleCellExperiment(data.filt.neuron5k)
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
ggsave(("new_analysis/plotUMAP_3backgroundmethods_neuron5k.afterpcaoutlier.pdf"), fig, width = 7,height=6, scale = 1.5)


p00 <- plotPCA(subset0,colour_by = "ident")+ 
  scale_fill_manual(values = c("grey", "red"))

p11 <- plotPCA(subset1,colour_by = "ident")+ 
  scale_fill_manual(values = c("grey", "green"))

p22 <- plotPCA(subset2,colour_by = "ident")+ 
  scale_fill_manual(values = c("grey", "blue"))

fig <- plot_grid(p00, p11, p22,
                 ncol = 2, labels = "AUTO")


ggsave(("new_analysis/plotPCA_3backgroundmethods_neuron5k.afterpcaoutlier.pdf"), fig, width = 7,height=6, scale = 1.5)

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

ggsave(("new_analysis/plottsne_3backgroundmethods_neuron5k.afterpcaoutlier.pdf"), fig, width = 7,height=6, scale = 1.5)

u1 <- plotUMAP(sce, colour_by = "log10_total_counts")+ ggtitle("total counts")
u2 <- plotUMAP(sce, colour_by = "log10_total_features_by_counts")+ ggtitle("total features")
u3 <- plotUMAP(sce, colour_by = "pct_counts_Mt")+ ggtitle("% mitochondrial genes")

t1 <- plotTSNE(sce, colour_by = "log10_total_counts")+ ggtitle("total counts")
t2 <- plotTSNE(sce, colour_by = "log10_total_features_by_counts")+ ggtitle("total features")
t3 <- plotTSNE(sce, colour_by = "pct_counts_Mt")+ ggtitle("% mitochondrial genes")

c1 <- plotPCA(sce, colour_by = "log10_total_counts")+ ggtitle("total counts")
c2 <- plotPCA(sce, colour_by = "log10_total_features_by_counts")+ ggtitle("total features")
c3 <- plotPCA(sce, colour_by = "pct_counts_Mt")+ ggtitle("% mitochondrial genes")

fig <- plot_grid(p000, p111, p222, u1, u2, u3, ncol = 3)
title <- ggdraw() + draw_label("After PCA based cell filtering", fontface='bold')
fig <- plot_grid(title, fig, ncol=1, rel_heights=c(0.1, 1))
ggsave(("new_analysis/plotumap_allQC_neuron5k.afterpcaoutlier.pdf"), fig, width = 7,height=6, scale = 1.5)

fig <- plot_grid(p00, p11, p22, c1, c2, c3, ncol = 3)
title <- ggdraw() + draw_label("After PCA based cell filtering", fontface='bold')
fig <- plot_grid(title, fig, ncol=1, rel_heights=c(0.1, 1))
ggsave(("new_analysis/plotpca_allQC_neuron5k.afterpcaoutlier.pdf"), fig, width = 7,height=6, scale = 1.5)

fig <- plot_grid(p0, p1, p2, t1, t2, t3, ncol = 3)
title <- ggdraw() + draw_label("After PCA based cell filtering", fontface='bold')
fig <- plot_grid(title, fig, ncol=1, rel_heights=c(0.1, 1))
ggsave(("new_analysis/plottsne_allQC_neuron5k.afterpcaoutlier.pdf"), fig, width = 7,height=6, scale = 1.5)


###############visualization manual threshold (univariate outlier) before and after filtering
#***before*** applying 2.5 SD
merge.methods <- readRDS("new_analysis/merge.methods.CISE13.rds")
merge.methods$log10GenesPerUMI <- log10(merge.methods$nFeature_RNA) / log10(merge.methods$nCount_RNA)
merge.methods$UMIperGenes <- merge.methods$nCount_RNA/merge.methods$nFeature_RNA
merge.methods$log10nUMI <- log10(merge.methods$nCount_RNA)
merge.methods$log10nGene <- log10(merge.methods$nFeature_RNA)
merge.methods$mitoRatio <- PercentageFeatureSet(object = merge.methods, pattern = "^mt-")
merge.methods$mitoRatio <- merge.methods@meta.data$mitoRatio / 100
merge.methods$pct_mt <- merge.methods@meta.data$mitoRatio * 100

sce <- as.SingleCellExperiment(merge.methods)
set.seed(1)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters=clusters)
sce <- normalize(sce)

sce <- runPCA(sce)
sce <- runTSNE(sce, dimred="PCA")
sce <- runUMAP(sce)

pdf(file='new_analysis/QC.beforeSD_CISE.pdf', width = 15, height = 10)
u1 <- plotUMAP(sce, colour_by = "log10nUMI")+ ggtitle("total counts")
u2 <- plotUMAP(sce, colour_by = "log10nGene")+ ggtitle("total features")
u3 <- plotUMAP(sce, colour_by = "pct_mt")+ ggtitle("% mitochondrial genes")
fig <- plot_grid(u1,u2,u3, ncol = 2)
title <- ggdraw() + draw_label("Before manual cell filtering", fontface='bold')
plot_grid(title, fig, ncol=1, rel_heights=c(0.1, 1))

t1 <- plotTSNE(sce, colour_by = "log10nUMI")+ ggtitle("total counts")
t2 <- plotTSNE(sce, colour_by = "log10nGene")+ ggtitle("total features")
t3 <- plotTSNE(sce, colour_by = "pct_mt")+ ggtitle("% mitochondrial genes")
fig <- plot_grid(t1,t2,t3,ncol = 2)
title <- ggdraw() + draw_label("Before manual cell filtering", fontface='bold')
plot_grid(title, fig, ncol=1, rel_heights=c(0.1, 1))

c1 <- plotPCA(sce, colour_by = "log10nUMI")+ ggtitle("total counts")
c2 <- plotPCA(sce, colour_by = "log10nGene")+ ggtitle("total features")
c3 <- plotPCA(sce, colour_by = "pct_mt")+ ggtitle("% mitochondrial genes")
fig <- plot_grid(c1,c2,c3,ncol = 2)
title <- ggdraw() + draw_label("Before manual cell filtering", fontface='bold')
plot_grid(title, fig, ncol=1, rel_heights=c(0.1, 1))
dev.off()

#***after*** applying 2.5 SD
data.filt <- readRDS("new_analysis/data.filt.CISE.rds")
sce <- as.SingleCellExperiment(data.filt)
set.seed(1)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters=clusters)
sce <- normalize(sce)

sce <- runPCA(sce)
sce <- runTSNE(sce, dimred="PCA")
sce <- runUMAP(sce)

# colData(sce)$methods <- "All three"
# colData(sce)$methods[sce$ident == "seurat200"] <- "seurat200"
# colData(sce)$methods[sce$ident == "emptydrops"] <- "emptydrops"
# colData(sce)$methods[sce$ident == "cellranger"] <- "cellranger"
# 
# methods <- rep("grey", nrow(reducedDim(sce, "TSNE")))
# methods[sce$ident=="emptydrops"] <- "green"
# methods[sce$ident=="cellranger"] <- "red"
# methods[sce$ident=="seurat200 "] <- "blue"

#subset0 <- subset(sce, ,methods==c("cellranger", "All three"))
p000 <- plotUMAP(sce, colour_by="ident")+ 
  scale_fill_manual(values = c("red", "green", "blue"))+
  xlab("UMAP1")+ ylab("UMAP2")+
  ggtitle("")



p00 <- plotPCA(sce, colour_by="ident")+ 
  scale_fill_manual(values = c("red", "green", "blue"))



p0 <- plotTSNE(sce, colour_by="ident")+ 
  scale_fill_manual(values =  c("red", "green", "blue"))+
  xlab("t-SNE1")+ ylab("t-SNE2")+
  ggtitle("")



fig <- plot_grid(p000, p00, p0,
                 ncol = 2, labels = "AUTO")

pdf(file='new_analysis/QC.afterSD_CISE.pdf', width = 15, height = 10)
u1 <- plotUMAP(sce, colour_by = "log10nUMI")+ ggtitle("total counts")
u2 <- plotUMAP(sce, colour_by = "log10nGene")+ ggtitle("total features")
u3 <- plotUMAP(sce, colour_by = "pct_mt")+ ggtitle("% mitochondrial genes")
fig <- plot_grid(p000, u1, u2, u3, ncol = 2)
title <- ggdraw() + draw_label("After manual cell filtering", fontface='bold')
plot_grid(title, fig, ncol=1, rel_heights=c(0.1, 1))

t1 <- plotTSNE(sce, colour_by = "log10nUMI")+ ggtitle("total counts")
t2 <- plotTSNE(sce, colour_by = "log10nGene")+ ggtitle("total features")
t3 <- plotTSNE(sce, colour_by = "pct_mt")+ ggtitle("% mitochondrial genes")
fig <- plot_grid(p0, t1, t2, t3, ncol = 2)
title <- ggdraw() + draw_label("After manual cell filtering", fontface='bold')
plot_grid(title, fig, ncol=1, rel_heights=c(0.1, 1))

c1 <- plotPCA(sce, colour_by = "log10nUMI")+ ggtitle("total counts")
c2 <- plotPCA(sce, colour_by = "log10nGene")+ ggtitle("total features")
c3 <- plotPCA(sce, colour_by = "pct_mt")+ ggtitle("% mitochondrial genes")
fig <- plot_grid(p00, c1, c2, c3, ncol = 2)
title <- ggdraw() + draw_label("After manual cell filtering", fontface='bold')
plot_grid(title, fig, ncol=1, rel_heights=c(0.1, 1))
dev.off()

########testing outlier detection 
library(scater)
cellranger <- subset(merge.methods, idents = "cellranger")
x <- SingleCellExperiment(list(counts=cellranger[["RNA"]]@counts))

sc.df<- perCellQCMetrics(x, 
                    subsets=list(Mito=grep("MT-", rownames(x))))

colData(x) <- cbind(colData(x), sc.df)
x <- runColDataPCA(x, variables=list(
  "sum", "detected", "subsets_Mito_percent"
), outliers = TRUE)

filtered <- x[,!x$outlier == TRUE]

#saveRDS(filtered, "/volumes/SD Card/Thesis/filtered_neuron5k_pca.rds")

filtered <- CreateSeuratObject(as(counts(filtered),"dgCMatrix"))
filtered <- NormalizeData(filtered)
filtered <- ScaleData(filtered)
saveRDS(filtered, "/volumes/SD Card/Thesis/filtered_neuron5k_pca_scaled.rds")
       