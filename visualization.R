library("SingleCellExperiment")
library("scater")
library("scran")

# RNA-seq
library("edgeR")

# Plotting
library("cowplot")
library("gridExtra")
library("mixtools")

# Tidyverse
library("tidyverse")
merge.methods <- readRDS("/Volumes/SD Card/Thesis/new_analysis/merge.methods.CISE13.rds")
my_object <- merge.methods
my_object <- SingleCellExperiment(list(counts=my_object[["RNA"]]@counts))
patterns <- c("^MT-", "^mt-")
is.mito<- grepl(paste(patterns, collapse="|"), rownames(my_object))
my_object<- calculateQCMetrics(my_object,feature_controls=list(Mt=is.mito))
cell_data <- as.data.frame(colData(my_object))
########
cellranger <- subset(my_object, methods =="cellranger")
emptydrops <- subset(my_object, methods =="emptydrops")
seurat200 <- subset(my_object, methods =="seurat200")

outlierHistogram <- function(data, x, mads = 1:3, bins = 100,
                             show_zero = FALSE) {
  x_data <- data[, x]
  med    <- median(x_data)
  MAD    <- mad(x_data, center = med, na.rm = TRUE)
  
  outs <- lapply(mads, function(n) {
    lower  <- med - n * MAD
    higher <- med + n * MAD
    n_low  <- sum(x_data < lower)
    n_high <- sum(x_data > higher)

    list(n = n, lower = lower, higher = higher,
         n_low = n_low, n_high = n_high)
    
  })
    gg <- ggplot(data, aes(x = !!ensym(x))) +
      geom_histogram(bins = 100, color = "black", fill = "gray"
                     )+ #aes(y=..density.., fill=..count..) inside of geom_histogram
      geom_vline(xintercept = med, colour = "blue", linetype = "dashed") +
      annotate("text", x = med, y = 0, label = "Median", colour = "blue",
               angle = 90, hjust = -0.1, vjust = -0.5) +
      annotate("text", x = med, y = Inf, label = round(med, 3), colour = "blue",
               angle = 90, hjust = 1, vjust = -0.5) +
      theme_minimal() +
      ylab("frequency")+
      theme(legend.position = "none")+
      ggtitle("")
      # stat_function(fun=dnorm, color="red", args=list(mean=mean(x_data),
      #                                                     sd=sd(x_data)))
      
      
    cols <- scales::brewer_pal(palette = "Reds",
                               direction = -1)(length(mads) + 1)
    
    for (i in seq_along(outs)) {
      out <- outs[[i]]
      
      if (show_zero) {
        show_low  <- TRUE
        show_high <- TRUE
      } else {
        show_low  <- out$n_low > 0
        show_high <- out$n_high > 0
      }
      
      if (show_low) {
        gg <- gg +
          geom_vline(xintercept = out$lower, colour = cols[i], linetype = "dashed") +
          annotate("text", x = out$lower, y = 0,
                   label = paste0(out$n, " MADs"), colour = cols[i],
                   angle = 90, hjust = -0.1, vjust = -0.5) +
          annotate("text", x = out$lower, y = Inf,
                   label = paste0(round(out$lower, 3),
                                  " (", out$n_low, " lower)"),
                   colour = cols[i], angle = 90, hjust = 1, vjust = -0.5)
      }
      
      if (show_high) {
        gg <- gg +
          geom_vline(xintercept = out$higher, colour = cols[i], linetype = "dashed") +
          annotate("text", x = out$higher, y = 0,
                   label = paste0(out$n, " MADs"), colour = cols[i],
                   angle = 90, hjust = -0.1, vjust = -0.5) +
          annotate("text", x = out$higher, y = Inf,
                   label = paste0(round(out$higher, 3),
                                  " (", out$n_high, " higher)"),
                   colour = cols[i], angle = 90, hjust = 1, vjust = -0.5)
      }
    }
    
    return(gg)
  }
  
#quality control
#outlierHistogram(cellranger, "log10nUMI", mads = c(2, 3, 4))

#pdf(file='new_analysis/fit_outliers_pbmc5.pdf', width = 15, height = 10)
outlierHistogram(cellranger, "log10nUMI", mads = 1:4)
# p <- plot(fitd)
# dev.off()

##-------------------------------------------------------------------------------------------------------------------------------------------------------
#histogram for mixture gaussian outlier detection:
cellranger <- subset(metadata, methods =="cellranger")
emptydrops <- subset(metadata, methods =="emptydrops")
seurat200 <- subset(metadata, methods =="seurat200")

outlier_Histogram <- function(data, x, sigmas = 1:3, bins = 100,
                             show_zero = FALSE) {
  x_data <- data[, x]
  
  mixmdl = normalmixEM(x_data, k=2,maxit=10000,epsilon=1e-08, maxrestarts = 1000)
  mu <- mixmdl$mu[2]
  sigma <- mixmdl$sigma[2]
  
  outs <- lapply(sigmas, function(n) {
    lower  <- mu - n * sigma
    higher <- mu + n * sigma
    n_low  <- sum(x_data < lower)
    n_high <- sum(x_data > higher)
    
    list(n = n, lower = lower, higher = higher,
         n_low = n_low, n_high = n_high)
  })
  sdnorm =
    function(x, mean=0, sd=1, lambda=1){lambda*dnorm(x, mean=mean, sd=sd)}
  
  gg <- ggplot(data, aes(x = !!ensym(x))) +
    geom_histogram(bins = 100, color = "black", fill = "gray",  
                   aes(y=..density.., fill=..count..)) +
    geom_vline(xintercept = mu, colour = "blue", linetype = "dashed") +
    annotate("text", x = mu, y = 0, label = "Mean", colour = "blue",
             angle = 90, hjust = -0.1, vjust = -0.5) +
    annotate("text", x = mu, y = Inf, label = round(mu, 3), colour = "blue",
             angle = 90, hjust = 1, vjust = -0.5)+
    stat_function(fun=sdnorm,args=list(mean=mixmdl$mu[1],sd=mixmdl$sigma[1], lambda=mixmdl$lambda[1]), color ="black")+
    stat_function(fun=sdnorm,args=list(mean=mixmdl$mu[2],sd=mixmdl$sigma[2], lambda=mixmdl$lambda[2]),color =  "black")+
    theme_minimal() +
    ylab("density")+
    theme(legend.position = "none")+
    ggtitle("cellranger CISE13")
  
  
  for (i in seq_along(outs)) {
    out <- outs[[i]]
    
    if (show_zero) {
      show_low  <- TRUE
      show_high <- TRUE
    } else {
      show_low  <- out$n_low > 0
      show_high <- out$n_high > 0
    }
    
    if (show_low) {
      gg <- gg +
        geom_vline(xintercept = out$lower, colour = "red", linetype = "dashed") +
        annotate("text", x = out$lower, y = 0,
                 label = paste0(out$n, " SD"), colour = "red",
                 angle = 90, hjust = -0.1, vjust = -0.5) +
        annotate("text", x = out$lower, y = Inf,
                 label = paste0(round(out$lower, 3),
                                " (", out$n_low, " lower)"),
                 colour = "red", angle = 90, hjust = 1, vjust = -0.5)
    }
    
    
    if (show_high) {
      gg <- gg +
        geom_vline(xintercept = out$higher, colour = "red", linetype = "dashed") +
        annotate("text", x = out$higher, y = 0,
                 label = paste0(out$n, "SD"), colour = "red",
                 angle = 90, hjust = -0.1, vjust = -0.5) +
        annotate("text", x = out$higher, y = Inf,
                 label = paste0(round(out$higher, 3),
                                " (", out$n_high, " higher)"),
                 colour = "red", angle = 90, hjust = 1, vjust = -0.5)
    }
  }
  
  
  return(gg)
}

outlier_Histogram(cellranger, "log10nGene", sigmas = c(2,2.5, 3, 4))

##-------------------------------------------------------------------------------------------------------------------------------------------------------
#nmads for feature and count
threshold = median(log(combined.object.neuron5@meta.data$nFeature_RNA)) - 2*mad(log(combined.object.neuron5@meta.data$nFeature_RNA))
threshold1 = median(log(combined.object.neuron5@meta.data$nFeature_RNA)) - 3*mad(log(combined.object.neuron5@meta.data$nFeature_RNA))
threshold2 = median(log(combined.object.neuron5@meta.data$nFeature_RNA)) - 4*mad(log(combined.object.neuron5@meta.data$nFeature_RNA))
thresholda = median(log(combined.object.neuron5@meta.data$nCount_RNA)) - 2*mad(log(combined.object.neuron5@meta.data$nCount_RNA))
thresholdb = median(log(combined.object.neuron5@meta.data$nCount_RNA)) - 3*mad(log(combined.object.neuron5@meta.data$nCount_RNA))
thresholdc = median(log(combined.object.neuron5@meta.data$nCount_RNA)) - 4*mad(log(combined.object.neuron5@meta.data$nCount_RNA))
#nmads for percent.mt 
thresholdi = median(log(combined.object.neuron5@meta.data$percent.mt)) + 2*mad(log(combined.object.neuron5@meta.data$percent.mt))
thresholdii = median(log(combined.object.neuron5@meta.data$percent.mt)) + 3*mad(log(combined.object.neuron5@meta.data$percent.mt))
thresholdiii = median(log(combined.object.neuron5@meta.data$percent.mt)) + 4*mad(log(combined.object.neuron5@meta.data$percent.mt))



plot1 <- FeatureScatter(object = combined.object.neuron5, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident",
                        pt.size = 2, cols = c("#8DC63F", "#7A52C7", "#00ADEF", "#EC008C","#7A52C7", "#F47920", "#00ECA6", "#E04F4A", "#4A80E0", "#99FFFF"))+
  geom_point(alpha = 1/500)+
  geom_hline(yintercept = c(exp(thresholdi), exp(thresholdii), exp(thresholdiii)), linetype="dashed", color = "black")+
  annotate("text", x = Inf, y = c(exp(thresholdi), exp(thresholdii), exp(thresholdiii)), label = c("2mads", "3mads","4mads"), colour = "black",
           angle = 0, hjust = 1, vjust = -0.5)+
  theme_minimal() +
  theme(legend.position = "bottom")


plot2 <- FeatureScatter(object = combined.object.neuron5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", 
                        pt.size = 2, cols = c("#8DC63F", "#7A52C7", "#00ADEF", "#EC008C","#7A52C7", "#F47920", "#00ECA6", "#E04F4A", "#4A80E0", "#99FFFF"))+
  geom_point(alpha = 1/550)+
  geom_hline(yintercept = c(exp(threshold), exp(threshold1), exp(threshold2)), linetype="dashed", color = "black")+
  annotate("text", x = Inf, y = c(exp(threshold), exp(threshold1), exp(threshold2)), label = c("2mads", "3mads","4mads"), colour = "black",
           angle = 0, hjust = 1, vjust = -0.5)+ geom_vline(xintercept = c(exp(thresholda), exp(thresholdb), exp(thresholdc)), linetype="dashed", color = "black")+
  annotate("text", y = 6300, x = c(exp(thresholda), exp(thresholdb), exp(thresholdc)), label = c("2mads", "3mads","4mads"), colour = "black",
           angle = 90, hjust =  -0.1, vjust = -0.5, size = 3.5 )+
  theme_minimal() +
  theme(legend.position = "bottom")
#plot2thresholdboth

plot5 <- FeatureScatter(object = combined.object.neuron5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", 
                        pt.size = 2, cols = c("#8DC63F", "#7A52C7", "#00ADEF", "#EC008C","#7A52C7", "#F47920", "#00ECA6", "#E04F4A", "#4A80E0","#99FFFF" ))+
  geom_point(alpha = 0.1)+
geom_vline(xintercept = c(exp(thresholda), exp(thresholdb), exp(thresholdc)), linetype="dashed", color = "black")+
  annotate("text", y = 6300, x = c(exp(thresholda), exp(thresholdb), exp(thresholdc)), label = c("2mads", "3mads","4mads"), colour = "black",
           angle = 90, hjust =  -0.1, vjust = -0.5, size = 3.5 )+
  theme_minimal() +
  theme(legend.position = "bottom")




#########for plot
#plot the effect of only pca on filtering cells on combined object with lebeld cells based on venn diagram
test_ <- test
test_meta <- test
# colData(test_)$all_cell <- "FALSE"
# test_.pca <- runPCA(test_, use_coldata=TRUE, detect_outliers=TRUE)
# kept_pca <- !colData(test_.pca)$outlier
# colData(test_.pca)$KeptPCA <- kept_pca
# colData(test_)$Kept[kept_pca]    <- "PCA"
# plot_test_ <- test_ %>% 
#   mutate(info_pca = "removed") %>%
#   mutate(info_pca = if_else(kept_pca, "PCA", info_pca)) %>%
#   mutate(info_pca = factor(info_pca, levels = c("removed","PCA"),
#                           labels = c("removed","PCA")))
# 
# plot_test_ <- as.data.frame(colData(plot_test_))
# 
# test_ <- as(counts(test_), "dgCMatrix")
# test_ <- CreateSeuratObject(test_)
# 
# test_ <- as.data.frame(test_@meta.data)
# test_$nCount_RNA <- NULL
# test_$nFeature_RNA <- NULL
# test_$pca_info <- "removed"
# test_$pca_info <- plot_test_$info_pca

test_meta <- as(counts(test_meta), "dgCMatrix")
test_meta <- CreateSeuratObject(test_meta)
for (i in rownames(test_meta@meta.data)){
  if (!(i %in% colnames(test.pca))) {
    test_meta@meta.data[i,"orig.ident"] = test_meta@meta.data["filtered","orig.ident"]
  } else {
    next
  }
}

test_meta@meta.data$orig.ident <- as.character(test_meta@meta.data$orig.ident)
test_meta@meta.data$orig.ident[is.na(test_meta@meta.data$orig.ident)] <- "Removed PCA"

test_meta@meta.data$orig.ident
#only pca filtered out
plot3 <- FeatureScatter(object = test_meta, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", 
                        pt.size = 2, cols= c("#8DC63F","#00ADEF","grey", "#7A52C7", "#EC008C", "#4A80E0","#E3635C","#F47920","yellow"))+
  geom_point(alpha = 1/700)+
  theme_minimal() +
  theme(legend.position = "bottom")

#all cells of regions without cells filtering
plot4 <- FeatureScatter(object = combined.object.neuron5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", 
                        pt.size = 2, cols= c("#8DC63F", "#7A52C7", "#00ADEF", "#EC008C","#7A52C7", "#F47920", "#00ECA6", "#E04F4A"))+
  geom_point(alpha = 0.1)+
  theme_minimal() +
  theme(legend.position = "bottom")

my_fig <- plot_grid(plot4, plot3,
                        ncol = 1, labels = "AUTO")
ggsave(("output_PBMC_1k_v3/qc_mine-thresholds_mads.pdf"), my_fig, width = 20,height=15, scale = 1.5, limitsize = FALSE)

res.dir <- paste0("/Volumes/SD Card/Thesis/output_PBMC_1k_v3/")
dir.create(res.dir) 
ggsave(
  grid.arrange(plot4, plot5, plot2,
               plot3, plot1,
               ncol=2),
  file=paste0(res.dir,"qc_mine-thresholds_mads.pdf"), width = 30,height=20)

###########################################################################
#UMI and genes plot with mito ratio gradients for joint upper and lower threshold:
#neuron5k:
ll.cd = 10^3.408
ll.ng = 10^3.262
ul.cd = 10^4.592
ul.ng = 10^3.836

ll.cd2 = 10^3.498
ll.ng2 = 10^3.254
ul.cd2 = 10^4.547
ul.ng2 = 10^3.838

ll.cd3 = 10^3.359
ll.ng3 = 10^3.212
ul.cd3 = 10^4.615
ul.ng3 = 10^3.858

cellranger <- subset(metadata, methods =="cellranger")
emptydrops <- subset(metadata, methods =="emptydrops")
seurat200 <- subset(metadata, methods =="seurat200")
pdf(file='new_analysis/new_cutoff/SD_cutoff/all_joint_cutoff_neuron5k.pdf', width = 15, height = 10)
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
  ggtitle("neuron_5k_v3")+
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
  ggtitle("neuron_5k_v3")+
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
  ggtitle("neuron_5k_v3")+
  theme_bw()
dev.off()


#neuron1k:
ll.cd = 10^3.38
ll.ng = 10^3.21
ul.cd = 10^4.697
ul.ng = 10^3.906

ll.cd2 = 10^3.421
ll.ng2 = 10^3.22
ul.cd2 = 10^4.687
ul.ng2 = 10^3.904

ll.cd3 = 10^3.367
ll.ng3 = 10^3.205
ul.cd3 = 10^4.705
ul.ng3 = 10^3.909

cellranger <- subset(metadata, methods =="cellranger")
emptydrops <- subset(metadata, methods =="emptydrops")
seurat200 <- subset(metadata, methods =="seurat200")
pdf(file='new_analysis/new_cutoff/SD_cutoff/all_joint_cutoff_neuron1k.pdf', width = 15, height = 10)
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
  ggtitle("neuron_1k_v3")+
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
  ggtitle("neuron_1k_v3")+
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
  ggtitle("neuron_1k_v3")+
  theme_bw()
dev.off()


#PBMC_1k:
ll.cd = 10^3.417
ll.ng = 10^2.938
ul.cd = 10^4.321
ul.ng = 10^3.692

ll.cd2 = 10^3.402
ll.ng2 = 10^2.947
ul.cd2 = 10^4.327
ul.ng2 = 10^3.691

ll.cd3 = 10^3.413
ll.ng3 = 10^2.939
ul.cd3 = 10^4.32
ul.ng3 = 10^3.693

cellranger <- subset(metadata, methods =="cellranger")
emptydrops <- subset(metadata, methods =="emptydrops")
seurat200 <- subset(metadata, methods =="seurat200")
pdf(file='new_analysis/new_cutoff/SD_cutoff/all_joint_cutoff_pbmc1k.pdf', width = 15, height = 10)
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
  ggtitle("PBMC_1k_v3")+
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
  ggtitle("PBMC_1k_v3")+
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
  ggtitle("PBMC_1k_v3")+
  theme_bw()
dev.off()

###############test:
library("ggplot2")
library("cowplot")
sdnorm =
  function(x, mean=0, sd=1, lambda=1){lambda*dnorm(x, mean=mean, sd=sd)}
g1 = c(sdnorm(cellranger$log10nUMI, mean=mixmdl$mu[1], sd=mixmdl$sigma[1], lambda=mixmdl$lambda[1]), sdnorm(cellranger$log10nUMI, mean=mixmdl$mu[2], sd=mixmdl$sigma[2], lambda=mixmdl$lambda[2]))
g2 = c(sdnorm(cellranger$log10nGene, mean=mixmdl2$mu[1], sd=mixmdl2$sigma[1], lambda=mixmdl2$lambda[1]), sdnorm(cellranger$log10nGene, mean=mixmdl2$mu[2], sd=mixmdl2$sigma[2], lambda=mixmdl2$lambda[2]))
group = as.factor(rep(c(1,2), each=1222))
df_exp = data.frame(G1=log2(g1 + 1) , G2=log2(g2 + 1), GROUP=group)

ll.cd = 3.417
ul.cd = 4.321
ll.ng = 2.938
ul.ng = 3.692

gg <- ggplot(cellranger, aes(x = log10nUMI,y=log10nGene,colour=mitoRatio)) +
  geom_point(size=0.7) +
  scale_colour_gradientn(colours = c("darkblue","darkgreen","green","darkorange2")) +
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
  ggtitle("cellranger PBMC")

gg_dist_g1 = ggplot(df_exp, aes(G1, fill=group)) + geom_density(alpha=.5)
gg_dist_g2 = ggplot(df_exp, aes(G2, fill=group)) + geom_density(alpha=.5)

#####################
threshold = median(log(cellranger$nUMI)) - 2*mad(log(cellranger$nUMI))

g1 = ggplot(cellranger) 
g1 = g1 + geom_histogram(aes(x=nUMI),bins=200,fill = "grey")
g1 = g1 + geom_density(aes(x=nUMI)) 
g1 + geom_vline(xintercept = exp(threshold), color="red")


#############################################################
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

merge.methods <- merge(cellranger, y = c(emptydrops, s_raw_200), add.cell.ids = c("cellranger", "emptydrops", "seurat200"))

methods<- rep("cellranger",ncol(merge.methods))
methods[Idents(merge.methods) == "emptydrops"] <- "emptydrops"
methods[Idents(merge.methods) == "seurat200"] <- "seurat200"


merge.methods <- AddMetaData(merge.methods, methods, col.name = "methods")
table(Idents(merge.methods))
table(merge.methods$orig.ident)

#stats <- barcodeRanks(merge.methods@assays$RNA)

merge.methods$log10nUMI <- log10(merge.methods$nCount_RNA)
merge.methods@meta.data %>%
  ggplot(aes(color=methods, x=log10nUMI, fill= methods)) +
  geom_density(alpha = 0.0001) +
  scale_x_log10() +
  theme_classic() +
  ylab("frequency")

# it should be colnames from each object not merged objects!!!
library("UpSetR")
my_list <- list(
  "CellRanger3" = colnames(s_filtered_cellranger3),
  "EmptyDrops" = colnames(s_raw_emptydrops),
  "Seurat200" = colnames(s_raw_200))
upset(fromList(my_list),sets.bar.color = c("red","blue","green"),
      main.bar.color = "black")


venn.plot <- venn.diagram(
  x = list(
    "CellRanger3" = colnames(s_filtered_cellranger3),
    "EmptyDrops" = colnames(s_raw_emptydrops),
    "Seurat200" = colnames(s_raw_200)
),
  euler.d = TRUE,
  filename=NULL,
  fil= c("red","blue","green"),
  cex = 2,
  cat.cex = 0.7,
  reverse = TRUE
)
