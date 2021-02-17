library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)


lymphgland <- readRDS('lymphgland.Rds')
DimPlot(lymphgland, reduction = 'tsne')
lymphgland@meta.data$Subclustering <- NULL
summary(lymphgland@meta.data$new_subclustering)
head(lymphgland@meta.data)

targetObj <- readRDS('tmp/dropseq.mtcut.norm.Rds')

### Transfer ###
transfer.anchors <- FindTransferAnchors(reference = lymphgland, query = targetObj, dims = 1:30)
saveRDS(transfer.anchors, 'tmp/transfer.anchors.Rds')
predictions <- TransferData(anchorset = transfer.anchors, refdata = lymphgland@meta.data$new_subclustering, dims = 1:30)
saveRDS(predictions, 'tmp/predictions.Rds')
head(predictions); summary(predictions$prediction.score.max)

targetObj@meta.data$newsubclustering <- as.character(predictions$predicted.id)
targetObj@meta.data$predscore <- predictions$prediction.score.max
targetObj <- subset(targetObj, predscore >= 0.5)
targetObj@meta.data$predscore <- NULL

targetObj@meta.data$newsubclustering <- factor(targetObj@meta.data$newsubclustering,
                                               levels = c("PSC", "PH 1", "PH 2", "PH 3", "PH 4", "PH 5", "PM 1", "PM 2",         
                                                          "LM 1", "LM 2", "CC 1", "CC 2", "GST-rich", "Adipohemocyte", "DV", "RG", "Neurons"))
targetObj@meta.data$celltype <- mapvalues(targetObj@meta.data$newsubclustering, 
                                          from = levels(targetObj@meta.data$newsubclustering),
                                          to = c("PSC", "PH", "PH", "PH", "PH", "PH", "PM", "PM",         
                                                 "LM", "LM", "CC", "CC", "GST-rich", "Adipohemocyte", "DV", "RG", "Neurons"))
Idents(targetObj) <- 'newsubclustering'
head(targetObj@meta.data)
summary(targetObj@meta.data$newsubclustering)
nrow(targetObj@meta.data)*0.001 # 5.748

summary(targetObj@meta.data$newsubclustering)[summary(targetObj@meta.data$newsubclustering) < nrow(targetObj@meta.data)*0.001]
targetObj <- subset(targetObj, cells = rownames(subset(targetObj@meta.data, !newsubclustering %in% c('Adipohemocyte', 'DV'))) ) # minor populations
targetObj <- subset(targetObj, cells = rownames(subset(targetObj@meta.data, !newsubclustering %in% c('RG', 'Neurons', 'PSC'))) ) # non-hemocytes
targetObj@meta.data <- droplevels(targetObj@meta.data)

DimPlot(targetObj, reduction = 'umap', label = T)
summary(targetObj@meta.data$newsubclustering)

saveRDS(targetObj, 'tmp/dropseq.mtcut.norm.predicted.Rds')


DimPlot(targetObj, reduction = 'tsne', cols = c('#f15fa6', '#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#1e659b',
                                                '#ecd5a5', '#f6ae72', 'orchid2', 'orchid4', '#71c0b0', '#177e7d', '#a4a4a4'))
ggsave('tsne/tsne.mitoCut.1_2.subclustering.seed210116302.pdf', units = 'cm', width = 13, height = 10)
DimPlot(targetObj, reduction = 'umap', cols = c('#f15fa6', '#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#1e659b',
                                                '#ecd5a5', '#f6ae72', 'orchid2', 'orchid4', '#71c0b0', '#177e7d', '#a4a4a4'))
ggsave('umap/umap.mitoCut.1_2.subclustering.seed210116363.pdf', units = 'cm', width = 13, height = 10)


