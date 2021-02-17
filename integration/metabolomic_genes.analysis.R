library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)
library(harmony)

dir.create('tmp')
dir.create('stats')
dir.create('tsne')
dir.create('tsne/compare')
dir.create('umap')
dir.create('umap/compare')

exprs <- 'merged.count.Rds'
labels <- 'merged.label.Rds'

data <- readRDS(exprs)
labelData <- readRDS(labels)
colnames(labelData)[1] <- 'dataset'

dropseq <- CreateSeuratObject(counts = data, project = "Dmel_harmony")
dropseq <- AddMetaData(object = dropseq, metadata = labelData)

plt <- ggplot(circulation@meta.data, aes(Library, nFeature_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5));plt
ggsave('stats/2-2.stats.nFeature_RNA.byLibrary.pdf', units = 'cm', width = 30, height = 8)
plt <- ggplot(circulation@meta.data, aes(Library, nCount_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5));plt
ggsave('stats/2-2.stats.nCount_RNA.byLibrary.pdf', units = 'cm', width = 30, height = 8)

plt <- ggplot(circulation@meta.data, aes(Library, percent.mt)) + geom_jitter(size = 0.25) + 
  theme_bw() +
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5));plt
ggsave('stats/2-2.stats.percent.mt.byLibrary.pdf', units = 'cm', width = 30, height = 8)

plt <- ggplot(circulation@meta.data, aes(percent.mt, col = Library)) + geom_density() +
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-2.stats.percent.mt.density_byLibrary.pdf', units = 'cm', width = 24, height = 12)
plt <- ggplot(circulation@meta.data, aes(nCount_RNA, col = Library)) + geom_density() +
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-2.stats.nCount_RNA.density_byLibrary.pdf', units = 'cm', width = 12, height = 6)
plt <- ggplot(circulation@meta.data, aes(nFeature_RNA, col = Library)) + geom_density() +
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-2.stats.nFeature_RNA.density_byLibrary.pdf', units = 'cm', width = 12, height = 6)

plot1 <- FeatureScatter(circulation, feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by = 'Library')
plot2 <- FeatureScatter(circulation, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'Library')
CombinePlots(plots = list(plot1, plot2))
ggsave('stats/2-2.stats.percent.mt.compare.pdf', units = 'cm', width = 40, height = 10)


### The standard workflow ###
circulation <- NormalizeData(circulation, normalization.method = "LogNormalize", scale.factor = 10000)
circulation <- FindVariableFeatures(object = circulation, selection.method = "vst", nfeatures = 2000)
circulation <- ScaleData(object = circulation, vars.to.regress = c('Library', 'nCount_RNA'))
circulation <- RunPCA(object = circulation, npcs = 60)
circulation <- JackStraw(object = circulation, num.replicate = 100, dims = 60)
circulation <- ScoreJackStraw(object = circulation, dims = 1:60)
JackStrawPlot(circulation, dims = 1:60) # 40 PCs 
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 25, height = 12)

circulation <- RunHarmony(circulation, group.by.vars = "dataset")


genes <- c("nw", "NimB3", "IM18", "Ance", "Hml", "NimC1", "Ama", "mthl4", "PPO1")
### UMAP ###
for (seed in c(210216150:210216169)){
  circulation <- RunUMAP(circulation, dims=1:40, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.2, seed.use = seed)
  DimPlot(object=circulation, group.by = "celltype",  reduction = 'umap') + scale_color_manual(values = c('#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4'))
  ggsave(paste0(c('umap/compare/umap.1_2.celltype.seed', seed, '.png'), collapse = ''), units = 'cm', width = 10, height = 7)
  
  DimPlot(object=circulation, group.by = "celltype",  reduction = 'umap', dims = c(1,3)) + scale_color_manual(values = c('#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4'))
  ggsave(paste0(c('umap/compare/umap.1_3.celltype.seed', seed, '.png'), collapse = ''), units = 'cm', width = 10, height = 7)
  
  DimPlot(object=circulation, group.by = "celltype",  reduction = 'umap', dims = c(2,3)) + scale_color_manual(values = c('#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4'))
  ggsave(paste0(c('umap/compare/umap.2_3.celltype.seed', seed, '.png'), collapse = ''), units = 'cm', width = 10, height = 7)
} # 210216157 (umap2-3); 

### tSNE ###
for (seed in c(210216100:210216139)){
  circulation <- RunTSNE(circulation, dims=1:40, reduction = "harmony", reduction.key='tSNE', dim.embed=3, seed.use = seed)
  DimPlot(object=circulation, group.by = "celltype",  reduction = 'tsne') + scale_color_manual(values = c('#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4'))
  ggsave(paste0(c('tsne/compare/tsne.1_2.celltype.seed', seed, '.png'), collapse = ''), units = 'cm', width = 10, height = 7)
} # 210216105; 210216135


### UMAP fix ###
circulation <- RunUMAP(circulation, dims=1:40, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.2, seed.use = 210216157)
FeaturePlot(circulation, features = genes, cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.combined.1_2.markers.seed210216157.png', units = 'cm', width = 35, height = 30)
FeaturePlot(circulation, features = genes, cols = c("grey","red"), reduction = "umap", dims = c(1,3))
ggsave('umap/umap.combined.1_3.markers.seed210216157.png', units = 'cm', width = 35, height = 30)
FeaturePlot(circulation, features = genes, cols = c("grey","red"), reduction = "umap", dims = c(2,3))
ggsave('umap/umap.combined.2_3.markers.seed210216157.png', units = 'cm', width = 35, height = 30)

DimPlot(object=circulation, group.by = "dataset",  reduction = 'umap') + scale_color_manual(values = c('steelblue1', 'steelblue3', 'red1', 'red3'))
ggsave('umap/umap.combined.1_2.dataset.seed210216157.pdf', units = 'cm', width = 14, height = 10)
DimPlot(object=circulation, group.by = "dataset",  reduction = 'umap', dims = c(1,3)) + scale_color_manual(values = c('steelblue1', 'steelblue3', 'red1', 'red3'))
ggsave('umap/umap.combined.1_3.dataset.seed210216157.pdf', units = 'cm', width = 14, height = 10)
DimPlot(object=circulation, group.by = "dataset",  reduction = 'umap', dims = c(2,3)) + scale_color_manual(values = c('steelblue1', 'steelblue3', 'red1', 'red3'))
ggsave('umap/umap.combined.2_3.dataset.seed210216157.pdf', units = 'cm', width = 11, height = 7)
DimPlot(object=circulation, group.by = "dataset",  reduction = 'umap', dims = c(2,3)) + 
  scale_color_manual(values = c('steelblue1', 'steelblue3', 'red1', 'red3')) + 
  facet_wrap(~dataset, ncol = 2) + theme(legend.position = 'none')
ggsave('umap/umap.combined.2_3.dataset_sep.seed210216157.pdf', units = 'cm', width = 12, height = 13)

FeaturePlot(circulation, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.combined.1_2.nFeature_RNA.seed210216157.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(circulation, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "umap", dims = c(1,3))
ggsave('umap/umap.combined.1_3.nFeature_RNA.seed210216157.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(circulation, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "umap", dims = c(2,3))
ggsave('umap/umap.combined.2_3.nFeature_RNA.seed210216157.pdf', units = 'cm', width = 12, height = 10)

FeaturePlot(circulation, features = 'nCount_RNA', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.combined.1_2.nCount_RNA.seed210216157.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(circulation, features = 'nCount_RNA', cols = c("grey","red"), reduction = "umap", dims = c(1,3))
ggsave('umap/umap.combined.1_3.nCount_RNA.seed210216157.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(circulation, features = 'nCount_RNA', cols = c("grey","red"), reduction = "umap", dims = c(2,3))
ggsave('umap/umap.combined.2_3.nCount_RNA.seed210216157.pdf', units = 'cm', width = 12, height = 10)

FeaturePlot(circulation, features = 'percent.mt', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.combined.1_2.percent.mt.seed210216157.pdf', units = 'cm', width = 11, height = 10)
FeaturePlot(circulation, features = 'percent.mt', cols = c("grey","red"), reduction = "umap", dims = c(1,3))
ggsave('umap/umap.combined.1_3.percent.mt.seed210216157.pdf', units = 'cm', width = 11, height = 10)
FeaturePlot(circulation, features = 'percent.mt', cols = c("grey","red"), reduction = "umap", dims = c(2,3))
ggsave('umap/umap.combined.2_3.percent.mt.seed210216157.pdf', units = 'cm', width = 11, height = 10)

DimPlot(object=circulation, group.by = "condition",  reduction = 'umap') + scale_color_manual(values = c('#3B3BEA', '#ED0100'))
ggsave('umap/umap.combined.1_2.condition.seed210216157.pdf', units = 'cm', width = 13.5, height = 10)
DimPlot(object=circulation, group.by = "condition",  reduction = 'umap', dims = c(1,3)) + scale_color_manual(values = c('#3B3BEA', '#ED0100'))
ggsave('umap/umap.combined.1_3.condition.seed210216157.pdf', units = 'cm', width = 13.5, height = 10)
DimPlot(object=circulation, group.by = "condition",  reduction = 'umap', dims = c(2,3)) + scale_color_manual(values = c('#3B3BEA', '#ED0100'))
ggsave('umap/umap.combined.2_3.condition.seed210216157.pdf', units = 'cm', width = 10.5, height = 7)
DimPlot(object=circulation, group.by = "condition",  reduction = 'umap', dims = c(2,3)) + 
  scale_color_manual(values = c('#3B3BEA', '#ED0100')) + facet_wrap(~condition) + theme(legend.position = 'none')
ggsave('umap/umap.combined.2_3.condition_sep.seed210216157.pdf', units = 'cm', width = 11, height = 7)


DimPlot(object=circulation, group.by = "celltype",  reduction = 'umap') + scale_color_manual(values = c('#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4'))
ggsave('umap/umap.combined.1_2.celltype.seed210216157.pdf', units = 'cm', width = 10, height = 7)
DimPlot(object=circulation, group.by = "celltype",  reduction = 'umap', dims = c(1,3)) + scale_color_manual(values = c('#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4'))
ggsave('umap/umap.combined.1_3.celltype.seed210216157.pdf', units = 'cm', width = 10, height = 7)
DimPlot(object=circulation, group.by = "celltype",  reduction = 'umap', dims = c(2,3)) + scale_color_manual(values = c('#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4'))
ggsave('umap/umap.combined.2_3.celltype.seed210216157.pdf', units = 'cm', width = 10, height = 7)

subcols <-  c('#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#1e659b', #PH
              '#ecd5a5', '#f6ae72', #PM
              'orchid2', 'orchid4', '#71c0b0', '#177e7d', # LM CC
              '#a4a4a4')
DimPlot(object=circulation, group.by = "subclustering",  reduction = 'umap') + 
  scale_color_manual(values = subcols)
ggsave('umap/umap.combined.1_2.subclustering.seed210216157.pdf', units = 'cm', width = 10, height = 7)
DimPlot(object=circulation, group.by = "subclustering",  reduction = 'umap', dims = c(1,3)) + 
  scale_color_manual(values = subcols)
ggsave('umap/umap.combined.1_3.subclustering.seed210216157.pdf', units = 'cm', width = 10, height = 7)
DimPlot(object=circulation, group.by = "subclustering",  reduction = 'umap', dims = c(2,3)) + 
  scale_color_manual(values = subcols)
ggsave('umap/umap.combined.2_3.subclustering.seed210216157.pdf', units = 'cm', width = 10, height = 7)

head(circulation@meta.data)

### tSNE fix ###
circulation <- RunTSNE(circulation, dims=1:40, reduction = "harmony", reduction.key='tSNE', dim.embed=3, seed.use = 210216135)
FeaturePlot(circulation, features = genes, cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.combined.1_2.markers.seed210216135.png', units = 'cm', width = 35, height = 30)
FeaturePlot(circulation, features = genes, cols = c("grey","red"), reduction = "tsne", dims = c(1,3))
ggsave('tsne/tsne.combined.1_3.markers.seed210216135.png', units = 'cm', width = 35, height = 30)
FeaturePlot(circulation, features = genes, cols = c("grey","red"), reduction = "tsne", dims = c(2,3))
ggsave('tsne/tsne.combined.2_3.markers.seed210216135.png', units = 'cm', width = 35, height = 30)

DimPlot(object=circulation, group.by = "dataset", reduction = 'tsne') + scale_color_manual(values = c('steelblue1', 'steelblue3', 'red1', 'red3'))
ggsave('tsne/tsne.combined.1_2.dataset.seed210216135.pdf', units = 'cm', width = 14, height = 10)
DimPlot(object=circulation, group.by = "dataset", reduction = 'tsne', dims = c(1,3)) + scale_color_manual(values = c('steelblue1', 'steelblue3', 'red1', 'red3'))
ggsave('tsne/tsne.combined.1_3.dataset.seed210216135.pdf', units = 'cm', width = 14, height = 10)
DimPlot(object=circulation, group.by = "dataset", reduction = 'tsne', dims = c(2,3)) + scale_color_manual(values = c('steelblue1', 'steelblue3', 'red1', 'red3'))
ggsave('tsne/tsne.combined.2_3.dataset.seed210216135.pdf', units = 'cm', width = 14, height = 10)

FeaturePlot(circulation, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.combined.1_2.nFeature_RNA.seed210216135.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(circulation, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "tsne", dims = c(1,3))
ggsave('tsne/tsne.combined.1_3.nFeature_RNA.seed210216135.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(circulation, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "tsne", dims = c(2,3))
ggsave('tsne/tsne.combined.2_3.nFeature_RNA.seed210216135.pdf', units = 'cm', width = 12, height = 10)

FeaturePlot(circulation, features = 'nCount_RNA', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.combined.1_2.nCount_RNA.seed210216135.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(circulation, features = 'nCount_RNA', cols = c("grey","red"), reduction = "tsne", dims = c(1,3))
ggsave('tsne/tsne.combined.1_3.nCount_RNA.seed210216135.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(circulation, features = 'nCount_RNA', cols = c("grey","red"), reduction = "tsne", dims = c(2,3))
ggsave('tsne/tsne.combined.2_3.nCount_RNA.seed210216135.pdf', units = 'cm', width = 12, height = 10)

FeaturePlot(circulation, features = 'percent.mt', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.combined.1_2.percent.mt.seed210216135.pdf', units = 'cm', width = 11, height = 10)
FeaturePlot(circulation, features = 'percent.mt', cols = c("grey","red"), reduction = "tsne", dims = c(1,3))
ggsave('tsne/tsne.combined.1_3.percent.mt.seed210216135.pdf', units = 'cm', width = 11, height = 10)
FeaturePlot(circulation, features = 'percent.mt', cols = c("grey","red"), reduction = "tsne", dims = c(2,3))
ggsave('tsne/tsne.combined.2_3.percent.mt.seed210216135.pdf', units = 'cm', width = 11, height = 10)

DimPlot(object=circulation, group.by = "condition",  reduction = 'tsne') + scale_color_manual(values = c('#3B3BEA', '#ED0100'))
ggsave('tsne/tsne.combined.1_2.condition.seed210216135.pdf', units = 'cm', width = 14, height = 10)
DimPlot(object=circulation, group.by = "condition",  reduction = 'tsne', dims = c(1,3)) + scale_color_manual(values = c('#3B3BEA', '#ED0100'))
ggsave('tsne/tsne.combined.1_3.condition.seed210216135.pdf', units = 'cm', width = 14, height = 10)
DimPlot(object=circulation, group.by = "condition",  reduction = 'tsne', dims = c(2,3)) + scale_color_manual(values = c('#3B3BEA', '#ED0100'))
ggsave('tsne/tsne.combined.2_3.condition.seed210216135.pdf', units = 'cm', width = 14, height = 10)


DimPlot(object=circulation, group.by = "celltype",  reduction = 'tsne') + scale_color_manual(values = c('#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4'))
ggsave('tsne/tsne.combined.1_2.celltype.seed210216135.pdf', units = 'cm', width = 10, height = 7)
DimPlot(object=circulation, group.by = "celltype",  reduction = 'tsne', dims = c(1,3)) + scale_color_manual(values = c('#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4'))
ggsave('tsne/tsne.combined.1_3.celltype.seed210216135.pdf', units = 'cm', width = 10, height = 7)
DimPlot(object=circulation, group.by = "celltype",  reduction = 'tsne', dims = c(2,3)) + scale_color_manual(values = c('#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4'))
ggsave('tsne/tsne.combined.2_3.celltype.seed210216135.pdf', units = 'cm', width = 10, height = 7)

subcols <-  c('#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#1e659b', #PH
              '#ecd5a5', '#f6ae72', #PM
              'orchid2', 'orchid4', '#71c0b0', '#177e7d', # LM CC
              '#a4a4a4')
DimPlot(object=circulation, group.by = "subclustering",  reduction = 'tsne') + 
  scale_color_manual(values = subcols)
ggsave('tsne/tsne.combined.1_2.subclustering.seed210216135.pdf', units = 'cm', width = 10, height = 7)
DimPlot(object=circulation, group.by = "subclustering",  reduction = 'tsne', dims = c(1,3)) + 
  scale_color_manual(values = subcols)
ggsave('tsne/tsne.combined.1_3.subclustering.seed210216135.pdf', units = 'cm', width = 10, height = 7)
DimPlot(object=circulation, group.by = "subclustering",  reduction = 'tsne', dims = c(2,3)) + 
  scale_color_manual(values = subcols)
ggsave('tsne/tsne.combined.2_3.subclustering.seed210216135.pdf', units = 'cm', width = 10, height = 7)


### markers ###
avg_logFC <- 2

Idents(circulation) <- 'celltype_condition'
markers <- FindAllMarkers(object = circulation, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.celltype_condition.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = circulation, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('degs/markers.celltype_condition.pdf', units = 'cm', width = 30, height = 20)

Idents(circulation) <- 'celltype'
markers <- FindAllMarkers(object = circulation, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.celltype.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = circulation, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('degs/markers.celltype.pdf', units = 'cm', width = 30, height = 20)

Idents(circulation) <- 'subclustering'
markers <- FindAllMarkers(object = circulation, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.subclustering.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = circulation, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('degs/markers.subclustering.pdf', units = 'cm', width = 40, height = 30)


### known markers - figure ###
#dir.create('knownmarkers')
genes <- c(#"Dl", "dome", "shg", "IM18", 'CecA2', 
           "Tep4", "Ance", "Nplp2", 
           "Hml", "Pxn", "NimC1", "eater", "vkg", 
           "lz", "peb", "PPO1", "PPO2", 
           "mthl4", "atilla", "PPO3", "CG18547", "CG3397", "dysf")


circulation@meta.data$celltype_condition <- factor(circulation@meta.data$celltype_condition, levels = rev(levels(circulation@meta.data$celltype_condition)))
DotPlot(circulation, group.by = 'celltype_condition', features = rev(genes)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('knownmarkers/dotplot.celltype_condition.knownmarkers_1.pdf', units = 'cm', width = 20, height = 10) 

circulation@meta.data$celltype <- factor(circulation@meta.data$celltype, 
                                         levels = rev(c("PH", "PM", "CC", "LM", "GST-rich")))
DotPlot(circulation, group.by = 'celltype', features = rev(genes)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('knownmarkers/dotplot.celltype.knownmarkers_1.pdf', units = 'cm', width = 20, height = 10) 


