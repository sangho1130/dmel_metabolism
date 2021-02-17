library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

dir.create('stats')
dir.create('tsne')
dir.create('tsne/compare')
dir.create('umap')
dir.create('umap/compare')
dir.create('tmp')

data <- read.delim('../tables/merged.expr.allGenes.txt', header = T, sep = '\t', check.names = F, row.names = 1)
labelData <- read.delim('../tables/merged.label.allGenes.txt', header = T, sep = '\t', check.names = F, row.names = 1)
identical(colnames(data), rownames(labelData))
head(labelData); nrow(labelData) # 6328

dropseq <- CreateSeuratObject(counts = data, project = "PI24Circ")
dropseq <- AddMetaData(object = dropseq, metadata = labelData, col.name = "Library")

plt <- ggplot(dropseq@meta.data, aes(Library, nCount_RNA)) + 
  geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-1.stats.nCount_RNA.byLibrary.pdf', units = 'cm', width = 6, height = 6)

plt <- ggplot(dropseq@meta.data, aes(Library, nFeature_RNA)) + 
  geom_jitter(size = 0.25) + 
  geom_hline(yintercept = c(200, 300), col = c('steelblue2', 'red2')) +
  theme_bw() + 
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-1.stats.nFeature_RNA.byLibrary.pdf', units = 'cm', width = 6, height = 8)

dropseq <- subset(x = dropseq, subset = nCount_RNA < 70000)
plt <- ggplot(dropseq@meta.data, aes(Library, nCount_RNA)) + 
  geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-1.stats.nCount_RNA.byLibrary_v2.pdf', units = 'cm', width = 6, height = 6)

saveRDS(dropseq, 'tmp/dropseq.Rds')


### Count threshold ###
split_dropseq <- SplitObject(dropseq, split.by='Library')
newdf <- data.frame(row.names = rownames(dropseq))
newlab <- data.frame()
for (tmpobj in split_dropseq) {
  head(tmpobj@meta.data)
  tmpobj <- subset(tmpobj, subset = nCount_RNA < as.integer(mean(tmpobj@meta.data$nCount_RNA) + 2*sd(tmpobj@meta.data$nCount_RNA)))
  tmpobj <- subset(tmpobj, subset = nFeature_RNA >= 300)
  
  newdf <- cbind(newdf, as.matrix(GetAssayData(tmpobj, slot = 'counts')))
  newlab <- rbind(newlab, tmpobj@meta.data)
}

dropseq <- CreateSeuratObject(counts = newdf, project = "PI24CIRC") #
dropseq <- AddMetaData(object = dropseq, metadata = newlab)
nrow(dropseq@meta.data) # 5996 cells
head(dropseq@meta.data)


plt <- ggplot(dropseq@meta.data, aes(Library, nCount_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-2.stats.nCount_RNA.byLibrary.pdf', units = 'cm', width = 6, height = 6)
plt <- ggplot(dropseq@meta.data, aes(Library, nFeature_RNA)) + geom_jitter(size = 0.25) +
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-2.stats.nFeature_RNA.byLibrary.pdf', units = 'cm', width = 6, height = 6)

####
dropseq[["percent.mt"]] <- PercentageFeatureSet(object = dropseq, pattern = "^mt:")
head(dropseq@meta.data)

plt <- ggplot(dropseq@meta.data, aes(Library, percent.mt)) + geom_jitter(size = 0.25) + 
  geom_hline(yintercept = 10, col = 'red2') +
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-2.stats.percent.mt.byLibrary.pdf', units = 'cm', width = 6, height = 6)


plot1 <- FeatureScatter(object = dropseq, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'Library') + geom_abline(intercept = 10, col = 'red2', slope = 0)
plot2 <- FeatureScatter(object = dropseq, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'Library')
CombinePlots(plots = list(plot1, plot2))
ggsave('stats/2-2.stats.percent.mt_nFeature_bynCount_RNA.byLibrary.pdf', units = 'cm', width = 26, height = 10)

dropseq.mtcut <- subset(dropseq, subset = percent.mt < 20)
nrow(dropseq.mtcut@meta.data) # 5986 
plot1 <- FeatureScatter(object = dropseq.mtcut, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'Library')
plot2 <- FeatureScatter(object = dropseq.mtcut, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'Library')
CombinePlots(plots = list(plot1, plot2))
ggsave('stats/2-2.stats.percent.mt_nFeature_bynCount_RNA.byLibrary.mitoCut.pdf', units = 'cm', width = 26, height = 10)

dropseq.mtcut.norm <- NormalizeData(object = dropseq.mtcut, normalization.method = "LogNormalize", scale.factor = 10000)
dropseq.mtcut.norm <- FindVariableFeatures(object = dropseq.mtcut.norm, selection.method = "vst", nfeatures = 2000)

top10 <- head(x = VariableFeatures(object = dropseq.mtcut.norm), 10)
plot1 <- VariableFeaturePlot(object = dropseq.mtcut.norm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE); plot2
ggsave('stats/2-3.FindVariableFeatures.mitoCut.pdf', units = 'cm', width = 17, height = 10)

head(dropseq.mtcut.norm@meta.data)
dropseq.mtcut.norm <- ScaleData(object = dropseq.mtcut.norm, vars.to.regress = c('Library', 'nCount_RNA'))
dropseq.mtcut.norm <- RunPCA(object = dropseq.mtcut.norm, features = VariableFeatures(dropseq.mtcut.norm), npcs = 50)
dropseq.mtcut.norm <- JackStraw(object = dropseq.mtcut.norm, num.replicate = 100, dims = 50)
dropseq.mtcut.norm <- ScoreJackStraw(object = dropseq.mtcut.norm, dims = 1:50)

JackStrawPlot(dropseq.mtcut.norm, dims = 1:50) # 33 PC
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 25, height = 12)
saveRDS(dropseq.mtcut.norm, 'tmp/dropseq.mtcut.norm.Rds')


#dropseq.mtcut.norm <- readRDS('tmp/dropseq.mtcut.norm.Rds')
genes <- c("E(spl)m3-HLH", "NimB3" , "IM18", "Ance", "Hml", "Nplp2", "NimC1", "PPO1" , "atilla")

### UMAP ###
for (seed in c(210116350:210116369)){
  dropseq.mtcut.norm <- RunUMAP(dropseq.mtcut.norm, dims=1:33, reduction = "pca", reduction.key='UMAP', n.components=3, min.dist=0.2, seed.use = seed)
  FeaturePlot(dropseq.mtcut.norm, features = genes, cols = c("grey","red"), reduction = "umap")
  ggsave(paste0(c('umap/compare/test.mitoCut.1_2.markers.seed', seed, '.png'), collapse = ''), units = 'cm', width = 35, height = 30)
} # 210116363

### tSNE ###
for (seed in c(210116300:210116319)){
  dropseq.mtcut.norm <- RunTSNE(dropseq.mtcut.norm, dims=1:33, reduction = "pca", reduction.key='tSNE', dim.embed=3, seed.use = seed)
  FeaturePlot(dropseq.mtcut.norm, features = genes, cols = c("grey","red"), reduction = "tsne")
  ggsave(paste0(c('tsne/compare/test.mitoCut.1_2.markers.seed', seed, '.png'), collapse = ''), units = 'cm', width = 35, height = 30)
} # 210116302; 210116305

dropseq.mtcut.norm <- RunUMAP(dropseq.mtcut.norm, dims=1:33, reduction = "pca", reduction.key='UMAP', n.components=3, min.dist=0.2, seed.use = 210116363)
FeaturePlot(dropseq.mtcut.norm, features = genes, cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.mitoCut.1_2.markers.seed210116363.png', units = 'cm', width = 35, height = 30)

DimPlot(object=dropseq.mtcut.norm, group.by = "Library", reduction = 'umap')
ggsave('umap/umap.mitoCut.1_2.libraries.seed210116363.pdf', units = 'cm', width = 14, height = 10)
FeaturePlot(dropseq.mtcut.norm, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.mitoCut.1_2.nFeature_RNA.seed210116363.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.mtcut.norm, features = 'percent.mt', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.mitoCut.1_2.percent.mt.seed210116363.pdf', units = 'cm', width = 11, height = 10)


dropseq.mtcut.norm <- RunTSNE(dropseq.mtcut.norm, dims=1:33, reduction = "pca", reduction.key='tSNE', dim.embed=3, seed.use = 210116302)
FeaturePlot(dropseq.mtcut.norm, features = genes, cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.mitoCut.1_2.markers.seed210116302.png', units = 'cm', width = 35, height = 30)

DimPlot(object=dropseq.mtcut.norm, group.by = "Library", reduction = 'tsne')
ggsave('tsne/tsne.mitoCut.1_2.libraries.seed210116302.pdf', units = 'cm', width = 14, height = 10)
FeaturePlot(dropseq.mtcut.norm, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.mitoCut.1_2.nFeature_RNA.seed210116302.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.mtcut.norm, features = 'percent.mt', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.mitoCut.1_2.percent.mt.seed210116302.pdf', units = 'cm', width = 11, height = 10)

saveRDS(dropseq.mtcut.norm, 'tmp/dropseq.mtcut.norm.Rds')

writeLabel <- as.matrix(dropseq.mtcut.norm@meta.data)
writeLabel <- data.frame(Barcode = rownames(writeLabel), writeLabel, check.names = F)
head(writeLabel)
write.table(writeLabel, 'pi24.circ.label', sep = '\t', quote = F, row.names = F, col.names = T)

