#####
rm(list=ls() )
options(stringsAsFactors = F)


####qc
rawdata_cs11 <- Read10X('/data3/dd/CS11_YS/hg19_v3/')
rawdata_cs15 <- Read10X('/data3/dd/cs15_10x/ys_cs15/outs/filtered_gene_bc_matrices/hg19/')

colnames(rawdata_cs11) <- paste(colnames(rawdata_cs11),'_cs11',sep="")
colnames(rawdata_cs15) <- paste(colnames(rawdata_cs15),'_cs15',sep="")

co_gene  <- intersect(rownames(rawdata_cs11), rownames(rawdata_cs15))
all_rawdata <- cbind(as.matrix(rawdata_cs11[co_gene,]), as.matrix(rawdata_cs15[co_gene,]))


all_anno <- data.frame(nGene = colSums(all_rawdata>0), nUMI = colSums(all_rawdata), row.names = colnames(all_rawdata) )

a <- ggplot(all_anno, aes(x=nGene, y=nUMI))+geom_point(col = 'grey')+stat_density_2d()+theme_bw()+theme(axis.title = element_text(size=15), axis.text = element_text(size=15))+geom_vline(xintercept = 1000, lty =2)

b <- ggplot(all_anno)+geom_histogram(aes(x=nGene, y= ..density..), fill='lightblue', col='black')+geom_density(aes(x=nGene, y= ..density..))+theme_bw()+theme(axis.title = element_text(size=15), axis.text = element_text(size=15))
c <- ggplot(all_anno)+geom_histogram(aes(x=nUMI, y= ..density..), fill='lightblue', col='black')+geom_density(aes(x=nUMI, y= ..density..))+theme_bw()+theme(axis.title = element_text(size=15), axis.text = element_text(size=15))

ggpubr::ggarrange(a,b,c, ncol = 3)

library(Seurat)
options(future.globals.maxSize = 8000 * 1024^2)

ana <- ReadH5AD('/data3/dd/project/cs11_cs15_combine_analysis/cs11_cs15_ys_combine2.h5ad')

# ys.list <- SplitObject(ana, split.by = "batch")
# #ys.list <- ys.list[c('0','1')]
# 
# for (i in 1:length(ys.list)) {
#   ys.list[[i]] <- SCTransform(ys.list[[i]], verbose = FALSE)
# }
# 
# ys.features <- SelectIntegrationFeatures(object.list = ys.list, nfeatures = 3000)
# ys.list <- PrepSCTIntegration(object.list = ys.list, anchor.features = ys.features, 
#                                     verbose = FALSE)
# 
# ys.anchors <- FindIntegrationAnchors(object.list = ys.list, normalization.method = "SCT", 
#                                            anchor.features = ys.features, verbose = FALSE)
# ys.integrated <- IntegrateData(anchorset = ys.anchors, normalization.method = "SCT", 
#                                      verbose = FALSE)
# 
# 
# ys.integrated <- RunPCA(ys.integrated, verbose = FALSE)
# ElbowPlot(ys.integrated, ndims = 50)
# 
# save(file='tmp.Rdata', list=ls())
# load('tmp.Rdata')
# rm(ys.anchors, ys.list)
# 
# library(umap)
# umap_data <- umap(ys.integrated@reductions$pca@cell.embeddings[,1:30])
# umap_coor <- umap_data$layout
# pca_data <- ys.integrated@reductions$pca@cell.embeddings
# ys.integrated@reductions$pca@cell.embeddings[] <- as.matrix(umap_coor)
# 
# names(ys.integrated@meta.data)
# a <- DimPlot(ys.integrated, group.by = 'batch')
# 
# FeaturePlot(ys.integrated, features = c('EPCAM','CDH5','PDGFRA','CD14','CD34','GATA1','GYPA','PF4','HBB'),ncol = 3,
#             cols = c('#E6E6E6','brown'),slot = 'counts')+NoLegend()
# 
# ys.integrated <- FindNeighbors(ys.integrated, dims = 1:10)
# ys.integrated <- FindClusters(ys.integrated, resolution = seq(0.1,2,0.1))             
# 
# b <- DimPlot(ys.integrated, group.by = 'integrated_snn_res.0.1', label = T)
# cowplot::plot_grid(a,b)
# #####
# FeaturePlot(ys.integrated, features = c('RGS1'))
# FeaturePlot(ys.integrated, features = c("EPCAM",'LUM','CD14','CD34','CDH5','GYPA','PF4','MYB'),ncol = 3)
# 
# 
# ####使用seurat的旧版的标准化的方法，进行校正
# ys2.list <- SplitObject(ana, split.by = "batch")
# 
# for (i in 1:length(ys2.list)) {
#   ys2.list[[i]] <- NormalizeData(ys2.list[[i]], verbose = FALSE)
#   ys2.list[[i]] <- FindVariableFeatures(ys2.list[[i]], selection.method = "vst", 
#                                              nfeatures = 3000, verbose = FALSE)
# }
# 
# 
# reference.list <- ys2.list[c("0", "1")]
# ys2.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
# ys2.integrated <- IntegrateData(anchorset = ys2.anchors, dims = 1:30)
# 
# DefaultAssay(ys2.integrated) <- "integrated"
# # Run the standard workflow for visualization and clustering
# ys2.integrated <- ScaleData(ys2.integrated, verbose = FALSE)
# ys2.integrated <- RunPCA(ys2.integrated, npcs = 30, verbose = FALSE)
# ElbowPlot(ys2.integrated)
# 
# library(umap)
# umap_data <- umap(ys2.integrated@reductions$pca@cell.embeddings[,1:10])
# umap_coor <- umap_data$layout
# pca_data <- ys2.integrated@reductions$pca@cell.embeddings[]
# ys2.integrated@reductions$pca@cell.embeddings[] <- as.matrix(umap_coor)
# 
# DimPlot(ys2.integrated, group.by = 'batch')
# FeaturePlot(ys2.integrated, features = c("EPCAM",'LUM','CD14','CD34','CDH5','GYPA','PF4','MYB'),ncol = 3)
# 
# 
####不去除批次效应
ana2 <- ana
ana2 <- SCTransform(ana2, return.only.var.genes = F)
ana2 <- RunPCA(ana2, features = VariableFeatures(ana2), npcs = 50)
ElbowPlot(ana2, ndims = 50)

library(umap)
umap_data <- umap(ana2@reductions$pca@cell.embeddings[,1:10])
umap_coor <- umap_data$layout
pca_data <- ana2@reductions$pca@cell.embeddings[]
ana2@reductions$pca@cell.embeddings[] <- as.matrix(umap_coor)

DimPlot(ana2, group.by = 'batch')
FeaturePlot(ana2, features = c("EPCAM",'LUM','CD14','CD34','CDH5','GYPA','PF4','MYB'),ncol = 3)

###去掉潜在的doublets
p <- DimPlot(ana2, group.by = 'batch')
doublets <- CellSelector(p)

ana_cells <- setdiff(colnames(ana2),doublets)
ana3 <- subset(ana, cells = ana_cells)

ys3.list <- SplitObject(ana3, split.by = "batch")

for (i in 1:length(ys3.list)) {
  ys3.list[[i]] <- NormalizeData(ys3.list[[i]], verbose = FALSE)
  ys3.list[[i]] <- FindVariableFeatures(ys3.list[[i]], selection.method = "vst", 
                                        nfeatures = 3000, verbose = FALSE)
}


reference.list <- ys3.list[c("0", "1")]
ys3.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
ys3.integrated <- IntegrateData(anchorset = ys3.anchors, dims = 1:30)

DefaultAssay(ys3.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
ys3.integrated <- ScaleData(ys3.integrated, verbose = FALSE)
ys3.integrated <- RunPCA(ys3.integrated, verbose = FALSE)
ElbowPlot(ys3.integrated)


ys3.integrated <- FindNeighbors(ys3.integrated, dims = 1:10)
ys3.integrated <- FindClusters(ys3.integrated, resolution = seq(0.1,2,0.1))


library(umap)
umap_data <- umap(ys3.integrated@reductions$pca@cell.embeddings[,1:10])
umap_coor <- umap_data$layout
pca_data <- ys3.integrated@reductions$pca@cell.embeddings[]
ys3.integrated@reductions$pca@cell.embeddings[] <- as.matrix(umap_coor)

DimPlot(ys3.integrated, group.by = 'batch')
FeaturePlot(ys3.integrated, features = c("rna_GJA5",'rna_CD14','rna_CD34','rna_CDH5','rna_GYPA','rna_PF4','rna_RUNX1','rna_CSF1R'),ncol = 3)

a <- DimPlot(ys3.integrated, group.by = 'batch')
b <- DimPlot(ys3.integrated, group.by = 'integrated_snn_res.0.2', label = T)
cowplot::plot_grid(a,b)

ys3.integrated@meta.data$cluster <- as.character(ys3.integrated@meta.data$integrated_snn_res.0.2)
ys3.integrated@meta.data$cluster <- plyr::mapvalues(ys3.integrated@meta.data$cluster,
                                                    as.character(0:7),
                                                    c('Mes1','EMP','Epi','Mac','ErP/MkP','Endo1','Endo2','Mes2'))

# cluster_order <- c('EMP','ErP/MkP','Mac', 'Endo1','Endo2','Epi','Mes1',"Mes2")
# col_cluster <- setNames(c('red3','orange3','green3','purple3','blue3','magenta','#C16200','#71226e'),
#                         cluster_order)
# 
# 
# library(ggplot2)
# a <- DimPlot(ys3.integrated, pt.size = 1, group.by = 'cluster', label = T, label.size = 8)+scale_color_manual(values = col_cluster)+
#   theme_bw()+theme(legend.text = element_text(size=15), axis.text = element_text(size = 15),axis.title = element_text(size = 15))+labs(x='UMAP1',y="UMAP2")
# 
# ys3.integrated@meta.data$batch <- plyr::mapvalues(ys3.integrated@meta.data$batch, c('0','1'), c('CS 11','CS 15'))
# 
# col_stage <- setNames(c('green3','red3'), c('CS 11','CS 15'))
# 
# b <- DimPlot(ys3.integrated, pt.size = 1, group.by = 'batch')+scale_color_manual(values = col_stage)+
#   theme_bw()+theme(legend.text = element_text(size=15), axis.text = element_text(size = 15),axis.title = element_text(size = 15))+labs(x='UMAP1',y="UMAP2")
# 
# ggpubr::ggarrange(a,b)
# 
# FeaturePlot(ys3.integrated, features = c("rna_LUM",'rna_EPCAM','rna_CDH5','rna_CD34','rna_MYB','rna_CD14','rna_GATA1','rna_GYPA','rna_PF4'),ncol = 3,
#             cols = c("#E6E6E6", 'brown'), order = T)
# 
unique(ys3.integrated@meta.data$cluster)

ys3.integrated <- SetIdent(ys3.integrated, value = 'cluster')
degs_all <- FindAllMarkers(ys3.integrated, logfc.threshold = log(1.5), only.pos = T)
degs_all <- degs_all[degs_all$p_val_adj<0.05,]
write.csv(file = '/home/dd/10x_cs11_cs15/DEGs_all_clusters_10x.csv', degs_all)
# 
###plot feature genes matrix
source('/data3/dd/code/plot_gene_violin_matrix.R')
gene <- c("RUNX1","CD34", "MYB", "GATA1","PF4", "KLF1", "CD14", "CD163", "CDH5", "SOX7", "EPCAM", "PDGFRA")

cluster_order <- c('EMP1','EMP2','ErP/MkP','Mac', 'Endo','Epi','Mes1',"Mes2")
col_cluster <- setNames(c('red3','orange3','green3','purple3','blue3','magenta','#C16200','#71226e'),
                        cluster_order)

anno <- ys3.integrated@meta.data[,'cluster',drop=F]
anno$cluster[anno$cluster=='EMP'] <- 'EMP2'
anno$cluster[anno$cluster=='Endo2'] <- 'EMP1'
anno$cluster[anno$cluster=='Endo1'] <- 'Endo'

library(reshape2)
anno$cluster <- factor(anno$cluster, levels = cluster_order)
result <- plot_gene_matrix(data = as.matrix(ys3.integrated@assays$RNA@data), gene = gene, anno = anno, col = col_cluster)
result+theme_bw()+theme(axis.text = element_text(size=12),axis.text.x =  element_text(angle = 60, hjust = 1), axis.title = element_blank())


ys3.integrated <- AddMetaData(ys3.integrated, metadata = anno)
DimPlot(ys3.integrated, pt.size = 1, group.by = 'cluster', label = T, label.size = 8)+scale_color_manual(values = col_cluster)+
  theme_bw()+theme(legend.text = element_text(size=15), axis.text = element_text(size = 15),axis.title = element_text(size = 15))+labs(x='UMAP1',y="UMAP2")

###EMP markers, use rawdata
anno$Cluster <- as.character(anno$cluster)
anno$Cluster[anno$Cluster %in% c('EMP1','EMP2')]<-'EMP'


ys3.integrated2 <- subset(ana2, cells = rownames(anno))

ys3.integrated2 <- AddMetaData(ys3.integrated, metadata = anno[colnames(ys3.integrated2),])
ys3.integrated2 <- SetIdent(ys3.integrated2, value = 'Cluster')
ys3.integrated2 <- subset(ys3.integrated2, idents = c('EMP','ErP/MkP','Mac'))
ys3.integrated2@meta.data$Cluster[!ys3.integrated2@meta.data$Cluster=='EMP']<- 'others'
ys3.integrated2 <- SetIdent(ys3.integrated2, value = 'Cluster')
###
deg_emp <- FindAllMarkers(ys3.integrated2, logfc.threshold = log(1.5), only.pos = T)
deg_emp <- deg_emp[deg_emp$p_val_adj<0.05,]

markers <- read.csv('/data3/dd/data/human_gene/human_cd_moculer.csv')
deg_emp <- deg_emp[rownames(deg_emp) %in% markers$symbols & deg_emp$cluster=='EMP',]
deg_emp$gene <- rownames(deg_emp)

ggplot(deg_emp)+geom_bar(stat = 'identity', aes(x = rev(1:nrow(deg_emp)), y= -log10(p_val_adj), fill = avg_logFC))+coord_flip()+
  scale_fill_gradientn(colours = colorRampPalette(c('#E6E6E6',' brown'))(100))+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(size=15))+
  geom_text(aes(x = rev(1:nrow(deg_emp)), y= -log10(p_val_adj)+15, label = gene))

write.csv(file = '/home/dd/10x_cs11_cs15/emp_marker_10x.csv', deg_emp)
###mac
main_cell <- colnames(ys3.integrated)[ys3.integrated@meta.data$cluster =='Mac']
mac <- subset(ana2, cells = main_cell)

save(file = 'ys_10x.Rdata', ana, ys3.integrated, mac)
mac <- FindVariableFeatures(mac, selection.method = 'vst', nfeatures = 2000)
mac <- ScaleData(mac)
mac <- RunPCA(mac, features = VariableFeatures(mac), npcs = 50)
ElbowPlot(mac, ndims = 50)
mac <- RunTSNE(mac, dims = 1:10)

mac <- FindNeighbors(mac, dims = 1:10)
mac <- FindClusters(mac, resolution = seq(0.1,2,0.1))

# umap_data2 <- umap(mac@reductions$pca@cell.embeddings[,1:10])
# umap_coor2 <- umap_data2$layout
# pca_data2 <- mac@reductions$pca@cell.embeddings[]
# mac@reductions$pca@cell.embeddings[] <- as.matrix(umap_coor2)
mac@meta.data$batch <- plyr::mapvalues(mac@meta.data$batch, c('0','1'),
                                       c('CS 11','CS 15'))


a <- DimPlot(mac, group.by = 'batch', reduction = 'pca')
b <- DimPlot(mac, group.by = 'SCT_snn_res.0.1', reduction = 'pca')
cowplot::plot_grid(a,b)

FeaturePlot(mac, features = c('BNIP3','LSP1','RGS1','ANXA2','MMP9','FABP3'),ncol = 3)
FeaturePlot(mac, features = c('CD34','MYB','GATA1','GYPA','PF4','CD14','RUNX1','nGene','CDH5'),ncol = 3)

####smartseq degs
VlnPlot(mac, features =degs_smart$gene[degs_smart$cluster=='Mac_1'][1:5],group.by = 'SCT_snn_res.0.1', ncol = 5, pt.size = 0.5)+

####
mac@meta.data$cluster <- as.character(mac@meta.data$SCT_snn_res.0.1)
mac@meta.data$cluster[mac@meta.data$cluster %in% c('1','2')]<-'Mac'
mac@meta.data$cluster[mac@meta.data$cluster %in% c('0','3')]<-'EMP'

mac <- SetIdent(mac, value = 'cluster')

####
mac <- SetIdent(mac, value = 'SCT_snn_res.0.1')
degs_mac <- FindAllMarkers(mac, logfc.threshold = log(1.5), only.pos = T)
degs_mac <- degs_mac[degs_mac$p_val_adj<0.05,]

degs_smart <- read.csv('/home/dd/deg_smartseq.csv', header = T)
degs_smart2 <- read.csv('/home/dd/degs_only_cs11_cs15.csv', header = T)

intersect(degs_smart2$gene[degs_smart2$cluster=='Mac_1'],degs_mac$gene[as.character(degs_mac$cluster)=='1'])
intersect(degs_smart2$gene[degs_smart2$cluster=='Mac_2'],degs_mac$gene[as.character(degs_mac$cluster)=='2'])


FeaturePlot(mac, features = degs_smart$gene[degs_smart$cluster=='Mac_1'][1:5],ncol = 5,
            cols = c('blue','green','yellow','red'), reduction = 'pca')




load('/data3/dd/STRT_Mac_DEG.Rdata')

mac_strt <- UpdateSeuratObject(mac_strt)

mac@meta.data$Cluster <- factor(mac@meta.data$SCT_snn_res.0.1, levels = c('1','0')) 
mac_strt@meta.data$Cluster <- factor(mac_strt@meta.data$cluster, levels = c('Mac_2','Mac_1')) 


a <- DotPlot(mac,col.min = -1, col.max = 1,features = c(result2$gene[result2$cluster=='Mac_2'][10:1],result2$gene[result2$cluster=='Mac_1'][10:1]), group.by = 'Cluster')+
  theme_bw()+theme(axis.text.y = element_text(size = 15),axis.text.x = element_text(size = 15, color = 'black', angle = 60, hjust = 1), axis.title = element_text(size= 15),legend.text = element_text(size = 15), legend.title = element_text(size=15))+scale_y_discrete(labels = c('Mac_2','Mac_1'))+
  labs(y="", title = '10X Genomics', x="")+scale_color_gradientn(colours = colorRampPalette(c('#E6E6E6','blue','green','yellow','red'))(100))

b <- DotPlot(mac_strt,  features = c(result2$gene[result2$cluster=='Mac_2'][10:1],result2$gene[result2$cluster=='Mac_1'][10:1]), group.by = 'Cluster')+
  theme_bw()+theme(axis.text.y = element_text(size = 15),axis.text.x = element_text(size = 15, angle = 60, hjust = 1),legend.text = element_text(size = 15), legend.title = element_text(size=15), axis.title = element_text(size= 15))+scale_y_discrete(labels = c('Mac_2','Mac_1'))+
  labs(y="", title = 'STRT-seq',x="")+scale_color_gradientn(colours = colorRampPalette(c('#E6E6E6','blue','green','yellow','red'))(100))

ggpubr::ggarrange(b,a, nrow = 2, ncol = 1, common.legend = T, legend = 'right')

####


VlnPlot(mac, features = degs_smart$gene[degs_smart$cluster=='Mac_1'][1:10], group.by = 'SCT_snn_res.0.1', ncol=5, pt.size = 0.8)+scale_x_discrete(labels = c('Mac_1','Mac_2'))
VlnPlot(mac, features = degs_smart$gene[degs_smart$cluster=='Mac_2'][1:10], group.by = 'SCT_snn_res.0.1', ncol=5, pt.size = 0.8)+scale_x_discrete(labels = c('Mac_1','Mac_2'))


####save data 
save(file = '/data3/dd/project/cs11_cs15_combine_analysis/cs11_cs15_ys.Rdata', ana, ys3.integrated, mac, col_cluster, col_stage,)

######plot dotplot

#######ggplot 
library(ggplot2)

genes <- c(result2$gene[result2$cluster=='Mac_1'][10:1],result2$gene[result2$cluster=='Mac_2'][10:1])

data <- mac_strt@assays$RNA@data
phmat <- t(scale(t(data)))
gene_data <-  data.frame(as.matrix(t(phmat[genes, ])))
gene_data$Cluster <- mac_strt@meta.data[rownames(gene_data), 'Cluster']

aggr_data <- aggregate(.~Cluster, gene_data, mean)
melt_data <- melt(aggr_data, id='Cluster')
names(melt_data) <- c('Cluster','gene','expression')
for (i in 1:nrow(melt_data)){
  
  celltype <- as.character(melt_data$Cluster)[i]
  cell <- rownames(mac_strt@meta.data)[as.character(mac_strt@meta.data$Cluster)==celltype]
  gene <- as.character(melt_data$gene)[i]
  
  
  melt_data$proportion[i] <- sum(data[gene, cell]>0)/length(cell)
  
}



ggplot(melt_data, aes(x = gene, y=Cluster, fill= expression, size =proportion))+geom_point(pch =21, color = "grey30")+
  scale_fill_gradientn(colours = colorRampPalette(c('#E6E6E6','blue','green','yellow','red'))(100))+theme_bw()+theme(axis.text.x = element_text(angle = 60, size = 10, hjust = 1))

##
data <- mac@assays$SCT@data
phmat <- t(scale(t(data)))
gene_data <-  data.frame(as.matrix(t(phmat[genes, ])))
gene_data$Cluster <- mac@meta.data[rownames(gene_data), 'Cluster']

aggr_data <- aggregate(.~Cluster, gene_data, mean)
melt_data <- melt(aggr_data, id='Cluster')
names(melt_data) <- c('Cluster','gene','expression')
for (i in 1:nrow(melt_data)){
  
  celltype <- as.character(melt_data$Cluster)[i]
  cell <- rownames(mac@meta.data)[as.character(mac@meta.data$Cluster)==celltype]
  gene <- as.character(melt_data$gene)[i]
  
  
  melt_data$proportion[i] <- sum(data[gene, cell]>0)/length(cell)
  
}



ggplot(melt_data, aes(x = gene, y=Cluster, fill= expression, size =proportion))+geom_point(pch =21, color = "grey30")+
  scale_fill_gradientn(colours = colorRampPalette(c('#E6E6E6','blue','green','yellow','red'))(100))+theme_bw()+theme(axis.text.x = element_text(angle = 60, size = 10, hjust = 1))
