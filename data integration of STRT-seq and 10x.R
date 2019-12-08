#####10x YS and STRT-seq YS
rm(list=ls() )
options(stringsAsFactors = F)


####qc
rawdata_cs11 <- Read10X('/data3/dd/CS11_YS/hg19_v3/')
rawdata_cs15 <- Read10X('/data3/dd/cs15_10x/ys_cs15/outs/filtered_gene_bc_matrices/hg19/')

colnames(rawdata_cs11) <- paste('ys_cs11:',colnames(rawdata_cs11),'x-0',sep="")
colnames(rawdata_cs15) <- paste('ys_cs15:',colnames(rawdata_cs15),'x-1',sep="")


load('/data3/dd/project/cs11_cs15_combine_analysis/fig1_v5.rds')
ys_anno <- all_anno[all_anno$Site=='YS', ]
ys_rawdata <- rawdata[, rownames(ys_anno)]

#####quality control
co_gene  <- intersect(rownames(rawdata_cs11), rownames(rawdata_cs15))
all_rawdata <- cbind(as.matrix(rawdata_cs11[co_gene,]), as.matrix(rawdata_cs15[co_gene,]))


all_anno <- data.frame(nGene = colSums(all_rawdata>0), nUMI = colSums(all_rawdata), row.names = colnames(all_rawdata) )

a <- ggplot(all_anno, aes(x=nGene, y=nUMI))+geom_point(col = 'grey')+stat_density_2d()+theme_bw()+theme(axis.title = element_text(size=15), axis.text = element_text(size=15))+geom_vline(xintercept = 1000, lty =2)

b <- ggplot(all_anno)+geom_histogram(aes(x=nGene, y= ..density..), fill='lightblue', col='black')+geom_density(aes(x=nGene, y= ..density..))+theme_bw()+theme(axis.title = element_text(size=15), axis.text = element_text(size=15))
c <- ggplot(all_anno)+geom_histogram(aes(x=nUMI, y= ..density..), fill='lightblue', col='black')+geom_density(aes(x=nUMI, y= ..density..))+theme_bw()+theme(axis.title = element_text(size=15), axis.text = element_text(size=15))

ggpubr::ggarrange(a,b,c, ncol = 3)

####
library(Seurat)
options(future.globals.maxSize = 8000 * 1024^2)

ana <- ReadH5AD('/data3/dd/project/cs11_cs15_combine_analysis/cs11_cs15_ys_combine2.h5ad')

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

#####
library(reshape2)
anno_10x <- colsplit(ana_cells, ':', names = c('batch', 'cell'))[,1,drop=F]
rownames(anno_10x) <- ana_cells
anno_10x$Stage <- toupper(colsplit(anno_10x$batch, '_', names = c('site', 'Stage'))[,2])

ys_anno$batch <- 'ys_STRT'
all_anno <- data.frame(rbind(anno_10x[, c('Stage','batch')],
                             ys_anno[,c('Stage','batch')]))


co_gene <- intersect(rownames(all_rawdata), rownames(ys_rawdata))
all_rawdata <- cbind(all_rawdata[co_gene, ana_cells],
                     ys_rawdata[co_gene,])

all_rawdata <- all_rawdata[, rownames(all_anno)]

####
ana3 <- CreateSeuratObject(all_rawdata, meta.data = all_anno)
ys3.list <- SplitObject(ana3, split.by = "batch")

for (i in 1:length(ys3.list)) {
  ys3.list[[i]] <- NormalizeData(ys3.list[[i]], verbose = FALSE)
  ys3.list[[i]] <- FindVariableFeatures(ys3.list[[i]], selection.method = "vst", 
                                        nfeatures = 3000, verbose = FALSE)
}


reference.list <- ys3.list[c("ys_cs11", "ys_cs15",'ys_STRT')]
ys3.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
ys3.integrated <- IntegrateData(anchorset = ys3.anchors, dims = 1:30)

DefaultAssay(ys3.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
ys3.integrated <- ScaleData(ys3.integrated, verbose = FALSE)
ys3.integrated <- RunPCA(ys3.integrated, verbose = FALSE)
ElbowPlot(ys3.integrated)


ys3.integrated <- FindNeighbors(ys3.integrated, dims = 1:10, k.param = 100)
ys3.integrated <- FindClusters(ys3.integrated, resolution = seq(0.1,2,0.1))


library(umap)
umap_data <- umap(ys3.integrated@reductions$pca@cell.embeddings[,1:10])
umap_coor <- umap_data$layout
pca_data <- ys3.integrated@reductions$pca@cell.embeddings[]
ys3.integrated@reductions$pca@cell.embeddings[] <- as.matrix(umap_coor)

a <- DimPlot(ys3.integrated, group.by = 'batch', pt.size = 1)
b <- DimPlot(ys3.integrated, group.by = 'Stage', pt.size = 1)+scale_color_manual(values = col_stage)
library(ggpubr)
ggarrange(a,b, ncol=2)


a <- DimPlot(ys3.integrated, group.by = 'batch')
b <- DimPlot(ys3.integrated, group.by = 'integrated_snn_res.0.1', label = T)
cowplot::plot_grid(a,b)

ys3.integrated@meta.data$cluster <- as.character(ys3.integrated@meta.data$integrated_snn_res.0.1)
FeaturePlot(ys3.integrated, reduction = 'pca', features = c('EPCAM','CDH5'))
ys3.integrated@meta.data$cluster <- plyr::mapvalues(ys3.integrated@meta.data$cluster,
                                                    as.character(0:7),
                                                    c('EMP2','Mes1','Mac','Epi','ErP/MkP','Endo','Mes2','EMP1'))

cluster_order <- c('EMP1','EMP2','ErP/MkP','Mac','Endo','Epi','Mes1',"Mes2")
col_cluster <- setNames(c('red3','orange3','green3','purple3','blue3','magenta','#C16200','#71226e'),
                         cluster_order)




library(ggplot2)
col_batch <- c('ys_cs11' = 'yellow3', 'ys_cs15' = 'pink', 'ys_STRT'= 'magenta')

DimPlot(ys3.integrated, pt.size = 1, group.by = 'cluster', label = T, label.size = 8)+scale_color_manual(values = col_cluster)+
   theme_bw()+theme(legend.text = element_text(size=15), axis.text = element_text(size = 15),axis.title = element_text(size = 15))+labs(x='UMAP1',y="UMAP2")

DimPlot(ys3.integrated, pt.size = 1, group.by = 'Stage', label = F, label.size = 8)+scale_color_manual(values = col_stage)+
   theme_bw()+theme(legend.text = element_text(size=15), axis.text = element_text(size = 15),axis.title = element_text(size = 15))+labs(x='UMAP1',y="UMAP2")

DimPlot(ys3.integrated, pt.size = 1, group.by = 'batch', label = F, label.size = 8)+scale_color_manual(values = col_batch)+
   theme_bw()+theme(legend.text = element_text(size=15), axis.text = element_text(size = 15),axis.title = element_text(size = 15))+labs(x='UMAP1',y="UMAP2")
 





ys3.integrated <- SetIdent(ys3.integrated, value = 'cluster')
degs_all <- FindAllMarkers(ys3.integrated, logfc.threshold = log(1.5), only.pos = T)
degs_all <- degs_all[degs_all$p_val_adj<0.05,]
 
top5 <- c()
for (i in cluster_order){
  
  tmp <- degs_all$gene[as.character(degs_all$cluster) ==i][1:5]
  top5 <- c(top5, tmp)
}

ys3.integrated@meta.data$cluster <- factor(ys3.integrated@meta.data$cluster,
                                           levels = cluster_order)
DotPlot(ys3.integrated, features = unique(top5), group.by = 'cluster', cols = c('#E6E6E6', 'brown'))+coord_flip()+theme_bw()+
  theme(legend.text = element_text(size=15), axis.text = element_text(size = 15),axis.title = element_text(size = 15), axis.title.y = element_text(angle = 60))+labs(x = "", y="")

###plot feature genes matrix
source('/data3/dd/code/plot_gene_violin_matrix.R')
gene <- c("RUNX1","CD34", "MYB", "GATA1","PF4", "KLF1", "CD14", "CD163", "CDH5", "SOX7", "EPCAM", "PDGFRA")

anno <- ys3.integrated@meta.data[,'cluster',drop=F]
library(reshape2)
anno$cluster <- factor(anno$cluster, levels = cluster_order)
result <- plot_gene_matrix(data = as.matrix(ys3.integrated@assays$RNA@data), gene = gene, anno = anno, col = col_cluster)
result+theme_bw()+theme(axis.text = element_text(size=12),axis.text.x =  element_text(angle = 60, hjust = 1), axis.title = element_blank(), legend.text = element_text(size = 15))


###EMP markers, use rawdata
ys3.integrated@meta.data$Cluster <- ys3.integrated@meta.data$cluster
ys3.integrated@meta.data$Cluster[ys3.integrated@meta.data$Cluster %in% c('EMP1','EMP2')]<- 'EMP'
ys3.integrated <- SetIdent(ys3.integrated, value = 'Cluster')

ys3.integrated2 <- subset(ys3.integrated, idents = c('EMP','ErP/MkP','Mac'))
ys3.integrated2@meta.data$Cluster[!ys3.integrated2@meta.data$Cluster=='EMP']<- 'others'
ys3.integrated2 <- SetIdent(ys3.integrated2, value = 'Cluster')

deg_emp <- FindAllMarkers(ys3.integrated2, logfc.threshold = log(1.5), only.pos = T)
deg_emp <- deg_emp[deg_emp$p_val_adj<0.05,]

markers <- read.csv('/data3/dd/data/human_gene/human_cd_moculer.csv')
deg_emp <- deg_emp[rownames(deg_emp) %in% markers$symbols & deg_emp$cluster=='EMP',]
deg_emp$gene <- rownames(deg_emp)

ggplot(deg_emp)+geom_bar(stat = 'identity', aes(x = rev(1:nrow(deg_emp)), y= -log10(p_val_adj), fill = avg_logFC))+coord_flip()+
  scale_fill_gradientn(colours = colorRampPalette(c('#E6E6E6',' brown'))(100))+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(size=15))+
  geom_text(aes(x = rev(1:nrow(deg_emp)), y= -log10(p_val_adj)+15, label = gene))

getwd()
write.csv(file = 'cs11_cs15_result/emp_marker.csv', deg_emp)
write.csv(file = 'cs11_cs15_result/all_deg.csv', degs_all )
save(file = 'cs11_cs15_result/cs11_cs15_STRT.Rdata', list = ls())
####
emp_cell <- rownames(ys3.integrated@meta.data)[ys3.integrated@meta.data$Cluster=='EMP' ]
emp <- subset(ys3.integrated, cells = emp_cell)
emp_data <- as.matrix(emp@assays$RNA@counts)
emp_data <- apply(emp_data,2,function(x) log2(1e4*x/sum(x)+1))
emp_anno <- emp@meta.data
emp_anno$PTPRC <- emp_data['PTPRC',rownames(emp_anno)]

source('/data3/dd/code/get_smoothdata.R')
smooth_data <- get_smooth_data(data = emp_data[, rownames(emp_anno)[order(emp_anno$PTPRC, decreasing = F)]],k=50)  

library(pheatmap)
pheatmap(emp_data[c('KIT','ITGA2B','MYB','PTPRC','GATA1','PF4','CSF1R'), rownames(emp_anno)[order(emp_anno$PTPRC)]],cluster_cols = F, cluster_rows = F,show_colnames = F,color = colorRampPalette(c('#E6E6E6','darkorange','red','brown'))(100))

####
library(ggpubr)
emp_anno$CD45 <- ifelse(emp_anno$PTPRC>0, 'pos','neg')
emp_anno <- data.frame(cbind(emp_anno, t(emp_data[c('KIT','ITGA2B','MYB','PTPRC','GATA1','PF4','CSF1R'),rownames(emp_anno)])))

a <- ggboxplot(emp_anno, x='CD45',y='KIT',fill='CD45',add = 'jitter')+stat_compare_means(method = 't.test')
b <- ggboxplot(emp_anno, x='CD45',y='ITGA2B',fill='CD45',add = 'jitter')+stat_compare_means(method = 't.test')
c <- ggboxplot(emp_anno, x='CD45',y='FCGR3A',fill='CD45',add = 'jitter')+stat_compare_means(method = 't.test')
d <- ggboxplot(emp_anno, x='CD45',y='GATA1',fill='CD45',add = 'jitter')+stat_compare_means(method = 't.test')
e <- ggboxplot(emp_anno, x='CD45',y='PF4',fill='CD45',add = 'jitter')+stat_compare_means(method = 't.test')
f <- ggboxplot(emp_anno, x='CD45',y='CSF1R',fill='CD45',add = 'jitter')+stat_compare_means(method = 't.test')

ggarrange(a,b,c,d,e,f,nrow = 6,ncol=1, legend = 'none')
