rm(list=ls())
options(stringsAsFactors = F)


###cluster detection in all cells passing qc
load('human_macro_development_data.rds')

####
library(Seurat)
tmp <- CreateSeuratObject(rawdata, meta.data = anno)
tmp <- NormalizeData(tmp, scale.factor = 1e5)
tmp <- ScaleData(tmp,vars.to.regress = c('nUMI','nGene'))

#find out hvg
tmp <- FindVariableGenes(tmp, x.low.cutoff = 0.5, y.cutoff = 0.5)
tmp <- RunPCA(tmp, pcs.compute = 50, pc.genes = tmp@var.genes )
PCElbowPlot(tmp, num.pc = 50)


tmp <- RunUMAP (tmp,dims.use = 1:20)
a <- DimPlot(tmp, reduction.use = 'umap',group.by = 'Stage',do.return = T)+scale_color_manual(values = col_stage)
b <- DimPlot(tmp, reduction.use = 'umap',group.by = 'Site',do.return = T)+scale_color_manual(values = col_site)
plot_grid(a,b)
FeaturePlot(tmp, features.plot = c('PTPRC','EPCAM','KRT8','COL1A1'),reduction.use = 'umap',max.cutoff = 3, no.legend = F, cols.use = c('#E6E6E6','red'))

epcam_cell <- DimPlot(tmp, reduction.use = 'umap',group.by = 'Stage',do.return = T,do.identify = T)
mes_cell <- DimPlot(tmp, reduction.use = 'umap',group.by = 'Stage',do.return = T,do.identify = T)

####filter out epi and mes
anno <- anno[!rownames(anno) %in% c(epcam_cell,mes_cell),]
rawdata <- rawdata[,rownames(anno)]

####mitigate batch effects between cs20 and cs23

anno$batch <- anno$Stage
anno$batch[!anno$batch %in% c('CS20','CS23')]<-'others'
main <- CreateSeuratObject(rawdata, meta.data = anno)
main <- NormalizeData(main, scale.factor = 1e5)
main <- ScaleData(main, vars.to.regress = c('nUMI','nGene','batch'))

#find hvg
main <- FindVariableGenes(main, x.low.cutoff = 1, y.cutoff = 1)
main <- RunPCA(main, pcs.compute = 50, pc.genes = main@var.genes )
PCElbowPlot(main, num.pc = 50)

main <- RunUMAP (main,dims.use = 1:20,min_dist = 0.3,spread=0.5)
main <- FindClusters(main, dims.use = 1:20, resolution =seq(0.1,2,0.1),
                     force.recalc = T)


DimPlot(main, reduction.use = 'umap',group.by = 'res.1.5',do.return = T,label.size = 7,pt.size = 1.5,do.label = F)+scale_color_manual(values = c(mclust::mclust.options()$classPlotColors,'pink','green3'))
DimPlot(main, reduction.use = 'umap',group.by = 'Stage',do.return = T,label.size = 7,pt.size = 1.5,do.label = F)+scale_color_manual(values = col_stage)
FeaturePlot(main, reduction.use = 'umap',features.plot = c('HBB','KLF1','SPI1','PF4','CPA3','RUNX3'))

#find sub-clusters in ErP-MKP
ana_cell <- rownames(main@meta.data)[main@meta.data$res.1.5=='2']
ana <- SubsetData(main, cells.use = ana_cell )
ana <- FindVariableGenes(ana, x.low.cutoff = 1, y.cutoff = 1)
ana <- RunPCA(ana, pc.genes = ana@var.genes, pcs.compute = 50)
PCElbowPlot(ana, num.pc = 50)

ana <- FindClusters(ana, dims.use = 1:10, resolution = seq(0.2,2,0.2))
DimPlot(ana, reduction.use = 'pca',group.by = 'res.0.8')
FeaturePlot(ana, reduction.use = 'pca', features.plot = c('PF4','KLF1','CPA3','GATA1'))

ana_anno <- ana@meta.data[,'res.0.8',drop=F]
names(ana_anno) <- 'cluster'
ana_anno$cluster <- paste('c',ana_anno$cluster,sep="")

main@meta.data$cluster <- main@meta.data$res.1.5
main@meta.data[rownames(main@meta.data) %in% rownames(ana_anno),]$cluster <-as.character(na.omit(ana_anno$cluster[match(rownames(main@meta.data),rownames(ana_anno))]))

DimPlot(main, reduction.use = 'umap',group.by = 'cluster',do.return = T,label.size = 7,pt.size = 1.5,do.label = F)+scale_color_manual(values = c(mclust::mclust.options()$classPlotColors,'pink','green3'))

#redefine the cluster
main@meta.data$cluster <- plyr::mapvalues(main@meta.data$cluster,
                                          unique(main@meta.data$cluster),
                                          c('Mac_2','MkP','ErP','GMP','Monocyte','Lymphoblast','LMP','EMP','Mast cell','Myeloblast','Mac_4','Mac_1','ILC','Mac_3','HSPC'))

cluster_order <- c('EMP','ErP','MkP','GMP','Myeloblast','Monocyte','Mac_1','Mac_2','Mac_3','Mac_4','LMP','Lymphoblast','ILC','Mast cell','HSPC')
all_anno <- main@meta.data[,c('Site','Stage','cluster')]
all_anno <- data.frame(all_anno, main@dr$umap@cell.embeddings[])
all_anno$Site[all_anno$Site=="FL"]<-  "Liver"
names(col_site)[3] <- 'Liver'
all_anno$Stage[all_anno$Stage=='CS11-early']<-'CS11'
names(col_stage)[3] <- 'CS11'


library(ggsci)
col_all <- pal_d3(palette = 'category20')(16)[-8]
names(col_all) <- sample(cluster_order,15)
p4 <- ggplot(all_anno,aes(x=UMAP1,y=UMAP2,col=factor(cluster,levels = cluster_order)))+geom_point()+scale_color_manual("",values = col_all)
p5 <- ggplot(all_anno,aes(x=UMAP1,y=UMAP2,col=Stage))+geom_point()+scale_color_manual("",values = col_stage)+theme(legend.position = c(0.8,0.2))
p6 <- ggplot(all_anno,aes(x=UMAP1,y=UMAP2,col=Site))+geom_point()+scale_color_manual("",values = col_site)+theme(legend.position = c(0.8,0.2))

######progenitor comparasion

progenitor_cell <- rownames(all_anno)[all_anno$cluster %in% c('EMP','ErP','MkP','GMP','Myeloblast','LMP','Lymphoblast','HSPC')]
progenitor_rawdata <- rawdata[,progenitor_cell]
progenitor_anno <- all_anno[progenitor_cell,]

pro <- CreateSeuratObject(progenitor_rawdata, meta.data = progenitor_anno)
pro <- NormalizeData(pro, scale.factor = 1e5)
pro <- ScaleData(pro)
pro <- SetAllIdent(pro, id='cluster')

deg_pro <- FindAllMarkers(pro, logfc.threshold = 0.5, only.pos = T, test.use = 'wilcox')
deg_pro <- deg_pro[deg_pro$p_val_adj<0.05, ]

top5 <- c()
for(i in c('EMP','ErP','MkP','GMP','Myeloblast','LMP','Lymphoblast','HSPC')){
  tmp_gene <- deg_pro$gene[deg_pro$cluster==i][1:5]
  top5 <- c(top5, tmp_gene)
}


gene_data <- data.frame(t(as.matrix(pro@data)),cluster=progenitor_anno$cluster,check.names = F)
average_data <- aggregate(.~cluster, gene_data, mean)
cluster_name <- average_data[,1]
average_data <- apply(average_data[,2:ncol(average_data)],2,as.numeric)
rownames(average_data) <- cluster_name
average_data <- t(average_data)
phmat1 <- t(scale(t(average_data)))
phmat1[phmat1> 1.5] <- 1.5
phmat1[phmat1 < -1.5] <- -1.5


#degs
library(pheatmap)
pheatmap(phmat1[top5,c('EMP','ErP','MkP','GMP','Myeloblast','LMP','Lymphoblast','HSPC')],cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("#440154" ,"#21908C", "#FDE725"))(100))


#scenic of progenitor
scenic <-  readRDS('F:/analysis/project/human_macro_development_project/macro_result_version3/4.2_binaryRegulonActivity_nonDupl.Rds')
cor_data <- cor(t(scenic),method = 'spearman')
filter_regulon <- rownames(cor_data)[rowSums(cor_data>0.3)>1]

#filter: activate 30% in at least one cluster
scenic_data <- data.frame(t(scenic[filter_regulon,rownames(progenitor_anno)]),check.names = F)
scenic_data$cluster <- progenitor_anno$cluster
scenic_aggr <- aggregate(.~cluster,scenic_data,sum)

cell_num <- table(progenitor_anno$cluster)
for (i in 1:nrow(scenic_aggr)){
  cl <- scenic_aggr$cluster[i]
  num <- cell_num[cl]
  
  scenic_aggr[i,2:ncol(scenic_aggr)]<- scenic_aggr[i,2:ncol(scenic_aggr)]/num
  
}

rownames(scenic_aggr) <- scenic_aggr$cluster
scenic_aggr <- scenic_aggr[,2:ncol(scenic_aggr)]
scenic_aggr <- scenic_aggr[,colSums(scenic_aggr>0.3)>=1]

scenic_filter <- scenic[colnames(scenic_aggr),]

cell<- c()
for(i in c('EMP','ErP','MkP','GMP','Myeloblast','LMP','Lymphoblast','HSPC')){
  tmp <- rownames(progenitor_anno)[progenitor_anno$cluster==i]
  cell <- c(cell, tmp)
}

#
ph <- pheatmap::pheatmap(scenic_filter[,cell],show_rownames = T,show_colnames = F,annotation_col = progenitor_anno[,c('cluster','Site'),drop=F],cluster_cols = F,
                   annotation_colors = list(cluster=col_all[c('EMP','ErP','MkP','GMP','Myeloblast','LMP','Lymphoblast','HSPC')],Site=col_site),clustering_method = 'ward.D2',cutree_rows = 7)


###change and difference between HSPC and EMP

hpc <- rownames(all_anno)[all_anno$cluster %in% c('HSPC','EMP')]
hpc_rawdata <- rawdata[,hpc]
hpc_anno <- all_anno[hpc,]

hpc <- CreateSeuratObject(hpc_rawdata, meta.data = hpc_anno)
hpc <- NormalizeData(hpc, scale.factor = 1e5)
hpc <- ScaleData(hpc)

hpc <- SetAllIdent(hpc, id='cluster')
deg_hpc <- FindAllMarkers(hpc, logfc.threshold = 0, only.pos = T, test.use = 'wilcox')
deg_hpc_filter <- deg_hpc[deg_hpc$p_val_adj<0.05 & deg_hpc$avg_logFC>0.25,]

deg_hpc$avg_logFC[deg_hpc$cluster=='HSPC'] <- 0-deg_hpc$avg_logFC[deg_hpc$cluster=='HSPC']


#show top10 genes on vocalno plots

library(dplyr)
deg_hpc_filter_order <- deg_hpc_filter %>% group_by(cluster) %>% arrange(-avg_logFC)
top10 <- c(deg_hpc_filter_order$gene[deg_hpc_filter_order$cluster=='EMP'][1:10],
           deg_hpc_filter_order$gene[deg_hpc_filter_order$cluster=='HSPC'][1:10],'HOXA6','HOXA10')

deg_hpc_filter_order$avg_logFC[deg_hpc_filter_order$cluster=='HSPC'] <- 0-deg_hpc_filter_order$avg_logFC[deg_hpc_filter_order$cluster=='HSPC']


deg_hpc_filter_order$label <- deg_hpc_filter_order$gene
deg_hpc_filter_order$label[!deg_hpc_filter_order$label %in% top10] <- NA
library(ggrepel)

deg_hpc$deg <- ifelse(deg_hpc$gene %in% deg_hpc_filter$gene,'DEG','non-DEG')
ggplot()+geom_point(data=deg_hpc, aes(x=avg_logFC,y=-log10(p_val),col=deg))+scale_color_manual(values = c('red','grey'))+
  geom_text_repel(data=as.data.frame(deg_hpc_filter_order),aes(x=avg_logFC,y=-log10(p_val),label=label))

###change of proportion
progenitor_proportion <- table(progenitor_anno$Stage, progenitor_anno$cluster)
progenitor_proportion <- data.frame(apply(progenitor_proportion,1,function(x) x/sum(x)))
progenitor_proportion$cluster <- rownames(progenitor_proportion)
progenitor_proportion_melt <- reshape::melt(progenitor_proportion, id='cluster')
names(progenitor_proportion_melt) <- c('cluster','stage','proportion')

ggplot(progenitor_proportion_melt[progenitor_proportion_melt$cluster %in% c('EMP','HSPC'),], aes(x=stage, y=proportion, col=cluster, group=cluster, shape=cluster))+geom_line(size=1.5)+
  scale_color_manual(values = col_all[c('EMP','HSPC')])+geom_point(size=3,col='#808080')+theme_light()

####Fig3 specification

rm(list=ls())
load('F:/analysis/project/human_macro_development_project/mcrophage_result_version5/Fig1/fig1_v5.rds')
cell.use <- rownames(all_anno)[all_anno$Site=='Liver' & (!all_anno$Stage %in% c('CS20','CS23')) & (!all_anno$cluster %in% c('ErP','MkP','Mast cell','LMP','Lymphoblast','Mac_1','Mac_2','Mac_3','Mac_4'))]
ana_rawdata <- rawdata[,cell.use]
ana_anno <- all_anno[cell.use,]

library(Seurat)
mono <- CreateSeuratObject(ana_rawdata, meta.data = ana_anno)
mono <- NormalizeData(mono, scale.factor = 1e5)
mono <- ScaleData(mono, vars.to.regress = c('nGene','nUMI'))
mono <- FindVariableGenes(mono, x.low.cutoff = 0.5, y.cutoff = 0.5)

#
mono <- RunPCA(mono, pc.genes = mono@var.genes,pcs.compute = 50)
PCElbowPlot(mono, num.pc = 50)

#DimPlot(mono, reduction.use = 'pca',group.by='Site',dim.1 = 1,dim.2 = 3)
a <- DimPlot(mono, reduction.use = 'pca',group.by='Stage',dim.1 = 1,dim.2 = 2,pt.size = 3,do.return = T)+scale_color_manual(values = col_stage)
b <- DimPlot(mono, reduction.use = 'pca',group.by='cluster',dim.1 = 1,dim.2 = 2,do.return = T,pt.size = 4)+scale_fill_manual(values = col_all)


FeaturePlot(mono, reduction.use = 'pca', features.plot = 'CEACAM8')

#filter out outlier(one cell)
cell_select <- DimPlot(mono, reduction.use = 'pca',group.by='cluster',dim.1 = 1,dim.2 = 3,do.return = T,pt.size = 4,do.identify = T)

#
pca_all<-FactoMineR::PCA(t(as.matrix(mono@data)[mono@var.genes,]),graph = F)

library(factoextra)
fviz_pca_ind(pca_all,
             label = "none", # hide individual labels
             habillage = factor(mono@meta.data$cluster),
             palette = col_all,
             addEllipses = F)+labs(main=NULL)

fviz_pca_ind(pca_all,axes = c(1,3),
             label = "none", # hide individual labels
             habillage = factor(mono@meta.data$cluster),
             palette = col_all,
             addEllipses = F
)+labs(main=NULL)



#mono@dr$pca@cell.embeddings[,1] <- 0-mono@dr$pca@cell.embeddings[,1]
mono_anno <- data.frame(mono@meta.data, mono@dr$pca@cell.embeddings[])[cell_select,]
ggplot()+geom_point(data=mono_anno, aes(x=PC1, y=PC3, fill=cluster),pch=21, col=scales::alpha('grey',alpha = 0.2),size=3)+
  scale_fill_manual(values = col_all)

ggplot()+geom_point(data=mono_anno, aes(x=PC1, y=PC2, fill=cluster),pch=21, col=scales::alpha('grey',alpha = 0.2),size=3)+
  scale_fill_manual(values = col_all)

col_map<-col_all[match(mono_anno$cluster,names(col_all))]
pairs(~PC1+PC2+PC3,data=mono_anno,main="PC plot matrix",bg=col_map,pch=21,cex=1.8)

#important gene
mono <- SubsetData(mono, cells.use = cell_select)
FeaturePlot(mono, reduction.use = 'pca', features.plot = c('CD34','MYB','MPO','CEACAM8','CCR2','HLA-DRA'),nCol = 3,dim.1 = 1,dim.2 = 3,cols.use = c('grey','red'),pt.size = 1.5)
FeaturePlot(mono, reduction.use = 'pca', features.plot = c('CD34','MYB','MPO','CEACAM8','CCR2','HLA-DRA'),nCol = 3,dim.1 = 1,dim.2 = 3,cols.use = c('grey','red'),pt.size = 1.5, no.legend = F)

#do monocle
do_monocle<-function(cell_cycle,genelow,rawdata,anno,col,gene=NULL,disper,reverse=NULL,reduce_brach=T,size,define_root=F,root_state,only_gamma,regress.batch=F,gamma.use){
  
  count_matrix <- as.matrix(rawdata)[,rownames(anno)]
  pd <- new("AnnotatedDataFrame", data =anno)
  HSMM <- newCellDataSet(count_matrix, phenoData = pd,expressionFamily=negbinomial())
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)
  #
  disp_table <- dispersionTable(HSMM)
  if(is.null(gene)){
    ordering_genes <- subset(disp_table,
                             mean_expression >= genelow &
                               dispersion_empirical >= disper* dispersion_fit)$gene_id
    
    ordering_genes <- setdiff(ordering_genes, cell_cycle)
  }else{
    ordering_genes <- gene
  }
  
  HSMM<- setOrderingFilter(HSMM, ordering_genes)
  plot_ordering_genes(HSMM)
  if(regress.batch==T){
    if(reduce_brach){
      if(only_gamma){
        HSMM <- reduceDimension(HSMM, max_components=2,ncenter=nrow(anno)-1, param.gamma=100,residualModelFormulaStr='~batch')
      }else{
        HSMM <- reduceDimension(HSMM, max_components=2, param.gamma=100,residualModelFormulaStr='~batch')
      }
    }else{
      HSMM <- reduceDimension(HSMM, max_components=2,residualModelFormulaStr='~batch')
    }}else{
      if(reduce_brach){
        if(only_gamma){
          HSMM <- reduceDimension(HSMM, max_components=2,ncenter=nrow(anno)-1, param.gamma=gamma.use)
        }else{
          HSMM <- reduceDimension(HSMM, max_components=2, param.gamma=gamma.use)
        }
      }else{
        HSMM <- reduceDimension(HSMM, max_components=2)
      } 
    }
  
  
  HSMM <- orderCells(HSMM, reverse=reverse)
  
  if(define_root){
    HSMM <- HSMM <- orderCells(HSMM,root_state = root_state)
  }
  
  
  print(plot_cell_trajectory(HSMM,color_by = "cluster",cell_size = size)+scale_color_manual(values=col))
  return(HSMM)
  
  
  
}

library(monocle)
monocle1 <- do_monocle(cell_cycle = NULL, rawdata = ana_rawdata[,rownames(mono_anno)],anno = data.frame(ana_anno[rownames(mono_anno),],t(as.matrix(mono@data))),col = col_all,size = 4,reduce_brach = T,only_gamma = T,gene = mono@var.genes,gamma.use = 100)

plot_cell_trajectory(monocle1,show_branch_points = T,color_by = "State",cell_size = 2)
a <- plot_cell_trajectory(monocle1,show_branch_points = F,color_by = "Pseudotime",cell_size = 3)+scale_color_gradientn(colours = colorRampPalette(c('blue','green','yellow','red'))(100))
b <- plot_cell_trajectory(monocle1,show_branch_points = F,color_by = "cluster",cell_size = 3)+scale_color_manual(values = col_all)
c <- plot_cell_trajectory(monocle1,show_branch_points = F,color_by = "Stage",cell_size = 3)+scale_color_manual(values = col_stage)

plot_grid(a,b)


p1 <- plot_cell_trajectory(monocle1,show_branch_points = F,color_by = "CD34",cell_size = 2)+scale_color_gradientn(colours = colorRampPalette(c('blue','green','yellow','red'))(100))
p2 <- plot_cell_trajectory(monocle1,show_branch_points = F,color_by = "MYB",cell_size = 2)+scale_color_gradientn(colours = colorRampPalette(c('blue','green','yellow','red'))(100))
p3 <- plot_cell_trajectory(monocle1,show_branch_points = F,color_by = "S100A8",cell_size = 2)+scale_color_gradientn(colours = colorRampPalette(c('blue','green','yellow','red'))(100))
p4 <- plot_cell_trajectory(monocle1,show_branch_points = F,color_by = "CEACAM8",cell_size = 2)+scale_color_gradientn(colours = colorRampPalette(c('blue','green','yellow','red'))(100))
p5 <- plot_cell_trajectory(monocle1,show_branch_points = F,color_by = "CCR2",cell_size = 2)+scale_color_gradientn(colours = colorRampPalette(c('blue','green','yellow','red'))(100))
p6 <- plot_cell_trajectory(monocle1,show_branch_points = F,color_by = "MPO",cell_size = 2)+scale_color_gradientn(colours = colorRampPalette(c('blue','green','yellow','red'))(100))
p7 <- plot_cell_trajectory(monocle1,show_branch_points = F,color_by = "HLA.DRA",cell_size = 2)+scale_color_gradientn(colours = colorRampPalette(c('blue','green','yellow','red'))(100))

plot_grid(p1,p2,p6,p3,p4,p7,nrow = 2,ncol = 3)

###########Fig4 
all_anno$select <- ifelse(all_anno$cluster %in% c('ErP','MkP','Mast cell','LMP','Lymphoblast','ILC','HSPC'),'others','select')
ggplot(all_anno, aes(x=UMAP1, y=UMAP2, col=select))+geom_point()+scale_color_manual(values = c('grey','darkred'))+theme_classic()


ana_cell <- rownames(all_anno)[!all_anno$cluster %in% c('ErP','MkP','Mast cell','LMP','Lymphoblast','ILC','HSPC')]
ana_rawdata <- rawdata[,ana_cell]
ana_anno <- all_anno[ana_cell,]

ana_anno$batch <- ana_anno$Stage
ana_anno$batch[!ana_anno$batch %in% c('CS20','CS23')]<-'others'

library(Seurat)
ana <- CreateSeuratObject(ana_rawdata, meta.data = ana_anno)
ana <- NormalizeData(ana, scale.factor = 1e5)
ana <- ScaleData(ana, vars.to.regress = c('nUMI','nGene','batch'))


#### find hvgs
ana <- FindVariableGenes(ana, x.low.cutoff = 0.3, y.cutoff = 0.75)
ana <- RunPCA(ana, pcs.compute = 50, pc.genes = ana@var.genes)
PCElbowPlot(ana, num.pc = 50)

ana <- FindClusters(ana, dims.use = 1:30, resolution =seq(0.1,2,0.1),
                    force.recalc = T,k.param = 10)

ana <- RunUMAP (ana,dims.use = 1:30,n_neighbors = 40,min_dist = 0.3, spread=0.5)

a <- DimPlot(ana, reduction.use = 'umap',group.by = 'res.1.6',do.return = T,label.size = 7,pt.size = 1.5,do.label = T)+scale_color_manual(values = c(mclust::mclust.options()$classPlotColors,'pink','green3'))
b <- DimPlot(ana, reduction.use = 'umap',group.by = 'Site',do.return = T,label.size = 7,pt.size = 1.5,do.label = T)+scale_color_manual(values = col_site)
c <- DimPlot(ana, reduction.use = 'umap',group.by = 'Stage',do.return = T,label.size = 7,pt.size = 1.5,do.label = T)+scale_color_manual(values = col_stage)
plot_grid(a,b,c,ncol = 3)

# ana@dr$umap@cell.embeddings[,1] <- 0-ana@dr$umap@cell.embeddings[,1]
ana@dr$umap@cell.embeddings[,2] <- 0-ana@dr$umap@cell.embeddings[,2]

###
ana@meta.data$cluster <- ana@meta.data$res.1.6
ana@meta.data$cluster <- plyr::mapvalues(ana@meta.data$cluster, 
                                         as.character(0:15),
                                         c('Myeloblast','EMP','Head_Mac1','Lung_Mac',"Monocyte",
                                           "YS_Mac1",'Blood_Mac','Skin_Mac','Head_Mac3','Head_Mac2',
                                           "EMP",'GMP','Head_Mac4','YS_Mac2','EMP','Liver_Mac'))


col_mac <- setNames(c("#98DF8AFF","#FFBB78FF","darkorange","darkgreen", 
                      "purple3","red3","pink","cyan","#1F77B4FF", 
                      "#E377C2FF","green","#8C564BFF","gold",'magenta'),
                    c('EMP',"GMP","Myeloblast","Monocyte",
                      "Blood_Mac","Liver_Mac","Skin_Mac","Lung_Mac","YS_Mac1",
                      "YS_Mac2","Head_Mac1","Head_Mac2","Head_Mac3",'Head_Mac4'))



cluster_order <- c('EMP','GMP',"Myeloblast","Monocyte",
                   "Blood_Mac","Liver_Mac","Skin_Mac","Lung_Mac","YS_Mac1",
                   "YS_Mac2","Head_Mac1","Head_Mac2","Head_Mac3",'Head_Mac4')


###
mac_anno <- data.frame(ana@meta.data[,c('Stage','Site','cluster')], ana@dr$umap@cell.embeddings[,1:2])

ggplot(mac_anno, aes(x=UMAP1, y=UMAP2, col=factor(cluster, levels = cluster_order)))+geom_point(size=1.5)+
  scale_color_manual("",values = col_mac)

ggplot(mac_anno, aes(x=UMAP1, y=UMAP2, col=Stage))+geom_point()+
  scale_color_manual("",values = col_stage)

ggplot(mac_anno, aes(x=UMAP1, y=UMAP2, col=Site))+geom_point()+
  scale_color_manual("",values = col_site)


#macrophage hetergenity and similarity and important gene
mac_cell <- rownames(mac_anno)[mac_anno$cluster %in% c('Liver_Mac','Blood_Mac','Lung_Mac','Skin_Mac','YS_Mac1','YS_Mac2','Head_Mac1','Head_Mac2',
                                                       "Head_Mac3",'Head_Mac4')]

mac_rawdata <- rawdata[,mac_cell]
mac_annot <- mac_anno[mac_cell,]

mac <- CreateSeuratObject(mac_rawdata, meta.data = mac_annot)
mac <- NormalizeData(mac, scale.factor = 1e5)
mac <- ScaleData(mac)

mac <- SetAllIdent(mac, id='cluster')
deg_mac <- FindAllMarkers(mac, only.pos = T, test.use = 'wilcox', logfc.threshold = 0.25)
deg_mac <- deg_mac[deg_mac$p_val_adj<0.05, ]

deg_mac <- filter_deg(deg_data = deg_mac, average_data = average_data)

top5_mac <- c()
for(i in c('Liver_Mac','Blood_Mac','Lung_Mac','Skin_Mac','YS_Mac1','YS_Mac2','Head_Mac1','Head_Mac2',"Head_Mac3",'Head_Mac4')){
  genetmp <- deg_mac$gene[deg_mac$cluster==i][1:5]
  top5_mac <- c(top5_mac, genetmp)
}

pheatmap(phmat2[as.character(na.omit(top5_mac)) ,c('Liver_Mac','Blood_Mac','Lung_Mac','Skin_Mac','YS_Mac1','YS_Mac2','Head_Mac1','Head_Mac2',"Head_Mac3",'Head_Mac4')], cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c("#440154" ,"#21908C", "#FDE725"))(100))


mac <- FindVariableGenes(mac, x.low.cutoff = 0.5, y.cutoff = 0.5)
hc <- hclust(dist(t(average_data[mac@var.genes, c('Liver_Mac','Blood_Mac','Lung_Mac','Skin_Mac','YS_Mac1','YS_Mac2','Head_Mac1','Head_Mac2',"Head_Mac3",'Head_Mac4')])),method = 'ward.D2')
dend3 <- as.dendrogram(hc) 

library(dendextend)
labels_colors(dend3) <- col_mac[labels(dend3)]
plot(dend3, horiz=T)
circlize_dendrogram(dend3, dend_track_height = 0.7)

###important genes of tissue-resided macrphage in science
head_gene <- c('SALL1','SALL3','MLXIPL','ZEB1','ETV1','MEIS3')
liver_gene <- c('ID3','SPIC','NR1H3','IRF7','ID1','NFE2','TCF12')
skin_gene <- c('ZBTB42','RUNX3','AHR')
lung_gene <- c('PPARG','KLF4','ATF5','CEBPB')

FeaturePlot(ana , reduction.use = 'umap', features.plot = head_gene)
FeaturePlot(ana , reduction.use = 'umap', features.plot = liver_gene)
FeaturePlot(ana , reduction.use = 'umap', features.plot = skin_gene)
FeaturePlot(ana , reduction.use = 'umap', features.plot = lung_gene)

gene_anno2 <- data.frame(gene=c(
                           rep('Liver',length(liver_gene)),
                           rep('Lung',length(lung_gene)),
                           rep('Skin', length(skin_gene)),
                           rep('Head',length(head_gene))),
                         row.names = c(liver_gene, lung_gene,skin_gene,head_gene))


pheatmap(phmat2[rownames(gene_anno2) ,c('Liver_Mac','Blood_Mac','Lung_Mac','Skin_Mac','YS_Mac1','YS_Mac2','Head_Mac1','Head_Mac2',"Head_Mac3",'Head_Mac4')], cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c("#1CB4E1" ,"black",'yellow'))(100),annotation_row =  gene_anno2)


#####head/liver/skin macrophage spefication 
rm( list = ls() )
options(stringsAsFactors = F)

load('F:/analysis/project/human_macro_development_project/mcrophage_result_version5/Fig4/fig4_v5.rds')
data <- apply(rawdata, 2, function(x) log2(1e5*x/sum(x)+1))

head_macro <- rownames(mac_anno)[as.character(mac_anno$Site)=='Head' &  (mac_anno$cluster %in% c('YS_Mac1','Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4'))]
mac_anno$head <- mac_anno$Stage
mac_anno$head[!rownames(mac_anno) %in% head_macro] <- 'others'
library(princurve)

head_dm <- as.matrix(mac_anno[head_macro,c('UMAP1','UMAP2')])
fit <- principal_curve(head_dm)
plot(fit,xlim=c(min(mac_anno$UMAP1),max(mac_anno$UMAP1)),
     ylim=c(min(mac_anno$UMAP2),max(mac_anno$UMAP2)),lwd=3)

col.use1 <- c(col_stage,'others'= scales::alpha('grey', alpha = 0.3))
points(as.matrix(mac_anno[,c('UMAP1','UMAP2')]),col=col.use1[match(mac_anno$head,names(col.use1))],pch=20,cex=1.5)
whiskers(head_dm, fit$s, col = "#808080",lwd=1)
legend('bottomleft', c('CS11','CS12','CS13','CS15','CS17','CS20','CS23'), fill = col_stage[c('CS11','CS12','CS13','CS15','CS17','CS20','CS23')])

mac_anno$Cluster <- mac_anno$cluster
mac_anno$Cluster[!rownames(mac_anno) %in% head_macro]<- 'others'
mac_anno$Cluster[mac_anno$Cluster=='YS_Mac1']<- "Head_Mac0"

plot(fit,xlim=c(min(mac_anno$UMAP1),max(mac_anno$UMAP1)),
     ylim=c(min(mac_anno$UMAP2),max(mac_anno$UMAP2)),lwd=3)

col.use2 <- c(col_mac,'Head_Mac0'='purple3','others'=scales::alpha('grey',alpha = 0.3))

points(as.matrix(mac_anno[,c('UMAP1','UMAP2')]),col=col.use2[match(mac_anno$Cluster,names(col.use2))],pch=20,cex=1.5)
whiskers(head_dm, fit$s, col = "#808080",lwd=1)
legend('bottomleft', c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4','others'), fill = col.use2[c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4','others')])

#gene to plot
head_anno <- mac_anno[mac_anno$Cluster %in% c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4'),'Cluster',drop=F]
names(head_anno) <- 'cluster'


head_gene <- c('CX3CR1','C1QB','SALL1','SPP1','C3','TMEM119','LYVE1', 'P2RY12')
gene_data <- data.frame(head_anno, t(as.matrix(data)[head_gene, rownames(head_anno)]))

library(ggplot2)

gene_data <- data.frame(cluster = mac_anno[head_macro, 'Cluster'], 
                         t(data[rowMeans(data[,head_macro])>0, head_macro]), check.names = F)


aggr_data <- aggregate(.~cluster, data = gene_data, FUN = mean)
average_data <- t(apply(aggr_data[,2:ncol(aggr_data)], 2, function(x) as.numeric(x)))
colnames(average_data) <- aggr_data[,1]

p1 <- ggplot()+geom_violin(data=gene_data, aes(x=cluster, y=SALL1, fill=cluster),scale = 'width')+
  geom_point(data=gene_data,aes(x=cluster, y=SALL1),size=0.5)+scale_fill_manual(values = col.use2[c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4')])+
  geom_segment(aes(x=1:4, y=as.numeric(average_data['SALL1',c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3')]),
                   xend=2:5, yend=as.numeric(average_data['SALL1',c('Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4')])),size=1)+
  theme(axis.text.x = element_blank())+labs(x="")

p2 <- ggplot()+geom_violin(data=gene_data, aes(x=cluster, y=CX3CR1, fill=cluster),scale = 'width')+
  geom_point(data=gene_data,aes(x=cluster, y=CX3CR1),size=0.5)+scale_fill_manual(values = col.use2[c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4')])+
  geom_segment(aes(x=1:4, y=as.numeric(average_data['CX3CR1',c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3')]),
                   xend=2:5, yend=as.numeric(average_data['CX3CR1',c('Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4')])),size=1)+
  theme(axis.text.x = element_blank())+labs(x="")


p3 <- ggplot()+geom_violin(data=gene_data, aes(x=cluster, y=C1QB, fill=cluster),scale = 'width')+
  geom_point(data=gene_data,aes(x=cluster, y=C1QB),size=0.5)+scale_fill_manual(values = col.use2[c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4')])+
  geom_segment(aes(x=1:4, y=as.numeric(average_data['C1QB',c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3')]),
                   xend=2:5, yend=as.numeric(average_data['C1QB',c('Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4')])),size=1)+
  theme(axis.text.x = element_blank())+labs(x="")

p4 <- ggplot()+geom_violin(data=gene_data, aes(x=cluster, y=SPP1, fill=cluster),scale = 'width')+
  geom_point(data=gene_data,aes(x=cluster, y=SPP1),size=0.5)+scale_fill_manual(values = col.use2[c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4')])+
  geom_segment(aes(x=1:4, y=as.numeric(average_data['SPP1',c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3')]),
                   xend=2:5, yend=as.numeric(average_data['SPP1',c('Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4')])),size=1)+
  theme(axis.text.x = element_blank())+labs(x="")

p5 <- ggplot()+geom_violin(data=gene_data, aes(x=cluster, y=C3, fill=cluster),scale = 'width')+
  geom_point(data=gene_data,aes(x=cluster, y=C3),size=0.5)+scale_fill_manual(values = col.use2[c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4')])+
  geom_segment(aes(x=1:4, y=as.numeric(average_data['C3',c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3')]),
                   xend=2:5, yend=as.numeric(average_data['C3',c('Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4')])),size=1)+
  theme(axis.text.x = element_blank())+labs(x="")

p6 <- ggplot()+geom_violin(data=gene_data, aes(x=cluster, y=LYVE1, fill=cluster),scale = 'width')+
  geom_point(data=gene_data,aes(x=cluster, y=LYVE1),size=0.5)+scale_fill_manual(values = col.use2[c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4')])+
  geom_segment(aes(x=1:4, y=as.numeric(average_data['LYVE1',c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3')]),
                   xend=2:5, yend=as.numeric(average_data['LYVE1',c('Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4')])),size=1)+
  theme(axis.text.x = element_blank())+labs(x="")

p7 <- ggplot()+geom_violin(data=gene_data, aes(x=cluster, y=P2RY12, fill=cluster),scale = 'width')+
  geom_point(data=gene_data,aes(x=cluster, y=P2RY12),size=0.5)+scale_fill_manual(values = col.use2[c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4')])+
  geom_segment(aes(x=1:4, y=as.numeric(average_data['P2RY12',c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3')]),
                   xend=2:5, yend=as.numeric(average_data['P2RY12',c('Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4')])),size=1)+
  theme(axis.text.x = element_blank())+labs(x="")

p8 <- ggplot()+geom_violin(data=gene_data, aes(x=cluster, y=TMEM119, fill=cluster),scale = 'width')+
  geom_point(data=gene_data,aes(x=cluster, y=TMEM119),size=0.5)+scale_fill_manual(values = col.use2[c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4')])+
  geom_segment(aes(x=1:4, y=as.numeric(average_data['TMEM119',c('Head_Mac0','Head_Mac1','Head_Mac2','Head_Mac3')]),
                   xend=2:5, yend=as.numeric(average_data['TMEM119',c('Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4')])),size=1)+
  theme(axis.text.x = element_blank())+labs(x="")

ggpubr::ggarrange(p2,p3,p4,p1,p5,p7,ncol = 2,nrow = 3, common.legend = T, legend = 'right')

####plot gene against stage
hp_genes <- read.csv('F:/data/human_gene/housekeeping_gene.csv')

mac_anno$nGene <- colSums(rawdata[, rownames(mac_anno)] > 0)
mac_anno$nUMI <- colSums(rawdata[, rownames(mac_anno)])
mac_anno$hp <- colMeans(data[intersect(hp_genes$gene, rownames(data)), rownames(mac_anno)])


mpath_data <- data[, head_macro]

gene_data2 <- data.frame(mac_anno[head_macro, c('Stage','nGene','nUMI','hp')], 
                        t(mpath_data[rowMeans(mpath_data)>0, head_macro]), check.names = F)


aggr_data2 <- aggregate(.~Stage, data = gene_data2, FUN = mean)
average_data2 <- t(apply(aggr_data2[,2:ncol(aggr_data2)], 2, function(x) as.numeric(x)))
colnames(average_data2) <- aggr_data2[,1]

head_gene2 <- c('GAPDH','ACTB','CD163','CX3CR1','C3','SALL1')
gene_data2 <- data.frame(mac_anno[rownames(head_anno), ], t(as.matrix(data)[head_gene2, rownames(head_anno)]))
gene_data2$nGene <- mac_anno[rownames(gene_data2) ,'nGene']
gene_data2$nUMI <- mac_anno[rownames(gene_data2) ,'nUMI']
gene_data2$hp <- mac_anno[rownames(gene_data2) ,'hp']

library(ggplot2)

a <- ggplot()+geom_violin(data=gene_data2, aes(x=Stage, y=nGene, fill=Stage),scale = 'width')+
  geom_point(data=gene_data2,aes(x=Stage, y=nGene),size=0.5)+scale_fill_manual(values = col_stage)+
  geom_segment(aes(x=1:6, y=as.numeric(average_data2['nGene',1:6]),
                   xend=2:7, yend=as.numeric(average_data2['nGene',2:7])),size=1)+
  labs(x="")



b <- ggplot()+geom_violin(data=gene_data2, aes(x=Stage, y=hp, fill=Stage),scale = 'width')+
  geom_point(data=gene_data2,aes(x=Stage, y=hp),size=0.5)+scale_fill_manual(values = col_stage)+
  geom_segment(aes(x=1:6, y=as.numeric(average_data2['hp',1:6]),
                   xend=2:7, yend=as.numeric(average_data2['hp',2:7])),size=1)+
  labs(x="", y='Housekeeping genes')


c <- ggplot()+geom_violin(data=gene_data2, aes(x=Stage, y=SALL1, fill=Stage),scale = 'width')+
  geom_point(data=gene_data2,aes(x=Stage, y=SALL1),size=0.5)+scale_fill_manual(values = col_stage)+
  geom_segment(aes(x=1:6, y=as.numeric(average_data2['SALL1',1:6]),
                   xend=2:7, yend=as.numeric(average_data2['SALL1',2:7])),size=1)+
  labs(x="")

d <- ggplot()+geom_violin(data=gene_data2, aes(x=Stage, y=C3, fill=Stage),scale = 'width')+
  geom_point(data=gene_data2,aes(x=Stage, y=C3),size=0.5)+scale_fill_manual(values = col_stage)+
  geom_segment(aes(x=1:6, y=as.numeric(average_data2['C3',1:6]),
                   xend=2:7, yend=as.numeric(average_data2['C3',2:7])),size=1)+
  labs(x="")

e <- ggplot()+geom_violin(data=gene_data2, aes(x=Stage, y=CD163, fill=Stage),scale = 'width')+
  geom_point(data=gene_data2,aes(x=Stage, y=CD163),size=0.5)+scale_fill_manual(values = col_stage)+
  geom_segment(aes(x=1:6, y=as.numeric(average_data2['CD163',1:6]),
                   xend=2:7, yend=as.numeric(average_data2['CD163',2:7])),size=1)+
  labs(x="")

f <- ggplot()+geom_violin(data=gene_data2, aes(x=Stage, y=CX3CR1, fill=Stage),scale = 'width')+
  geom_point(data=gene_data2,aes(x=Stage, y=CX3CR1),size=0.5)+scale_fill_manual(values = col_stage)+
  geom_segment(aes(x=1:6, y=as.numeric(average_data2['CX3CR1',1:6]),
                   xend=2:7, yend=as.numeric(average_data2['CX3CR1',2:7])),size=1)+
  labs(x="")


ggpubr::ggarrange(a,b,c,d,f,e, ncol = 3, nrow = 2, legend = F)

#使用Mpath, slingshot, TSCAN 推断5个群体的时序
# library(Mpath)
# library(ggplot2)
# source('F:/code/function/plot_mpath.R')
# mpath <- plot_mpath(anno = mac_anno[head_macro,'Cluster',drop=F],data = mpath_data,
#            col = col.use2[11:16], dm = mac_anno[, c('UMAP1','UMAP2')], dim.use = 'UMAP', alpha = 0.3)
# 
# mpath + theme_bw()+theme(panel.grid = element_blank(), legend.position = 'none')
# 
# #slingshot
# library(slingshot)
# pca_data <- prcomp(t(mpath_data[rowMeans(mpath_data)>0,]),scale. = T)$x[,1:30]
# 
# l <- getLineages(pca_data, as.character(mac_anno[head_macro,'Cluster']))
# slt <- slingshot(data = mac_anno[head_macro, c('UMAP1','UMAP2')], clusterLabels = as.character(mac_anno[head_macro,'Cluster']),
#                  lineages = l)
# 
# col_map <- col.use2[match(mac_anno$Cluster, names(col.use2))]
# plot(mac_anno[, c('UMAP1','UMAP2')], col = col_map, pch=16)
# lines(slt, lwd=2, type = 'lineages',col='black')

###使用anova 寻找差异表达基因
gene_data <- data.frame(cluster = mac_anno[head_macro,'Cluster'], 
                        t(mpath_data[rowMeans(mpath_data)>0, head_macro]), check.names = F)

aggr_data <- aggregate(.~cluster, data = gene_data, FUN = mean)
average_data <- t(apply(aggr_data[,2:ncol(aggr_data)], 2, function(x) as.numeric(x)))
colnames(average_data) <- aggr_data[,1]


anova.test <- function(gene, data.use, anno.use){
  
  groups <- as.character(anno.use[,1])
  expr <- data.use[gene, rownames(anno.use)]
  
  result <- aov(expr~group, data = data.frame(group=groups, expr= expr))
  
  pvalue <- summary.aov(result)[[1]]$`Pr(>F)`[1]
  Fvalue <- summary.aov(result)[[1]]$`F value`[1]
  
  return(c(pvalue, Fvalue))
}


gene.test <- rownames(average_data)[rowSums(average_data>=1)>=1]
anova_result <- data.frame(t(sapply(gene.test, function(gene) anova.test(gene = gene, data.use = mpath_data, anno.use = mac_anno[head_macro, 'Cluster',drop=F]))))
names(anova_result) <- c('pvalue','F value')
anova_result$padj <- p.adjust(anova_result$pvalue, method = 'fdr')

markers <- read.csv('F:/data/human_gene/human_cd_moculer.csv', header = T)
tfs <- read.delim('F:/data/human_gene/TF_human_db.txt')


intersect(rownames(anova_result)[anova_result$padj< 1e-2], tfs$Symbol)
intersect(rownames(anova_result)[anova_result$padj< 1e-2], markers$symbols)

library(pheatmap)

library(fpc)
phmat <- t(scale(t(average_data)))

library(factoextra)
a <- fviz_nbclust(phmat[rownames(anova_result)[anova_result$padj< 1e-2],], kmeans, method = "wss",k.max=10)
b <- fviz_nbclust(phmat[rownames(anova_result)[anova_result$padj< 1e-2],], FUNcluster = kmeans, method = "silhouette")

library(cluster)
clusterings <- pam(phmat[rownames(anova_result)[anova_result$padj< 1e-2],], k=5)
clusterings <- data.frame(clusterings$clustering)
names(clusterings) <- 'cluster'
clusterings$cluster <- as.character(clusterings$cluster)

gene_order <- c()
for(i in c('2','5','4','1','3')){
  
  tmpgene <- rownames(clusterings)[clusterings$cluster==i]
  gene_order <- c(gene_order, tmpgene)
}

pheatmap(phmat[gene_order,],
         cluster_cols = F, show_rownames = F,show_colnames = T,cluster_rows = F, annotation_row = clusterings,
         color = colorRampPalette(c("navy" ,"white", "red"))(100),annotation_col = data.frame(Cluster=colnames(phmat), row.names = colnames(phmat)),
         annotation_colors = list(Cluster=col.use2[11:15], cluster = setNames(c('#004e66','blue3','red3','darkgreen','#71226e'),as.character(1:5))))

write.csv(file = 'F:/analysis/project/human_macro_development_project/revised_for_190822/tables_revised_need/head_pattern.csv',clusterings)
####
gene_anno <- data.frame(average_data[rownames(clusterings), ],gene=rownames(clusterings), pattern = clusterings[,1], check.names = F)
gene_anno <- reshape2::melt(gene_anno, id= c('gene','pattern'))
names(gene_anno)[3:4] <- c('cluster','expression')

# pattern_anno <- aggregate(.~pattern, 
#                           data.frame(average_data[rownames(clusterings), ], pattern = clusterings$cluster, check.names = F),
#                           mean)


ggplot(gene_anno)+geom_line(aes(x=cluster, y=expression, group=gene, col=pattern), size = 0.1)+facet_wrap(.~pattern, nrow = 2, scales = 'free_y')+
  scale_color_manual(values = setNames(c('#004e66','blue3','red3','darkgreen','#71226e'),as.character(1:5)))+
  theme_bw()+theme(axis.text.x = element_text(angle = 60, hjust = 1))+stat_summary(aes(x=cluster,y=expression,group=1),fun.data="mean_cl_boot",color="black",fill="black",alpha=0.2,size=1.1,geom="smooth")

source('F:/code/function/go_analysis.R')

go1 <- go_analysis(gene = rownames(clusterings)[clusterings$cluster=='2'][1:100],
                    org = 'hsa', title = 'GO term', col = '#004e66', length = 10, on = 'BP')

go2 <- go_analysis(gene = rownames(clusterings)[clusterings$cluster=='5'][1:100],
                    org = 'hsa', title = 'GO term', col = 'blue3', length = 10, on = 'BP')

go3 <- go_analysis(gene = rownames(clusterings)[clusterings$cluster=='4'][1:100],
                    org = 'hsa', title = 'GO term', col = 'red3', length = 10, on = 'BP')

go4 <- go_analysis(gene = rownames(clusterings)[clusterings$cluster=='1'][1:100],
                    org = 'hsa', title = 'GO term', col = 'darkgreen', length = 10, on = 'BP')

go5 <- go_analysis(gene = rownames(clusterings)[clusterings$cluster=='3'][1:100],
                    org = 'hsa', title = 'GO term', col = '#71226e', length = 10, on = 'BP')


####plot gene,
cluster_data <- data.frame(t(mpath_data), cluster= mac_anno[colnames(mpath_data),'Cluster'], check.names = F)
cluster_data <- reshape::melt(cluster_data, id='cluster')
names(cluster_data)[2:3] <- c('gene','expression')

library(ggpubr)
cluster_data$gene <- as.character(cluster_data$gene)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
cluster_sum <- summarySE(cluster_data, measurevar="expression", groupvars=c("gene","cluster"))

pd <- position_dodge(0.1)

ggplot(cluster_sum[cluster_sum$gene %in% c('CD163','JUND','IGF1','TMSB4X','CD1C'),], aes(x=cluster, y=expression, col= cluster))+ 
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se), width=0.5, position=pd, size=1)+
  geom_line(position=pd, group=1, col='black',size=0.5)+theme_bw()+geom_point()+
  scale_color_manual(values = col.use2)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = 'none',panel.grid = element_blank())+
  labs(x="")+facet_wrap('gene', scales = 'free_y',ncol = 1)



#top5 tfs
clusterings$p_val_adj <- anova_result[rownames(clusterings),'padj']
clusterings <- clusterings[order(clusterings$p_val_adj,decreasing = F),]

top5tf <- c()
for (i in c('2','5','4','1','3')){
  
  tmpgene <- intersect(rownames(clusterings)[clusterings$cluster==i], tfs$Symbol)[1:15]
  top5tf <- c(top5tf, tmpgene)
}

pheatmap(phmat[as.character(na.omit(top5tf)),],
         cluster_cols = F, show_rownames = T,show_colnames = T,cluster_rows = F, annotation_row = clusterings[,'cluster',drop=F],
         color = colorRampPalette(c("navy" ,"white", "red"))(100),annotation_col = data.frame(Cluster=colnames(phmat), row.names = colnames(phmat)),
         annotation_colors = list(Cluster=col.use2[11:15], cluster = setNames(c('#004e66','blue3','red3','darkgreen','#71226e'),as.character(1:5))),
         border_color = NA, fontsize_row = 8)


############################################################################
#########################Liver macrophage development#######################
############################################################################
mac_cell <- rownames(mac_anno)[mac_anno$cluster %in% c('Liver_Mac','Blood_Mac','Lung_Mac','Skin_Mac','YS_Mac1','YS_Mac2','Head_Mac1','Head_Mac2',"Head_Mac3",'Head_Mac4')]
mac_anno$fl <- as.character(mac_anno$Stage)
mac_anno$fl[!mac_anno$Site=='Liver']<- 'others'
mac_anno$fl[!rownames(mac_anno) %in% mac_cell] <- 'others'
#mac_anno$fl[mac_anno$UMAP1 <5] <- 'others'


cell_order <- c(rownames(mac_anno)[mac_anno$fl=='others'],
                rownames(mac_anno)[mac_anno$fl=='CS12'],
                rownames(mac_anno)[mac_anno$fl=='CS15'],
                rownames(mac_anno)[mac_anno$fl=='CS17'],
                rownames(mac_anno)[mac_anno$fl=='CS23'])

mac_anno <- mac_anno[cell_order,] 
ggplot(mac_anno,aes(x=UMAP1,y=UMAP2,col=fl))+geom_point()+scale_color_manual('Stage',values = c(col_stage,'others'=scales::alpha('grey',alpha = 1)))+
  theme_bw()+theme(panel.grid = element_blank(),legend.key = element_rect(size=0.3,color='#808080'),legend.title = element_text(face='bold'))


####DPT
library(destiny)
fl.use <- rownames(mac_anno)[!mac_anno$fl=='others']
fl_data <- data[, fl.use]
fl_data <- fl_data[rowMeans(fl_data)>0, ]

pca <- prcomp(t(fl_data), scale. = T)
dm <- DiffusionMap(pca$x[,1:10])
dpt <- DPT(dm, tips = 7)
plot(eigenvalues(dm), ylim = 0:1, pch = 20)


dm_data <- data.frame(DC1 = dm$DC1, DC2 = 0-dm$DC2, dptval = dpt$dpt,
                      stage = mac_anno[fl.use, 'Stage'], ID= 1:length(fl.use))

dm_data <- data.frame(dm_data, t(fl_data[c('CDH5','CD5L','SPIC','VCAM1','C1QB','CCR2','S100A4','IL17RA','CD163','TIMD4','VSIG4'),]))

p1 <- ggplot(dm_data, aes(x=DC1, y= DC2))+geom_point(aes(fill=stage),size=3, pch=21, col='black')+scale_fill_manual(values = col_stage)+
  theme_bw()

p2 <- ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, fill=dptval),size=3, pch=21, col='black')+scale_fill_gradientn('DPT',colours = colorRampPalette(c('blue','green','yellow','red'))(100))+
  theme_bw()

ggpubr::ggarrange(p1,p2, nrow = 2)

p3 <- ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, col=C1QB),size=2)+scale_color_gradientn("",colours = colorRampPalette(c('blue','green','yellow','red'))(100))+
  theme_bw()+labs(title = 'C1QB')

p4 <- ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, col=CD5L),size=2)+scale_color_gradientn("",colours = colorRampPalette(c('blue','green','yellow','red'))(100))+
  theme_bw()+labs(title = 'CD5L')

p5 <- ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, col=SPIC),size=2)+scale_color_gradientn("", colours = colorRampPalette(c('blue','green','yellow','red'))(100))+
  theme_bw()+labs(title = 'SPIC')

p6 <- ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, col=VCAM1),size=2)+scale_color_gradientn("",colours = colorRampPalette(c('blue','green','yellow','red'))(100))+
  theme_bw()+labs(title = 'VCAM1')

p7 <- ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, col=CCR2),size=2)+scale_color_gradientn("",colours = colorRampPalette(c('blue','green','yellow','red'))(100))+
  theme_bw()+labs(title = 'CCR2')

p8 <- ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, col=S100A4),size=2)+scale_color_gradientn("",colours = colorRampPalette(c('blue','green','yellow','red'))(100))+
  theme_bw()+labs(title = 'S100A4')

p9 <- ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, col=IL17RA),size=2)+scale_color_gradientn("",colours = colorRampPalette(c('blue','green','yellow','red'))(100))+
  theme_bw()+labs(title = 'IL17RA')

p10 <- ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, col=CD163),size=2)+scale_color_gradientn("",colours = colorRampPalette(c('blue','green','yellow','red'))(100))+
  theme_bw()+labs(title = 'CD163')

p11 <- ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, col=TIMD4),size=2)+scale_color_gradientn("",colours = colorRampPalette(c('blue','green','yellow','red'))(100))+
  theme_bw()+labs(title = 'TIMD4')

p12 <- ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, col=VSIG4),size=2)+scale_color_gradientn("",colours = colorRampPalette(c('blue','green','yellow','red'))(100))+
  theme_bw()+labs(title = 'VSIG4')

p13 <- ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, col=CDH5),size=2)+scale_color_gradientn("",colours = colorRampPalette(c('blue','green','yellow','red'))(100))+
  theme_bw()+labs(title = 'CDH5')

plot_grid(p1,p13)


ggpubr::ggarrange(p3,p4,p5,p6,ncol = 2, nrow = 2)
ggpubr::ggarrange(p7,p8,p9,p10,p11,p12,ncol = 2, nrow = 3)

library(ggbeeswarm)
library(ggthemes)

p7 <- ggplot(dm_data, aes(x=DC1, y= stage, col=stage))+
  scale_color_manual(values = col_stage)+geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() 

p8 <- ggplot(dm_data, aes(x=dptval, y= stage, col=stage))+
  scale_color_manual(values = col_stage)+geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() 

ggarrange(p7,p8)

#############find pattern
gene_data2 <- data.frame(cluster = mac_anno[fl.use,'Stage'], 
                        t(fl_data[rowMeans(fl_data)>0, fl.use]), check.names = F)

aggr_data2 <- aggregate(.~cluster, data = gene_data2, FUN = mean)
average_data2 <- t(apply(aggr_data2[,2:ncol(aggr_data2)], 2, function(x) as.numeric(x)))
colnames(average_data2) <- aggr_data2[,1]

gene.test2 <- rownames(average_data2)[rowSums(average_data2 >=1 ) >= 1]
anova_result2 <- data.frame(t(sapply(gene.test2, function(gene) anova.test(gene = gene, data.use = fl_data, anno.use = mac_anno[fl.use, 'Stage',drop=F]))))
names(anova_result2) <- c('pvalue','F value')
anova_result2$padj <- p.adjust(anova_result2$pvalue, method = 'fdr')

intersect(rownames(anova_result2)[anova_result2$padj< 1e-2], tfs$Symbol)
intersect(rownames(anova_result2)[anova_result2$padj< 1e-2], markers$symbols)

phmat2 <- t(scale(t(average_data2)))

c <- fviz_nbclust(phmat2[rownames(anova_result2)[anova_result2$padj< 1e-2],], kmeans, method = "wss",k.max=10)
d <- fviz_nbclust(phmat2[rownames(anova_result2)[anova_result2$padj< 1e-2],], FUNcluster = kmeans, method = "silhouette")

clusterings2 <- pam(phmat2[rownames(anova_result2)[anova_result2$padj< 1e-2],], k=3)
clusterings2 <- data.frame(clusterings2$clustering)
names(clusterings2) <- 'cluster'
clusterings2$cluster <- as.character(clusterings2$cluster)

gene_order2 <- c()
for(i in c('2','3','1')){
  
  tmpgene <- rownames(clusterings2)[clusterings2$cluster==i]
  gene_order2 <- c(gene_order2, tmpgene)
}

write.csv(file = 'F:/analysis/project/human_macro_development_project/revised_for_190822/tables_revised_need/liver_pattern.csv',clusterings2)

pheatmap(phmat2[gene_order2,],
         cluster_cols = F, show_rownames = F,show_colnames = T,cluster_rows = F, annotation_row = clusterings2,
         color = colorRampPalette(c("navy" ,"white", "red"))(100),annotation_col = data.frame(Cluster=colnames(phmat2), row.names = colnames(phmat2)),
         annotation_colors = list(Cluster=col_stage[c('CS12','CS15','CS17','CS23')], cluster = setNames(c('red3','green3','blue3','purple3','magenta'),c('1','2','3','4','5'))),border_color = NA)


pheatmap(phmat2[intersect(gene_order2, c(tfs$Symbol, markers$symbols)),],
         cluster_cols = F, show_rownames = T,show_colnames = T,cluster_rows = F, annotation_row = clusterings2,
         color = colorRampPalette(c("navy" ,"white", "red"))(100),annotation_col = data.frame(Cluster=colnames(phmat2), row.names = colnames(phmat2)),
         annotation_colors = list(Cluster=col_stage[c('CS12','CS15','CS17','CS23')], cluster = setNames(c('red3','green3','blue3','purple3','magenta'),c('1','2','3','4','5'))),border_color = NA)

write.csv(file = 'F:/analysis/project/human_macro_development_project/revised_for_190822/head/anova_result.csv',anova_result)
write.csv(file = 'F:/analysis/project/human_macro_development_project/revised_for_190822/head/head_pattern.csv',clusterings)

write.csv(file = 'F:/analysis/project/human_macro_development_project/revised_for_190822/liver/anova_result.csv',anova_result2)
write.csv(file = 'F:/analysis/project/human_macro_development_project/revised_for_190822/liver/liver_pattern.csv',clusterings2)

####skin macrophage deduce using DPT
library(destiny)
load('F:/analysis/project/human_macro_development_project/mcrophage_result_version5/Fig4/fig4_v5.rds')
skin.use <- rownames(mac_anno)[mac_anno$Site=='Skin' &
                                 mac_anno$cluster %in% c('YS_Mac1','YS_Mac2',
                                                         "Blood_Mac",'Liver_Mac','Skin_Mac','Lung_Mac',
                                                         'Head_Mac1','Head_Mac2','Head_Mac3','Head_Mac4')]

mac_data <- apply(rawdata, 2,function(x) log2(1e5*x/sum(x)+1))
skin_data <- mac_data[, skin.use]
skin_data <- skin_data[rowMeans(skin_data)>0, ]

pca <- prcomp(t(skin_data), scale. = T)
dm <- DiffusionMap(pca$x[,1:15])
dpt <- DPT(dm)
plot(eigenvalues(dm), ylim = 0:1, pch = 20)


dm_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],PC3 = pca$x[,3],DC1 = dm$DC1, DC2 = dm$DC2, dptval = dpt$dpt,
                      stage = mac_anno[skin.use, 'Stage'], ID= 1:length(skin.use))

rownames(dm_data) <- colnames(skin_data)

dm_data <- data.frame(dm_data, t(mac_data[c('CD14','CD1A','MRC1','CD1C','CD207','HLA-DRA') , rownames(dm_data)]), check.names = F)

a <- ggplot(dm_data, aes(x=DC1, y= DC2))+geom_point(aes(fill=stage),size=3, pch=21, col='black')+scale_fill_manual(values = col_stage)+
  theme_bw()

b <- ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, fill= CD1A),size=3, pch=21, col='black')+scale_fill_gradientn('DPT',colours = colorRampPalette(c('blue','green','yellow','red'))(100))+
  theme_bw()

c <- ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, fill=CD1C),size=3, pch=21, col='black')+scale_fill_gradientn(colours = colorRampPalette(c('blue','green','yellow','red'))(100))+
  theme_bw()

d <- ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, fill=CD207),size=3, pch=21, col='black')+scale_fill_gradientn(colours = colorRampPalette(c('blue','green','yellow','red'))(100))+
  theme_bw()

ggpubr::ggarrange(a,b,c,d, ncol = 4)


ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, fill = `HLA-DRA`),size=3, pch=21, col='black')+scale_fill_gradientn(colours = colorRampPalette(c('#E6E6E6','blue','green','yellow','red'))(100))+
  theme_bw()

ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, fill = `MRC1`),size=3, pch=21, col='black')+scale_fill_gradientn(colours = colorRampPalette(c('#E6E6E6','blue','green','yellow','red'))(100))+
  theme_bw()

ggplot(dm_data)+geom_point(aes(x=DC1, y= DC2, fill = `CD14`),size=3, pch=21, col='black')+scale_fill_gradientn(colours = colorRampPalette(c('#E6E6E6','blue','green','yellow','red'))(100))+
  theme_bw()

####cell proportions 
load('F:/analysis/project/human_macro_development_project/mcrophage_result_version5/Fig1/fig1_v5.rds')
scales::show_col(col_all)

plot_single_stage <- function(time, anno){
  
  anno1 <- anno
  anno1$cluster[!anno1$Stage==time]<- 'others'
  
  anno1 <- anno1[c(rownames(anno1)[anno1$cluster=='others'],
                   rownames(anno1)[!anno1$cluster=='others']),]
  g1 <- ggplot(anno1, aes(x=UMAP1, y=UMAP2, col=cluster))+geom_point()+scale_color_manual(values = c(col_all, 'others'='#E6E6E6'))+theme_bw()
  
  num2 <- data.frame(table(anno1[anno1$Stage==time, ]$cluster))
  names(num2) <- c('cluster','number')
  num2$number <- round(100*num2$number/sum(num2$number),0)
  
  
  g2 <- ggpie(num2, x="number", fill='cluster')+scale_fill_manual(values = col_all)+
    theme(legend.position = 'none')+
    ggrepel::geom_label_repel(label = paste(num2$cluster, ' (',num2$number,'%',')',sep=""))
  
  ggarrange(g1,g2)
}

plot_single_stage(time = 'CS11', anno = all_anno)
plot_single_stage(time = 'CS12', anno = all_anno)
plot_single_stage(time = 'CS13', anno = all_anno)
plot_single_stage(time = 'CS15', anno = all_anno)
plot_single_stage(time = 'CS17', anno = all_anno)
plot_single_stage(time = 'CS23', anno = all_anno)


####Fig7 macrophage comparasion
####mainly: Liver, Lung, Skin, Head
# Head_Mac4, Skin_Mac, Liver_Mac, Lung_Mac mainly in cs23
our_anno <- mac_anno[mac_anno$cluster %in% c('Head_Mac4','Liver_Mac','Lung_Mac','Skin_Mac'),]
our_anno$cluster <- plyr::mapvalues(our_anno$cluster,c('Head_Mac4','Liver_Mac','Lung_Mac','Skin_Mac'),
                                    c('Head_Mac','Liver_Mac','Lung_Mac','Skin_Mac'))
our_rawdata <- rawdata[,rownames(our_anno)]




#other datas
load('F:/analysis/project/human_macro_development_project/rawdata/fl_10x.rds')
load('F:/analysis/project/human_macro_development_project/rawdata/lung_10x.rds')
load('F:/analysis/project/human_macro_development_project/rawdata/lc_mono.rds')
load('F:/analysis/project/human_macro_development_project/rawdata/microglia_head_rawdata.rds')



#arrange data
fl_anno$Site <- 'Liver'
fl_anno$Stage <- 'Adult'
#fl_anno <- fl_anno[sample(rownames(fl_anno),30),]
fl_rawdata <- fl_rawdata[,rownames(fl_anno)]
rownames(fl_rawdata) <- reshape2::colsplit(rownames(fl_rawdata), pattern = "_",names = c('id1','id2'))[,1]




lc_mono_anno <- lc_mono_anno[lc_mono_anno$Site=='Skin',c('Site','Stage')]
lc_mono_anno$Stage <- 'Child'
lc_mono_rawdata <- lc_mono_rawdata[, rownames(lc_mono_anno)]

###
lung_anno <- data.frame(Site=rep('Lung',ncol(lung_rawdata)),
                        Stage= rep('Adult',ncol(lung_rawdata)),
                        row.names = colnames(lung_rawdata))

#lung_anno <- lung_anno[sample(rownames(lung_anno),30),]
lung_rawdata <- lung_rawdata[,rownames(lung_anno)]

###
microglia_anno <- data.frame(Site=rep('Head',ncol(microglia_data)),
                             Stage= rep('Adult',ncol(microglia_data)),
                             row.names = colnames(microglia_data))


###combine data
co_gene <- intersect(rownames(our_rawdata), intersect(rownames(fl_rawdata),
                                            intersect(rownames(lung_rawdata),
                                            intersect(rownames(lc_mono_rawdata),rownames(microglia_data)))))

all_rawdata <- cbind(our_rawdata[co_gene,],
                     fl_rawdata[co_gene,],
                     lung_rawdata[co_gene,],
                     lc_mono_rawdata[co_gene,],
                     microglia_data[co_gene,])

all_anno <- data.frame(rbind(our_anno[,c('Site','Stage')],
                             fl_anno[,c('Site','Stage')],
                             lung_anno[,c('Site','Stage')],
                             lc_mono_anno[,c('Site','Stage')],
                             microglia_anno[,c('Site','Stage')]))

all_anno <- all_anno[all_anno$Site %in% c('Head','Lung','Skin','Liver'),]
all_rawdata <- all_rawdata[,rownames(all_anno)]



all_anno$Stage[!all_anno$Stage %in% c('Adult','Child')]<- 'embryo'
all_anno$group <- paste(all_anno$Site," ( ",all_anno$Stage,' )',sep="")
all_anno$batch <- plyr::mapvalues(all_anno$group, unique(all_anno$group),
                                  c(rep('batch1',4),"batch2",'batch3','batch4','batch5'))

all_anno$method <- ifelse(all_anno$group=='Head ( Adult )','ss2','ss' )

table(all_anno$group)
####
all_anno$Group <- all_anno$group
all_anno$Group[!all_anno$Group %in% c("Head ( Adult )",'Head ( embryo )')]<- 'others'

library(Seurat)
ana <- CreateSeuratObject(all_rawdata, meta.data = all_anno)
ana <- NormalizeData(ana ,scale.factor = 1e4)
ana <- ScaleData(ana, vars.to.regress = c("nUMI"))
ana <- FindVariableGenes(ana, x.low.cutoff = 0.5, y.cutoff = 1)
ana <- RunPCA(ana, pc.genes = ana@var.genes, pcs.compute = 50)
PCElbowPlot(ana, num.pc = 50)




a <- DimPlot(ana, reduction.use = 'pca', group.by = 'Site',do.return = T,dim.1 = 1,dim.2 = 3)+scale_color_manual(values = col_site)
b <- DimPlot(ana, reduction.use = 'pca', group.by = 'Stage',do.return = T,dim.1 = 1, dim.2 = 3)
c <- DimPlot(ana, reduction.use = 'pca', group.by = 'group',do.return = T,dim.1 = 1, dim.2 = 2)+scale_color_manual(values = mclust::mclust.options()$classPlotColors)
ggpubr::ggarrange(a,b,c,ncol = 3,nrow = 1)

ana <-  RunTSNE(ana, dims.use = 1:6)
ana <-  RunUMAP(ana, dims.use = 1:6,n_neighbors = 30,min_dist = 0.5,spread=1)

a <- DimPlot(ana, reduction.use = 'umap', group.by = 'Site',do.return = T)+scale_color_manual(values = col_site)
b <- DimPlot(ana, reduction.use = 'umap', group.by = 'Stage',do.return = T)
c <- DimPlot(ana,pt.size=2, reduction.use = 'umap', group.by = 'group',do.return = T)+scale_color_manual(values = mclust::mclust.options()$classPlotColors)

ggpubr::ggarrange(a,b,c,ncol = 3)


####umap+kmeans
km <- kmeans(ana@dr$umap@cell.embeddings,centers = 5)
km_anno <- data.frame(km$cluster)
names(km_anno) <- 'cluster'

ana@meta.data$km_cluster <- km_anno[rownames(ana@meta.data),,drop=F]$cluster
d <- DimPlot(ana,pt.size=2, reduction.use = 'umap', group.by = 'km_cluster',do.return = T)+scale_color_manual(values = mclust::mclust.options()$classPlotColors)
plot_grid(c,d)

ana@meta.data$km_cluster <- plyr::mapvalues(ana@meta.data$km_cluster,
                                            c(1:5),c('Lung macrophage',
                                                     'Unspecified macrophage',
                                                     'Skin macrophage',
                                                     'Liver macrophage',
                                                     'Head macrophage'))






