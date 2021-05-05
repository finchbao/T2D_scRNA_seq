##############################################
'
cell cluster identification
'
##############################################

library(Cairo)
library(Seurat)
library(dplyr)
library(magrittr)
library(patchwork)
library(monocle)
library(reshape2)
library(ggplot2)
library(DESeq2)
library(BiocParallel)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(magick)
############################################## data input
setwd('/data/phenomics/cell_report/')

folders <- list.files(pattern="*txt.gz")

SceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = read.table(folder,row.names = 1,header = T),
                     min.cells =1 , 
                     min.features = 100)
})

gse101207_list<-SceList

############################################## label
for (i in 1: 5 ) {
  gse101207_list[[i]]$group<-'normal'
}

gse101207_list[[8]]$group<-'normal'
gse101207_list[[6]]$group<-'T2D'
gse101207_list[[7]]$group<-'T2D'
gse101207_list[[9]]$group<-'T2D'

for (i in 1: 5 ) {
  gse101207_list[[i]]$group2<-paste('normal',c(i))
}

gse101207_list[[8]]$group2<-'normal6'
gse101207_list[[6]]$group2<-'T2D1'
gse101207_list[[7]]$group2<-'T2D2'
gse101207_list[[9]]$group2<-'T2D3'

############################################## merge data
for (i in 1:length(gse101207_list)) {
  gse101207_list[[i]] <- NormalizeData(gse101207_list[[i]], verbose = FALSE)
  gse101207_list[[i]] <- FindVariableFeatures(gse101207_list[[i]], selection.method = "vst", nfeatures = 2000)
}

for (i in 1:length(gse101207_list)) {
  gse101207_list[[i]][["percent.mt"]] <- PercentageFeatureSet(gse101207_list[[i]], pattern = "^MT-")
}#### fliter MT rna 

gse101207_list_TC<- FindIntegrationAnchors(object.list = gse101207_list, dims = 1:20)
gse101207_list_TC1 <- IntegrateData(anchorset = gse101207_list_TC, dims = 1:20)

############################################## dimensionality reduction
gse101207_list_TC2 <- ScaleData(gse101207_list_TC1)
gse101207_list_TC2 <- RunPCA(gse101207_list_TC2)
tsne1<-RunTSNE(object = gse101207_list_TC2,dims.use = 1:20,do.fast = TRUE,check_duplicates = FALSE)
tsne1 <-RunUMAP(tsne1,dims=1:20)

############################################## cluster
tsne1<- FindNeighbors(tsne1, reduction = 'pca', dims = 1:20)
tsne1 <- FindClusters(object = tsne1,resolution =0.15)

FeaturePlot(tsne1, features = c("INS", "GCG","SST","PPY",'GHRL',"KRT19","COL1A2", "REG1A",'FABP4','ANGPT2'),
            pt.size = 0.00001,label =F,ncol = 2,cols = c('#FFCCCC','#FF6666','#FF0033'))& NoLegend() &theme (title = element_text (size=10))

new.cluster.ids <- c("beta",'alpha','alpha','sigma',
                     'beta','PDCs','PCS','gamma',
                     'Endothelial Cells','Acinar','alpha')
new.cluster.ids <- levels(tsne1)
names(new.cluster.ids) <- levels(tsne1)
tsne1 <- RenameIdents(tsne1, new.cluster.ids)
DimPlot(tsne1, reduction = "umap", label = F, pt.size = 0.5) +theme_bw()+theme(panel.grid=element_blank())+NoLegend()

saveRDS(tsne1,'/data/phenomics/cell_report/final_merge_gse101207.rds')##### save result

markers <- FindAllMarkers(object = tsne1, only.pos = TRUE)
beta_markers = row.names(subset(markers,cluster='beta'))
##############################################
'
beta cell trajectory construction in pesudotime
'
##############################################

gse101207_beta_cell =subset(gse101207, idents = c('beta')) ####subset beta cell
data <- as(as.matrix(gse101207_beta_cell@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data =gse101207_beta_cell@meta.data) 
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData) 
HSMM <- newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())

##############################################estimate data
sce <- estimateSizeFactors(HSMM)
sce <- estimateDispersions(sce)
disp_table <- dispersionTable(sce)

##############################################filter
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.01& dispersion_empirical >= 1 * dispersion_fit)
sce <- setOrderingFilter(sce, unsup_clustering_genes$gene_id) 
plot_ordering_genes(sce)

##############################################trajectory construction
sce <- reduceDimension(sce, max_components = 2,method = 'DDRTree')
sce <- orderCells(sce)
sce<- orderCells(sce, root_state = 4)
plot_cell_trajectory(sce, color_by = "Pseudotime",
                     cell_size = 0.5)+guides(colour = guide_legend(override.aes = list(size=1)))
plot_cell_trajectory(sce, color_by = "group",
                     cell_size = 0.5,show_state_number = F)+guides(colour = guide_legend(override.aes = list(size=1)))
plot_cell_trajectory(sce, color_by = "group2",
                     cell_size = 0.5,show_state_number = F)+guides(colour = guide_legend(override.aes = list(size=1)))
plot_cell_trajectory(sce, color_by = "State",
                     cell_size = 0.5,show_state_number = T)+guides(colour = guide_legend(override.aes = list(size=1)))

##############################################marker gene
sce_gene = row.names(fData(sce))
b = c('INS','EIF2AK3','ATG5','ATG7','CASP3','CASP6','CASP7','GPD2','CLEC16A','IRS2','MAFA','SKAP2',
      'STRA6','VMAT2','HRD1','SIN3A','SIN3B','PTEN')
b<-intersect(sce_gene,b)
plot_genes_branched_heatmap(sce[b,],branch_point = 2,show_rownames = T,cluster_rows = F,num_clusters = 3)
dev.off()


de_cluster_one <- differentialGeneTest(sce[disp_table$gene_id,],
                                       fullModelFormulaStr = '~group',
                                       cores = 32)
##############################################beam
Beam_res = BEAM(sce,branch_point = 2,cores = 12)
Beam_res = Beam_res[order(Beam_res$qval),]
beam_deg= row.names(subset(Beam_res, qval<0.001))
beam_deg_0.001 = plot_genes_branched_heatmap(sce[beam_deg,], branch_point = 2,num_clusters = 4,cores = 12, show_rownames = F,return_heatmap = T)
plot_genes_branched_heatmap(sce[beam_deg,], branch_point = 2,num_clusters = 4,cores = 12, show_rownames = F,return_heatmap = F,)+theme(legend.box = "vertical")

write.csv(beam_deg_0.001$annotation_row,'/data/phenomics/cell_report/beam_DEG_cluster_0.001.csv')
##############################################
'
DEseq2 for differential expression
'
##############################################

register(MulticoreParam(workers=32))#####set cores number

##############################################preprocess
mtx =as.matrix(sce@assayData[["exprs"]])
countData <- mtx[apply(mtx, 1, sum) > 1 , ] ####fliter
countData=countData+1

##############################################label set
label_bystate = data.frame(cell = sce@assayData[["exprs"]]@Dimnames[[2]], state = as.character(sce@phenoData@data[["State"]]))
label_bystate$state[label_bystate$state=="1"] = 'B'
label_bystate$state[label_bystate$state=="2"] = 'B'
label_bystate$state[label_bystate$state=="3"] = 'B'
label_bystate$state[label_bystate$state=="4"] = 'A'
label_bystate$state[label_bystate$state=="5"] = 'C'
label_bystate$state[label_bystate$state=="6"] = 'C'
label_bystate$state[label_bystate$state=="7"] = 'C'

label_bygroup = data.frame(cell = sce@assayData[["exprs"]]@Dimnames[[2]], group = sce@phenoData@data[["group"]])
##############################################dds
dds_bystate<- DESeqDataSetFromMatrix(countData = countData, colData = label_bystate, design= ~ state) 
dds_bygroup <- DESeqDataSetFromMatrix(countData = countData, colData = label_bygroup, design= ~ group) 
dds_res_bystate<- DESeq(dds_bystate,parallel = T)
dds_res_bygroup <- DESeq(dds_bygroup,parallel = T)

############################################## result
res_bystate_ab<- results(dds_res_bystate, contrast=c("state","A","B"))
res_bystate_bc<- results(dds_res_bystate, contrast=c("state","B","C"))
res_bystate_ac<- results(dds_res_bystate, contrast=c("state","A","C"))
res_bygroup <- results(dds_res_bygroup)

mtx_bystate_ab = as.data.frame(res_bystate_ab)
mtx_bystate_bc = as.data.frame(res_bystate_bc)
mtx_bystate_ac = as.data.frame(res_bystate_ac)
mtx_bygroup = as.data.frame(res_bygroup)

mtx_bystate_ab[,7] = rownames(mtx_bystate_ab)
mtx_bystate_ac[,7] = rownames(mtx_bystate_ac)

##############################################fliter
obv.gene_bygroup = row.names(subset(mtx_bygroup, padj<0.001))
obv.gene_bystate_ab = row.names(subset(mtx_bystate_ab, padj<0.05))
obv.gene_bystate_bc = row.names(subset(mtx_bystate_bc,padj<0.05))
obv.gene_bystate_cb = row.names(subset(mtx_bystate_bc,padj<0.05 & log2FoldChange<= 0))
obv.gene_bystate_ac = row.names(subset(mtx_bystate_ac,  padj<0.05))


paper_t2d = read.csv('/data/phenomics/cell_report/2.csv',header = F)
paper_t2d = paper_t2d[,1]

paper_ob = read.csv('/data/phenomics/cell_report/1.csv',header = F)
paper_ob = paper_ob[,1]

intersect(obv.gene_bystate_cb,intersect(obv.gene_bystate_bc,
                                        intersect(obv.gene_bystate_ab,obv.gene_bystate_ac)))


bc = intersect(paper_t2d,obv.gene_bystate_bc)
intersect(paper_ob,obv.gene_bystate_ab)
intersect(obv.gene_bystate_ab,obv.gene_bystate_cb)

diff_t2d = setdiff(obv.gene_bystate_bc,paper_t2d)
diff_ob = setdiff(obv.gene_bystate_ab,paper_ob)
diff_3 = setdiff(diff_t2d,diff_ob)

write.csv(bc,'/data/phenomics/cell_report/bc.csv')
write.csv(diff_t2d,'/data/phenomics/cell_report/diff_t2d.csv')
write.csv(diff_ob,'/data/phenomics/cell_report/diff_ob.csv')
write.csv(diff_3,'/data/phenomics/cell_report/diff_3.csv')
##############################################
'
visualization
'
##############################################
##############################################fig1
p1 = FeaturePlot(tsne1, features = c("INS", "GCG","SST","PPY",'GHRL',"KRT19","COL1A2", "REG1A",'FABP4','ANGPT2'),
                 pt.size = 0.00001,label =F,ncol = 2,cols = c('#FFCCCC','#FF6666','#FF0033')) &theme_void()& theme(title = element_text (size=8),plot.title = element_text(hjust = 0.5))& NoLegend()
p2 = DimPlot(tsne1, reduction = "umap",label = F,pt.size = 0.5,repel = T)+NoLegend()
p3 = DimPlot(tsne1, reduction = "umap",repel = TRUE,group.by = 'group',pt.size = 0.5)+labs(title = "")+
  theme(legend.position = c(0.05, 0.1))+scale_colour_hue(labels=c("Normal","T2D"))

((((p2+plot_spacer()+plot_layout(heights = c(1, 0.0025)))/(plot_spacer()+p3+plot_layout(heights = c(0.0025,1)))))|p1)+
  plot_layout(widths = c(1.25, 1))& 
  theme(plot.tag = element_text(size = 20,hjust = 0,vjust = 0.5))
##############################################fig2
p4 = plot_cell_trajectory(sce, color_by = "Pseudotime",
                          cell_size = 0.5)+guides(colour = guide_legend(override.aes = list(size=1)))+ theme(legend.direction = "horizontal")+
  scale_color_gradient2(mid = "#009966",high='#99CC33')+NoLegend()

p5 = plot_cell_trajectory(sce, color_by = "group",
                          cell_size = 0.5,show_state_number = F)+guides(colour = guide_legend(override.aes = list(size=1)))+ 
  theme(legend.position = c(0.8, 0.8),legend.text = element_text(size=15),legend.title =element_text(size = 0),legend.key.width=unit(1,'cm'))+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_colour_hue(labels=c("Normal","T2D"))

sce_gene = row.names(fData(sce))
b = c('INS','EIF2AK3','ATG5','ATG7','CASP3','CASP6','CASP7','GPD2','CLEC16A','IRS2','MAFA','SKAP2',
      'STRA6','VMAT2','HRD1','SIN3A','SIN3B','PTEN')
b<-intersect(sce_gene,b)

g = plot_genes_branched_heatmap(sce[b,],branch_point = 2,show_rownames = T,cluster_rows = F,num_clusters = 3,
                                return_heatmap = T)
dev.off()
p9 = as.ggplot(g)

table(pData(sce)$State)
unique = sce

state_sample = table(pData(unique)$group2,pData(unique)$State)
colnames(state_sample) = c('B','B','B','A','C','C','C')
a_group<-table(pData(sce)$group2)
for (i in 1:9) {
  for (j in 1:7) {
    state_sample[i,j] <-state_sample[i,j]/a_group[i]
  }
}

library(reshape2)
plot_a1<-melt(state_sample)

mycolors = brewer.pal(14, "RdYlGn")
cols<-brewer.pal(8, "RdBu")
pal<-colorRampPalette(cols)
mycolors<-pal(12)

p7 = ggplot(plot_a1,aes(x=as.character(Var2),weight=value)) + 
  geom_bar(aes(fill=factor(Var1))) + 
  theme_classic() +labs(x = "Branch", y = "Distribution of Cells")+  
  scale_x_discrete(breaks=c("A", "B", "C"),labels=c("Normal", "Obesity-like",'T2D-like'))+
  scale_fill_manual(values=mycolors, name="Group",labels=c("Normal 1", "Normal 2", "Normal 3","Normal 4", "Normal 5", "Normal 6","T2D 1","T2D 2","T2D 3"))


plot_a1$Var2 <- factor(plot_a1$Var2, levels=c('A', 'B', 'C'))
p8 = ggplot(plot_a1,aes(x=Var1,weight=value)) + geom_bar(aes(fill=factor(Var2)),position = 'fill') + 
  theme_classic()+scale_fill_manual(name="Branch",
                                    breaks=c("A", "B", "C"),
                                    labels=c("Normal", "Obesity-like", "T2D-like"),
                                    values=c("#FF6666", "#FFFF66", "#99CC66")) +  
  scale_x_discrete(breaks=c("normal 1", "normal 2", "normal 3","normal 4", "normal 5", "normal6",'T2D1','T2D2','T2D3'),
                   labels=c("N1", "N2",'N3',"N4", "N5",'N6','T2D1','T2D2','T2D3'))+
  labs(x = "Group", y = "Distribution of Cells")


p6 <- image_read("/data/phenomics/cell_report/Rplot01.pdf")



(((p4/p5)|plot_spacer())+plot_layout(widths = c(1, 1.5)))/((p7|p8)+plot_layout(widths = c(1.6, 2)))+plot_layout(heights = c(3,1))& 
  theme(plot.tag = element_text(size = 20,hjust = 0,vjust = 0.5))


##############################################fig3
first = read.csv('first_up_to_down.csv',row.names = 1)
second = read.csv('second_up_to_down.csv',row.names = 1)
third = read.csv('third_up_to_down.csv',row.names = 1)
four = read.csv('four_up_to_down.csv',row.names = 1)

p1 = ggplot(first,aes(Ratio,Description)) +
  geom_point() +
  geom_point(aes(size=Counts,color=-Log.q.value.)) +
  scale_color_gradient2(low="#3399CC",mid='#99CC99',high = "#FF9900")+
  labs(title = "Cluster1")+ylab('')+scale_y_discrete(position = "right")+
  theme(plot.title = element_text(hjust = 0.5))

p2 = ggplot(second,aes(Ratio,Description)) +
  geom_point() +
  geom_point(aes(size=Counts,color=-Log.q.value.)) +
  scale_color_gradient2(low="#3399CC",mid='#99CC99',high = "#FF9900")+
  labs(title = "Cluster2")+ylab('')+scale_y_discrete(position = "right")+
  theme(plot.title = element_text(hjust = 0.5))

p3 = ggplot(third,aes(Ratio,Description)) +
  geom_point() +
  geom_point(aes(size=Counts,color=-Log.q.value.)) +
  scale_color_gradient2(low="#3399CC",mid='#99CC99',high = "#FF9900")+
  labs(title = "Cluster3")+ylab('')+scale_y_discrete(position = "right")+
  theme(plot.title = element_text(hjust = 0.5))

p4= ggplot(four,aes(Ratio,Description)) +
  geom_point() +
  geom_point(aes(size=Counts,color=-Log.q.value.)) +
  scale_color_gradient2(low="#3399CC",mid='#99CC99',high = "#FF9900")+
  labs(title = "Cluster4")+ylab('')+scale_y_discrete(position = "right")+
  theme(plot.title = element_text(hjust = 0.5))

(p1/p2/p3/p4)
###############################fig4
x1 <- list(This=obv.gene_bystate_ab, Fang=paper_ob)
venn(x1,zcolor='style')
x2 <- list(This=obv.gene_bystate_ac, Fang=paper_t2d)
venn(x2,zcolor='style')
x3 <- list(This=obv.gene_bystate_bc, Fang=paper_t2d)
venn(x3,zcolor='style')

ob = read.csv('ob.csv')
t2d = read.csv('t2d.csv')


p5 = ggplot(ob,aes(Ratio,Description)) +
  geom_point() +
  geom_point(aes(size=Counts,color=-Log.q.value.)) +
  scale_color_gradient2(low="#3399CC",mid='#99CC99',high = "#FF9900")+
  labs(title = "OB")+ylab('')+scale_y_discrete(position = "right")+
  theme(plot.title = element_text(hjust = 0.5))

p6= ggplot(t2d,aes(Ratio,Description)) +
  geom_point() +
  geom_point(aes(size=Counts,color=-Log.q.value.)) +
  scale_color_gradient2(low="#3399CC",mid='#99CC99',high = "#FF9900")+
  labs(title = "T2D")+ylab('')+scale_y_discrete(position = "right")+
  theme(plot.title = element_text(hjust = 0.5))
p5/p6