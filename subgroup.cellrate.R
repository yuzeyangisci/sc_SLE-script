library('SeuratObject')
library('Seurat')
#library(DoubletFinder)
library(tidyverse)
library(patchwork)
library('dplyr')
library('ggplot2')
library(clustree)
library(ggpubr)
library(cowplot)
library(reshape2)
library(RColorBrewer)



load('F:/10X_SLE/pbmc.12sample.res2.integrated.seurat.object.RData')
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9,  
                         10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                         20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                         30, 31, 32, 33, 34, 35, 36, 37)

new.cluster.ids <- c(
  "CD4+ T-cells",  "CD4+ T-cells",  "CD4+ T-cells",  "CD4+ T-cells",  "NK cells",
  "CD4+ T-cells",  "CD4+ T-cells",  "CD4+ T-cells",  "CD8+ T-cells",  "CD8+ T-cells",
  "B-cells",  "B-cells",  "CD4+ T-cells",  "CD8+ T-cells",  "CD8+ T-cells",
  "B-cells",  "CD8+ T-cells",  "B-cells",  "B-cells",  "B-cells",
  "NK cells",  "CD8+ T-cells",  "CD4+ T-cells",  "B-cells",  "CD8+ T-cells",
  "Monocytes",  "B-cells",  "Monocytes",  "Monocytes",  "Monocytes",
  "Monocytes",  "Monocytes",  "platelet",  "platelet",  "erythrocytes",
  "pDC",  "B-cells",  "CD4+ T-cells")

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

pbmc$cell.type=Idents(pbmc)
DimPlot(pbmc, reduction = "umap", label = TRUE, repel = T, pt.size = 0.5, label.size = 10) # + NoLegend()

DimPlot(pbmc, reduction = "umap",label = T)

#######
table(pbmc$orig.ident)#?鿴????ϸ????
prop.table(table(Idents(pbmc)))
table(Idents(pbmc), pbmc$orig.ident)#???鲻ͬϸ??Ⱥϸ????
Cellratio <- prop.table(table(Idents(pbmc), pbmc$orig.ident), margin = 2)#??????????????ͬϸ??Ⱥ????
Cellratio <- data.frame(Cellratio)
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#??????תΪ??????
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]


###???ӷ?????Ϣ
sample <- c("AT1.filter", "AT2.filter", "AT3.filter",
            "NC1.filter", "NC2.filter", "NC3.filter",
            "NP1.filter", "NP2.filter", "NP3.filter",
            "SLE_C.filter", "SLE_PE.filter", "SLE_PE31.filter")
group <- c("AT","AT","AT",
           "NC","NC","NC",
           "NP","NP","NP",
           "C","C","C")

#sample <- c("NC1.filter", "NC2.filter", "NC3.filter",
#            "SLE_C.filter", "SLE_PE.filter", "SLE_PE31.filter")
#group <- c("NC","NC","NC",
#           "C","C","C")


samples <- data.frame(sample, group)#???????ݿ?

rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R??????
cellper$group <- samples[rownames(cellper),'group']#R??????

###??ͼչʾ
pplist = list()
sce_groups = c("CD4+ T-cells", "NK cells", "CD8+ T-cells",
               "B-cells", "Monocytes", "platelet", "erythrocytes","pDC")


for(group_ in sce_groups){
  cellper_  = cellper %>% select(one_of(c('sample','group',group_)))#ѡ??һ??????
  colnames(cellper_) = c('sample','group','percent')#??ѡ????????????
  cellper_$percent = as.numeric(cellper_$percent)#??ֵ??????
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))#???·?λ??
  #print(group_)
  #print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot??ͼ
    geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  ###????t????????
  labely = max(cellper_$percent)
  cc<-compare_means(percent ~ group,  data = cellper_,method = "t.test")#,  p.adjust.method = "fdr")
  print(group_)
  print(cc)
  write.table(group_,"F:/10X_SLE/cellratio/all.cellratio.txt", row.names=FALSE,append = TRUE)
  write.table(cc, "F:/10X_SLE/cellratio/all.cellratio.txt", row.names=FALSE,append = TRUE)
  my_comparisons <- list( c("NC", "C") ,c("AT","NP"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test",label = "p.adj")#, p.adjust.method = "bonferroni")
  pplist[[group_]] = pp1
}
pp.all <- pplist



## subgroup
pbmc.T=subset(x = pbmc, subset = seurat_clusters == c("0",	"1",	"2",	"3",	"5",	"6",	"7",	"8",	"9",	"12",	"13",	"14",	"16",	"21",	"22",	"24",	"37"))
pbmc.B=subset(x = pbmc, subset = seurat_clusters == c("10",	"11",	"15",	"17",	"18",	"19",	"23",	"26",	"36"))
pbmc.Monocytes=subset(x = pbmc, subset = seurat_clusters == c("25",	"27",	"28",	"29",	"30",	"31"))
pbmc.NK=subset(x = pbmc, subset = seurat_clusters == c("4","20"))
pbmc.erythrocytes=subset(x = pbmc, subset = seurat_clusters == c("34"))
pbmc.pDC=subset(x = pbmc, subset = seurat_clusters == c("35"))
pbmc.platelet=subset(x = pbmc, subset = seurat_clusters == c("32","33"))

pbmc.TNK <- subset(x = pbmc, subset = seurat_clusters == c("0",	"1",	"2",	"3",	"4", "5",	"6",	"7",	"8",	"9",	"12",	"13",	"14",	"16", "20",	"21",	"22",	"24",	"37"))

## T
all.genes <- rownames(pbmc.T)
pbmc.T <- ScaleData(pbmc.T, features = all.genes)
pbmc.T <- FindVariableFeatures(object = pbmc.T)
pbmc.T <- RunPCA(pbmc.T, features = VariableFeatures(object = pbmc.T))

pbmc.T <- RunUMAP(pbmc.T, reduction = "pca", dims = 1:20)
pbmc.T <- FindNeighbors(pbmc.T, reduction = "pca", dims = 1:10)
pbmc.T <- FindClusters(pbmc.T, resolution = 0.1)

p2 <- DimPlot(pbmc.T, reduction = "umap",label = T)
p2


## TNK
all.genes <- rownames(pbmc.TNK)
pbmc.TNK <- ScaleData(pbmc.TNK, features = all.genes)
pbmc.TNK <- FindVariableFeatures(object = pbmc.TNK)
pbmc.TNK <- RunPCA(pbmc.TNK, features = VariableFeatures(object = pbmc.TNK))

pbmc.TNK <- RunUMAP(pbmc.TNK, reduction = "pca", dims = 1:20)
pbmc.TNK <- FindNeighbors(pbmc.TNK, reduction = "pca", dims = 1:10)
pbmc.TNK <- FindClusters(pbmc.TNK, resolution = 0.8)

#p1 <- DimPlot(pbmc.TNK, reduction = "umap",group.by = "new.ident")
p2 <- DimPlot(pbmc.TNK, reduction = "umap",label = T)
p2


pbmc.TNK.markers <- FindAllMarkers(pbmc.TNK, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.TNK.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(x=pbmc.TNK.markers,
          "pbmc.TNK.FindAllMarkers.markers.csv")

current.cluster.ids <- c(  0, 1, 2, 3, 4, 5, 6, 7)

new.cluster.ids <- c(   "memory CD4 T",   "naive CD8 T",  "naive CD4 T",  "CD16 NK",
                        "CD56 NK",  "memory CD8 T",  "regulatory CD4 T",  "gdT")
  
names(new.cluster.ids) <- levels(pbmc.TNK)
pbmc.TNK <- RenameIdents(pbmc.TNK, new.cluster.ids)

pbmc.TNK$cell.type=Idents(pbmc.TNK)
DimPlot(pbmc.TNK, reduction = "umap", label = TRUE, repel = T, pt.size = 0.5, label.size = 10) # + NoLegend()

DimPlot(pbmc.TNK, reduction = "umap",label = T)

  

table(pbmc.TNK$orig.ident)#?鿴????ϸ????
prop.table(table(Idents(pbmc.TNK)))
table(Idents(pbmc.TNK), pbmc.TNK$orig.ident)#???鲻ͬϸ??Ⱥϸ????
Cellratio <- prop.table(table(Idents(pbmc.TNK), pbmc.TNK$orig.ident), margin = 2)#??????????????ͬϸ??Ⱥ????
Cellratio <- data.frame(Cellratio)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#??????תΪ??????
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]


###???ӷ?????Ϣ
sample <- c("AT1.filter", "AT2.filter", "AT3.filter",
            "NC1.filter", "NC2.filter", "NC3.filter",
            "NP1.filter", "NP2.filter", "NP3.filter",
            "SLE_C.filter", "SLE_PE.filter", "SLE_PE31.filter")
group <- c("AT","AT","AT",
           "NC","NC","NC",
           "NP","NP","NP",
           "C","C","C")

samples <- data.frame(sample, group)#???????ݿ?

rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R??????
cellper$group <- samples[rownames(cellper),'group']#R??????

###??ͼչʾ
pplist = list()
sce_groups = c( "memory CD4 T",   "naive CD8 T",  "naive CD4 T",  "CD16 NK",
                "CD56 NK",  "memory CD8 T",  "regulatory CD4 T",  "gdT")

for(group_ in sce_groups){
  cellper_  = cellper %>% select(one_of(c('sample','group',group_)))#ѡ??һ??????
  colnames(cellper_) = c('sample','group','percent')#??ѡ????????????
  cellper_$percent = as.numeric(cellper_$percent)#??ֵ??????
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))#???·?λ??
  #print(group_)
  #print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot??ͼ
    geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  ###????t????????
  labely = max(cellper_$percent)
  cc<-compare_means(percent ~ group,  data = cellper_,method = "t.test")#,  p.adjust.method = "fdr")
  print(group_)
  print(cc)
  write.table(group_,"F:/10X_SLE/cellratio/TNK.cellratio.txt", row.names=FALSE,append = TRUE)
  write.table(cc, "F:/10X_SLE/cellratio/TNK.cellratio.txt", row.names=FALSE,append = TRUE)
  my_comparisons <- list( c("NC", "C") ,c("AT","NP"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test",label = "p.adj")#, p.adjust.method = "bonferroni")
  pplist[[group_]] = pp1
}
pp.TNK <- pplist





##############
## B
all.genes <- rownames(pbmc.B)
pbmc.B <- ScaleData(pbmc.B, features = all.genes)
pbmc.B <- FindVariableFeatures(object = pbmc.B)
pbmc.B <- RunPCA(pbmc.B, features = VariableFeatures(object = pbmc.B))

pbmc.B <- RunUMAP(pbmc.B, reduction = "pca", dims = 1:20)
pbmc.B <- FindNeighbors(pbmc.B, reduction = "pca", dims = 1:10)
pbmc.B <- FindClusters(pbmc.B, resolution = 0.1)

p2 <- DimPlot(pbmc.B, reduction = "umap",label = T)
p2

p1<-DotPlot(pbmc.B, features = c( "MS4A1",	 "TCL1A",		
                               "CD27",	 "IGHA1",	 "IGHG1",
                             	 "CD38"),assay = 'integrated')
p1+p2
pbmc.B.markers <- FindAllMarkers(pbmc.B, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.B.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(x=pbmc.B.markers,
          "pbmc.B.FindAllMarkers.markers.csv")


current.cluster.ids <- c(  0, 1, 2, 3)

new.cluster.ids <- c(  "memory B", "naïve B", "plasma cells", "naïve B")

names(new.cluster.ids) <- levels(pbmc.B)
pbmc.B <- RenameIdents(pbmc.B, new.cluster.ids)

pbmc.B$cell.type=Idents(pbmc.B)
DimPlot(pbmc.B, reduction = "umap", label = TRUE, repel = T, pt.size = 0.5, label.size = 10) # + NoLegend()

DimPlot(pbmc.B, reduction = "umap",label = T)



table(pbmc.B$orig.ident)#?鿴????ϸ????
prop.table(table(Idents(pbmc.B)))
table(Idents(pbmc.B), pbmc.B$orig.ident)#???鲻ͬϸ??Ⱥϸ????
Cellratio <- prop.table(table(Idents(pbmc.B), pbmc.B$orig.ident), margin = 2)#??????????????ͬϸ??Ⱥ????
Cellratio <- data.frame(Cellratio)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#??????תΪ??????
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]


###???ӷ?????Ϣ
sample <- c("AT1.filter", "AT2.filter", "AT3.filter",
            "NC1.filter", "NC2.filter", "NC3.filter",
            "NP1.filter", "NP2.filter", "NP3.filter",
            "SLE_C.filter", "SLE_PE.filter", "SLE_PE31.filter")
group <- c("AT","AT","AT",
           "NC","NC","NC",
           "NP","NP","NP",
           "C","C","C")

samples <- data.frame(sample, group)#???????ݿ?

rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R??????
cellper$group <- samples[rownames(cellper),'group']#R??????

###??ͼչʾ
pplist = list()
sce_groups = c( "memory B", "naïve B", "plasma cells", "naïve B")

for(group_ in sce_groups){
  cellper_  = cellper %>% select(one_of(c('sample','group',group_)))#ѡ??һ??????
  colnames(cellper_) = c('sample','group','percent')#??ѡ????????????
  cellper_$percent = as.numeric(cellper_$percent)#??ֵ??????
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))#???·?λ??
  #print(group_)
  #print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot??ͼ
    geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  ###????t????????
  labely = max(cellper_$percent)
  cc<-compare_means(percent ~ group,  data = cellper_,method = "t.test")#,  p.adjust.method = "fdr")
  print(group_)
  print(cc)
  write.table(group_,"F:/10X_SLE/cellratio/B.cellratio.txt", row.names=FALSE,append = TRUE)
  write.table(cc, "F:/10X_SLE/cellratio/B.cellratio.txt", row.names=FALSE,append = TRUE)
  my_comparisons <- list( c("NC", "C") ,c("AT","NP"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test",label = "p.adj")#, p.adjust.method = "bonferroni")
  pplist[[group_]] = pp1
}
pp.B <- pplist



## Monocytes
all.genes <- rownames(pbmc.Monocytes)
pbmc.Monocytes <- ScaleData(pbmc.Monocytes, features = all.genes)
pbmc.Monocytes <- FindVariableFeatures(object = pbmc.Monocytes)
pbmc.Monocytes <- RunPCA(pbmc.Monocytes, features = VariableFeatures(object = pbmc.Monocytes))
 
pbmc.Monocytes <- RunUMAP(pbmc.Monocytes, reduction = "pca", dims = 1:20)
pbmc.Monocytes <- FindNeighbors(pbmc.Monocytes, reduction = "pca", dims = 1:10)
pbmc.Monocytes <- FindClusters(pbmc.Monocytes, resolution = 0.1)

p2 <- DimPlot(pbmc.Monocytes, reduction = "umap",label = T)
p2

p1<-DotPlot(pbmc.Monocytes, features = c( "NKG7",	 "GNLY",	 "MYL9",	 "PF4",	 "GPX1", 	"LZY",
                                  "CD74",	 "S100A4",	 "FTL",	 "VIM",	 "SA100A10", 	"IFI27"),assay = 'integrated')
p1+p2

pbmc.Monocytes.markers <- FindAllMarkers(pbmc.Monocytes, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.Monocytes.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(x=pbmc.Monocytes.markers,
          "pbmc.Monocytes.FindAllMarkers.markers.csv")

current.cluster.ids <- c(  0, 1)

new.cluster.ids <- c(  "CD16 mono",  "CD14 mono")

names(new.cluster.ids) <- levels(pbmc.Monocytes)
pbmc.Monocytes <- RenameIdents(pbmc.Monocytes, new.cluster.ids)

pbmc.Monocytes$cell.type=Idents(pbmc.Monocytes)
DimPlot(pbmc.Monocytes, reduction = "umap", label = TRUE, repel = T, pt.size = 0.5, label.size = 10) # + NoLegend()

DimPlot(pbmc.Monocytes, reduction = "umap",label = T)



table(pbmc.Monocytes$orig.ident)#?鿴????ϸ????
prop.table(table(Idents(pbmc.Monocytes)))
table(Idents(pbmc.Monocytes), pbmc.Monocytes$orig.ident)#???鲻ͬϸ??Ⱥϸ????
Cellratio <- prop.table(table(Idents(pbmc.Monocytes), pbmc.Monocytes$orig.ident), margin = 2)#??????????????ͬϸ??Ⱥ????
Cellratio <- data.frame(Cellratio)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#??????תΪ??????
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]


###???ӷ?????Ϣ
sample <- c("AT1.filter", "AT2.filter", "AT3.filter",
            "NC1.filter", "NC2.filter", "NC3.filter",
            "NP1.filter", "NP2.filter", "NP3.filter",
            "SLE_C.filter", "SLE_PE.filter", "SLE_PE31.filter")
group <- c("AT","AT","AT",
           "NC","NC","NC",
           "NP","NP","NP",
           "C","C","C")

samples <- data.frame(sample, group)#???????ݿ?

rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R??????
cellper$group <- samples[rownames(cellper),'group']#R??????

###??ͼչʾ
pplist = list()
sce_groups = c( "CD16 mono",  "CD14 mono")

for(group_ in sce_groups){
  cellper_  = cellper %>% select(one_of(c('sample','group',group_)))#ѡ??һ??????
  colnames(cellper_) = c('sample','group','percent')#??ѡ????????????
  cellper_$percent = as.numeric(cellper_$percent)#??ֵ??????
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))#???·?λ??
  #print(group_)
  #print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot??ͼ
    geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  ###????t????????
  labely = max(cellper_$percent)
  cc<-compare_means(percent ~ group,  data = cellper_,method = "t.test")#,  p.adjust.method = "fdr")
  print(group_)
  print(cc)
  write.table(group_,"F:/10X_SLE/cellratio/Monocytes.cellratio.txt", row.names=FALSE,append = TRUE)
  write.table(cc, "F:/10X_SLE/cellratio/Monocytes.cellratio.txt", row.names=FALSE,append = TRUE)
  my_comparisons <- list( c("NC", "C") ,c("AT","NP"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test",label = "p.adj")#, p.adjust.method = "bonferroni")
  pplist[[group_]] = pp1
}
pp.Monocytes <- pplist



###
## all
# "CD4+ T-cells", "NK cells", "CD8+ T-cells",
# "B-cells", "Monocytes", "platelet", "platelet","pDC"

## TNK
#"memory CD4 T",   "naive CD8 T",  "naive CD4 T",  "CD16 NK",
#"CD56 NK",  "memory CD8 T",  "regulatory CD4 T",  "gdT"

## B
# "memory B", "naïve B", "plasma cells"

## Mono
#"CD16 mono",  "CD14 mono"

plot_grid(pp.TNK[['naive CD8 T']],
          pp.TNK[['memory CD8 T']],
          pp.TNK[['naive CD4 T']],
          pp.TNK[['memory CD4 T']],
          pp.TNK[['CD16 NK']],
          pp.TNK[['gdT']],
          pp.TNK[['CD56 NK']],
          pp.TNK[['regulatory CD4 T']],
          pp.B[['memory B']],
          pp.B[['naïve B']],
          pp.B[['plasma cells']],
          pp.Monocytes[['CD16 mono']],
          pp.Monocytes[['CD14 mono']],
          pp.all[['platelet']],
          pp.all[['pDC']],
          ncol=4)




