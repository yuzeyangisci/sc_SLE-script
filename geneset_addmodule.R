library(patchwork)
library(ggplot2)
library(Seurat)
library(dplyr)

library(ggpubr)
library(ggplot2)
library(cowplot)
#library(ggthemr)     #???????ڽ??ܵ????????ð?
#library(ggsignif)    #????ggsignif


geneset <-  read.csv('F:/E-bak/WORK/PB/10X/gene_set.csv',header = T)
geneset.matrix <- as.matrix(geneset)


########
########
load('F:/E-bak/WORK/PB/10X/SLE/sub/out_sub_0426_TNK/pbmc.TNK.integrated.seurat.object.RData')
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
new.cluster.ids <- c("naive CD8 T",
                     "memory CD8 T",
                     "naive CD4 T",
                     "memory CD4 T",
                     "CD16 NK",
                     "memory CD8 T",
                     "naive CD4 T",
                     "memory CD4 T",
                     "gdT",
                     "memory CD8 T",
                     "naive CD4 T",
                     "CD56 NK",
                     "regulatory CD4 T")
names(new.cluster.ids) <- levels(pbmc.TNK)
pbmc.TNK <- RenameIdents(pbmc.TNK, new.cluster.ids)

pbmc.TNK$cell.type=Idents(pbmc.TNK)
DimPlot(pbmc.TNK, reduction = "umap", label = TRUE, repel = T, pt.size = 0.5, label.size = 10) # + NoLegend()



############
## TNK: M5.6,M3.4,M1.2,M3.3,M4.2,M4.7, M4.15, M3.6, M4.1
##
geneset.TNK <- c("M5.6","M3.4","M1.2","M3.3","M4.2","M4.7", "M4.15", "M3.6", "M4.1")
pplist = list()
for(gene.set in geneset.TNK){

  gene.TNK <- list()
  gene.TNK <- unlist(as.list(geneset.matrix[,gene.set]))
  gene.TNK <- gene.TNK[gene.TNK!=""]  
  gene.TNK <- list(gene.TNK)
  DefaultAssay(pbmc.TNK)<- "RNA"
  pbmc.TNK <- AddModuleScore(object = pbmc.TNK,
                             features = gene.TNK,
                             name = "TNK.AddModuleScore")
  VlnPlot(pbmc.TNK,features = "TNK.AddModuleScore1",split.by = "new.ident") #+ geom_boxplot()
  
  mydata=FetchData(pbmc.TNK,
                   vars = c("UMAP_1","UMAP_2","new.ident","cell.type","TNK.AddModuleScore1"))
  head(mydata)
  compare_means(TNK.AddModuleScore1 ~ new.ident,  data = mydata, method = "t.test")
  pp1 <- mydata %>%
    mutate(new.ident = factor(new.ident, levels = c("C","NC","AT","NP"))) %>%
    ggplot(aes(x=new.ident, y=TNK.AddModuleScore1,colour=new.ident,fill=new.ident)) +
    geom_boxplot(outlier.shape = NA)+
    facet_grid(~cell.type) +RotatedAxis()
  my_comparisons <- list( c("NC", "C") ,c("AT","NP"))
  
  head(mydata)
  pp1 <- pp1 + stat_compare_means( comparisons = my_comparisons, label = "p.signif")
  
  pp1 <- pp1 + labs( y = "T_NK cellular status", title = gene.set) 
  pplist[[gene.set]] = pp1
  
}

# M5.6,M3.4,M1.2,M3.3,M4.2,M4.7, M4.15, M3.6, M4.1

pdf('TNK.pdf',height = 6, width = 15)

plot(pplist[['M5.6']])
plot(pplist[['M3.4']])
plot(pplist[['M1.2']])
plot(pplist[['M3.3']])
plot(pplist[['M4.2']])
plot(pplist[['M4.7']])
plot(pplist[['M4.15']])
plot(pplist[['M3.6']])
plot(pplist[['M4.1']])

dev.off()






##B
########
########
load('F:/E-bak/WORK/PB/10X/SLE/sub/pbmc.B.integrated.seurat.object.RData')

current.cluster.ids <- c(0, 1, 2, 3, 4)
new.cluster.ids <- c("naive_B", "memory_B", "naive_B", "memory_B", "plasma_cell_B")
names(new.cluster.ids) <- levels(pbmc.B)
pbmc.B <- RenameIdents(pbmc.B, new.cluster.ids)
pbmc.B$cell.type=Idents(pbmc.B)

DimPlot(pbmc.B, reduction = "umap", label = TRUE, repel = T, pt.size = 0.5, label.size = 10) # + NoLegend()



############
## B:M5.6,M3.4,M1.2,M3.3,M4.2,M4.7, M4.10, M4.11
## TNK: M5.6,M3.4,M1.2,M3.3,M4.2,M4.7, M4.15, M3.6, M4.1
##
geneset.B <- c("M5.6","M3.4","M1.2","M3.3","M4.2","M4.7", "M4.10", "M4.11")
pplist = list()
for(gene.set in geneset.B){
  
  gene.B <- list()
  gene.B <- unlist(as.list(geneset.matrix[,gene.set]))
  gene.B <- gene.B[gene.B!=""]  
  gene.B <- list(gene.B)
  DefaultAssay(pbmc.B)<- "RNA"
  pbmc.B <- AddModuleScore(object = pbmc.B,
                             features = gene.B,
                             name = "B.AddModuleScore")
  VlnPlot(pbmc.B,features = "B.AddModuleScore1",split.by = "new.ident") #+ geom_boxplot()
  
  mydata=FetchData(pbmc.B,
                   vars = c("UMAP_1","UMAP_2","new.ident","cell.type","B.AddModuleScore1"))
  head(mydata)
  cc<-compare_means(B.AddModuleScore1 ~ new.ident,  data = mydata, group.by = "cell.type")#method = "t.test")
  pp1 <- mydata %>%
    mutate(new.ident = factor(new.ident, levels = c("C","NC","AT","NP"))) %>%
    ggplot(aes(x=new.ident, y=B.AddModuleScore1,colour=new.ident,fill=new.ident)) +
    geom_boxplot(outlier.shape = NA)+
    facet_grid(~cell.type) +RotatedAxis()
  my_comparisons <- list( c("NC", "C") ,c("AT","NP"))
  
  head(mydata)
  pp1 <- pp1 + stat_compare_means( comparisons = my_comparisons, label = "p.signif")
  
  pp1 <- pp1 + labs( y = "B cellular status", title = gene.set) 
  pplist[[gene.set]] = pp1

  print(gene.set)
  print(cc)
  write.table(gene.set,"F:/10X_SLE/geneset/B.geneset.txt", row.names=FALSE,append = TRUE)
  write.table(cc, "F:/10X_SLE/geneset/B.geneset.txt", row.names=FALSE,append = TRUE)
  
}



pdf('B.pdf',height = 6, width = 15)

# M5.6,M3.4,M1.2,M3.3,M4.2,M4.7, M4.10, M4.11

plot(pplist[['M5.6']])
plot(pplist[['M3.4']])
plot(pplist[['M1.2']])
plot(pplist[['M3.3']])
plot(pplist[['M4.2']])
plot(pplist[['M4.7']])
plot(pplist[['M4.10']])
plot(pplist[['M4.11']])

dev.off()






## MOno
load('F:/E-bak/WORK/PB/10X/SLE/sub/pbmc.Monocytes.integrated.seurat.object.RData')

current.cluster.ids <- c(0, 1, 2)
new.cluster.ids <- c("CD16_monocyte", "CD14_monocyte", "CD14_monocyte")
names(new.cluster.ids) <- levels(pbmc.Monocytes)
pbmc.Monocytes <- RenameIdents(pbmc.Monocytes, new.cluster.ids)
pbmc.Monocytes$cell.type=Idents(pbmc.Monocytes)

DimPlot(pbmc.Monocytes, reduction = "umap", label = TRUE, repel = T, pt.size = 0.5, label.size = 10) # + NoLegend()


############
## Mono:M5.6,M3.4,M1.2,M3.3,M4.2,M3.2,M7.27
geneset.Mono <- c("M5.6","M3.4","M1.2","M3.3","M4.2","M3.2","M7.27")

pplist = list()
for(gene.set in geneset.Mono){
  
  gene.Mono <- list()
  gene.Mono <- unlist(as.list(geneset.matrix[,gene.set]))
  gene.Mono <- gene.Mono[gene.Mono!=""]  
  gene.Mono <- list(gene.Mono)
  DefaultAssay(pbmc.Monocytes)<- "RNA"
  pbmc.Monocytes <- AddModuleScore(object = pbmc.Monocytes,
                             features = gene.Mono,
                             name = "Mono.AddModuleScore")
  VlnPlot(pbmc.Monocytes,features = "Mono.AddModuleScore1",split.by = "new.ident") #+ geom_boxplot()
  
  mydata=FetchData(pbmc.Monocytes,
                   vars = c("UMAP_1","UMAP_2","new.ident","cell.type","Mono.AddModuleScore1"))
  head(mydata)
  cc<-compare_means(Mono.AddModuleScore1 ~ new.ident,  data = mydata, group.by = "cell.type")# method = "t.test")
  pp1 <- mydata %>%
    mutate(new.ident = factor(new.ident, levels = c("C","NC","AT","NP"))) %>%
    ggplot(aes(x=new.ident, y=Mono.AddModuleScore1,colour=new.ident,fill=new.ident)) +
    geom_boxplot(outlier.shape = NA)+
    facet_grid(~cell.type) +RotatedAxis()
  my_comparisons <- list( c("NC", "C") ,c("AT","NP"))
  
  head(mydata)
  pp1 <- pp1 + stat_compare_means( comparisons = my_comparisons, label = "p.signif")
  
  pp1 <- pp1 + labs( y = "Monocytes cellular status", title = gene.set) 
  pplist[[gene.set]] = pp1

  print(gene.set)
  print(cc)
  write.table(gene.set,"F:/10X_SLE/geneset/Mono.geneset.txt", row.names=FALSE,append = TRUE)
  write.table(cc, "F:/10X_SLE/geneset/Mono.geneset.txt", row.names=FALSE,append = TRUE)
  
}

pdf('Mono.pdf',height = 6, width = 15)

# M5.6,M3.4,M1.2,M3.3,M4.2,M3.2,M7.27
plot(pplist[['M5.6']])
plot(pplist[['M3.4']])
plot(pplist[['M1.2']])
plot(pplist[['M3.3']])
plot(pplist[['M4.2']])
plot(pplist[['M3.2']])
plot(pplist[['M7.27']])

dev.off()
