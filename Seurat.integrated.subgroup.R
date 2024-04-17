library('SeuratObject')
library('Seurat')
library(tidyverse)
library(patchwork)
library('dplyr')
library('ggplot2')


# load('./out_4group/pbmc.12sample.RNA.seurat.object.RData')
load('./out_4group/pbmc.12sample.integrated.seurat.object.RData')



# # DefaultAssay(pbmc) <- "integrated"
# # # DefaultAssay(pbmc) <- "RNA"
# # 
# # all.genes <- rownames(pbmc)
# # pbmc <- ScaleData(pbmc, features = all.genes)
# # pbmc <- FindVariableFeatures(object = pbmc)
# # pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# # 
# # pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:20)
# # pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:10)
# # pbmc <- FindClusters(pbmc, resolution = 0.2)
# #####################################################################
# 
# # p1 <- DimPlot(pbmc, reduction = "umap",group.by = "orig.ident")
# # p2 <- DimPlot(pbmc, reduction = "umap",label = T)
# 
# 
# ## new.cluster.ids
# new.cluster.ids <- c("CD4+ T-cells", "CD8+ T-cells","CD4+ T-cells", "B-cells",
#                      "B-cells", "NK cells", "Monocytes", "Monocytes",
#                      "Platelet", "Erythrocytes", "B-cells")
# names(new.cluster.ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new.cluster.ids)
# DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) # + NoLegend()
# 
# ggsave(filename = "./out_sub_0426/pbmc.integrated.rename.umap.png",
#        device = "png",width = 15, height = 10, units = "cm")
# 
# 
# 
# ################################################################
# ## subgroup
# 
## cell
pbmc.T <- subset(x = pbmc, subset = seurat_clusters == c("0","1","2"))
pbmc.B <- subset(x = pbmc, subset = seurat_clusters == c("3","4","10"))
pbmc.Monocytes <- subset(x = pbmc, subset = seurat_clusters == c("6","7"))
pbmc.NK <- subset(x = pbmc, subset = seurat_clusters == "5")
pbmc.TNK <- subset(x = pbmc, subset = seurat_clusters == c("0","1","2","5"))
# 
# 
# 
# 
# ## T
# all.genes <- rownames(pbmc.T)
# pbmc.T <- ScaleData(pbmc.T, features = all.genes)
# pbmc.T <- FindVariableFeatures(object = pbmc.T)
# pbmc.T <- RunPCA(pbmc.T, features = VariableFeatures(object = pbmc.T))
# 
# pbmc.T <- RunUMAP(pbmc.T, reduction = "pca", dims = 1:20)
# pbmc.T <- FindNeighbors(pbmc.T, reduction = "pca", dims = 1:10)
# pbmc.T <- FindClusters(pbmc.T, resolution = 0.1)
# 
# p1 <- DimPlot(pbmc.T, reduction = "umap",group.by = "new.ident")
# p2 <- DimPlot(pbmc.T, reduction = "umap",label = T)
# 
# p1+p2
# ggsave(filename = "./out_sub_0426/pbmc.T.umap.png",
#        device = "png",width = 15, height = 10, units = "cm")
# 
# 
# pbmc.T.markers <- FindAllMarkers(pbmc.T, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# write.csv(x=pbmc.T.markers, "./out_sub_0426/pbmc.T.markers.csv")
# pbmc.T.markers %>%
#     group_by(cluster) %>%
#     slice_max(n = 2, order_by = avg_log2FC)
# 
# pbmc.T.markers %>%
#     group_by(cluster) %>%
#     top_n(n = 10, wt = avg_log2FC) -> top10
# 
# pdf("./out_sub_0426/pbmc.T.integrated.heatmap.pdf",width=15,height=15)
# DoHeatmap(pbmc.T, features = top10$gene) # + NoLegend()
# dev.off()
# 
# 
# 
# ## B
# all.genes <- rownames(pbmc.B)
# pbmc.B <- ScaleData(pbmc.B, features = all.genes)
# pbmc.B <- FindVariableFeatures(object = pbmc.B)
# pbmc.B <- RunPCA(pbmc.B, features = VariableFeatures(object = pbmc.B))
# 
# pbmc.B <- RunUMAP(pbmc.B, reduction = "pca", dims = 1:20)
# pbmc.B <- FindNeighbors(pbmc.B, reduction = "pca", dims = 1:10)
# pbmc.B <- FindClusters(pbmc.B, resolution = 0.1)
# 
# p1 <- DimPlot(pbmc.B, reduction = "umap",group.by = "new.ident")
# p2 <- DimPlot(pbmc.B, reduction = "umap",label = T)
# 
# p1+p2
# ggsave(filename = "./out_sub_0426/pbmc.B.umap.png",
#        device = "png",width = 15, height = 10, units = "cm")
# 
# pbmc.B.markers <- FindAllMarkers(pbmc.B, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# write.csv(x=pbmc.B.markers, "./out_sub_0426/pbmc.B.markers.csv")
# pbmc.B.markers %>%
#     group_by(cluster) %>%
#     slice_max(n = 2, order_by = avg_log2FC)
# 
# pbmc.B.markers %>%
#     group_by(cluster) %>%
#     top_n(n = 10, wt = avg_log2FC) -> top10
# 
# pdf("./out_sub_0426/pbmc.B.integrated.heatmap.pdf",width=15,height=15)
# DoHeatmap(pbmc.B, features = top10$gene) # + NoLegend()
# dev.off()
# 
# 
# 
# ## Monocytes
# all.genes <- rownames(pbmc.Monocytes)
# pbmc.Monocytes <- ScaleData(pbmc.Monocytes, features = all.genes)
# pbmc.Monocytes <- FindVariableFeatures(object = pbmc.Monocytes)
# pbmc.Monocytes <- RunPCA(pbmc.Monocytes, features = VariableFeatures(object = pbmc.Monocytes))
# 
# pbmc.Monocytes <- RunUMAP(pbmc.Monocytes, reduction = "pca", dims = 1:20)
# pbmc.Monocytes <- FindNeighbors(pbmc.Monocytes, reduction = "pca", dims = 1:10)
# pbmc.Monocytes <- FindClusters(pbmc.Monocytes, resolution = 0.1)
# 
# p1 <- DimPlot(pbmc.Monocytes, reduction = "umap",group.by = "new.ident")
# p2 <- DimPlot(pbmc.Monocytes, reduction = "umap",label = T)
# 
# p1+p2
# ggsave(filename = "./out_sub_0426/pbmc.Monocytes.umap.png",
#        device = "png",width = 15, height = 10, units = "cm")
# 
# pbmc.Monocytes.markers <- FindAllMarkers(pbmc.Monocytes, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# write.csv(x=pbmc.Monocytes.markers, "./out_sub_0426/pbmc.Monocytes.markers.csv")
# pbmc.Monocytes.markers %>%
#     group_by(cluster) %>%
#     slice_max(n = 2, order_by = avg_log2FC)
# 
# pbmc.Monocytes.markers %>%
#     group_by(cluster) %>%
#     top_n(n = 10, wt = avg_log2FC) -> top10
# 
# pdf("./out_sub_0426/pbmc.Monocytes.integrated.heatmap.pdf",width=15,height=15)
# DoHeatmap(pbmc.Monocytes, features = top10$gene) # + NoLegend()
# dev.off()
# 
# 
# ## NK
# all.genes <- rownames(pbmc.NK)
# pbmc.NK <- ScaleData(pbmc.NK, features = all.genes)
# pbmc.NK <- FindVariableFeatures(object = pbmc.NK)
# pbmc.NK <- RunPCA(pbmc.NK, features = VariableFeatures(object = pbmc.NK))
# 
# pbmc.NK <- RunUMAP(pbmc.NK, reduction = "pca", dims = 1:20)
# pbmc.NK <- FindNeighbors(pbmc.NK, reduction = "pca", dims = 1:10)
# pbmc.NK <- FindClusters(pbmc.NK, resolution = 0.1)
# 
# p1 <- DimPlot(pbmc.NK, reduction = "umap",group.by = "new.ident")
# p2 <- DimPlot(pbmc.NK, reduction = "umap",label = T)
# 
# p1+p2
# ggsave(filename = "./out_sub_0426/pbmc.NK.umap.png",
#        device = "png",width = 15, height = 10, units = "cm")
# 
# pbmc.NK.markers <- FindAllMarkers(pbmc.NK, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# write.csv(x=pbmc.NK.markers, "./out_sub_0426/pbmc.NK.markers.csv")
# pbmc.NK.markers %>%
#     group_by(cluster) %>%
#     slice_max(n = 2, order_by = avg_log2FC)
# 
# pbmc.NK.markers %>%
#     group_by(cluster) %>%
#     top_n(n = 10, wt = avg_log2FC) -> top10
# 
# pdf("./out_sub_0426/pbmc.NK.integrated.heatmap.pdf",width=15,height=15)
# DoHeatmap(pbmc.NK, features = top10$gene) # + NoLegend()
# dev.off()


## TNK
all.genes <- rownames(pbmc.TNK)
pbmc.TNK <- ScaleData(pbmc.TNK, features = all.genes)
pbmc.TNK <- FindVariableFeatures(object = pbmc.TNK)
pbmc.TNK <- RunPCA(pbmc.TNK, features = VariableFeatures(object = pbmc.TNK))

pbmc.TNK <- RunUMAP(pbmc.TNK, reduction = "pca", dims = 1:20)
pbmc.TNK <- FindNeighbors(pbmc.TNK, reduction = "pca", dims = 1:10)
pbmc.TNK <- FindClusters(pbmc.TNK, resolution = 0.8)

p1 <- DimPlot(pbmc.TNK, reduction = "umap",group.by = "new.ident")
p2 <- DimPlot(pbmc.TNK, reduction = "umap",label = T)

p1+p2
ggsave(filename = "./out_sub_0426_TNK/pbmc.TNK.umap.png",
       device = "png",width = 15, height = 10, units = "cm")

pbmc.TNK.markers <- FindAllMarkers(pbmc.TNK, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(x=pbmc.TNK.markers, "./out_sub_0426_TNK/pbmc.TNK.markers.csv")
pbmc.TNK.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

pbmc.TNK.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

pdf("./out_sub_0426_TNK/pbmc.TNK.integrated.heatmap.pdf",width=15,height=15)
DoHeatmap(pbmc.TNK, features = top10$gene) # + NoLegend()
dev.off()





# save(pbmc.T, file='./out_sub_0426_TNK/pbmc.T.integrated.seurat.object.RData')
# save(pbmc.B, file='./out_sub_0426_TNK/pbmc.B.integrated.seurat.object.RData')
# save(pbmc.Monocytes, file='./out_sub_0426_TNK/pbmc.Monocytes.integrated.seurat.object.RData')
# save(pbmc.NK, file='./out_sub_0426_TNK/pbmc.NK.integrated.seurat.object.RData')
save(pbmc.TNK, file='./out_sub_0426_TNK/pbmc.TNK.integrated.seurat.object.RData')



