library('SeuratObject')
library('Seurat')
library(tidyverse)
library(patchwork)
library('dplyr')
library('ggplot2')



# load('out_4group/pbmc.pbmc_NP.integrated.seurat.object.RData')
# load('out_4group/pbmc.pbmc_NC.integrated.seurat.object.RData')
# load('out_4group/pbmc.pbmc_AT.integrated.seurat.object.RData')
# load('out_4group/pbmc.pbmc_C.integrated.seurat.object.RData')


load('/share/RD/PB_RD/10X/SLE/SLE_seurat/out/NP123.RP.seurat.object.RData')
load('/share/RD/PB_RD/10X/SLE/SLE_seurat/out/NC123.RP.seurat.object.RData')
load('/share/RD/PB_RD/10X/SLE/SLE_seurat/out/AT123.RP.seurat.object.RData')
load('/share/RD/PB_RD/10X/SLE/SLE_seurat/out/SLE_C123.RP.seurat.object.RData')


load('/share/RD/PB_RD/10X/SLE/SLE_seurat/4group/out_4group/pbmc.12sample.integrated.seurat.object.RData')

#####################################################################


###################################################################
###################################################################
pbmc.list <- list(pbmc_NP,pbmc_NC,pbmc_AT,pbmc_C)
pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, dims = 1:20)
pbmc <- IntegrateData(anchorset=pbmc.anchors, dims = 1:20)

DefaultAssay(pbmc) <- "integrated"
# DefaultAssay(pbmc) <- "RNA"
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:20)
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
#####################################################################

p1 <- DimPlot(pbmc, reduction = "umap",group.by = "orig.ident")
p2 <- DimPlot(pbmc, reduction = "umap",label = T)

p1+p2
ggsave(filename = "./pbmc.12sample.integrated.umap.png",
       device = "png",width = 50,height = 30,units = "cm")

save(pbmc, file='./pbmc.12sample.integrated.seurat.object.RData')



# DotPlot(pbmc, features = c( "CD19", "CD79A", "CD79B","MS4A1", 
#                            "CD3D", "CD3E", "CD3G",
#                            "IL7R", "CD8A", "CD8B", "CD4",
#                            "HLA-DQA1", "HLA-DQB1", "CLEC10A",
#                            "CD14", "VCAN", "FCN1",
#                            "GP9","PF4"))
# 
# ggsave(filename = "./out_4group/pbmc.12sample.integrated.dotplot.png",
#        device = "png",width = 50,height = 15,units = "cm")


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

write.csv(x=pbmc.markers,
          "./pbmc.12sample.integrated.FindAllMarkers.markers.csv")

 
