library(dplyr)
library(Seurat)
library(patchwork)

library(DoubletFinder)


suppressMessages(require(DoubletFinder))



NC1.data <- Read10X(data.dir = "/share/RD/PB_RD/10X/SLE/SLE_matrix/NC-1/")
NC2.data <- Read10X(data.dir = "/share/RD/PB_RD/10X/SLE/SLE_matrix/NC-2/")
NC3.data <- Read10X(data.dir = "/share/RD/PB_RD/10X/SLE/SLE_matrix/NC-3/")

NC1 <- CreateSeuratObject(counts = NC1.data, project = "NC1", min.cells = 3, min.features = 200) 
NC2 <- CreateSeuratObject(counts = NC2.data, project = "NC2", min.cells = 3, min.features = 200) 
NC3 <- CreateSeuratObject(counts = NC3.data, project = "NC3", min.cells = 3, min.features = 200) 

##doublet
nExp1 <- round(ncol(NC1) * 0.046)
nExp2 <- round(ncol(NC2) * 0.050)
nExp3 <- round(ncol(NC3) * 0.060)


### QC and normalize, preprocessing
################################################################################
# for NC1 data
################################################################################
NC1[["percent.mt"]] <- PercentageFeatureSet(NC1, pattern = "^MT-")
VlnPlot(NC1, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
NC1 <- subset(NC1, subset = nFeature_RNA >200 & nCount_RNA >2000 & nCount_RNA<60000 & percent.mt < 10 )
NC1 <- NormalizeData(NC1, normalization.method = "LogNormalize", scale.factor = 10000)
NC1 <- FindVariableFeatures(NC1, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(NC1)
NC1 <- ScaleData(NC1, features = all.genes)

NC1 <- RunPCA(NC1, features = VariableFeatures(object = NC1))
VizDimLoadings(NC1, dims = 1:2, reduction = "pca")
DimPlot(NC1, reduction = "pca")

NC1 <- FindNeighbors(NC1, dims = 1:20)
NC1 <- FindClusters(NC1, resolution = 0.5)

NC1 <- RunUMAP(NC1, dims = 1:20)
DimPlot(NC1, reduction = "umap")

NC1 <- doubletFinder_v3(NC1, pN = 0.25, pK = 0.09, nExp = nExp1, PCs = 1:10)
DF1.name = colnames(NC1@meta.data)[grepl("DF.classification", colnames(NC1@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(NC1, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(NC1, group.by = DF1.name) + NoAxes())
VlnPlot(NC1, features = "nFeature_RNA", group.by = DF1.name, pt.size = 0.1)

NC1 = NC1[, NC1@meta.data[, DF1.name] == "Singlet"]
dim(NC1)



################################################################################
# for NC2 data
################################################################################
NC2[["percent.mt"]] <- PercentageFeatureSet(NC2, pattern = "^MT-")
VlnPlot(NC2, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
# NC2 <- subset(NC2, subset = nFeature_RNA < 3000 & percent.mt < 10 & percent.HB < 10)
NC2 <- subset(NC2, subset = nFeature_RNA >200 & nCount_RNA >2000 & nCount_RNA<60000 & percent.mt < 10 )

NC2 <- NormalizeData(NC2, normalization.method = "LogNormalize", scale.factor = 10000)
NC2 <- FindVariableFeatures(NC2, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(NC2)
NC2 <- ScaleData(NC2, features = all.genes)

NC2 <- RunPCA(NC2, features = VariableFeatures(object = NC2))
VizDimLoadings(NC2, dims = 1:2, reduction = "pca")
DimPlot(NC2, reduction = "pca")

NC2 <- FindNeighbors(NC2, dims = 1:20)
NC2 <- FindClusters(NC2, resolution = 0.5)

NC2 <- RunUMAP(NC2, dims = 1:20)
DimPlot(NC2, reduction = "umap")

NC2 <- doubletFinder_v3(NC2, pN = 0.25, pK = 0.09, nExp = nExp2, PCs = 1:10)

DF2.name = colnames(NC2@meta.data)[grepl("DF.classification", colnames(NC2@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(NC2, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(NC2, group.by = DF2.name) + NoAxes())

VlnPlot(NC2, features = "nFeature_RNA", group.by = DF2.name, pt.size = 0.1)

NC2 = NC2[, NC2@meta.data[, DF2.name] == "Singlet"]
dim(NC2)


# for NC3 data
NC3[["percent.mt"]] <- PercentageFeatureSet(NC3, pattern = "^MT-")
VlnPlot(NC3, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
NC3 <- subset(NC3, subset = nFeature_RNA >200 & nCount_RNA >2000 & nCount_RNA<60000 & percent.mt < 10 )

NC3 <- NormalizeData(NC3, normalization.method = "LogNormalize", scale.factor = 10000)
NC3 <- FindVariableFeatures(NC3, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(NC3)
NC3 <- ScaleData(NC3, features = all.genes)

NC3 <- RunPCA(NC3, features = VariableFeatures(object = NC3))
VizDimLoadings(NC3, dims = 1:2, reduction = "pca")
DimPlot(NC3, reduction = "pca")

NC3 <- FindNeighbors(NC3, dims = 1:20)
NC3 <- FindClusters(NC3, resolution = 0.5)

NC3 <- RunUMAP(NC3, dims = 1:20)
DimPlot(NC3, reduction = "umap")

# define the expected number of doublet cellscells.
NC3 <- doubletFinder_v3(NC3, pN = 0.25, pK = 0.09, nExp = nExp3, PCs = 1:10)

# name of the DF prediction can change, so extract the correct column name.
DF3.name = colnames(NC3@meta.data)[grepl("DF.classification", colnames(NC3@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(NC3, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(NC3, group.by = DF3.name) + NoAxes())

VlnPlot(NC3, features = "nFeature_RNA", group.by = DF3.name, pt.size = 0.1)

NC3 = NC3[, NC3@meta.data[, DF3.name] == "Singlet"]
dim(NC3)
####
####

save(NC1,NC2,NC3, file='/share/RD/PB_RD/10X/SLE/SLE_seurat/out/NC123.seurat.object.RData')



