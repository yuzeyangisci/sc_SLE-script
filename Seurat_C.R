library(dplyr)
library(Seurat)
library(patchwork)

library(DoubletFinder)


suppressMessages(require(DoubletFinder))



SLE_C1.data <- Read10X(data.dir = "/share/RD/PB_RD/10X/SLE/SLE_matrix/SLE_C")
SLE_C2.data <- Read10X(data.dir = "/share/RD/PB_RD/10X/SLE/SLE_matrix/SLE-PE")
SLE_C3.data <- Read10X(data.dir = "/share/RD/PB_RD/10X/SLE/SLE_matrix/PE31")

SLE_C1 <- CreateSeuratObject(counts = SLE_C1.data, project = "SLE-C", min.cells = 3, min.features = 200) 
SLE_C2 <- CreateSeuratObject(counts = SLE_C2.data, project = "SLE-PE", min.cells = 3, min.features = 200) 
SLE_C3 <- CreateSeuratObject(counts = SLE_C3.data, project = "SLE-PE31", min.cells = 3, min.features = 200) 

##doublet
nExp1 <- round(ncol(SLE_C1) * 0.059)
nExp2 <- round(ncol(SLE_C2) * 0.036)
nExp3 <- round(ncol(SLE_C3) * 0.079)


### QC and normalize, preprocessing
################################################################################
# for SLE_C1 data
################################################################################
SLE_C1[["percent.mt"]] <- PercentageFeatureSet(SLE_C1, pattern = "^MT-")
VlnPlot(SLE_C1, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
SLE_C1 <- subset(SLE_C1, subset = nFeature_RNA >200 & nCount_RNA >2000 & nCount_RNA<60000 & percent.mt < 10 )
SLE_C1 <- NormalizeData(SLE_C1, normalization.method = "LogNormalize", scale.factor = 10000)
SLE_C1 <- FindVariableFeatures(SLE_C1, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(SLE_C1)
SLE_C1 <- ScaleData(SLE_C1, features = all.genes)

SLE_C1 <- RunPCA(SLE_C1, features = VariableFeatures(object = SLE_C1))
VizDimLoadings(SLE_C1, dims = 1:2, reduction = "pca")
DimPlot(SLE_C1, reduction = "pca")

SLE_C1 <- FindNeighbors(SLE_C1, dims = 1:20)
SLE_C1 <- FindClusters(SLE_C1, resolution = 0.5)

SLE_C1 <- RunUMAP(SLE_C1, dims = 1:20)
DimPlot(SLE_C1, reduction = "umap")

SLE_C1 <- doubletFinder_v3(SLE_C1, pN = 0.25, pK = 0.09, nExp = nExp1, PCs = 1:10)
DF1.name = colnames(SLE_C1@meta.data)[grepl("DF.classification", colnames(SLE_C1@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(SLE_C1, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(SLE_C1, group.by = DF1.name) + NoAxes())
VlnPlot(SLE_C1, features = "nFeature_RNA", group.by = DF1.name, pt.size = 0.1)

SLE_C1 = SLE_C1[, SLE_C1@meta.data[, DF1.name] == "Singlet"]
dim(SLE_C1)



################################################################################
# for SLE_C2 data
################################################################################
SLE_C2[["percent.mt"]] <- PercentageFeatureSet(SLE_C2, pattern = "^MT-")
VlnPlot(SLE_C2, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
# SLE_C2 <- subset(SLE_C2, subset = nFeature_RNA < 3000 & percent.mt < 10 & percent.HB < 10)
SLE_C2 <- subset(SLE_C2, subset = nFeature_RNA >200 & nCount_RNA >2000 & nCount_RNA<60000 & percent.mt < 10 )

SLE_C2 <- NormalizeData(SLE_C2, normalization.method = "LogNormalize", scale.factor = 10000)
SLE_C2 <- FindVariableFeatures(SLE_C2, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(SLE_C2)
SLE_C2 <- ScaleData(SLE_C2, features = all.genes)

SLE_C2 <- RunPCA(SLE_C2, features = VariableFeatures(object = SLE_C2))
VizDimLoadings(SLE_C2, dims = 1:2, reduction = "pca")
DimPlot(SLE_C2, reduction = "pca")

SLE_C2 <- FindNeighbors(SLE_C2, dims = 1:20)
SLE_C2 <- FindClusters(SLE_C2, resolution = 0.5)

SLE_C2 <- RunUMAP(SLE_C2, dims = 1:20)
DimPlot(SLE_C2, reduction = "umap")

SLE_C2 <- doubletFinder_v3(SLE_C2, pN = 0.25, pK = 0.09, nExp = nExp2, PCs = 1:10)

DF2.name = colnames(SLE_C2@meta.data)[grepl("DF.classification", colnames(SLE_C2@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(SLE_C2, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(SLE_C2, group.by = DF2.name) + NoAxes())

VlnPlot(SLE_C2, features = "nFeature_RNA", group.by = DF2.name, pt.size = 0.1)

SLE_C2 = SLE_C2[, SLE_C2@meta.data[, DF2.name] == "Singlet"]
dim(SLE_C2)


# for SLE_C3 data
SLE_C3[["percent.mt"]] <- PercentageFeatureSet(SLE_C3, pattern = "^MT-")
VlnPlot(SLE_C3, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
SLE_C3 <- subset(SLE_C3, subset = nFeature_RNA >200 & nCount_RNA >2000 & nCount_RNA<60000 & percent.mt < 10 )

SLE_C3 <- NormalizeData(SLE_C3, normalization.method = "LogNormalize", scale.factor = 10000)
SLE_C3 <- FindVariableFeatures(SLE_C3, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(SLE_C3)
SLE_C3 <- ScaleData(SLE_C3, features = all.genes)

SLE_C3 <- RunPCA(SLE_C3, features = VariableFeatures(object = SLE_C3))
VizDimLoadings(SLE_C3, dims = 1:2, reduction = "pca")
DimPlot(SLE_C3, reduction = "pca")

SLE_C3 <- FindNeighbors(SLE_C3, dims = 1:20)
SLE_C3 <- FindClusters(SLE_C3, resolution = 0.5)

SLE_C3 <- RunUMAP(SLE_C3, dims = 1:20)
DimPlot(SLE_C3, reduction = "umap")

# define the expected number of doublet cellscells.
SLE_C3 <- doubletFinder_v3(SLE_C3, pN = 0.25, pK = 0.09, nExp = nExp3, PCs = 1:10)

# name of the DF prediction can change, so extract the correct column name.
DF3.name = colnames(SLE_C3@meta.data)[grepl("DF.classification", colnames(SLE_C3@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(SLE_C3, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(SLE_C3, group.by = DF3.name) + NoAxes())

VlnPlot(SLE_C3, features = "nFeature_RNA", group.by = DF3.name, pt.size = 0.1)

SLE_C3 = SLE_C3[, SLE_C3@meta.data[, DF3.name] == "Singlet"]
dim(SLE_C3)
####
####

save(SLE_C1,SLE_C2,SLE_C3, file='/share/RD/PB_RD/10X/SLE/SLE_seurat/out/SLE_C123.seurat.object.RData')



