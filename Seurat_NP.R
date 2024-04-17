library(dplyr)
library(Seurat)
library(patchwork)

library(DoubletFinder)


suppressMessages(require(DoubletFinder))



NP1.data <- Read10X(data.dir = "/share/RD/PB_RD/10X/SLE/SLE_matrix/NP-1/")
NP2.data <- Read10X(data.dir = "/share/RD/PB_RD/10X/SLE/SLE_matrix/NP-2/")
NP3.data <- Read10X(data.dir = "/share/RD/PB_RD/10X/SLE/SLE_matrix/NP-3/")

NP1 <- CreateSeuratObject(counts = NP1.data, project = "NP1", min.cells = 3, min.features = 200) 
NP2 <- CreateSeuratObject(counts = NP2.data, project = "NP2", min.cells = 3, min.features = 200) 
NP3 <- CreateSeuratObject(counts = NP3.data, project = "NP3", min.cells = 3, min.features = 200) 

##doublet
nExp1 <- round(ncol(NP1) * 0.038)
nExp2 <- round(ncol(NP2) * 0.051)
nExp3 <- round(ncol(NP3) * 0.052)



### QC and normalize, preprocessing
################################################################################
# for NP1 data
################################################################################
NP1[["percent.mt"]] <- PercentageFeatureSet(NP1, pattern = "^MT-")
VlnPlot(NP1, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
NP1 <- subset(NP1, subset = nFeature_RNA >200 & nCount_RNA >2000 & nCount_RNA<60000 & percent.mt < 10 )
NP1 <- NormalizeData(NP1, normalization.method = "LogNormalize", scale.factor = 10000)
NP1 <- FindVariableFeatures(NP1, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(NP1)
NP1 <- ScaleData(NP1, features = all.genes)

NP1 <- RunPCA(NP1, features = VariableFeatures(object = NP1))
VizDimLoadings(NP1, dims = 1:2, reduction = "pca")
DimPlot(NP1, reduction = "pca")

NP1 <- FindNeighbors(NP1, dims = 1:20)
NP1 <- FindClusters(NP1, resolution = 0.5)

NP1 <- RunUMAP(NP1, dims = 1:20)
DimPlot(NP1, reduction = "umap")

NP1 <- doubletFinder_v3(NP1, pN = 0.25, pK = 0.09, nExp = nExp1, PCs = 1:10)
DF1.name = colnames(NP1@meta.data)[grepl("DF.classification", colnames(NP1@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(NP1, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(NP1, group.by = DF1.name) + NoAxes())
VlnPlot(NP1, features = "nFeature_RNA", group.by = DF1.name, pt.size = 0.1)

NP1 = NP1[, NP1@meta.data[, DF1.name] == "Singlet"]
dim(NP1)



################################################################################
# for NP2 data
################################################################################
NP2[["percent.mt"]] <- PercentageFeatureSet(NP2, pattern = "^MT-")
VlnPlot(NP2, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
# NP2 <- subset(NP2, subset = nFeature_RNA < 3000 & percent.mt < 10 & percent.HB < 10)
NP2 <- subset(NP2, subset = nFeature_RNA >200 & nCount_RNA >2000 & nCount_RNA<60000 & percent.mt < 10 )

NP2 <- NormalizeData(NP2, normalization.method = "LogNormalize", scale.factor = 10000)
NP2 <- FindVariableFeatures(NP2, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(NP2)
NP2 <- ScaleData(NP2, features = all.genes)

NP2 <- RunPCA(NP2, features = VariableFeatures(object = NP2))
VizDimLoadings(NP2, dims = 1:2, reduction = "pca")
DimPlot(NP2, reduction = "pca")

NP2 <- FindNeighbors(NP2, dims = 1:20)
NP2 <- FindClusters(NP2, resolution = 0.5)

NP2 <- RunUMAP(NP2, dims = 1:20)
DimPlot(NP2, reduction = "umap")

NP2 <- doubletFinder_v3(NP2, pN = 0.25, pK = 0.09, nExp = nExp2, PCs = 1:10)

DF2.name = colnames(NP2@meta.data)[grepl("DF.classification", colnames(NP2@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(NP2, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(NP2, group.by = DF2.name) + NoAxes())

VlnPlot(NP2, features = "nFeature_RNA", group.by = DF2.name, pt.size = 0.1)

NP2 = NP2[, NP2@meta.data[, DF2.name] == "Singlet"]
dim(NP2)


# for NP3 data
NP3[["percent.mt"]] <- PercentageFeatureSet(NP3, pattern = "^MT-")
VlnPlot(NP3, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
NP3 <- subset(NP3, subset = nFeature_RNA >200 & nCount_RNA >2000 & nCount_RNA<60000 & percent.mt < 10 )

NP3 <- NormalizeData(NP3, normalization.method = "LogNormalize", scale.factor = 10000)
NP3 <- FindVariableFeatures(NP3, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(NP3)
NP3 <- ScaleData(NP3, features = all.genes)

NP3 <- RunPCA(NP3, features = VariableFeatures(object = NP3))
VizDimLoadings(NP3, dims = 1:2, reduction = "pca")
DimPlot(NP3, reduction = "pca")

NP3 <- FindNeighbors(NP3, dims = 1:20)
NP3 <- FindClusters(NP3, resolution = 0.5)

NP3 <- RunUMAP(NP3, dims = 1:20)
DimPlot(NP3, reduction = "umap")

# define the expected number of doublet cellscells.
NP3 <- doubletFinder_v3(NP3, pN = 0.25, pK = 0.09, nExp = nExp3, PCs = 1:10)

# name of the DF prediction can change, so extract the correct column name.
DF3.name = colnames(NP3@meta.data)[grepl("DF.classification", colnames(NP3@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(NP3, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(NP3, group.by = DF3.name) + NoAxes())

VlnPlot(NP3, features = "nFeature_RNA", group.by = DF3.name, pt.size = 0.1)

NP3 = NP3[, NP3@meta.data[, DF3.name] == "Singlet"]
dim(NP3)
####
####

save(NP1,NP2,NP3, file='/share/RD/PB_RD/10X/SLE/SLE_seurat/out/NP123.seurat.object.RData')



