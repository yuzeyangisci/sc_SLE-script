library(dplyr)
library(Seurat)
library(patchwork)

library(DoubletFinder)


suppressMessages(require(DoubletFinder))



AT1.data <- Read10X(data.dir = "/share/RD/PB_RD/10X/SLE/SLE_matrix/AT-1/")
AT2.data <- Read10X(data.dir = "/share/RD/PB_RD/10X/SLE/SLE_matrix/AT-2/")
AT3.data <- Read10X(data.dir = "/share/RD/PB_RD/10X/SLE/SLE_matrix/AT-3/")

AT1 <- CreateSeuratObject(counts = AT1.data, project = "AT1", min.cells = 3, min.features = 200) 
AT2 <- CreateSeuratObject(counts = AT2.data, project = "AT2", min.cells = 3, min.features = 200) 
AT3 <- CreateSeuratObject(counts = AT3.data, project = "AT3", min.cells = 3, min.features = 200) 

##doublet
nExp1 <- round(ncol(AT1) * 0.054)
nExp2 <- round(ncol(AT2) * 0.079)
nExp3 <- round(ncol(AT3) * 0.113)


### QC and normalize, preprocessing
################################################################################
# for AT1 data
################################################################################
AT1[["percent.mt"]] <- PercentageFeatureSet(AT1, pattern = "^MT-")
VlnPlot(AT1, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
AT1 <- subset(AT1, subset = nFeature_RNA >200 & nCount_RNA >2000 & nCount_RNA<60000 & percent.mt < 10 )
AT1 <- NormalizeData(AT1, normalization.method = "LogNormalize", scale.factor = 10000)
AT1 <- FindVariableFeatures(AT1, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(AT1)
AT1 <- ScaleData(AT1, features = all.genes)

AT1 <- RunPCA(AT1, features = VariableFeatures(object = AT1))
VizDimLoadings(AT1, dims = 1:2, reduction = "pca")
DimPlot(AT1, reduction = "pca")

AT1 <- FindNeighbors(AT1, dims = 1:20)
AT1 <- FindClusters(AT1, resolution = 0.5)

AT1 <- RunUMAP(AT1, dims = 1:20)
DimPlot(AT1, reduction = "umap")

AT1 <- doubletFinder_v3(AT1, pN = 0.25, pK = 0.09, nExp = nExp1, PCs = 1:10)
DF1.name = colnames(AT1@meta.data)[grepl("DF.classification", colnames(AT1@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(AT1, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(AT1, group.by = DF1.name) + NoAxes())
VlnPlot(AT1, features = "nFeature_RNA", group.by = DF1.name, pt.size = 0.1)

AT1 = AT1[, AT1@meta.data[, DF1.name] == "Singlet"]
dim(AT1)



################################################################################
# for AT2 data
################################################################################
AT2[["percent.mt"]] <- PercentageFeatureSet(AT2, pattern = "^MT-")
VlnPlot(AT2, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
# AT2 <- subset(AT2, subset = nFeature_RNA < 3000 & percent.mt < 10 & percent.HB < 10)
AT2 <- subset(AT2, subset = nFeature_RNA >200 & nCount_RNA >2000 & nCount_RNA<60000 & percent.mt < 10 )

AT2 <- NormalizeData(AT2, normalization.method = "LogNormalize", scale.factor = 10000)
AT2 <- FindVariableFeatures(AT2, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(AT2)
AT2 <- ScaleData(AT2, features = all.genes)

AT2 <- RunPCA(AT2, features = VariableFeatures(object = AT2))
VizDimLoadings(AT2, dims = 1:2, reduction = "pca")
DimPlot(AT2, reduction = "pca")

AT2 <- FindNeighbors(AT2, dims = 1:20)
AT2 <- FindClusters(AT2, resolution = 0.5)

AT2 <- RunUMAP(AT2, dims = 1:20)
DimPlot(AT2, reduction = "umap")

AT2 <- doubletFinder_v3(AT2, pN = 0.25, pK = 0.09, nExp = nExp2, PCs = 1:10)

DF2.name = colnames(AT2@meta.data)[grepl("DF.classification", colnames(AT2@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(AT2, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(AT2, group.by = DF2.name) + NoAxes())

VlnPlot(AT2, features = "nFeature_RNA", group.by = DF2.name, pt.size = 0.1)

AT2 = AT2[, AT2@meta.data[, DF2.name] == "Singlet"]
dim(AT2)


# for AT3 data
AT3[["percent.mt"]] <- PercentageFeatureSet(AT3, pattern = "^MT-")
VlnPlot(AT3, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
AT3 <- subset(AT3, subset = nFeature_RNA >200 & nCount_RNA >2000 & nCount_RNA<60000 & percent.mt < 10 )

AT3 <- NormalizeData(AT3, normalization.method = "LogNormalize", scale.factor = 10000)
AT3 <- FindVariableFeatures(AT3, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(AT3)
AT3 <- ScaleData(AT3, features = all.genes)

AT3 <- RunPCA(AT3, features = VariableFeatures(object = AT3))
VizDimLoadings(AT3, dims = 1:2, reduction = "pca")
DimPlot(AT3, reduction = "pca")

AT3 <- FindNeighbors(AT3, dims = 1:20)
AT3 <- FindClusters(AT3, resolution = 0.5)

AT3 <- RunUMAP(AT3, dims = 1:20)
DimPlot(AT3, reduction = "umap")

# define the expected number of doublet cellscells.
AT3 <- doubletFinder_v3(AT3, pN = 0.25, pK = 0.09, nExp = nExp3, PCs = 1:10)

# name of the DF prediction can change, so extract the correct column name.
DF3.name = colnames(AT3@meta.data)[grepl("DF.classification", colnames(AT3@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(AT3, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(AT3, group.by = DF3.name) + NoAxes())

VlnPlot(AT3, features = "nFeature_RNA", group.by = DF3.name, pt.size = 0.1)

AT3 = AT3[, AT3@meta.data[, DF3.name] == "Singlet"]
dim(AT3)
####
####

save(AT1,AT2,AT3, file='/share/RD/PB_RD/10X/SLE/SLE_seurat/out/AT123.seurat.object.RData')



