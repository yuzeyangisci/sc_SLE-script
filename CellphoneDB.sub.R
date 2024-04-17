library(Seurat)

# # pbmc_small
# 
# #########################################################################################
# ## B   0   naive B
# ##     1   memory B
# ##     2   naive B
# ##     3   memory B
# ##     4   plasma cell B
# 
# load('/share/RD/PB_RD/10X/SLE/SLE_seurat/4group/out_sub/pbmc.B.integrated.seurat.object.RData')
# 
# current.cluster.ids <- c(0, 1, 2, 3, 4)
# new.cluster.ids <- c("naive B", "memory B", "naive B", "memory B", "plasma cell B")
# names(new.cluster.ids) <- levels(pbmc.B)
# pbmc.B <- RenameIdents(pbmc.B, new.cluster.ids)
# 
# 
# 
# 
# pbmc.B$cell_type <- plyr::mapvalues(x = pbmc.B$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
# 
# head(pbmc.B@meta.data)
# 
# write.table(as.matrix(pbmc.B@assays$RNA@data), 'cellphonedb_count.B.txt', sep='\t', quote=F)
# meta_data <- cbind(rownames(pbmc.B@meta.data), pbmc.B@meta.data[,'cell_type', drop=F])  
# meta_data <- as.matrix(meta_data)
# meta_data[is.na(meta_data)] = "Unkown"
# 
# write.table(meta_data, 'cellphonedb_meta.B.txt', sep='\t', quote=F, row.names=F)
# 
# 
# # #########################################################################################
# # ## Mono    0   CD16 monocyte
# # ##         1   CD14 monocyte
# # ##         2   CD14 monocyte
# # 
# # load('/share/RD/PB_RD/10X/SLE/SLE_seurat/4group/out_sub/pbmc.Monocytes.integrated.seurat.object.RData')
# # 
# # current.cluster.ids <- c(0, 1, 2)
# # new.cluster.ids <- c("CD16 monocyte", "CD14 monocyte", "CD14 monocyte")
# # names(new.cluster.ids) <- levels(pbmc.Monocytes)
# # pbmc.Monocytes <- RenameIdents(pbmc.Monocytes, new.cluster.ids)
# # 
# # pbmc.Monocytes$cell_type <- plyr::mapvalues(x = pbmc.Monocytes$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
# # 
# # head(pbmc.Monocytes@meta.data)
# # 
# # write.table(as.matrix(pbmc.Monocytes@assays$RNA@data), 'cellphonedb_count.Monocytes.txt', sep='\t', quote=F)
# # meta_data <- cbind(rownames(pbmc.Monocytes@meta.data), pbmc.Monocytes@meta.data[,'cell_type', drop=F])  
# # meta_data <- as.matrix(meta_data)
# # meta_data[is.na(meta_data)] = "Unkown"
# # 
# # write.table(meta_data, 'cellphonedb_meta.Monocytes.txt', sep='\t', quote=F, row.names=F)
# # 
# # 
# 
# 
# 
#########################################################################################
## T-NK    0   naive CD8 T
##         1   memory CD8 T
##         2   naive CD4 T
##         3   memory CD4 T
##         4   CD16 NK
##         5   memory CD8 T
##         6   naive CD4 T
##         7   memory CD4 T
##         8   gdT
##         9   memory CD8 T
##         10  naive CD4 T
##         11  CD56 NK
##         12  regulatory CD4 T


# load('/share/RD/PB_RD/10X/SLE/SLE_seurat/4group/out_sub/pbmc.NK.integrated.seurat.object.RData')
load('/share/RD/PB_RD/10X/SLE/SLE_seurat/4group/out_sub_0426_TNK/pbmc.TNK.integrated.seurat.object.RData')

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


pbmc.T  <- subset(x = pbmc.TNK, subset = seurat_clusters == c("1","2","3","5","6","7","8","9","10","12"))
pbmc.NK <- subset(x = pbmc.TNK, subset = seurat_clusters == c("4","11"))

pbmc.T$cell_type <- plyr::mapvalues(x = pbmc.T$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

head(pbmc.T@meta.data)

write.table(as.matrix(pbmc.T@assays$RNA@data), 'cellphonedb_count.T.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(pbmc.T@meta.data), pbmc.T@meta.data[,'cell_type', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown"

write.table(meta_data, 'cellphonedb_meta.T.txt', sep='\t', quote=F, row.names=F)



pbmc.NK$cell_type <- plyr::mapvalues(x = pbmc.NK$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

head(pbmc.NK@meta.data)

write.table(as.matrix(pbmc.NK@assays$RNA@data), 'cellphonedb_count.NK.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(pbmc.NK@meta.data), pbmc.NK@meta.data[,'cell_type', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown"

write.table(meta_data, 'cellphonedb_meta.NK.txt', sep='\t', quote=F, row.names=F)



# pbmc.TNK.C <- subset(x = pbmc, subset = new.ident == "C")
# pbmc.TNK.C1 <- subset(x = pbmc, subset = orig.ident == "SLE_C.filter")
# pbmc.TNK.C2 <- subset(x = pbmc, subset = orig.ident == "SLE_PE.filter")
# pbmc.TNK.C3 <- subset(x = pbmc, subset = orig.ident == "SLE_PE31.filter")
# pbmc.TNK.NC <- subset(x = pbmc, subset = new.ident == "NC")


# ###################################################################################################
# ## C
# pbmc.TNK.C$cell_type <- plyr::mapvalues(x = pbmc.TNK.C$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
# 
# head(pbmc.TNK.C@meta.data)
# 
# write.table(as.matrix(pbmc.TNK.C@assays$RNA@data), 'cellphonedb_count.TNK.txt', sep='\t', quote=F)
# meta_data <- cbind(rownames(pbmc.TNK.C@meta.data), pbmc.TNK.C@meta.data[,'cell_type', drop=F])  
# meta_data <- as.matrix(meta_data)
# meta_data[is.na(meta_data)] = "Unkown"
# 
# write.table(meta_data, 'cellphonedb_meta.pbmc.TNK.C.txt', sep='\t', quote=F, row.names=F)
# 
# 
# 
