library(dplyr)
library(Seurat)
library(patchwork)

library(DoubletFinder)



load('/share/RD/PB_RD/10X/SLE/SLE_seurat/out/SLE_C123.seurat.object.RData')

#####
counts <- GetAssayData(SLE_C1, slot="counts", assay="RNA")   

# cells
genes.use.SLE_C1 <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",
                  rownames(SLE_C1), value=TRUE, invert=TRUE) 

counts.sub.SLE_C1 <- counts[genes.use.SLE_C1,]
SLE_C1.filter <- CreateSeuratObject(counts=counts.sub.SLE_C1, project = "SLE_C.filter")

#####
counts <- GetAssayData(SLE_C2, slot="counts", assay="RNA")   

# cells
genes.use.SLE_C2 <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",
                  rownames(SLE_C2), value=TRUE, invert=TRUE) 

counts.sub.SLE_C2 <- counts[genes.use.SLE_C2,]
SLE_C2.filter <- CreateSeuratObject(counts=counts.sub.SLE_C2, project = "SLE_PE.filter")

###
counts <- GetAssayData(SLE_C3, slot="counts", assay="RNA")   

# cells
genes.use.SLE_C3 <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",
                  rownames(SLE_C3), value=TRUE, invert=TRUE) 

counts.sub.SLE_C3 <- counts[genes.use.SLE_C3,]
SLE_C3.filter <- CreateSeuratObject(counts=counts.sub.SLE_C3, project = "SLE_PE31.filter")

save(SLE_C1.filter,SLE_C2.filter,SLE_C3.filter, file='/share/RD/PB_RD/10X/SLE/SLE_seurat/out/SLE_C123.RP.seurat.object.RData')




