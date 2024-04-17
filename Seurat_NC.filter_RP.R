library(dplyr)
library(Seurat)
library(patchwork)

library(DoubletFinder)



load('/share/RD/PB_RD/10X/SLE/SLE_seurat/out/NC123.seurat.object.RData')

#####
counts <- GetAssayData(NC1, slot="counts", assay="RNA")   

# cells
genes.use.NC1 <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",
                  rownames(NC1), value=TRUE, invert=TRUE) 

counts.sub.NC1 <- counts[genes.use.NC1,]
NC1.filter <- CreateSeuratObject(counts=counts.sub.NC1, project = "NC1.filter")

#####
counts <- GetAssayData(NC2, slot="counts", assay="RNA")   

# cells
genes.use.NC2 <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",
                  rownames(NC2), value=TRUE, invert=TRUE) 

counts.sub.NC2 <- counts[genes.use.NC2,]
NC2.filter <- CreateSeuratObject(counts=counts.sub.NC2, project = "NC2.filter")

###
counts <- GetAssayData(NC3, slot="counts", assay="RNA")   

# cells
genes.use.NC3 <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",
                  rownames(NC3), value=TRUE, invert=TRUE) 

counts.sub.NC3 <- counts[genes.use.NC3,]
NC3.filter <- CreateSeuratObject(counts=counts.sub.NC3, project = "NC3.filter")

save(NC1.filter,NC2.filter,NC3.filter, file='/share/RD/PB_RD/10X/SLE/SLE_seurat/out/NC123.RP.seurat.object.RData')




