library(dplyr)
library(Seurat)
library(patchwork)

library(DoubletFinder)



load('/share/RD/PB_RD/10X/SLE/SLE_seurat/out/NP123.seurat.object.RData')

#####
counts <- GetAssayData(NP1, slot="counts", assay="RNA")   

# cells
genes.use.NP1 <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",
                  rownames(NP1), value=TRUE, invert=TRUE) 

counts.sub.NP1 <- counts[genes.use.NP1,]
NP1.filter <- CreateSeuratObject(counts=counts.sub.NP1, project = "NP1.filter")

#####
counts <- GetAssayData(NP2, slot="counts", assay="RNA")   

# cells
genes.use.NP2 <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",
                  rownames(NP2), value=TRUE, invert=TRUE) 

counts.sub.NP2 <- counts[genes.use.NP2,]
NP2.filter <- CreateSeuratObject(counts=counts.sub.NP2, project = "NP2.filter")

###
counts <- GetAssayData(NP3, slot="counts", assay="RNA")   

# cells
genes.use.NP3 <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",
                  rownames(NP3), value=TRUE, invert=TRUE) 

counts.sub.NP3 <- counts[genes.use.NP3,]
NP3.filter <- CreateSeuratObject(counts=counts.sub.NP3, project = "NP3.filter")

save(NP1.filter,NP2.filter,NP3.filter, file='/share/RD/PB_RD/10X/SLE/SLE_seurat/out/NP123.RP.seurat.object.RData')




