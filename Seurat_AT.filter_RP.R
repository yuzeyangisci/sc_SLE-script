library(dplyr)
library(Seurat)
library(patchwork)

library(DoubletFinder)



load('/share/RD/PB_RD/10X/SLE/SLE_seurat/out/AT123.seurat.object.RData')

#####
counts <- GetAssayData(AT1, slot="counts", assay="RNA")   

# cells
genes.use.AT1 <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",
                  rownames(AT1), value=TRUE, invert=TRUE) 

counts.sub.AT1 <- counts[genes.use.AT1,]
AT1.filter <- CreateSeuratObject(counts=counts.sub.AT1, project = "AT1.filter")

#####
counts <- GetAssayData(AT2, slot="counts", assay="RNA")   

# cells
genes.use.AT2 <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",
                  rownames(AT2), value=TRUE, invert=TRUE) 

counts.sub.AT2 <- counts[genes.use.AT2,]
AT2.filter <- CreateSeuratObject(counts=counts.sub.AT2, project = "AT2.filter")

###
counts <- GetAssayData(AT3, slot="counts", assay="RNA")   

# cells
genes.use.AT3 <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",
                  rownames(AT3), value=TRUE, invert=TRUE) 

counts.sub.AT3 <- counts[genes.use.AT3,]
AT3.filter <- CreateSeuratObject(counts=counts.sub.AT3, project = "AT3.filter")

save(AT1.filter,AT2.filter,AT3.filter, file='/share/RD/PB_RD/10X/SLE/SLE_seurat/out/AT123.RP.seurat.object.RData')




