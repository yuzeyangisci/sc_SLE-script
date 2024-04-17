##
library(clusterProfiler)
library(org.Hs.eg.db)


##################################
## T_C_NC
degenes_T_C_NC <- read.csv('F:/10X_SLE/TNK/DEsingle.results.classified.T.C_NC.csv',
                          header = T,stringsAsFactors = F,sep = ',')

degenes_T_C_NC_up <- subset(degenes_T_C_NC, State=="up")

genelist_T_C_NC_up <- degenes_T_C_NC_up$X
genelist_T_C_NC_up[duplicated(genelist_T_C_NC_up)]

gene_T_C_NC.df_up <- bitr(genelist_T_C_NC_up, fromType="SYMBOL",
                         toType="ENTREZID", 
                         OrgDb = "org.Hs.eg.db",
                         drop = T)

go_T_C_NC_up <- enrichGO(gene = gene_T_C_NC.df_up$ENTREZID, OrgDb = org.Hs.eg.db, 
                        ont='ALL',
                        pvalueCutoff = 0.01,
                        qvalueCutoff =0.9,
                        readable = T)

p1 <- dotplot(go_T_C_NC_up,showCategory=20)

write.csv(go_T_C_NC_up,"GO.T_C_NC.up.geneID.csv")


degenes_T_C_NC_down <- subset(degenes_T_C_NC, State=="down")

genelist_T_C_NC_down <- degenes_T_C_NC_down$X
genelist_T_C_NC_down[duplicated(genelist_T_C_NC_down)]

gene_T_C_NC.df_down <- bitr(genelist_T_C_NC_down, fromType="SYMBOL",
                           toType="ENTREZID", 
                           OrgDb = "org.Hs.eg.db",
                           drop = T)

go_T_C_NC_down <- enrichGO(gene = gene_T_C_NC.df_down$ENTREZID, OrgDb = org.Hs.eg.db, 
                          ont='ALL',
                          pvalueCutoff = 0.01,
                          qvalueCutoff =0.9,
                          readable = T)

p2 <- dotplot(go_T_C_NC_down,showCategory=20)

write.csv(go_T_C_NC_down,"GO.C_NC.down.geneID.csv")

p1+p2

##################################
## NK_C_NC
degenes_NK_C_NC <- read.csv('F:/10X_SLE/TNK/DEsingle.results.classified.NK.C_NC.csv',
                           header = T,stringsAsFactors = F,sep = ',')

degenes_NK_C_NC_up <- subset(degenes_NK_C_NC, State=="up")

genelist_NK_C_NC_up <- degenes_NK_C_NC_up$X
genelist_NK_C_NC_up[duplicated(genelist_NK_C_NC_up)]

gene_NK_C_NC.df_up <- bitr(genelist_NK_C_NC_up, fromType="SYMBOL",
                          toType="ENTREZID", 
                          OrgDb = "org.Hs.eg.db",
                          drop = T)

go_NK_C_NC_up <- enrichGO(gene = gene_NK_C_NC.df_up$ENTREZID, OrgDb = org.Hs.eg.db, 
                         ont='ALL',
                         pvalueCutoff = 0.01,
                         qvalueCutoff =0.9,
                         readable = T)

p1 <- dotplot(go_NK_C_NC_up,showCategory=20)

write.csv(go_NK_C_NC_up,"GO.NK_C_NC.up.geneID.csv")


degenes_NK_C_NC_down <- subset(degenes_NK_C_NC, State=="down")

genelist_NK_C_NC_down <- degenes_NK_C_NC_down$X
genelist_NK_C_NC_down[duplicated(genelist_NK_C_NC_down)]

gene_NK_C_NC.df_down <- bitr(genelist_NK_C_NC_down, fromType="SYMBOL",
                            toType="ENTREZID", 
                            OrgDb = "org.Hs.eg.db",
                            drop = T)

go_NK_C_NC_down <- enrichGO(gene = gene_NK_C_NC.df_down$ENTREZID, OrgDb = org.Hs.eg.db, 
                           ont='ALL',
                           pvalueCutoff = 0.01,
                           qvalueCutoff =0.9,
                           readable = T)

p2 <- dotplot(go_NK_C_NC_down,showCategory=20)

write.csv(go_NK_C_NC_down,"GO.C_NC.down.geneID.csv")

p1+p2
























##################################
## AT_NP
degenes_AT_NP <- read.csv('D:/WORK/PB/10X/DEsingle/DEsingle_SLE/DEsingle.results.classified.AT_NP.csv',
                         header = T,stringsAsFactors = F,sep = ',')

degenes_AT_NP_up <- subset(degenes_AT_NP, State=="up")

genelist_AT_NP_up <- degenes_AT_NP_up$X
genelist_AT_NP_up[duplicated(genelist_AT_NP_up)]

gene_AT_NP.df_up <- bitr(genelist_AT_NP_up, fromType="SYMBOL",
                        toType="ENTREZID", 
                        OrgDb = "org.Hs.eg.db",
                        drop = T)

go_AT_NP_up <- enrichGO(gene = gene_AT_NP.df_up$ENTREZID, OrgDb = org.Hs.eg.db, 
                       ont='ALL',
                       pvalueCutoff = 0.9,
                       qvalueCutoff =0.9,
                       readable = T)

p <- dotplot(go_AT_NP_up,showCategory=20)
p
write.csv(go_AT_NP_up,"GO.AT_NP.up.geneID.csv")


degenes_AT_NP_down <- subset(degenes_AT_NP, State=="down")

genelist_AT_NP_down <- degenes_AT_NP_down$X
genelist_AT_NP_down[duplicated(genelist_AT_NP_down)]

gene_AT_NP.df_down <- bitr(genelist_AT_NP_down, fromType="SYMBOL",
                          toType="ENTREZID", 
                          OrgDb = "org.Hs.eg.db",
                          drop = T)

go_AT_NP_down <- enrichGO(gene = gene_AT_NP.df_down$ENTREZID, OrgDb = org.Hs.eg.db, 
                         ont='ALL',
                         pvalueCutoff = 0.9,
                         qvalueCutoff =0.9,
                         readable = T)

p <- dotplot(go_AT_NP_down,showCategory=20)
p
write.csv(go_AT_NP_down,"GO.AT_NP.down.geneID.csv")





#################################
## C_PE
degenes_C_PE <- read.csv('D:/WORK/PB/10X/DEsingle/inner/DEsingle.results.classified.C_PE.csv',
                    header = T,stringsAsFactors = F,sep = ',')

degenes_C_PE_up <- subset(degenes_C_PE, State=="up")

genelist_C_PE_up <- degenes_C_PE_up$X
genelist_C_PE_up[duplicated(genelist_C_PE_up)]

gene_C_PE.df_up <- bitr(genelist_C_PE_up, fromType="SYMBOL",
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db",
                   drop = T)

go_C_PE_up <- enrichGO(gene = gene_C_PE.df_up$ENTREZID, OrgDb = org.Hs.eg.db, 
                  ont='ALL',
                  pvalueCutoff = 0.9,
                  qvalueCutoff =0.9,
                  readable = T)

p <- dotplot(go_C_PE_up,showCategory=20)
p
write.csv(go_C_PE_up,"GO.C_PE.up.geneID.csv")


degenes_C_PE_down <- subset(degenes_C_PE, State=="down")

genelist_C_PE_down <- degenes_C_PE_down$X
genelist_C_PE_down[duplicated(genelist_C_PE_down)]

gene_C_PE.df_down <- bitr(genelist_C_PE_down, fromType="SYMBOL",
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db",
                   drop = T)

go_C_PE_down <- enrichGO(gene = gene_C_PE.df_down$ENTREZID, OrgDb = org.Hs.eg.db, 
                  ont='ALL',
                  pvalueCutoff = 0.9,
                  qvalueCutoff =0.9,
                  readable = T)

p <- dotplot(go_C_PE_down,showCategory=20)
p
write.csv(go_C_PE_down,"GO.C_PE.down.geneID.csv")

## PE_PE31
degenes_PE_PE31 <- read.csv('D:/WORK/PB/10X/DEsingle/inner/DEsingle.results.classified.PE_PE31.csv',
                    header = T,stringsAsFactors = F,sep = ',')

degenes_PE_PE31_up <- subset(degenes_PE_PE31, State=="up")

genelist_PE_PE31_up <- degenes_PE_PE31_up$X
genelist_PE_PE31_up[duplicated(genelist_PE_PE31_up)]

gene_PE_PE31.df_up <- bitr(genelist_PE_PE31_up, fromType="SYMBOL",
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db",
                   drop = T)

go_PE_PE31_up <- enrichGO(gene = gene_PE_PE31.df_up$ENTREZID, OrgDb = org.Hs.eg.db, 
                  ont='ALL',
                  pvalueCutoff = 0.9,
                  qvalueCutoff =0.9,
                  readable = T)

p <- dotplot(go_PE_PE31_up,showCategory=20)
p
write.csv(go_PE_PE31_up,"GO.PE_PE31.up.geneID.csv")


degenes_PE_PE31_down <- subset(degenes_PE_PE31, State=="down")

genelist_PE_PE31_down <- degenes_PE_PE31_down$X
genelist_PE_PE31_down[duplicated(genelist_PE_PE31_down)]

gene_PE_PE31.df_down <- bitr(genelist_PE_PE31_down, fromType="SYMBOL",
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db",
                   drop = T)

go_PE_PE31_down <- enrichGO(gene = gene_PE_PE31.df_down$ENTREZID, OrgDb = org.Hs.eg.db, 
                  ont='ALL',
                  pvalueCutoff = 0.9,
                  qvalueCutoff =0.9,
                  readable = T)

p <- dotplot(go_PE_PE31_down,showCategory=20)
p
write.csv(go_PE_PE31_down,"GO.PE_PE31.down.geneID.csv")

## C_PE31
degenes_C_PE31 <- read.csv('D:/WORK/PB/10X/DEsingle/inner/DEsingle.results.classified.C_PE31.csv',
                    header = T,stringsAsFactors = F,sep = ',')

degenes_C_PE31_up <- subset(degenes_C_PE31, State=="up")

genelist_C_PE31_up <- degenes_C_PE31_up$X
genelist_C_PE31_up[duplicated(genelist_C_PE31_up)]

gene_C_PE31.df_up <- bitr(genelist_C_PE31_up, fromType="SYMBOL",
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db",
                   drop = T)

go_C_PE31_up <- enrichGO(gene = gene_C_PE31.df_up$ENTREZID, OrgDb = org.Hs.eg.db, 
                  ont='ALL',
                  pvalueCutoff = 0.9,
                  qvalueCutoff =0.9,
                  readable = T)

p <- dotplot(go_C_PE31_up,showCategory=20)
p
write.csv(go_C_PE31_up,"GO.C_PE31.up.geneID.csv")


degenes_C_PE31_down <- subset(degenes_C_PE31, State=="down")

genelist_C_PE31_down <- degenes_C_PE31_down$X
genelist_C_PE31_down[duplicated(genelist_C_PE31_down)]

gene_C_PE31.df_down <- bitr(genelist_C_PE31_down, fromType="SYMBOL",
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db",
                   drop = T)

go_C_PE31_down <- enrichGO(gene = gene_C_PE31.df_down$ENTREZID, OrgDb = org.Hs.eg.db, 
                  ont='ALL',
                  pvalueCutoff = 0.9,
                  qvalueCutoff =0.9,
                  readable = T)

p <- dotplot(go_C_PE31_down,showCategory=20)
p
write.csv(go_C_PE31_down,"GO.C_PE31.down.geneID.csv")

## AT1_AT2
degenes_AT1_AT2 <- read.csv('D:/WORK/PB/10X/DEsingle/inner/DEsingle.results.classified.AT1_AT2.csv',
                    header = T,stringsAsFactors = F,sep = ',')

degenes_AT1_AT2_up <- subset(degenes_AT1_AT2, State=="up")

genelist_AT1_AT2_up <- degenes_AT1_AT2_up$X
genelist_AT1_AT2_up[duplicated(genelist_AT1_AT2_up)]

gene_AT1_AT2.df_up <- bitr(genelist_AT1_AT2_up, fromType="SYMBOL",
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db",
                   drop = T)

go_AT1_AT2_up <- enrichGO(gene = gene_AT1_AT2.df_up$ENTREZID, OrgDb = org.Hs.eg.db, 
                  ont='ALL',
                  pvalueCutoff = 0.9,
                  qvalueCutoff =0.9,
                  readable = T)

p <- dotplot(go_AT1_AT2_up,showCategory=20)
p
write.csv(go_AT1_AT2_up,"GO.AT1_AT2.up.geneID.csv")


degenes_AT1_AT2_down <- subset(degenes_AT1_AT2, State=="down")

genelist_AT1_AT2_down <- degenes_AT1_AT2_down$X
genelist_AT1_AT2_down[duplicated(genelist_AT1_AT2_down)]

gene_AT1_AT2.df_down <- bitr(genelist_AT1_AT2_down, fromType="SYMBOL",
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db",
                   drop = T)

go_AT1_AT2_down <- enrichGO(gene = gene_AT1_AT2.df_down$ENTREZID, OrgDb = org.Hs.eg.db, 
                  ont='ALL',
                  pvalueCutoff = 0.9,
                  qvalueCutoff =0.9,
                  readable = T)

p <- dotplot(go_AT1_AT2_down,showCategory=20)
p
write.csv(go_AT1_AT2_down,"GO.AT1_AT2.down.geneID.csv")

## AT2_AT3
degenes_AT2_AT3 <- read.csv('D:/WORK/PB/10X/DEsingle/inner/DEsingle.results.classified.AT2_AT3.csv',
                    header = T,stringsAsFactors = F,sep = ',')

degenes_AT2_AT3_up <- subset(degenes_AT2_AT3, State=="up")

genelist_AT2_AT3_up <- degenes_AT2_AT3_up$X
genelist_AT2_AT3_up[duplicated(genelist_AT2_AT3_up)]

gene_AT2_AT3.df_up <- bitr(genelist_AT2_AT3_up, fromType="SYMBOL",
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db",
                   drop = T)

go_AT2_AT3_up <- enrichGO(gene = gene_AT2_AT3.df_up$ENTREZID, OrgDb = org.Hs.eg.db, 
                  ont='ALL',
                  pvalueCutoff = 0.9,
                  qvalueCutoff =0.9,
                  readable = T)

p <- dotplot(go_AT2_AT3_up,showCategory=20)
p
write.csv(go_AT2_AT3_up,"GO.AT2_AT3.up.geneID.csv")


degenes_AT2_AT3_down <- subset(degenes_AT2_AT3, State=="down")

genelist_AT2_AT3_down <- degenes_AT2_AT3_down$X
genelist_AT2_AT3_down[duplicated(genelist_AT2_AT3_down)]

gene_AT2_AT3.df_down <- bitr(genelist_AT2_AT3_down, fromType="SYMBOL",
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db",
                   drop = T)

go_AT2_AT3_down <- enrichGO(gene = gene_AT2_AT3.df_down$ENTREZID, OrgDb = org.Hs.eg.db, 
                  ont='ALL',
                  pvalueCutoff = 0.9,
                  qvalueCutoff =0.9,
                  readable = T)

p <- dotplot(go_AT2_AT3_down,showCategory=20)
p
write.csv(go_AT2_AT3_down,"GO.AT2_AT3.down.geneID.csv")

## AT1_AT3
degenes_AT1_AT3 <- read.csv('D:/WORK/PB/10X/DEsingle/inner/DEsingle.results.classified.AT1_AT3.csv',
                    header = T,stringsAsFactors = F,sep = ',')

degenes_AT1_AT3_up <- subset(degenes_AT1_AT3, State=="up")

genelist_AT1_AT3_up <- degenes_AT1_AT3_up$X
genelist_AT1_AT3_up[duplicated(genelist_AT1_AT3_up)]

gene_AT1_AT3.df_up <- bitr(genelist_AT1_AT3_up, fromType="SYMBOL",
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db",
                   drop = T)

go_AT1_AT3_up <- enrichGO(gene = gene_AT1_AT3.df_up$ENTREZID, OrgDb = org.Hs.eg.db, 
                  ont='ALL',
                  pvalueCutoff = 0.9,
                  qvalueCutoff =0.9,
                  readable = T)

p <- dotplot(go_AT1_AT3_up,showCategory=20)
p
write.csv(go_AT1_AT3_up,"GO.AT1_AT3.up.geneID.csv")


degenes_AT1_AT3_down <- subset(degenes_AT1_AT3, State=="down")

genelist_AT1_AT3_down <- degenes_AT1_AT3_down$X
genelist_AT1_AT3_down[duplicated(genelist_AT1_AT3_down)]

gene_AT1_AT3.df_down <- bitr(genelist_AT1_AT3_down, fromType="SYMBOL",
                   toType="ENTREZID", 
                   OrgDb = "org.Hs.eg.db",
                   drop = T)

go_AT1_AT3_down <- enrichGO(gene = gene_AT1_AT3.df_down$ENTREZID, OrgDb = org.Hs.eg.db, 
                  ont='ALL',
                  pvalueCutoff = 0.9,
                  qvalueCutoff =0.9,
                  readable = T)

p <- dotplot(go_AT1_AT3_down,showCategory=20)
p
write.csv(go_AT1_AT3_down,"GO.AT1_AT3.down.geneID.csv")


