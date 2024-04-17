# 安装{BioManager}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# 安装{EnhancedVolcano}
BiocManager::install("EnhancedVolcano")

# 载入
library(EnhancedVolcano)

res.T <- read.csv("F:/10X_SLE/TNK/DEsingle.results.classified.T.C_NC.csv",header=TRUE,row.names = 1)
res.NK <- read.csv("F:/10X_SLE/TNK/DEsingle.results.classified.NK.C_NC.csv",header=TRUE,row.names = 1)

res <- res.T
#res$threshold<-as.factor(ifelse(res$log2FoldChange >= 2,'Up',ifelse(res$pvalue<0.05 & res$log2FoldChange <= -2,'Down','Not')))
res<-data.frame(res)
res$log2FoldChange <- log2(res$foldChange)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'p_val_adj')
p1<-EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',pCutoff = 10e-2,
                xlim = c(-5, 5), #                pCutoff = 2,
                title = 'T_C/NC ')


p2<-EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue.adj.FDR',pCutoff = 10e-2,
                xlim = c(-5, 5), #                pCutoff = 2,
                title = 'T_C/NC ')
p1+p2
