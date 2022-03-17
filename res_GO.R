library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图

item = 'TYPE2'
fig_width = 1600
info = read.csv(paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/', item, '_limma_sig.csv'), header=T)
print(info)
GO_database <- 'org.Hs.eg.db'
gene <- bitr(info$SYMBOL, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

GO<-enrichGO(gene$ENTREZID,#GO富集分析
             OrgDb = GO_database,
             keyType = "ENTREZID",#设定读取的gene ID类型
             ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
             pvalueCutoff = 0.05,#设定p值阈值
             qvalueCutoff = 0.05,#设定q值阈值
             readable = T)


tiff(paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/', item, '_GO_barplot.tiff'),
     width=1.5*fig_width, height=fig_width, res=300)
barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
dev.off()

tiff(paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/', item, '_GO_dotplot.tiff'),
     width=1.5*fig_width, height=fig_width, res=300)
dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图和点状图
dev.off()

tiff(paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/', item, '_GO_relationship.tiff'),
     width=1.5*fig_width, height=1.5*fig_width, res=300)
enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE)#基因-通路关联网络图
dev.off()


GO2 <- pairwise_termsim(GO)

tiff(paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/', item, '_GO_pathway_relationship.tiff'),
     width=1.5*fig_width, height=1.5*fig_width, res=300)
enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk")#通路间关联网络图
dev.off()

info = read.csv(paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/', item, '_limma.csv'), header=T)
names(info) <- c('SYMBOL','Log2FoldChange','pvalue','padj')
info_merge <- merge(info,gene,by='SYMBOL')#合并转换后的基因ID和Log2FoldChange
GSEA_input <- info_merge$Log2FoldChange
names(GSEA_input) = info_merge$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)
GSEA_GO <- gseGO(GSEA_input, OrgDb=org.Hs.eg.db, pvalueCutoff = 0.5)#GSEA富集分析
tiff(paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/', item, '_GSEA_ridgeplot.tiff'),
     width=1.5*fig_width, height=fig_width, res=300)
ridgeplot(GSEA_GO) 
dev.off()

#tiff(paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/', item, '_GSEA_kegg.tiff'),
#     width=1600, height=1200, res=300)
#gseaplot2(GSEA_KEGG,1)
#dev.off()

#gseaplot2(GSEA_KEGG,1:30)#30是根据ridgeplot中有30个富集通路得到的
