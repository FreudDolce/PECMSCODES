library('pheatmap')

file = '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/gsva_result.csv'
pathEs <- read.csv(file, header=T, row.names=1)

tiff('~/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/ClusterHeatMap.tiff', width=2500, height=900, res=300)
gsva_fig <- pheatmap(pathEs, show_colnames = FALSE, show_rownames = TRUE)

dev.off()
