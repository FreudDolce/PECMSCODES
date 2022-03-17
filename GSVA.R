library(BiocManager)
library(GSVA)
library(GSEABase)
library(limma)
library(pheatmap)

# This part is used to generate the .csv file of gsva result

args <- commandArgs(T)
print (paste0('mRNA expression file: ', args[1]))
print (paste0('Gmtfile: ', args[2]))
print (paste0('Save path :', args[3]))

pathSet = getGmt(args[2])


save_pheatmap_pdf <- function(x, filename, width=10, height=4.5) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

stu_data = read.csv(args[1])
rownames(stu_data) = stu_data[, 1]
stu_data = stu_data[, -1]

pathEs = gsva(expr = as.matrix(stu_data), gset.idx.list = pathSet, parallel.sz = 4)#, kcdf = 'Poisson')

write.table(pathEs, paste0(args[3], 'gsva_result.csv'), row.names = TRUE, col.names = TRUE, sep = ',')

gsva_fig <- pheatmap(pathEs, show_colnames = FALSE, show_rownames = TRUE)

v <- 2: 6

for (i in v){
   col_cluster <- cutree(gsva_fig$tree_col, k = i)
   write.csv(col_cluster, paste0(args[3], 'ClusterResult/cluster_', i, '.csv'))
   row_cluster <- cutree(gsva_fig$tree_row, k = i)
   write.csv(row_cluster, paste0(args[3], 'ClusterResult/p_cluster_', i, '.csv'))
}

save_pheatmap_pdf(gsva_fig, paste0(args[3], 'gsva_heatmap.pdf'))
