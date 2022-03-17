library(BiocManager)
library(limma)
library(edgeR)


args <- commandArgs(T)
print (paste0('Expression data: ', args[1]))
print (paste0('Class data: ', args[2]))
print (paste0('Number of classes: ', args[3]))
print (paste0('Data save as: ', args[4]))

expre_data <- read.csv(args[1], header = T)
rownames(expre_data) = expre_data[, 1]
expre_data = expre_data[, -1]
dge <- DGEList(counts = expre_data)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE, prior.count = 3)


class_data <- read.csv(args[2], header = T)
group <- class_data[, 2]

design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
rownames(design)=colnames(data)

class_list <- unique(class_data$x)

contrast.matrix <- makeContrasts('COL', levels=design)
fit <- lmFit(logCPM, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
##step3
tT = topTable(fit2, adjust = 'BH', sort.by = 'logFC', n = Inf)
write.csv(tT, args[4])
