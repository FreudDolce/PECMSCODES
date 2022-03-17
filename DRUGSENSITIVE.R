library(pRRophetic)
data(cgp2016ExprRma)
dim(cgp2016ExprRma)

args <- commandArgs(T)
print (paste0('mRNA expression file: ', args[1]))
print (paste0('Drug: ', args[2]))
print (paste0('Save as: ', args[3]))

exprdata = read.csv(args[1], header=T)
rownames(exprdata) = exprdata[, 1]
exprdata = exprdata[, -1]

Testexprdata <- as.matrix(exprdata)

drugsensitive <- pRRopheticPredict(testMatrix = Testexprdata,
                                   drug = args[2],
                                   tissueType = "all", 
                                   batchCorrect = "eb",
                                   selection = 1,
                                   dataset = "cgp2016")

write.table(drugsensitive, args[3], row.names = TRUE, col.name = TRUE, sep = ',')
