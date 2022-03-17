library('maftools')

args <- commandArgs(T)
print (paste0('Mutationfile: ', args[1]))
print (paste0('Clinical annotation: ', args[2]))
print (paste0('Picture save in: ', args[3]))

setwd('~/Documents/MANU/BI/PANC_ECM/DATA/'

paad = read.maf(maf = args[1], clinicalData = args[2])#, isTCGA = T)
print(paad)

jpeg(filename = args[3], width = 1200, height = 800, res = 120)

oncop = oncoplot(maf = paad,
		 top = 20,
		 clinicalFeatures = c('define'),
		 sortByAnnotation = T)
dev.off()
