library('maftools')

COEF = '2'
item = 'TYPE2'
mutfile = '/Users/freud/Documents/MANU/BI/PANC_ECM/DATA/SomaticMutation/TCGA.PAAD.mutect.fea333b5-78e0-43c8-bf76-4c78dd3fac92.DR-10.0.somatic.maf.gz'

clinicalfile = paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_TCGA_', item, '_COEF_', COEF, '.csv')
savepath = paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/FALL_TCGA_', item, '_COEF_', COEF, '.tiff')

paad = read.maf(maf = mutfile, clinicalData = clinicalfile, isTCGA = T)
#print(cpt)

tiff(filename = savepath, width = 2000, height = 750, res = 300)

oncop = oncoplot(maf = paad,
		 top = 10,
		 clinicalFeatures = c('lasso_result'),
		 sortByAnnotation = T)
dev.off()

tmb_tcga <- tmb(maf=paad, logScale=F)
print(tmb_tcga)
write.csv(tmb_tcga, '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TMB_TCGA.csv')
