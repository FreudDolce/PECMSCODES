library('maftools')

COEF = '2'
item = 'TYPE2'
mutfile = '~/Documents/MANU/BI/PANC_ECM/DATA/SomaticMutation/cpt_comb.maf'
clinicalfile = paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_CPT_', item, '_COEF_', COEF, '.csv')
savepath = paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/CPT_', item, '_COEF_', COEF, '.tiff')

cpt = read.maf(maf = mutfile, clinicalData = clinicalfile)#, isTCGA = T)
#print(cpt)

tiff(filename = savepath, width = 2000, height = 750, res = 300)

oncop = oncoplot(maf = cpt,
		 top = 10,
		 clinicalFeatures = c('lasso_result'),
		 sortByAnnotation = T)
dev.off()

tmb_cpt <- tmb(maf=cpt, logScale=F)
print(tmb_cpt)
write.csv(tmb_cpt, '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TMB_CPT.csv')
