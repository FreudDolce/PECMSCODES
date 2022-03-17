library(survival)
library(survminer)

coef <- '2'
db <- 'CPT'
item <- 'TYPE2'
surv_file <- paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_', db, '_', item, '_COEF_', coef, '.csv')

print(paste0('Survival file: ', surv_file))
print(paste0('Analysis: ', item))

figure_savename <- paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/', db,'_', item, '_COEF_', coef, '.tiff')
print(paste0('Save as: ', figure_savename))

surv_data <- read.csv(surv_file)
rownames(surv_data) <- surv_data[, 1]
surv_data <- surv_data[, -1]

surv <- Surv(time = surv_data$days_to_death, surv_data$vital_status == 1)
sfit <- survfit(formula = surv~lasso_result, data = surv_data)

tiff(figure_savename, width=1800, height=1600, res=300)
surv_p <- ggsurvplot(sfit,
	   data = surv_data,
	   palette = 'npg',
	   conf.int = T, 
	   conf.int.alpha = 0.1, 
	   risk.table = T, 
	   tables.height = 0.35,
	   pval = T, pval.size = 4, 
	   surv.median.line = 'hv', 
	   xlab = 'Follow up time (d)', 
	   legend = c(0.8, 0.75), 
	   legend.title = '', 
	   break.x.by = 400, 
	   xlim = c(0, 2800)
)
#print(surv_p)
print(sfit)
dev.off()
