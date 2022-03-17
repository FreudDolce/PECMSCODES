library(survival)
library(survminer)

surv_file <- '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/IMV_clinical.csv'
figure_savename <- '~/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/IMV_surv.tiff'

surv_data <- read.csv(surv_file)
rownames(surv_data) <- surv_data[, 1]

surv <- Surv(time = surv_data$os, surv_data$censOS == 1)
sfit <- survfit(formula = surv~lasso_result, data = surv_data)
print(sfit)

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
	   break.x.by = 5,
	   xlim = c(0, 25)
)
print(surv_p)
dev.off()
