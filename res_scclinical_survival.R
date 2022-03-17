library(survival)
library(survminer)

surv_file <- paste0('~/Documents/MANU/BI/PANC_ECM/DATA/Clinical_data_exp/Clinical_data.csv')

figure_savename <- paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/Scclinical_Surv.tiff')


surv_data <- read.csv(surv_file)
rownames(surv_data) <- surv_data[, 1]
surv_data <- surv_data[, -1]

surv <- Surv(time = surv_data$OS, surv_data$STATUS== 1)
sfit <- survfit(formula = surv~PECMS_CLASS, data = surv_data)
print(sfit)

tiff(figure_savename, width=1150, height=1200, res=300)
surv_p <- ggsurvplot(sfit,
	   data = surv_data,
	   palette = 'npg',
	   conf.int = T, 
	   conf.int.alpha = 0.1, 
	   tables.height = 0.35,
	   pval = T, pval.size = 4, 
	   surv.median.line = 'hv', 
	   xlab = 'Follow up time (d)', 
	   legend = c(0.8, 0.75), 
	   legend.title = '', 
	   break.x.by = 6,
	   xlim = c(0, 30)
)
print(surv_p)
dev.off()
