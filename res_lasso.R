library(glmnet)
library(foreign)

item = 'TYPE2'

lassofile <- paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/', item, '_for_limma_norm.csv')
lassoinfo <- read.csv(lassofile, header=T, row.names=1)
n_cols <- dim(lassoinfo)[2]

X <- as.matrix(lassoinfo[, 1: n_cols - 1])
y <- as.matrix(lassoinfo[, n_cols])

f1 = glmnet(X, y, family="gaussian", nlambda=100, alpha=1)
tiff(paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/Lasso/', item, '_lasso_lambda.tiff'), width=1600, height=1600, res=300)

f1_plot <- plot(f1, xvar='lambda', label=T)

dev.off()

cvfit = cv.glmnet(X, y)
tiff(paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/Lasso/', item, '_lasso_cvfit.tiff'), width=1600, height=1600, res=300)

cv_plot <- plot(cvfit)

dev.off()

par_1 <- cvfit$lambda.min
par_2 <- cvfit$lambda.1se

l.coef1 <- coef(cvfit$glmnet.fit, s=par_1, exact = F)
l.coef2 <- coef(cvfit$glmnet.fit, s=par_2, exact = F)

coef1 = as.matrix(l.coef1)
coef2 = as.matrix(l.coef2)

write.csv(coef1, paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/', item, '_lasso_coef1.csv'))
write.csv(coef2, paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/', item, '_lasso_coef2.csv'))
