library(survival)
library(survminer)

COEF <- 1
db <- 'TCGA'
item <- 'TYPE2'

stu_file <- paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_', db, '_', item, '_COEF_', COEF, '.csv')
stu_data <- read.csv(stu_file, header=T, row.names=1)

head(stu_data)
