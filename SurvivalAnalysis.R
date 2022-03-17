library('survival')
library('survminer')

args <- commandArgs(TRUE)
print(paste('Reading class: ', args[1]))

survival_file <- args[1] 
survival_data <- read.csv(survival_file, header = TRUE)
rownames(survival_data) = survival_data[, 1]
survival_data = survival_data[, -1]
#survival_data$surv <- Surv(time = survival_data$days_to_death, survival_data$vital_status == 2)
surv <- Surv(time = survival_data$days_to_death, survival_data$vital_status == 2)

sfit <- survfit(formula = surv~class, data=survival_data)
summary(sfit)

survplot <- ggsurvplot( sfit, # 创建的拟合对象
           		data = survival_data,  # 指定变量数据来
           		conf.int = TRUE, # 显示置信区间
           		pval = TRUE, # 添加P值
           		surv.median.line = "hv",  # 添加中位生存时间线
           		risk.table = TRUE, # 添加风险表
           		xlab = "Follow up time(d)", # 指定x轴标签
           		legend = c(0.8,0.75), # 指定图例位置
           		legend.title = "", # 设置图例标题，这里设置不显示标题，用空格替代
           		break.x.by = 100) 
ggsave('survplot.pdf')
sdiff <- survdiff(formula = surv~class, data=survival_data)

schiq <- sdiff$chisq
p.value <- 1 - pchisq(sdiff$chisq, length(sdiff$n) -1)

plist <- matrix(c(schiq, p.value), nrow = 1)
write.csv(plist, 'survival_p.csv')
