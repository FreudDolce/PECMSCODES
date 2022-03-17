library(limma)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggrepel)

TEST = 'TYPE2'
exp_file = paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/', TEST, '_EXP_SYMBOL.csv')
#class_file = paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/', TEST, '_ori_class_for_limma.csv')
class_file = paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/', TEST, '_post_lasso_coef_2.csv')
exp_info = read.csv(exp_file, header=T, row.names=1, stringsAsFactors = F)

dge <- DGEList(counts=exp_info)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

group_file <- read.csv(class_file, header=T, row.names=1)
group_list <- factor(group_file$lasso_class, levels=unique(group_file$lasso_class))
design <- model.matrix(~0 + group_list)
print(design)
colnames(design) <- unique(group_file$lasso_class)
rownames(design) <- rownames(group_file)

cont.wt <- makeContrasts(STROMA_HIGH-STROMA_LOW, levels=design)
fit <- lmFit(logCPM, design)
fit2 <- contrasts.fit(fit, cont.wt)
fit2 <- eBayes(fit2)
tT = topTable(fit2, adjust="BH",sort.by="logFC",n=Inf)
write.csv(tT, paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/', TEST, '_post_limma.csv'))
diffsig = na.omit(tT)

foldChange = 1.5
padj = 0.0001

All_diffSig <- diffsig[(diffsig$adj.P.Val < padj & (diffsig$logFC>foldChange | diffsig$logFC < (-foldChange))),]
write.csv(All_diffSig, paste0('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/', TEST, '_post_limma_sig.csv'))

diffup <-  All_diffSig[(All_diffSig$P.Value < padj & (All_diffSig$logFC > foldChange)),]
diffdown <- All_diffSig[(All_diffSig$P.Value < padj & (All_diffSig < -foldChange)),]

logFC <- diffsig$logFC
deg.padj <- diffsig$P.Value
data <- data.frame(logFC = logFC, padj = deg.padj)
data$group[(data$padj > padj | data$padj == "NA") | (data$logFC < foldChange) & data$logFC > -foldChange] <- "Not"
data$group[(data$padj <= padj & data$logFC > 1)] <-  "Up"
data$group[(data$padj <= padj & data$logFC < -1)] <- "Down"

#pdf('volcano.pdf',width = 7,height = 6.5)
tiff(paste0('/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/', TEST, '_post_limma_volcano.tiff') ,width = 1800,height = 1600, res=300, pointsize = 8)
label = subset(diffsig,P.Value < padj & abs(logFC) > foldChange)
label1 = rownames(label)

colnames(diffsig)[1] = 'log2FC'
Significant=ifelse((diffsig$P.Value < padj & abs(diffsig$log2FC)> foldChange), ifelse(diffsig$log2FC > foldChange,"Up","Down"), "Not")

ggplot(diffsig, aes(log2FC, -log10(P.Value)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("#0072B5","grey","#BC3C28"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-foldChange, foldChange), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(padj),colour="black", linetype="dashed")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  labs(x="log2(FoldChange)",y="-log10(Pvalue)")+
  theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
  str(diffsig, max.level = c(-1, 1))+theme_bw()

dev.off()
