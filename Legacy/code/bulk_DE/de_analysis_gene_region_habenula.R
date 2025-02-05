library(SummarizedExperiment)
library(recount)
library(edgeR)
library(limma)
library(jaffelab)
library(RColorBrewer)
library(ggplot2)

load("/dcl01/ajaffe/data/lab/brain_swap/count_data/RNAseq_Collection_postQC_n5536_11dataset_2021-01_geneRSE.Rdata")

## filter to relevant studies and samples
rel_dataset = c("BrainSeq_Phase2_DLPFC",
	"BrainSeq_Phase2_HIPPO", "BrainSeq_Phase3_Caudate", "BrainSeq_Phase4and5",
	"Habenula", "Nicotine_NAc", "psychENCODE_BP", "psychENCODE_Mood",
	"VA_PTSD")

keepIndex = which(rse_gene$Dataset %in% rel_dataset &
	rse_gene$Sex == "M" & rse_gene$Age > 20 & rse_gene$Age < 69 & rse_gene$Dx =="Control")
rse_gene = rse_gene[,keepIndex]


colData(rse_gene)$Hab_ID<- factor(ifelse(colData(rse_gene)$Region=="Habenula","Habenula","Other"))
rpkm<- getRPKM(rse_gene,"Length")
Means<-rowMeans(rpkm[rowData(rse_gene)$Symbol=="GPR151"])
rpkm[rowData(rse_gene)$Symbol=="GPR151"]
Means[rowData(rse_gene)$Symbol=="POU4F1"]

palette(brewer.pal(11,"Set3"))
DT_GPR151<-as.data.frame(log2(rpkm[rowData(rse_gene)$Symbol=="GPR151"]+1))
DT_GPR151$Region<-rse_gene$Hab_ID
colnames(DT_GPR151)[1]<-"RPKM"

DT_POU4F1<-as.data.frame(log2(rpkm[rowData(rse_gene)$Symbol=="POU4F1"]+1))
DT_POU4F1$Region<-rse_gene$Hab_ID
colnames(DT_POU4F1)[1]<-"RPKM"

 summary(lm(log2(RPKM+1) ~ Region, data = DT_POU4F1))

p <- ggplot(DT_POU4F1, aes(x=Region, y=RPKM)) +
  geom_boxplot(outlier.shape = NA)
p<-p + geom_jitter(shape=16, size=1,alpha=.5, position=position_jitter(0.2),aes(colour=as.factor(rse_gene$Region)))
p<-p + ggtitle("POU4F1") + ylab("log2(RPKM)") + ylim(0, .75) +theme(text = element_text(size=20))
p$labels$colour<-"Region"
pdf("pdfs/POU4F1.pdf")
print(p)
dev.off()
t.test(DT_POU4F1$RPKM)
p <- ggplot(DT_GPR151, aes(x=Region, y=RPKM)) +
  geom_boxplot(outlier.shape = NA)
p<-p + geom_jitter(shape=16, size=1,alpha=.5, position=position_jitter(0.2),aes(colour=as.factor(rse_gene$Region)))
p<-p + ggtitle("GPR151") + ylab("log2(RPKM)") + ylim(0, 1)+theme(text = element_text(size=20))
p$labels$colour<-"Region"

pdf("pdfs/GPR151.pdf")
print(p)
dev.off()

rse_gene = rse_gene[rowMeans(getRPKM(rse_gene,"Length")) > 0.2,]
pd<- colData(rse_gene)

num_cols <- unlist(lapply(pd, is.numeric))
pd_num <- pd[ , num_cols]
pd_num<-as.data.frame(pd_num)

pdf(file="pdfs/check_confounders.pdf")
for(i in 1:ncol(pd_num)){
    boxplot(pd_num[,i] ~rse_gene$Region, xlab ="Region",ylab=colnames(pd_num)[i],las = 2)
}
dev.off()

dim(rse_gene)
#[1] 25911  1266

table(rse_gene$Region)
  # Amygdala   BasoAmyg    Caudate       dACC      DLPFC   Habenula      HIPPO
       # 124         74        147         74        213         34        150
# MedialAmyg       mPFC        NAc       sACC
        # 75         69        165        141
#### explore human
geneExprs = log2(getRPKM(rse_gene,"Length")+1)
pca = prcomp(t(geneExprs))
pca_vars = getPcaVars(pca)
pca_vars_lab = paste0("PC", seq(along=pca_vars), ": ",
	pca_vars, "% Var Expl")


##########
## pc1 vs pc2
pdf("pdfs/PCA_plots_gene_Exprs_control.pdf",w=9)
par(mar=c(8,6,2,2),cex.axis=1.8,cex.lab=1.8)
palette(brewer.pal(9,"Pastel1"))
plot(pca$x, pch=21, bg=factor(rse_gene$Dx),cex=1.2,
	xlab = pca_vars_lab[1], ylab = pca_vars_lab[2])
legend("topleft", levels(factor(rse_gene$Dx)), col=1:length(levels(factor(rse_gene$Dx))), pch=15,cex=2)

dev.off()

###
table(rse_gene$Dx,rse_gene$Region)




pdf("pdfs/PCA_plots_gene_Exprs_all_regions.pdf",w=9)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
palette(brewer.pal(11,"Set3"))
## pc1 vs pc2
plot(pca$x, pch=21, bg=factor(rse_gene$Region),cex=1.5,
	xlab = pca_vars_lab[1], ylab = pca_vars_lab[2])
legend("topright", inset=c(-0.2,0), levels(factor(rse_gene$Region)), col=1:length(levels(factor(rse_gene$Region))), pch=15,cex=1.5)
plot(pca$x[,3:4], pch=21, bg=factor(rse_gene$Region),cex=1.5,
	xlab = pca_vars_lab[3], ylab = pca_vars_lab[4])
legend("topright", inset=c(-0.2,0), levels(factor(rse_gene$Region)), col=1:length(levels(factor(rse_gene$Region))), pch=15,cex=1.5)
plot(pca$x[,5:6], pch=21, bg=factor(rse_gene$Region),cex=1.5,
	xlab = pca_vars_lab[5], ylab = pca_vars_lab[6])
legend("topright", inset=c(-0.2,0), levels(factor(rse_gene$Region)), col=1:length(levels(factor(rse_gene$Region))), pch=15,cex=1.5)
plot(pca$x[,7:8], pch=21, bg=factor(rse_gene$Region),cex.axis=2, cex.lab = 2,
	xlab = pca_vars_lab[7], ylab = pca_vars_lab[8])
legend("topright", inset=c(-0.23,0), levels(factor(rse_gene$Region)), col=1:length(levels(factor(rse_gene$Region))), pch=15,cex=1.2)

dev.off()

pdf("pdfs/PCA_plots_gene_Exprs_all_regions_7_8.pdf",w=9)
par(mar=c(5.1, 4.7, 4.1, 8.1), xpd=TRUE)
palette(brewer.pal(11,"Set3"))
plot(pca$x[,7:8], pch=21, bg=factor(rse_gene$Region),cex.axis=2, cex.lab = 1.8,
	xlab = pca_vars_lab[7], ylab = pca_vars_lab[8])
legend("topright", inset=c(-0.235,0), levels(factor(rse_gene$Region)), col=1:length(levels(factor(rse_gene$Region))), pch=15,cex=1.2)
dev.off()

dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
dge = calcNormFactors(dge)

## mean-variance
mod = model.matrix(~Hab_ID + Age + Hab_ID*mitoRate,
	data=colData(rse_gene))

pdf("pdfs/voom_gene_region.pdf")
vGene = voom(dge,mod,plot=TRUE)
dev.off()

fitGene = lmFit(vGene)

ebGene = eBayes(fitGene)
outGene = topTable(ebGene,coef=2,
	p.value = 1,number=nrow(rse_gene),sort="none")

pdf("pdfs/hist_pval_regions.pdf")
hist(outGene$P.Value)
dev.off()

# The basic scatter plot: x is "log2FoldChange", y is "pvalue"
pdf("pdfs/volcano_plot_all.pdf")
ggplot(data=outGene, aes(x=logFC, y=adj.P.Val)) + geom_point()
dev.off()

table(outGene$adj.P.Val == 0)
table(outGene$adj.P.Val < 0.05)
table(outGene$adj.P.Val < 0.1)


sigGene =  topTable(ebGene,coef=2,
	p.value = 0.1,number=nrow(rse_gene))

pdf("pdfs/volcano_plot_sig.pdf")
ggplot(data=sigGene, aes(x=logFC, y=adj.P.Val)) + geom_point()
dev.off()

head(sigGene[,c("Symbol","logFC", "P.Value","AveExpr")])
sigGene[sigGene$logFC > 0,c("Symbol","logFC", "P.Value")]
sigGene[sigGene$logFC <  0,c("Symbol","logFC", "P.Value")]

write.csv(outGene, file = "tables/de_stats_allExprs_region.csv")
write.csv(sigGene, file = "tables/de_stats_fdr10_sorted_region.csv")


exprs = vGene$E[rownames(sigGene),]
exprsClean = cleaningY(exprs, mod, 2)
#cleanGeneExprs = cleaningY(geneExprs_hs, mod[,!is.na(eBGene$p.value[1,])], P=3)
pdf("pdfs/DE_boxplots_byRegion.pdf",w=10)
par(mar=c(8,6,4,2),cex.axis=1.8,cex.lab=1.8, cex.main=1.8)
palette(brewer.pal(12,"Dark2"))
for(i in 1:nrow(sigGene)) {
	yy = exprs[i,]
	boxplot(yy ~ rse_gene$Hab_ID, las=3,outline=FALSE,
		ylim=range(yy), ylab="Normalized log2 Exprs", xlab="",
		main = paste(sigGene$Symbol[i], "-", sigGene$gencodeID[i]),
		pch = 21,cex=1.3)
	ll = ifelse(sigGene$logFC[i] > 0, "topleft", "topright")
	legend(ll, paste0("adj.p=", signif(sigGene$adj.P.Val[i],3)), cex=1.3)
}
dev.off()

pdf("pdfs/DE_boxplots_byRegion_clean.pdf",w=10)
par(mar=c(8,6,4,2),cex.axis=1.8,cex.lab=1.8, cex.main=1.8)
palette(brewer.pal(12,"Dark2"))
for(i in 1:nrow(sigGene)) {
	yy = exprsClean[i,]
	boxplot(yy ~ rse_gene$Hab_ID, las=3,outline=FALSE,
		ylim=range(yy), ylab="Normalized log2 Exprs", xlab="",
		main = paste(sigGene$Symbol[i], "-", sigGene$gencodeID[i]),
		pch = 21,cex=1.3)
	ll = ifelse(sigGene$logFC[i] > 0, "topleft", "topright")
	legend(ll, paste0("adj.p=", signif(sigGene$adj.P.Val[i],3)), cex=1.3)
}
dev.off()


pdf(file = "pdfs/DE_boxplots_split_by_region.pdf")
par(mar=c(8,6,4,2),cex.axis=1.8,cex.lab=1.8, cex.main=1.8)
palette(brewer.pal(12,"Paired"))
for(i in nrow(sigGene)) {
	yy = exprsClean[i,]
	boxplot(yy ~ rse_gene$Region, las=3,outline=FALSE,
		ylim=range(yy), ylab="Normalized log2 Exprs", xlab="",
		main = paste(sigGene$Symbol[i], "-", sigGene$gencodeID[i]),
		pch = 21,cex=1.3)
	stripchart(yy ~ rse_gene$Region, vertical = TRUE,
    method = "jitter", add = TRUE, pch = 20, alpha=.1, col = "blue")
	ll = ifelse(sigGene$logFC[i] > 0, "topleft", "topright")
	legend(ll, paste0("adj.p=", signif(sigGene$adj.P.Val[i],3)), cex=1.3)
}
dev.off()

library(clusterProfiler)
library(org.Hs.eg.db)

## get significant genes by sign
sigGene = outGene[outGene$adj.P.Val < 0.005,]
sigGeneList = split(as.character(sigGene$EntrezID), sign(sigGene$logFC))
sigGeneList = lapply(sigGeneList, function(x) x[!is.na(x)])
geneUniverse = as.character(outGene$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

## do GO and KEGG
goBP_Adj <- compareCluster(sigGeneList, fun = "enrichGO",
	universe = geneUniverse, OrgDb = org.Hs.eg.db,
	ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 1,
	qvalueCutoff  = 1,	readable= TRUE)

goMF_Adj <- compareCluster(sigGeneList, fun = "enrichGO",
	universe = geneUniverse, OrgDb = org.Hs.eg.db,
	ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 1,
	qvalueCutoff  = 1,	readable= TRUE)

goCC_Adj <- compareCluster(sigGeneList, fun = "enrichGO",
	universe = geneUniverse, OrgDb = org.Hs.eg.db,
	ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 1,
	qvalueCutoff  = 1,	readable= TRUE)

kegg_Adj <- compareCluster(sigGeneList, fun = "enrichKEGG",
	universe = geneUniverse,  pAdjustMethod = "BH",
	pvalueCutoff  = 1, qvalueCutoff  = 1)

dir.create("rdas")
save(goBP_Adj, goCC_Adj, goMF_Adj, kegg_Adj,
	file = "rdas/gene_set_objects_p005.rda")

goList = list(BP = goBP_Adj, MF = goMF_Adj, CC = goCC_Adj, KEGG = kegg_Adj)
goDf = dplyr::bind_rows(lapply(goList, as.data.frame), .id = "Ontology")
goDf = goDf[order(goDf$pvalue),]

write.csv(goDf, file = "tables/geneSet_output_regions.csv", row.names=FALSE)

options(width=130)
goDf[goDf$p.adjust < 0.05, c(1:5)]

#############################
## compare to postmortem ####
#############################



