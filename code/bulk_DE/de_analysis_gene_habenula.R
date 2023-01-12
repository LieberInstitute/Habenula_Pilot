###
library(SummarizedExperiment)
library(recount)
library(edgeR)
library(limma)
library(jaffelab)
library(RColorBrewer)

dir.create("pdfs")
dir.create("tables")

### load data
load("count_data/rse_gene_Roche_Habenula_qcAndAnnotated_n69.Rdata")

## keep subset of samples


## filter for expressed
rse_gene = rse_gene[rowMeans(getRPKM(rse_gene,"Length")) > 0.2,]

##############
## metrics ###
##############

pd<- colData(rse_gene)

## check if ratios of cell changed by batch
pdf("pdfs/variable_metrics_dx.pdf")
boxplot(rse_gene$RIN ~ rse_gene$PrimaryDx,xlab="Dx")
boxplot(rse_gene$AgeDeath ~ rse_gene$PrimaryDx,xlab="Dx")
boxplot(rse_gene$mitoRate ~ rse_gene$PrimaryDx,las=3,xlab="Dx")
boxplot(rse_gene$totalAssignedGene ~ rse_gene$PrimaryDx,las=3,xlab="Dx")
dev.off()


#### explore human
geneExprs = log2(getRPKM(rse_gene,"Length")+1)
pca = prcomp(t(geneExprs))
pca_vars = getPcaVars(pca)
pca_vars_lab = paste0("PC", seq(along=pca_vars), ": ",
	pca_vars, "% Var Expl")


##########
pdf("pdfs/PCA_plots_gene_Exprs_dx.pdf",w=9)
par(mar=c(8,6,2,2),cex.axis=1.8,cex.lab=1.8)
palette(brewer.pal(4,"Dark2"))

## pc1 vs pc2
plot(pca$x, pch=21, bg=factor(rse_gene$PrimaryDx),cex=1.2,
	xlab = pca_vars_lab[1], ylab = pca_vars_lab[2])
legend("bottomleft", levels(factor(rse_gene$PrimaryDx)), col=1:2, pch=15,cex=2)

dev.off()


####################
## modeling ########
####################

dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
dge = calcNormFactors(dge)

## mean-variance
mod = model.matrix(~PrimaryDx + mitoRate +AgeDeath,
	data=colData(rse_gene))

pdf("pdfs/voom_gene_dx.pdf")
vGene = voom(dge,mod,plot=TRUE)
dev.off()



fitGeneDupl = lmFit(vGene)

ebGeneDupl = eBayes(fitGeneDupl)
outGeneDupl = topTable(ebGeneDupl,coef=2,
	p.value = 1,number=nrow(rse_gene),sort="none")

pdf("pdfs/hist_pval_dx.pdf")
hist(outGeneDupl$P.Value)
dev.off()

table(outGeneDupl$adj.P.Val < 0.05)
table(outGeneDupl$adj.P.Val < 0.1)

sigGeneDupl =  topTable(ebGeneDupl,coef=2,
	p.value = 0.2,number=nrow(rse_gene))

sigGeneDupl[,c("Symbol","logFC", "P.Value","AveExpr")]
sigGeneDupl[sigGeneDupl$logFC > 0,c("Symbol","logFC", "P.Value")]
sigGeneDupl[sigGeneDupl$logFC <  0,c("Symbol","logFC", "P.Value")]

write.csv(outGeneDupl, file = "tables/de_stats_allExprs_dx.csv")
write.csv(sigGeneDupl, file = "tables/de_stats_fdr10_sorted_dx.csv")

###################
## check plots ####
###################

exprs = vGene$E[rownames(sigGeneDupl),]
#exprsClean = cleaningY(exprs, mod, 2)


### make boxplots
# cleanGeneExprs = cleaningY(geneExprs_hs, mod[,!is.na(eBGene$p.value[1,])], P=3)
pdf("pdfs/DE_boxplots_byDx.pdf",w=10)
par(mar=c(8,6,4,2),cex.axis=1.8,cex.lab=1.8, cex.main=1.8)
palette(brewer.pal(4,"Dark2"))
for(i in 1:nrow(sigGeneDupl)) {
	yy = exprs
	boxplot(yy ~ rse_gene$PrimaryDx, las=3,outline=FALSE,
		ylim=range(yy), ylab="Normalized log2 Exprs", xlab="",
		main = paste(sigGeneDupl$Symbol, "-", sigGeneDupl$gencodeID))
	ll = ifelse(sigGeneDupl$logFC > 0, "topleft", "topright")
	legend(ll, paste0("p=", signif(sigGeneDupl$P.Value[i],3)), cex=1.3)
}
dev.off()

### gene ontology
library(clusterProfiler)
library(org.Hs.eg.db)

## get significant genes by sign
sigGene = outGeneDupl[outGeneDupl$P.Value < 0.005,]
sigGeneList = split(as.character(sigGene$EntrezID), sign(sigGene$logFC))
sigGeneList = lapply(sigGeneList, function(x) x[!is.na(x)])
geneUniverse = as.character(outGeneDupl$EntrezID)
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
	file = "rdas/gene_set_objects_p005_dx.rda")

goList = list(BP = goBP_Adj, MF = goMF_Adj, CC = goCC_Adj, KEGG = kegg_Adj)
goDf = dplyr::bind_rows(lapply(goList, as.data.frame), .id = "Ontology")
goDf = goDf[order(goDf$pvalue),]

write.csv(goDf, file = "tables/geneSet_output_dx.csv", row.names=FALSE)

options(width=130)
goDf[goDf$p.adjust < 0.05, c(1:5)]



