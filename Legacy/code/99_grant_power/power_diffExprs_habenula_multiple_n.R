###

library(RNASeqPower)
library(lattice)
library(RColorBrewer)
library(matrixStats)
library(recount)
# load data for SDs
load("/dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/Data/rdas/General/RSEs/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata")
geneMap<-rowData(rse_gene)
pd<-colData(rse_gene)
## adults
#aIndex = which(pd$age > 13)
#pd2 = pd[aIndex,]
#geneRpkm2 = log2(as.matrix(geneRpkm[,aIndex])+1)
geneRpkm<-getRPKM(rse_gene,"Length")
# filter
keepIndex = which((2^rowMeans(geneRpkm) -1) > 0.1)
geneRpkm2 = geneRpkm[keepIndex,]
geneMap2 = geneMap[keepIndex,]


## get numbers
sds = rowSds(geneRpkm2)
diffs = seq(0, 4, 0.01)
N1 = c(60,60)
N2 = c(70,70)
N3 = c(80,80)
N <- list(N1,N2,N3)
#####
dat<-list()
powerMat<-list()

for(i in 1:length(N)){
    print(i)
dat[[i]] = expand.grid(theSdQuant = c(0.5,0.75, 0.9),
	alpha0 = c(1e-4, 1e-5, 1e-6),	meanDiff = diffs)
dat[[i]]$power = 1 - pchisq(qchisq(1-dat[[i]]$alpha0, df=1),
	(dat[[i]]$meanDiff/(quantile(sds, prob= dat[[i]]$theSdQuant)*
		sqrt(1/N[[i]][1] + 1/N[[i]][2])))^2)
dat[[i]]$group = paste0(dat[[i]]$theSdQuant, ":", dat[[i]]$alpha0)
powerMat[[i]] = do.call("cbind", split(dat[[i]]$power, dat[[i]]$group))
rownames(powerMat[[i]]) = diffs
write.csv(powerMat[[i]], file=paste0("tables/powerMatrix_n",N[[i]][1],"exprs.csv"))
}


# plot

pdf("pdfs/rnaseq_power_calc_multiple_ns.pdf")
for(i in 1:3){
    print(i)
par(mar=c(5,6,2,2))
palette(brewer.pal(8,"Dark2"))
matplot(diffs, powerMat[[i]],type="l",
	col = rep(1:3, times=3), xlim = c(0.05,4),
	lwd = 5, lty = c(rep(3:1,each=3)), main = paste0("N =",N[[i]][1]),
	xlab="Log2 Fold Change", ylab="Power",
	cex.axis = 2, cex.lab=2)
abline(h=0.8, lty=2)
legend("bottomright", c("50th", "75th", "90th"),
	lty = 3:1, col="black",cex=1.75,lwd=3)
legend("right", legend=paste0("p=", eval(10^seq(-4, -6))),
	pch=15, col=1:3,cex=1.5,pt.cex=1.5)
}
dev.off()

