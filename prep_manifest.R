###
library(jaffelab)

fqPath = "/dcl02/lieber/ajaffe/Roche_Habenula/preprocessed_data/FASTQ/"
# fqPath = "/dcl01/ajaffe/data/Nina/Roche_20181023/fastq/"
fqFiles = list.files(fqPath, recur=TRUE)

fqFiles = fqFiles[!grepl("^Und", fqFiles)]

leftReads = fqFiles[grep("_R1_", fqFiles)]
rightReads = fqFiles[grep("_R2_", fqFiles)]

man = data.frame(leftReads = paste0(fqPath, leftReads),
	leftMd5 = 0,rightReads = paste0(fqPath, rightReads), 
	rightMd5 = 0, SampleID = ss(leftReads, "_"),stringsAsFactors=FALSE) 
	
write.table(man, file="preprocessed_data/samples.manifest", 
	sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)
	