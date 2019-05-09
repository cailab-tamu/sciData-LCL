library(Matrix)
library(geuvPack)
library(preprocessCore)
library(limma)

data("geuFPKM")
geneLength <- read.csv("hg38_genes_length.tsv", sep = "\t", stringsAsFactors = FALSE)
gL <- geneLength[,4]
names(gL) <- geneLength[,1]
geneInfo <- cbind(geuFPKM@rowRanges$gene_id,geuFPKM@rowRanges$gene_name)
rownames(geneInfo) <- geneInfo[,1]
popInfo <- read.csv("populationInfo/GD660.GeneQuantCount.txt.gz", sep = "\t")
gNames <- popInfo[,2]
popInfo <- popInfo[gNames %in% geneInfo[,1],]
gNames <- gNames[gNames %in% geneInfo[,1]]
gNames <- geneInfo[as.vector(gNames),2]
colnames(popInfo) <- unlist(lapply(strsplit(colnames(popInfo), "\\."), function(X){X[1]}))
popInfo <- popInfo[,colnames(geuFPKM)]
rownames(popInfo) <- make.unique(gNames)

CEU <- popInfo[,geuFPKM@colData$popcode == "CEU",]
CEU <- CEU[rownames(CEU) %in% geneLength[,1],]
CEU <- CEU/gL[rownames(CEU)]
CEU <- t(t(CEU)/(colSums(CEU)/1e6))
CEU <- cbind(rowMeans(CEU))

YRI <- popInfo[,geuFPKM@colData$popcode == "YRI",]
YRI <- YRI[rownames(YRI) %in% geneLength[,1],]
YRI <- YRI/gL[rownames(YRI)]
YRI <- t(t(YRI)/(colSums(YRI)/1e6))
YRI <- cbind(rowMeans(YRI))

bulkRNAseq <- read.csv("bulkFQ/RNAseq.tsv", sep = "\t")
bulkRNAseq <- bulkRNAseq[rownames(bulkRNAseq) %in% geneLength[,1],]

GM12878 <- readMM("../GM12878/GRCh38/matrix.mtx")
rownames(GM12878) <- read.csv("../GM12878/GRCh38/genes.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)[,2]
GM12878 <- cbind(apply(GM12878,1,sum))

GM18502 <- readMM("../GM18502/GRCh38/matrix.mtx")
rownames(GM18502) <- read.csv("../GM18502/GRCh38/genes.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)[,2]
GM18502 <- cbind(apply(GM18502,1,sum))

sharedGenes <- table(c(make.unique(rownames(GM12878)),
                       make.unique(rownames(GM18502)),
                       make.unique(rownames(CEU)),
                       make.unique(rownames(YRI)),
                       make.unique(rownames(bulkRNAseq))))
sharedGenes <- names(sharedGenes[sharedGenes == 5])

bulkRNAseq <- bulkRNAseq[sharedGenes,]
GM18502 <- GM18502[sharedGenes,]
GM12878 <- GM12878[sharedGenes,]
CEU <- CEU[sharedGenes,]
YRI <- YRI[sharedGenes,]

all <- cbind(bulkRNAseq)
all <- all/gL[rownames(all)]
all <- t(t(all)/(apply(all,2,sum)/1e6))
all <- cbind(all,GM12878,GM18502,CEU,YRI)
all <- all[rowSums(all == 0) == 0,]

colnames(all) <- c("bGM12878", "bGM18502", "scGM12878", "scGM18502", "CEU", "YRI")
all <- log(all)

cor(all,method = "sp")
all <- as.data.frame(all)

pdf("Figure5.pdf", width = 3.3*3, height = 2.5*2)
layout(matrix(c(1,2,5,5,3,4,5,5), nrow = 2, ncol = 4, byrow = TRUE))
par(mar=c(3,3,2.5,1), mgp=c(1.5,0.5,0))
plot(all[,"scGM12878"],all[,"bGM12878"], pch = 16, 
     col=densCols(cbind(all[,"scGM12878"],all[,"bGM12878"])), 
     xlab = "", 
     ylab = "", las=1)
corV <- round(cor(all[,"scGM12878"],all[,"bGM12878"], method = 'sp'),2)
mtext(parse(text = "log[10] (paste('total UMI'))"), side = 1,  line = 2)
mtext(parse(text = "log[10] (paste('TPM'))"), side = 2,  line = 1.5)
mtext("SC vs. BULK", side = 3, line = 1, cex = 1,font = 2)
mtext("GM12878 vs. GM12878", side = 3, line = 0.1, cex = 0.7)
legend("topleft", legend = parse(text = paste('rho == ', corV)), col = "red", lty = 1, bty = "n")
abline(lm(all[,"bGM12878"]~all[,"scGM12878"]), col = "red")

plot(all[,"scGM18502"],all[,"bGM18502"], pch = 16, 
     col=densCols(cbind(all[,"scGM18502"],all[,"bGM18502"])), 
     xlab="", ylab="",las=1)
mtext(parse(text = "log[10] (paste('total UMI'))"), side = 1,  line = 2)
mtext(parse(text = "log[10] (paste('TPM'))"), side = 2,  line = 1.5)
mtext("SC vs. BULK", side = 3, line = 1, cex = 1,font = 2)
mtext("GM18502 vs. GM18502", side = 3, line = 0.1, cex = 0.7)
corV <- round(cor(all[,"scGM18502"],all[,"bGM18502"], method = 'sp'),2)
legend("topleft", legend = parse(text = paste('rho == ', corV)), col = "red", lty = 1, bty = "n")
abline(lm(all[,"bGM18502"]~all[,"scGM18502"]), col = "red")

plot(all[,"scGM12878"],all[,"CEU"], pch = 16, 
     col=densCols(cbind(all[,"scGM12878"],all[,"CEU"])), 
     xlab = "", 
     ylab = "", las=1)
corV <- round(cor(all[,"scGM12878"],all[,"CEU"], method = 'sp'),2)
mtext(parse(text = "log[10] (paste('total UMI'))"), side = 1,  line = 2)
mtext(parse(text = "log[10](paste('average TPM'))"), side = 2,  line = 1.5)
mtext("SC vs. BULK AVERAGE", side = 3, line = 1, cex = 1,font = 2)
mtext("GM12878 vs. CEU Population", side = 3, line = 0.1, cex = 0.7)
legend("topleft", legend = parse(text = paste('rho == ', corV)), col = "red", lty = 1, bty = "n")
abline(lm(all[,"CEU"]~all[,"scGM12878"]), col = "red")

plot(all[,"scGM18502"],all[,"YRI"], pch = 16, 
     col=densCols(cbind(all[,"scGM18502"],all[,"YRI"])),
     xlab = "", 
     ylab = "", las=1)
corV <- round(cor(all[,"scGM18502"],all[,"YRI"], method = 'sp'),2)
mtext(parse(text = "log[10] (paste('total UMI'))"), side = 1,  line = 2)
mtext(parse(text = "log[10](paste('average TPM'))"), side = 2, line = 1.5)
mtext("SC vs. BULK AVERAGE", side = 3, line = 1, cex = 1,font = 2)
mtext("GM18502 vs. YRI Population", side = 3, line = 0.1, cex = 0.7)
legend("topleft", legend = parse(text = paste('rho == ', corV)), col = "red", lty = 1, bty = "n")
abline(lm(all[,"YRI"]~all[,"scGM18502"]), col = "red")

geneLength <- read.csv("hg38_genes_length.tsv", sep = "\t", stringsAsFactors = FALSE)
gL <- geneLength[,4]
names(gL) <- geneLength[,1]
geneInfo <- cbind(geuFPKM@rowRanges$gene_id,geuFPKM@rowRanges$gene_name)
rownames(geneInfo) <- geneInfo[,1]
popInfo <- read.csv("populationInfo/GD660.GeneQuantCount.txt.gz", sep = "\t")
gNames <- popInfo[,2]
popInfo <- popInfo[gNames %in% geneInfo[,1],]
gNames <- gNames[gNames %in% geneInfo[,1]]
gNames <- geneInfo[as.vector(gNames),2]
colnames(popInfo) <- unlist(lapply(strsplit(colnames(popInfo), "\\."), function(X){X[1]}))
popInfo <- popInfo[,colnames(geuFPKM)]
rownames(popInfo) <- make.unique(gNames)

CEU <- popInfo[,geuFPKM@colData$popcode == "CEU",]
CEU <- CEU[rownames(CEU) %in% geneLength[,1],]

YRI <- popInfo[,geuFPKM@colData$popcode == "YRI",]
YRI <- YRI[rownames(YRI) %in% geneLength[,1],]

bulkRNAseq <- read.csv("bulkFQ/RNAseq.tsv", sep = "\t")
bulkRNAseq <- bulkRNAseq[rownames(bulkRNAseq) %in% geneLength[,1],]

GM12878 <- readMM("../GM12878/GRCh38/matrix.mtx")
rownames(GM12878) <- read.csv("../GM12878/GRCh38/genes.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)[,2]
GM12878 <- cbind(apply(GM12878,1,sum))

GM18502 <- readMM("../GM18502/GRCh38/matrix.mtx")
rownames(GM18502) <- read.csv("../GM18502/GRCh38/genes.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)[,2]
GM18502 <- cbind(apply(GM18502,1,sum))

sharedGenes <- table(c(make.unique(rownames(GM12878)),
                       make.unique(rownames(GM18502)),
                       make.unique(rownames(CEU)),
                       make.unique(rownames(YRI)),
                       make.unique(rownames(bulkRNAseq))))
sharedGenes <- names(sharedGenes[sharedGenes == 5])

bulkRNAseq <- bulkRNAseq[sharedGenes,]
GM18502 <- GM18502[sharedGenes,]
GM12878 <- GM12878[sharedGenes,]

CEU <- CEU[sharedGenes,]
YRI <- YRI[sharedGenes,]

colnames(CEU) <- paste0("CEU_",colnames(CEU))
colnames(YRI) <- paste0("YRI_", colnames(YRI))

all <- cbind(bulkRNAseq,GM12878,GM18502, CEU, YRI)
dType <- c(rep("b",2), rep("sc",2), rep("geu", ncol(CEU)), rep("geu", ncol(YRI)))
dType <- as.factor(dType)
dType <- relevel(dType,ref="sc")
all <- removeBatchEffect(all,batch = dType)

all <- all/gL[rownames(all)]
all <- round(t(t(all)/(apply(all,2,sum))) * 1e6)
all <- all[rowSums(all == 0) == 0,]

colnames(all)[1:4] <- c("bGM12878", "bGM18502", "scGM12878", "scGM18502")
allN <- normalize.quantiles(all)
colnames(allN) <- colnames(all)
rownames(allN) <- rownames(all)
all <- round(allN)

cList <- rep("gray", ncol(all))
cList[grepl("YRI", colnames(all))] <- rgb(1,0,0,0.5)
cList[grepl("scGM18502", colnames(all))] <- rgb(1,0,0,1)
cList[grepl("bGM18502", colnames(all))] <- rgb(1,0,0,1)
cList[grepl("CEU", colnames(all))] <- rgb(0,0,1,0.5)
cList[grepl("scGM12878", colnames(all))] <- rgb(0,0,1,1)
cList[grepl("bGM12878", colnames(all))] <- rgb(0,0,1,1)

par(mar=c(3,3,2.5,1), mgp=c(1.5,0.5,0))
outData <- prcomp(scale(t(all)))$x[,1:2]
colnames(outData) <- c("PC1","PC2")
plot(outData, pch=16, col=cList, las=1, xlab = "", ylab = "", cex = 1.5)
mtext(parse(text = "PC1"), side = 1,  line = 2)
mtext(parse(text = "PC2"), side = 2,  line = 1.5)
mtext("PCA", side = 3, line = 1, cex = 1,font = 2)
mtext("ALL SAMPLES", side = 3, line = 0.1, cex = 0.7)
text(outData[1,1],outData[1,2]+5,"BULK-GM12878")
text(outData[2,1],outData[2,2]-5,"BULK-GM18502")
text(outData[3,1],outData[3,2]+5,"SC-GM12878")
text(outData[4,1],outData[4,2]-5,"SC-GM18502")
abline(h=0, lty=2, col=rgb(0,0,0,0.1))
legend("bottomright", legend = c("CEU", "YRI"), bty = "n", horiz = TRUE, pch = 16, col = c(rgb(0,0,1,0.5),rgb(1,0,0,0.5)))
dev.off()
