library(Matrix)
library(geuvPack)

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

GM12878 <- readMM("GM12878/GRCh38/matrix.mtx")
rownames(GM12878) <- read.csv("GM12878/GRCh38/genes.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)[,2]
GM12878 <- cbind(apply(GM12878,1,sum))

GM18502 <- readMM("GM18502/GRCh38/matrix.mtx")
rownames(GM18502) <- read.csv("GM18502/GRCh38/genes.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)[,2]
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


pdf("Figure5.pdf", width = 3.5*4, height = 3.5)
par(mfrow=c(1,4))
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
dev.off()

# "log10(RNA-seq FPKM)\nCEU:SRR038449",
# "log10(RNA-seq FPKM)\nYRI:SRR6355954",
png("Figure5PCA.png", width = 900*4, height = 900, res = 300)
par(mfrow=c(1,4))
par(mar=c(3,3,2.5,1), mgp=c(1.5,0.5,0))
yLim = c(-10,10)
o <- cbind(all$scGM12878, all$bGM12878)
o <- prcomp(o)
summary(o)
plot(o$x, pch = 16, col=densCols(o$x), xlab="PC1 (89.65%)", ylab="PC2 (10.35%)", las = 1, ylim=yLim)
mtext("SC vs. BULK", side = 3, line = 1, cex = 1,font = 2)
mtext("GM12878 vs. GM12878", side = 3, line = 0.1, cex = 0.7)


o <- cbind(all$scGM18502, all$bGM18502)
o <- prcomp(o)
summary(o)
plot(o$x, pch = 16, col=densCols(o$x), xlab="PC1 (79.91%)", ylab="PC2 (20.09%)", las = 1, ylim=yLim)
mtext("SC vs. BULK", side = 3, line = 1, cex = 1,font = 2)
mtext("GM18502 vs. GM18502", side = 3, line = 0.1, cex = 0.7)



o <- cbind(all$scGM12878, all$CEU)
o <- prcomp(o)
summary(o)
plot(o$x, pch = 16, col=densCols(o$x), xlab="PC1 (89.06%)", ylab="PC2 (10.94%)", las = 1, ylim=yLim)
mtext("SC vs. BULK AVERAGE", side = 3, line = 1, cex = 1,font = 2)
mtext("GM12878 vs. CEU Population", side = 3, line = 0.1, cex = 0.7)

o <- cbind(all$scGM18502, all$YRI)
o <- prcomp(o)
summary(o)
plot(o$x, pch = 16, col=densCols(o$x), xlab="PC1 (89.58%)", ylab="PC2 (10.42%)", las = 1, ylim=yLim)
mtext("SC vs. BULK AVERAGE", side = 3, line = 1, cex = 1,font = 2)
mtext("GM18502 vs. YRI Population", side = 3, line = 0.1, cex = 0.7)
dev.off()