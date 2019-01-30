library(Matrix)

GM12878 <- readMM("GM12878/GRCh38/matrix.mtx")
rownames(GM12878) <- read.csv("GM12878/GRCh38/genes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1]
GM18502 <- readMM("GM18502/GRCh38/matrix.mtx")
rownames(GM18502) <- read.csv("GM18502/GRCh38/genes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1]
MIX <- readMM("mix/GRCh38/matrix.mtx")
rownames(MIX) <- read.csv("mix/GRCh38/genes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1]

MT <- read.csv("MT.txt", sep = "\t", stringsAsFactors = FALSE)

png("Figure2.png", width = 3000, height = 900, res = 300)
layout(matrix(c(1,2,3,1,2,3,1,2,3,4,5,6), 4, 3, byrow = TRUE))
# GM12878
MTp <- apply(GM12878[rownames(GM12878) %in% MT[,1],],2,sum)/apply(GM12878,2,sum)
Lsize1 <- apply(GM12878,2,sum)
par(mar=c(0.15,3,1,1), mgp=c(1.6,0.5,0))
d1 <- cbind(Lsize1,MTp)
plot(d1, ylim = c(0,0.4), ylab = "", main = "GM12878", pch = 20, col=densCols(d1),xaxt='n', xlab = "", las=1)
mtext("Mitochondrial Transcript Percent", side = 2, line = 2, cex = 0.8)
abline(h=0.1, col="red", lty=2)
minO <- min(boxplot.stats(Lsize1)$out)
lines(x = c(minO,minO), y = c(-1,0.33), lty = 3)
legend("topright", legend = c("Mitochondrial Proportion Threshold", "Library Size Threshold"), bty = "n", lty = c(2,3), col = c("red", "black"))

#Mixture
MTp <- apply(MIX[rownames(MIX) %in% MT[,1],],2,sum)/apply(MIX,2,sum)
Lsize3 <- apply(MIX,2,sum)
d3 <- cbind(Lsize3,MTp)
par(mar=c(0.15,3,1,1), mgp=c(1.6,0.5,0))
plot(d3, ylim = c(0,0.4),  ylab = "", main = "MIXTURE GM12878/GM18502", pch = 20, col=densCols(d3),xaxt='n', xlab = "", las=1)
mtext("Mitochondrial Transcript Percent", side = 2, line = 2, cex = 0.8)
abline(h=0.1, col="red", lty=2)
minO <- min(boxplot.stats(Lsize3)$out)
lines(x = c(minO,minO), y = c(-1,0.33), lty = 3)
legend("topright", legend = c("Mitochondrial Proportion Threshold", "Library Size Threshold"), bty = "n", lty = c(2,3), col = c("red", "black"))

#GM18502
MTp <- apply(GM18502[rownames(GM18502) %in% MT[,1],],2,sum)/apply(GM18502,2,sum)
Lsize2 <- apply(GM18502,2,sum)
par(mar=c(0.15,3,1,1), mgp=c(1.6,0.5,0))
d2 <- cbind(Lsize2,MTp)
plot(d2, ylim = c(0,0.4), ylab = "", main = "GM18502", pch = 20, col=densCols(d2),xaxt='n', xlab = "", las=1)
mtext("Mitochondrial Transcript Percent", side = 2, line = 2, cex = 0.8)
abline(h=0.1, col="red", lty=2)
minO <- min(boxplot.stats(Lsize2)$out)
lines(x = c(minO,minO), y = c(-1,0.33), lty = 3)
legend("topright", legend = c("Mitochondrial Proportion Threshold", "Library Size Threshold"), lty = c(2,3), col = c("red", "black"), bty="n",bg = "white")

#GM12878
par(mar=c(3,3,0,1), mgp=c(1.5,0.5,0))
boxplot(Lsize1, las = 1, horizontal = TRUE, xlab = "", pch =20)
mtext("Library Size", side = 1, line = 1.5, cex=0.8)

#Mixture
par(mar=c(3,3,0,1), mgp=c(1.5,0.5,0))
boxplot(Lsize3, las = 1, horizontal = TRUE, xlab = "", pch =20)
mtext("Library Size", side = 1, line = 1.5, cex=0.8)

#GM18502
par(mar=c(3,3,0,1), mgp=c(1.5,0.5,0))
boxplot(Lsize2, las = 1, horizontal = TRUE, xlab = "", pch =20)
mtext("Library Size", side = 1, line = 1.5, cex=0.8)

dev.off()