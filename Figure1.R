library(Matrix)
# Single Cell Matrices
GM12878 <- readMM("GM12878/GRCh38/matrix.mtx")
rownames(GM12878) <- read.csv("GM12878/GRCh38/genes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1]
GM18502 <- readMM("GM18502/GRCh38/matrix.mtx")
rownames(GM18502) <- read.csv("GM18502/GRCh38/genes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1]
MIX <- readMM("mix/GRCh38/matrix.mtx")
rownames(MIX) <- read.csv("mix/GRCh38/genes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1]

# Cell Assignation from NNMF
cellL <- read.csv("cellLineage.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)[,2]

# Growth Curve
allCC <- read.csv("CC.csv", header = FALSE)
CC <- allCC[,2:5]
n <- apply(!apply(CC, 1, is.na),2,sum)
sd <- apply(CC, 1, function(x){sd(x, na.rm=TRUE)})
se <- sd/sqrt(n)
avg <- apply(CC, 1, function(x){mean(x, na.rm=TRUE)})
names(avg) <- rep(c(0:4),2)

# Figure 1 Output file
png(filename = "Figure1.png", width = 3000, height = 900, res = 300)
par(mar=c(3,4,1,1), mgp=c(3,0.5,0), mfrow=c(1,3))
  # Growth Curve
plot(avg[1:5], type = "b", pch = 16, las=1, ylim = c(0,1e6), col= "black", ylab = "", xaxt="n", xlab="")
axis(side = 1, at = 1:5, labels = 0:4)
points(avg[6:10], type = "b", pch = 15, col = "red")
arrows(x0 = 1:5, y0 = (avg-se)[1:5], y1 = (avg+se)[1:5], x1 = 1:5, angle = 90, length = 0.05,code = 3)
arrows(x0 = 1:5, y0 = (avg-se)[6:10], y1 = (avg+se)[6:10], x1 = 1:5, angle = 90, length = 0.05,code = 3, col="red")
legend("bottomright", legend = c("GM12878", "GM18502"), bty = "n", pch = c(16,15), col = c("black", "red"))
mtext("Time (Days)", side=1, line=1.5, cex=.8,las=1)
mtext("Viable Cell Count / mL", side = 2, line = 3, cex=.8)

rNames <- table(c(rownames(GM12878),rownames(MIX)))
rNames <- names(rNames[rNames == 2])
A <- apply(GM12878[rNames,],1,mean)
B <- apply(MIX[rNames,cellL == "EUR"],1,mean)
C <- cbind(A,B)
C <- C[apply(C == 0, 1, sum) < 1,]
C <- log(C)
par(mar=c(3,3,1,1), mgp=c(2.2,0.5,0))
  # GM12878
plot(C, pch = 16, las = 1,ylab="", xlab="", col=densCols(C), xlim = c(-8.5,7), ylim = c(-8.5,7), main = "GM12878")
legend("bottomright", legend = c(expression(paste(rho ," = 0.993")), "n = 18410"), bty="n")
mtext("Expression Level from Pure Sample, log", side = 2, line = 1.5, cex=.8)
mtext("Expression Level from Mixture, log", side=1, line=1.5, cex=0.8)
abline(0,1, col = "red")

rNames <- table(c(rownames(GM18502),rownames(MIX)))
rNames <- names(rNames[rNames == 2])
A <- apply(GM18502[rNames,],1,mean)
B <- apply(MIX[rNames,cellL == "AFR"],1,mean)
C <- cbind(A,B)
C <- C[apply(C == 0, 1, sum) < 1,]
X <- C[,1]
W <- as.numeric(solve(t(X) %*% X) %*% X)
C <- log(C)
par(mar=c(3,3,1,1), mgp=c(2.2,0.5,0))
  # GM18502
plot(C, pch = 16, ylab = "", las = 1, xlab="", col=densCols(C), xlim = c(-8.5,7), ylim = c(-8.5,7), main="GM18502")
legend("bottomright", legend = c(expression(paste(rho ," = 0.993")), "n = 18888"), bty="n")
mtext("Expression Level from Pure Sample, log", side = 2, line = 1.5, cex=.8)
mtext("Expression Level from Mixture, log", side=1, line=1.5, cex=0.8)
abline(0,1, col = "red")
dev.off()
