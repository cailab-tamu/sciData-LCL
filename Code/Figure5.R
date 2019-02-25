# Image scale function from https://gist.githubusercontent.com/menugget/7689145/raw/2381856aaa12ea5db0791fc4e4ac1eed8189cdfe/image.scale.2.r
image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}

library(Seurat)

GM12878 <- Read10X("GM12878/GRCh38/")
GM18502 <- Read10X("GM18502/GRCh38/")

GM12878 <- CreateSeuratObject(GM12878)
GM12878 <- NormalizeData(GM12878)
GM12878 <- ScaleData(GM12878, genes.use = unlist(cc.genes))
GM18502 <- CreateSeuratObject(GM18502)
GM18502 <- NormalizeData(GM18502)
GM18502 <- ScaleData(GM18502, genes.use = unlist(cc.genes))

GM12878 <- CellCycleScoring(GM12878, g2m.genes = cc.genes$g2m.genes, s.genes = cc.genes$s.genes)
GM18502 <- CellCycleScoring(GM18502, g2m.genes = cc.genes$g2m.genes, s.genes = cc.genes$s.genes)

CEU <- as.matrix(GM12878@raw.data)
YRI <- as.matrix(GM18502@raw.data)

CEU <- CEU[rownames(CEU) %in% rownames(YRI),]
YRI <- YRI[rownames(YRI) %in% rownames(CEU),]

YRI <- YRI[rownames(CEU),]

CEU_S <- CEU[,GM12878@meta.data$Phase == "S"]
CEU_S <- CEU_S[order(rowSums(CEU_S),decreasing = TRUE),order(colSums(CEU_S),decreasing = FALSE)]

CEU_G1 <- CEU[,GM12878@meta.data$Phase == "G1"]
CEU_G1 <- CEU_G1[order(rowSums(CEU_G1),decreasing = TRUE),order(colSums(CEU_G1),decreasing = FALSE)]

CEU_G2 <- CEU[,GM12878@meta.data$Phase == "G2M"]
CEU_G2 <- CEU_G2[order(rowSums(CEU_G2),decreasing = TRUE),order(colSums(CEU_G2),decreasing = FALSE)]

sGenes <- sort(unique(c(rownames(CEU_S)[1:200],rownames(CEU_G1)[1:150],rownames(CEU_G2)[1:200]))[1:200])

CEU_S <- CEU_S[sGenes,]
CEU_S <- CEU_S[order(rowSums(CEU_S),decreasing = TRUE),order(colSums(CEU_S),decreasing = TRUE)]

CEU_G1 <- CEU_G1[sGenes,]
CEU_G1 <- CEU_G1[order(rowSums(CEU_G1),decreasing = TRUE),order(colSums(CEU_G1),decreasing = TRUE)]

CEU_G2 <- CEU_G2[sGenes,]
CEU_G2 <- CEU_G2[order(rowSums(CEU_G2),decreasing = TRUE),order(colSums(CEU_G2),decreasing = TRUE)]
####

YRI_S <- YRI[,GM18502@meta.data$Phase == "S"]
YRI_S <- YRI_S[order(rowSums(YRI_S),decreasing = TRUE),order(colSums(YRI_S),decreasing = FALSE)]

YRI_G1 <- YRI[,GM18502@meta.data$Phase == "G1"]
YRI_G1 <- YRI_G1[order(rowSums(YRI_G1),decreasing = TRUE),order(colSums(YRI_G1),decreasing = FALSE)]

YRI_G2 <- YRI[,GM18502@meta.data$Phase == "G2M"]
YRI_G2 <- YRI_G2[order(rowSums(YRI_G2),decreasing = TRUE),order(colSums(YRI_G2),decreasing = FALSE)]

sGenes <- sort(unique(c(rownames(YRI_S)[1:200],rownames(YRI_G1)[1:200],rownames(YRI_G2)[1:150]))[1:200])

YRI_S <- YRI_S[sGenes,]
YRI_S <- YRI_S[order(rowSums(YRI_S),decreasing = TRUE),order(colSums(YRI_S),decreasing = TRUE)]

YRI_G1 <- YRI_G1[sGenes,]
YRI_G1 <- YRI_G1[order(rowSums(YRI_G1),decreasing = TRUE),order(colSums(YRI_G1),decreasing = TRUE)]

YRI_G2 <- YRI_G2[sGenes,]
YRI_G2 <- YRI_G2[order(rowSums(YRI_G2),decreasing = TRUE),order(colSums(YRI_G2),decreasing = TRUE)]


CEU <- cbind(CEU_G1,CEU_S,CEU_G2)
YRI <- cbind(YRI_G1, YRI_S, YRI_G2)

O <- cbind(CEU,YRI)
O <- log10(O+1)
sGenes <- sort(rowSums(O), decreasing = TRUE)
O <- O[names(sGenes),]
O[O > 4] <- 4

sGenes <- sort(rowSums(O), decreasing = FALSE)
O <- O[names(sGenes),]
n <- nrow(O)
newCol <- colorRampPalette(c("blue4","blue3","blue2","blue1","blue", "white", "red", "red1","red2","red3","red4"))
  #colorRampPalette(c("red4","red3","red2","red1","red","white","blue","blue1","blue2","blue3","blue4"))
png("F5.png", width = 3000, height = 900, res = 300)
layout(matrix(c(rep(1,19),2), 1, 20, byrow = TRUE))
opar <- par()
par(mar=c(3,2,3,1),mgp=c(1.5,0.5,0))
image(x=t(as.matrix(O)), useRaster=TRUE, col=newCol(50), yaxt="n", xaxt="n", bty = "n")
CD <- ncol(CEU)/ncol(O)
abline(v = CD, lwd = 2, col = "white")
atV <- c(ncol(CEU_G1)/2,
         ncol(CEU_G1)+(ncol(CEU_S)/2),
         ncol(CEU_G1)+ncol(CEU_S)+(ncol(CEU_G2)/2),
         ncol(CEU_G1)+ncol(CEU_S)+ncol(CEU_G2)+(ncol(YRI_G1)/2),
         ncol(CEU_G1)+ncol(CEU_S)+ncol(CEU_G2)+ncol(YRI_G1)+(ncol(YRI_S)/2),
         ncol(CEU_G1)+ncol(CEU_S)+ncol(CEU_G2)+ncol(YRI_G1)+ncol(YRI_S)+(ncol(YRI_G2)/2)
         )/ncol(O)
axis(side = 3,at = atV, labels = c("G1", "S", "G2M", "G1", "S", "G2M"), font = 2, cex.axis=1)
atL <- c(4633,1113,1299,2525,1367,1297)
axis(side = 1, labels = atL, at = atV)
mtext("GM12878",side = 3, at = (ncol(CEU)/2)/ncol(O),line = 1.5, font = 2)
mtext("GM18502",side = 3, at = (ncol(CEU)+(ncol(YRI)/2))/ncol(O),line = 1.5, font = 2)
mtext("Top 200 Genes",side = 2, line = 0.5, cex = 0.8)
mtext("Number of Cells",side = 1, line = 1.5, cex = 0.8)
# axis(side=2, at= (1/(n-1))*seq(1,n,2), labels = sGenes[seq(1,n,2)], las=2, cex.axis = 0.7)
# axis(side=4, at= ((1/(n-1))*seq(0,(n-2),2)), labels = sGenes[seq(2,n,2)], las=2, cex.axis = 0.7)
box()
par(mar=c(10,2,3,1))
image.scale(O,col = newCol(50),horiz = FALSE, ylab = "", xlab="", las = 1, cex.axis =0.8)
mtext(parse(text = "log[10]('Expression + 1')"), side = 2, cex = 0.8, line = 1)
dev.off()

