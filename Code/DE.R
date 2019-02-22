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
GM12878 <- CreateSeuratObject(raw.data = GM12878, project = "GM12878")
GM12878 <- NormalizeData(GM12878)

GM18502 <- Read10X("GM18502/GRCh38/")
GM18502 <- CreateSeuratObject(raw.data = GM18502, project = "GM18502")
GM18502 <- NormalizeData(GM18502)

merged <- MergeSeurat(GM12878, GM18502, add.cell.id1 = "GM12878", add.cell.id2 = "GM18502")
rm(GM12878)
rm(GM18502)
DE <- FindMarkers(merged, ident.1 = "GM18502")
merged <- ScaleData(merged, genes.use = rownames(DE))
DE_values <- merged@scale.data[rownames(DE),]
DE_values <- t(scale(t(DE_values)))
M2 <- merged@scale.data
M2 <- M2[sort(rownames(DE)),]
M2[M2 > 2] <- 2
M2[M2 < -2] <- -2
colors <- c(
  rgb(25/255,25/255,255/255),
  rgb(70/255,70/255,255/255),
  rgb(115/255,115/255,255/255),
  rgb(160/255,160/255,255/255),
  #rgb(205/255,205/255,255/255),
  rgb(1,1,1),
  #rgb(255/255,205/255,205/255),
  rgb(255/255,160/255,160/255),
  rgb(255/255,115/255,115/255),
  rgb(255/255, 70/255,70/255),
  rgb(255/255,25/255,25/255)
)

png("Figure4.png", height = 2.9*300, width = 6*300, res = 300)
nCells <- table(merged@meta.data$orig.ident)
layout(matrix(c(rep(1,60),rep(2,60),rep(3,60),4,4,4,5,5,5,5,5,5,5), ncol = 19, byrow = FALSE))
par(mar=c(0.5,3.5,1.5,3.5), mgp = c(1,0.5,0))

gNames <- sort(rownames(M2)[1:80], decreasing = TRUE)
image(t(M2[gNames,]), col = colors, yaxt="n", xaxt="n", bty = "n")
axis(side = 3, at = c((nCells[1]/2), nCells[1]+(nCells[2]/2))/sum(nCells), labels = c("GM12878", "GM18502"), las=1, cex.axis=0.8)
abline(v = nCells[1]/sum(nCells), col= "white", lwd = 3)
axis(side=2, at= (1/79)*seq(1,80,2), labels = gNames[seq(1,80,2)], las=2, cex.axis = 0.6)
axis(side=4, at= ((1/79)*seq(0,78,2)), labels = gNames[seq(2,80,2)], las=2, cex.axis = 0.6)
box()

gNames <- sort(rownames(M2)[81:160], decreasing = TRUE)
image(t(M2[gNames,]), col = colors, yaxt="n", xaxt="n", bty = "n")
axis(side = 3, at = c((nCells[1]/2), nCells[1]+(nCells[2]/2))/sum(nCells), labels = c("GM12878", "GM18502"), las=1, cex.axis=0.8)
abline(v = nCells[1]/sum(nCells), col= "white", lwd = 3)
axis(side=2, at= (1/79)*seq(1,80,2), labels = gNames[seq(1,80,2)], las=2, cex.axis = 0.6)
axis(side=4, at= ((1/79)*seq(0,78,2)), labels = gNames[seq(2,80,2)], las=2, cex.axis = 0.6)
box()


gNames <- sort(rownames(M2)[161:240], decreasing = TRUE)
image(t(M2[gNames,]), col = colors, yaxt="n", xaxt="n", bty = "n")
axis(side = 3, at = c((nCells[1]/2), nCells[1]+(nCells[2]/2))/sum(nCells), labels = c("GM12878", "GM18502"), las=1, cex.axis=0.8)
abline(v = nCells[1]/sum(nCells), col= "white", lwd = 3)
axis(side=2, at= (1/79)*seq(1,80,2), labels = gNames[seq(1,80,2)], las=2, cex.axis = 0.6)
axis(side=4, at= ((1/79)*seq(0,78,2)), labels = gNames[seq(2,80,2)], las=2, cex.axis = 0.6)
box()

par(mar=c(1,1.5,2.3,0.5), mgp = c(1.5,0.5,0))
image.scale(M2, col = colors, ylab="", xlab = "", horiz = FALSE, las=2, cex.axis=0.7)
box()
dev.off()
