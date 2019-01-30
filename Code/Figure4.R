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
png("Figure4.png", height = 2*300, width = 8*300, res = 300)
layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2), ncol = 21, nrow = 1))
par(mar=c(4,1.5,4,0), mgp = c(1,0.5,0))
image((M2), col = colors, yaxt="n", xaxt="n", bty = "n")
nCells <- table(merged@meta.data$orig.ident)
axis(side = 2, at = c((nCells[1]/2), nCells[1]+(nCells[2]/2))/sum(nCells), labels = c("GM12878", "GM18502"), las=3, cex.axis=0.7)
abline(h = nCells[1]/sum(nCells), col= "white", lwd = 3)
axis(side = 1, at = seq(from=0,to=nrow(DE)-1,by=2)*(1/(nrow(DE))), labels = rownames(DE)[seq(2,nrow(DE),2)], las = 2, cex.axis = 0.65, lwd = 0.3)
axis(side = 3, at = seq(from=1,to=nrow(DE),by=2)*(1/(nrow(DE))), labels = rownames(DE)[seq(1,nrow(DE),2)], las = 2, cex.axis = 0.65, lwd = 0.3)
box()
par(mar=c(4,2,4,0.5), mgp = c(1.5,0.5,0))
image.scale(M2, col = colors, ylab="", xlab = "", horiz = FALSE, las=2, cex.axis=0.7)
box()
dev.off()

