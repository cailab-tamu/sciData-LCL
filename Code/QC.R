library(Seurat)

cellQC <- function(cellLine, outName){
  gE <- Read10X(paste0(cellLine,"/GRCh38/"))
  gE <- CreateSeuratObject(raw.data = gE)
  gE <- NormalizeData(gE)
  gE <- ScaleData(gE, genes.use = unlist(cc.genes))
  gE <- CellCycleScoring(gE, g2m.genes = cc.genes$g2m.genes, s.genes = cc.genes$s.genes)
  cellQC <- data.frame(
    UMI = gE@meta.data$nUMI, 
    geneNumber = gE@meta.data$nGene,
    mtProportion = round(colSums(gE@raw.data[grepl("^MT-", rownames(gE@raw.data)),])/colSums(gE@raw.data),3),
    cellPhase = gE@meta.data$Phase
  )
  write.table(cellQC, quote = FALSE, sep = "\t", file = paste0(outName, "_cellQC.tsv"))
}

cellQC("GM12878", "GM12878")
cellQC("GM18502", "GM18502")
cellQC("mix", "MIXTURE")