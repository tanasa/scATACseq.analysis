###################################################################################
###################################################################################
library(ArchR)
library(parallel)
library(hexbin)
library(presto)
library(harmony)

set.seed(1)
addArchRThreads(threads = 1)
addArchRGenome("hg38")  ### it is a masked version of hg38

###################################################################################
###################################################################################
# to change the FOLDER to : "/media/bogdan/easystore/ATAC_seq_scATAC_seq/DATASET_BMMC"
###################################################################################
###################################################################################

pathFragments="./"
# pathFragments="/media/bogdan/easystore/ATAC_seq_scATACse/"
inputFiles <- list.files(pathFragments, pattern = ".gz", full.names = TRUE)
# inputFiles <- "atac_v1_pbmc_10k_fragments.tsv.gz"
names(inputFiles) <- gsub("_fragments.tsv.gz", "", list.files(pathFragments, pattern = ".gz"))
inputFiles <- inputFiles[!grepl(".tbi", inputFiles)]

inputFiles

# inputFiles <- getTutorialData("Hematopoiesis")
# inputFiles="scATAC_BMMC_R1.fragments.tsv.gz"

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 1, # Dont set this too high because you can always increase later
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles

# "scATAC_BMMC_R1.fragments.tsv.gz.arrow"

dim(proj@cellColData$TSSEnrichment) ### 4985 cells

################################################################################### DOUBLETS
###################################################################################

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

# > doubScores
# [[1]]
# List of length 2
# names(2): doubletScore doubletEnrich

################################################################################### PROJECT
###################################################################################
###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
################################################################################### the STRUCTURE of the object PROJECT

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

getAvailableMatrices(proj)
getAvailableMatrices(proj)

#  “GeneScoreMatrix” “TileMatrix”

# class: ArchRProject
# outputDirectory: /media/bogdan/easystore/ATAC_seq_scATAC_seq/DATASET_BMMC/HemeTutorial
# samples(1): scATAC_BMMC_R1.fragments.tsv.gz
# sampleColData names(1): ArrowFiles
# cellColData names(15): Sample TSSEnrichment ... DoubletEnrichment
# BlacklistRatio
# numberOfCells(1): 4932
# medianTSS(1): 15.2575
# medianFrags(1): 2771

 # proj@projectMetadata
 # proj@projectSummary  
 # proj@sampleColData
 # proj@sampleMetadata
 # proj@cellColData  #### that is very useful !!!!
# colnames(proj@cellColData)
# [1] "Sample"            "TSSEnrichment"     "ReadsInTSS"      
# [4] "ReadsInPromoter"   "ReadsInBlacklist"  "PromoterRatio"    
# [7] "PassQC"            "NucleosomeRatio"   "nMultiFrags"      
#  [10] "nMonoFrags"        "nFrags"            "nDiFrags"        
# [13] "DoubletScore"      "DoubletEnrichment" "BlacklistRatio"
# proj@cellColData
 
#   proj@cellColData$Sample
#   proj@cellColData$ReadsInPromoter
#   proj@cellColData$PassQC
#   proj@cellColData$nMonoFrags
#   proj@cellColData$DoubletScore
#   proj@cellColData$TSSEnrichment
#   proj@cellColData$ReadsInBlacklist
#   proj@cellColData$NucleosomeRatio
#   proj@cellColData$nFrags
#   proj@cellColData$DoubletEnrichment
#   proj@cellColData$ReadsInTSS
#   proj@cellColData$PromoterRatio
#   proj@cellColData$nMultiFrags
#   proj@cellColData$nDiFrags
#   proj@cellColData$BlacklistRatio
     
rownames(proj@cellColData)
dim(rownames(proj@cellColData))
# length(rownames(proj@cellColData))
# [1] 4932

write.table("THE_NUMBER_of_CELLS_before_filtering", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t")
           
write.table(dim(proj@cellColData$TSSEnrichment), file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t")
           
           
# proj@cellMetadata
# proj@reducedDims
# proj@embeddings  
# proj@peakSet
# proj@peakAnnotation
# proj@geneAnnotation
# proj@genomeAnnotation
# proj@imputeWeights
# proj@genomeAnnotation
# proj@imputeWeights  

# IT MAY MAKE more SENSE to use the QUANTILE ????????

png(filename = "HISTOGRAM.cellColData$ReadsInPromoter.png", width = 480, height = 480, units = "px", pointsize = 12)
hist(proj@cellColData$ReadsInPromoter)
dev.off()

median(proj@cellColData$ReadsInPromoter)
quantile(proj@cellColData$ReadsInPromoter)

png(filename = "HISTOGRAM.cellColData$PassQC.png", width = 480, height = 480, units = "px", pointsize = 12)
hist(proj@cellColData$PassQC)
dev.off()

median(proj@cellColData$PassQC)

png(filename = "HISTOGRAM.proj@cellColData$nMonoFrags.png", width = 480, height = 480, units = "px", pointsize = 12)
hist(proj@cellColData$nMonoFrags)
dev.off()

median(proj@cellColData$nMonoFrags)

png(filename = "HISTOGRAM.cellColData$DoubletScore.png", width = 480, height = 480, units = "px", pointsize = 12)
hist(proj@cellColData$DoubletScore)
dev.off()

median(proj@cellColData$DoubletScore)

png(filename = "HISTOGRAM.cellColData$TSSEnrichment.png", width = 480, height = 480, units = "px", pointsize = 12)
hist(proj@cellColData$TSSEnrichment)
dev.off()

median(proj@cellColData$TSSEnrichment)

png(filename = "HISTOGRAM.cellColData$ReadsInBlacklist.png", width = 480, height = 480, units = "px", pointsize = 12)
hist(proj@cellColData$ReadsInBlacklist)
dev.off()

median(proj@cellColData$ReadsInBlacklist)

png(filename = "HISTOGRAM.cellColData$NucleosomeRatio.png", width = 480, height = 480, units = "px", pointsize = 12)
hist(proj@cellColData$NucleosomeRatio)
dev.off()

median(proj@cellColData$NucleosomeRatio)

png(filename = "HISTOGRAM.cellColData.nFrags.png", width = 480, height = 480, units = "px", pointsize = 12)
hist(proj@cellColData$nFrags)
dev.off()

median(proj@cellColData$nFrags)

png(filename = "HISTOGRAM.cellColData.DoubletEnrichment.png", width = 480, height = 480, units = "px", pointsize = 12)
hist(proj@cellColData$DoubletEnrichment)
dev.off()

median(proj@cellColData$DoubletEnrichment)

png(filename = "HISTOGRAM.cellColData.ReadsInTSS.png", width = 480, height = 480, units = "px", pointsize = 12)
hist(proj@cellColData$ReadsInTSS)
dev.off()

median(proj@cellColData$ReadsInTSS)

png(filename = "HISTOGRAM.cellColData.PromoterRatio.png", width = 480, height = 480, units = "px", pointsize = 12)
hist(proj@cellColData$PromoterRatio)
dev.off()

median(proj@cellColData$PromoterRatio)

png(filename = "HISTOGRAM.cellColData.nMultiFrags.png", width = 480, height = 480, units = "px", pointsize = 12)
hist(proj@cellColData$nMultiFrags)
dev.off()

median(proj@cellColData$nMultiFrags)

png(filename = "HISTOGRAM.cellColData.nDiFrags.png", width = 480, height = 480, units = "px", pointsize = 12)
hist(proj@cellColData$nDiFrags)
dev.off()

median(proj@cellColData$nDiFrags)

png(filename = "HISTOGRAM.cellColData.BlacklistRatio.png", width = 480, height = 480, units = "px", pointsize = 12)
hist(proj@cellColData$BlacklistRatio)
dev.off()

median(proj@cellColData$BlacklistRatio)

###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
###################################################################################
median(proj@cellColData$ReadsInPromoter)
median(proj@cellColData$PassQC)
median(proj@cellColData$nMonoFrags)
median(proj@cellColData$DoubletScore)
median(proj@cellColData$TSSEnrichment)
median(proj@cellColData$ReadsInBlacklist)
median(proj@cellColData$NucleosomeRatio)
median(proj@cellColData$nFrags)
median(proj@cellColData$DoubletEnrichment)
median(proj@cellColData$ReadsInTSS)
median(proj@cellColData$PromoterRatio)
median(proj@cellColData$nMultiFrags)
median(proj@cellColData$nDiFrags)
median(proj@cellColData$BlacklistRatio)

# median(proj@cellColData$)
# median(proj@cellColData$)

# DoubletScore      
# DoubletEnrichment

# write.table("", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
#            append = TRUE, quote = TRUE, sep = "\t")
           
# write.table(, file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
#            append = TRUE, quote = TRUE, sep = "\t")

write.table("median(proj@cellColData$ReadsInPromoter)", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table(median(proj@cellColData$ReadsInPromoter), file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table("median(proj@cellColData$PassQC)", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table(median(proj@cellColData$PassQC), file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table("median(proj@cellColData$nMonoFrags)", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)

write.table(median(proj@cellColData$nMonoFrags), file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)

write.table("median("proj@cellColData$DoubletScore", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table(median("proj@cellColData$DoubletScore", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)

write.table("median(proj@cellColData$ReadsInBlacklist)", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table(median(proj@cellColData$ReadsInBlacklist), file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table("median(proj@cellColData$NucleosomeRatio)", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table(median(proj@cellColData$NucleosomeRatio), file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table("median(proj@cellColData$nFrags)", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table(median(proj@cellColData$nFrags), file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table("median(proj@cellColData$DoubletEnrichment)", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table(median(proj@cellColData$DoubletEnrichment), file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table("median(proj@cellColData$ReadsInTSS)", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table(median(proj@cellColData$ReadsInTSS), file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)            
                       
write.table("median(proj@cellColData$PromoterRatio)", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table(median(proj@cellColData$PromoterRatio), file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)                                    
           

write.table("median(proj@cellColData$nMultiFrags)", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table(median(proj@cellColData$nMultiFrags), file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)            
                       
write.table("median(proj@cellColData$nDiFrags)", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table(median(proj@cellColData$nDiFrags), file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)                                    
           
           
write.table("median(proj@cellColData$BlacklistRatio)", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
write.table(median(proj@cellColData$BlacklistRatio), file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)  
           
# write.table("", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
#            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
# write.table(, file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
#            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)  

# write.table("", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
#            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
# write.table(, file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
#            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)  
           
# write.table("", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
#            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)
           
# write.table(, file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
#            append = TRUE, quote = TRUE, sep = "\t", row.names = FALSE)          

###################################################################################
######################################################################################################################################################################
###################################################################################

# shall we have 2 samples in the DATA :
# > str(proj@cellColData$Sample)
# Formal class 'Rle' [package "S4Vectors"] with 4 slots
#  ..@ values         : chr [1:2] "atac_v1_pbmc_10k" "scATAC_BMMC_R1.fragments.tsv.gz"

# > str(proj@cellColData$Sample@values)
# chr [1:2] "atac_v1_pbmc_10k" "scATAC_BMMC_R1.fragments.tsv.gz
 
###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
###################################################################################

proj <- filterDoublets(ArchRProj = proj)

# At this moment, we are left with 4689 CELLS .
# length(proj@cellColData$DoubletEnrichment)
# [1] 4985 cells

write.table("THE_NUMBER_of_CELLS_after_filtering_DOUBLETS", file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t")
           
write.table(dim(proj@cellColData$TSSEnrichment), file = "THE.SUMMARY.NUMBER.OF.CELLS.txt",
            append = TRUE, quote = TRUE, sep = "\t")

######################################################################################################################################################################
######################################################################################################################################################################
###################################################################################
################################################################################### displaying THE TSS ENRICHMENT vs FRAGMENTS
###################################################################################
###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################

df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
df

## TO ADJUST the CUTOFFs shall these CUTOFFS be different than the CUTOFFS that we have used initially, for FILTERING :

p <- ggPoint(
    x = df[,1],
    y = df[,2],
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

p

plotPDF(p, name = "PICTURE.after.FILTERING.DOUBLETS.TSS-vs-Frags.pdf", ArchRProj = proj, addDOC = FALSE)


proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")

table(proj$Clusters)


# CLUSTERING can be done by using two methods :

# method = "Seurat",
# resolution = 0.8

# method = "scran",
# k = 15

# in order to VISUALIZE the DATA that has been ADDED :

# head(proj@reducedDims$IterativeLSI$matSVD)
# head(proj@reducedDims$IterativeLSI)

# proj@cellMetadata
# proj@reducedDims
# proj@embeddings  
# proj@peakSet
# proj@peakAnnotation
# proj@geneAnnotation
# proj@genomeAnnotation
# proj@imputeWeights
# proj@genomeAnnotation
# proj@imputeWeights  

# the CLUSTERS are being added to the each CELL :

#                                                      Clusters
#                                                   <character>
#scATAC_BMMC_R1.fragments.tsv.gz#TTATGTCAGTGATTAG-1          C3
#scATAC_BMMC_R1.fragments.tsv.gz#GCATTGAAGATTCCGT-1          C7
#scATAC_BMMC_R1.fragments.tsv.gz#TATGTTCAGGGTTCCC-1          C7

###################################################################################
################################################################################### UMAP visualization
###################################################################################
################################################################################### VISUALIZATION

# proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

# there is another element that is introduced in the DATA STRUCTURE :

# @ embeddings      :Formal class 'SimpleList' [package "S4Vectors"] with 4 slots
#  .. .. ..@ listData       :List of 1
#  .. .. .. ..$ UMAP:
 
# proj@embeddings
 
# in order to understand the structure of this OBJECT :
 
# head(proj@embeddings@listData$UMAP$df)
                                                   IterativeLSI#UMAP_Dimension_1
# scATAC_BMMC_R1.fragments.tsv.gz#TTATGTCAGTGATTAG-1                     0.6728396
# scATAC_BMMC_R1.fragments.tsv.gz#GCATTGAAGATTCCGT-1                     2.1138837
# scATAC_BMMC_R1.fragments.tsv.gz#TATGTTCAGGGTTCCC-1                     0.6712529
# scATAC_BMMC_R1.fragments.tsv.gz#AGTTACGAGAACGTCG-1                     1.6911861
# scATAC_BMMC_R1.fragments.tsv.gz#TTACTCAGTTCGGGAA-1                     0.2246994
# scATAC_BMMC_R1.fragments.tsv.gz#GCACCTTAGACTAGCG-1                     2.2809383
                                                   IterativeLSI#UMAP_Dimension_2
# scATAC_BMMC_R1.fragments.tsv.gz#TTATGTCAGTGATTAG-1                   -0.08073505
# scATAC_BMMC_R1.fragments.tsv.gz#GCATTGAAGATTCCGT-1                   -5.35591357
# scATAC_BMMC_R1.fragments.tsv.gz#TATGTTCAGGGTTCCC-1                   -3.29497879
# scATAC_BMMC_R1.fragments.tsv.gz#AGTTACGAGAACGTCG-1                    1.46791155
# scATAC_BMMC_R1.fragments.tsv.gz#TTACTCAGTTCGGGAA-1                   -3.28392308
# scATAC_BMMC_R1.fragments.tsv.gz#GCACCTTAGACTAGCG-1                   -5.60076659

###################################################################################
###################################################################################
###################################################################################
###################################################################################

proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

# nNeighbors = 30,
# minDist = 0.5,
# metric = "cosine"

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP",
                    size = 2, baseSize = 10)
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP",
                   size = 2, baseSize = 10)

ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2, name = "FIGURE.after.FILTERING.Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

###################################################################################
################################################################################### TSNE visualization
###################################################################################
################################################################################### VISUALIZATION

proj <- addTSNE(ArchRProj = proj, reducedDims = "IterativeLSI")

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "TSNE",
                    size = 2, baseSize = 10)
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "TSNE",
                   size = 2, baseSize = 10)

ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2, name = "FIGURE.after.FILTERING.Plot-TSNE-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

###################################################################################
###################################################################################
###################################################################################
################################################################################### to impute the WEIGHTS :

proj <- addImputeWeights(proj)

unique(proj@cellColData$Clusters)

# str(proj@imputeWeights@listData)
# str(proj@imputeWeights@listData$Weights)
# str(proj@imputeWeights@listData$Params)

###################################################################################
###################################################################################
###################################################################################
################################################################################### before COMPUTING the CLUSTERS :

# another way to get the clusters is :
x = unique(proj@cellColData$Clusters)
# [1] "C7" "C8" "C4" "C1" "C2" "C6" "C3" "C5"

markersGS <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",                     ### to STRATIFY CELL GROUPS
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

# each CLUSTER has a LIST of GENES :
#> head(MDF)
#  group group_name seqnames     start       end strand      name  idx   Log2FC
#1     1         C1     chr2  45169037  45173216      1      SIX3  277 1.476506
#2     1         C1     chr1 246939315 246955685      1 LINC01341 2252 1.630452
#3     2         C2     chr3  12328984  12475855      1     PPARG   73 2.546852
#4     2         C2    chr18   9708228   9862553      1     RAB31   53 2.375811


M = markerList %>% as_tibble()
MDF = as.data.frame(M)
MDFL = dim(MDF)[1]

for (i in 1:MDFL) {
     
     CLUSTER = MDF$group_name[i] ;
     GENE = MDF$name[i] ;
     p = plotEmbedding(
         ArchRProj = proj,
         colorBy = "GeneScoreMatrix",
         name = GENE,
         embedding = "UMAP",
         imputeWeights = getImputeWeights(proj),
         size = 2, baseSize = 10)
       
       do.call(cowplot::plot_grid, c(list(ncol = 3), p))
       plotPDF(p, name = paste(GENE, "Plot.MARKER_genes", "cluster", CLUSTER, "gene", GENE, "pdf", sep="."),
                  ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

}

### IT WAS from the OLD TUTORIAL !!

#p <- plotEmbedding(
#    ArchRProj = proj,
#    colorBy = "GeneScoreMatrix",
#    name = markerGenes,
#    embedding = "UMAP",
#    imputeWeights = getImputeWeights(proj),
#    size = 2, baseSize = 10
#)


#do.call(cowplot::plot_grid, c(list(ncol = 3), p))
#plotPDF(p, name = "Plot.UMAP_set_genes.pdf",
#        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
       
#plotPDF(plotList = p,
#    name = "Plot.UMAP-Marker-Genes-W-Imputation.pdf",
#    ArchRProj = proj,
#    addDOC = FALSE, width = 5, height = 5)


# REARRANGE for grid plotting

# p2 <- lapply(p, function(x){
#    x + guides(color = FALSE, fill = FALSE) +
#    theme_ArchR(baseSize = 6.5) +
#    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
#    theme(
#        axis.text.x=element_blank(),
#        axis.ticks.x=element_blank(),
#        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank()
#    )
# })

#do.call(cowplot::plot_grid, c(list(ncol = 3), p2))
#plotPDF(p2, name = "Plot.UMAP_set_genes.pdf",
#        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

###################################################################################
###################################################################################
###################################################################################
################################################################################### in order to display only 1 GENE

# p$CD4
# plotPDF(p$CD4, name = "UMAP_set_genes.CD14.pdf",
#         ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# grid::grid.newpage()
# grid::grid.draw(p$CD4)

###################################################################################
###################################################################################
################################################################################### displayin the GENOME BROWSER
###################################################################################

M = markerList %>% as_tibble()
MDF = as.data.frame(M)
MDFL = dim(MDF)[1]

for (i in 1:MDFL) {
     
     CLUSTER = MDF$group_name[i] ;
     GENE = MDF$name[i] ;
     
     p <- plotBrowserTrack(
          ArchRProj = proj,
          groupBy = "Clusters",
          geneSymbol = GENE,
          upstream = 50000,
          downstream = 50000,
          plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
          sizes = c(10, 1.5, 3, 4)
      )

      do.call(cowplot::plot_grid, c(list(ncol = 3), p))
      plotPDF(plotList = p,
        name = paste(GENE, "Plot.MARKER_genes.in.BROWSER", "cluster", CLUSTER, "gene", GENE, "pdf", sep="."),
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)
}  
   
###################################################################################
###################################################################################

# Sys.Date()
# ArchRBrowser(ArchRProj = proj)

proj <- saveArchRProject(ArchRProj = proj)

###################################################################################
###################################################################################
################################################################################### 
################################################################################### 
###################################################################################
###################################################################################

# PROJ :

# class: ArchRProject
# outputDirectory: /media/bogdan/easystore/ATAC_seq_scATAC_seq/DATASET_PBMC/HemeTutorial
# samples(1): atac_v1_pbmc_10k
# sampleColData names(1): ArrowFiles
# cellColData names(16): Sample TSSEnrichment ... BlacklistRatio Clusters
# numberOfCells(1): 8383
# medianTSS(1): 17.208
# medianFrags(1): 9432

###################################################################################
###################################################################################
################################################################################### SAVING the ARCh PROJECT
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
################################################################################### creating ARCH PROJECT

### we have worked with a project that is called : "proj"

proj <- saveArchRProject(ArchRProj = proj)

#############################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
################################################################################### THE STRUCTURE of the ARCh PROJECT

#proj@projectMetadata
#proj@projectSummary
#proj@sampleColData
#proj@sampleMetadata
#proj@cellColData             ########## the most important ACCESSOR and MATRIX.
#proj@cellMetadata
#proj@reducedDims
#proj@embeddings
#proj@peakSet
#proj@peakAnnotation
#proj@geneAnnotation
#proj@genomeAnnotation
#proj@imputeWeights

# colnames(proj@cellColData)
# [1] "Sample"            "TSSEnrichment"     "ReadsInTSS"      
# [4] "ReadsInPromoter"   "ReadsInBlacklist"  "PromoterRatio"    
# [7] "PassQC"            "NucleosomeRatio"   "nMultiFrags"      
# [10] "nMonoFrags"        "nFrags"            "nDiFrags"        
# [13] "DoubletScore"      "DoubletEnrichment" "BlacklistRatio"  
# [16] "Clusters"  

################ THE MEMORY that we have USED :
# > paste0("Memory Size = ", round(object.size(proj) / 10^6, 3), " MB")
# [1] "Memory Size = 50.673 MB"

################ AVAILABLE MATRICES
# [1] "GeneScoreMatrix" "TileMatrix

projHeme2 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme2,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "BioClassification",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

projHeme3 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme2,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = TRUE,
    force= TRUE,
    # groupList = groupList,
    groupRNA = "BioClassification",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)

###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
########################################################################## USING the ACCESSORS : $$$$

head(proj$Sample)
head(proj$cellNames)
head(proj$DoubletScore)
head(proj$TSSEnrichment)

quantile(proj$TSSEnrichment)

################################################################################### SUBSETTING by CELLS : the first cells
########################################

proj[1:100, ]
proj[proj$cellNames[1:100]]

###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
###################################################################################
###################################################################################  SUBSETTING by CELLS : function of TSS ENRICHMENT
######################################## SUBSETTING by CELLS :

# starting with 8383 CELLS
#  
idxPass <- which(proj$TSSEnrichment >= 8)
cellsPass <- proj$cellNames[idxPass]
proj[cellsPass, ]

###################################################################################
###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################

######################################################################################################################################################################
################################################################################### OBTAINING COLUMNS and PERFORMING OPERATIONS :
######################################################################################################################################################################

df <- getCellColData(proj, select = "nFrags")
df

df <- getCellColData(proj, select = c("log10(nFrags)", "nFrags - 1"))
df

######################################################################################################################################################################
############################################################################################################################ PLOTTING QC metrics :

df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
df

## TO ADJUST the CUTOFFs shall these CUTOFFS be different than the CUTOFFS that we have used initially, for FILTERING :

p <- ggPoint(
    x = df[,1],
    y = df[,2],
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

p

plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = proj, addDOC = FALSE)

# ggPoint    
# This function is a wrapper around ggplot geom_point to allow for a
# more intuitive plotting of ArchR data.

###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
###################################################################################
###  Plotting Sample Statistics from an ArchRProject

### make a RIDGE PLOT :

###################################################################################

p1 <- plotGroups(
    ArchRProj = proj,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "ridges"
)

plotPDF(p1, name = "PICTURE.TSSenrichment.pdf", ArchRProj = proj, addDOC = FALSE)

p3 <- plotGroups(
    ArchRProj = proj,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "ridges"
)

plotPDF(p3, name = "PICTURE.log10_FRAGS.pdf", ArchRProj = proj, addDOC = FALSE)
   
### make a VIOLIN PLOT :

p2 <- plotGroups(
    ArchRProj = proj,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

p4 <- plotGroups(
    ArchRProj = proj,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj, addDOC = FALSE, width = 4, height = 4)

###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
###################################################################################

## PLOTTING SAMPLE FRAGMENT DISTRIBUTION :
## Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.

p1 <- plotFragmentSizes(ArchRProj = proj)
p1

p2 <- plotTSSEnrichment(ArchRProj = proj)
p2

plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
################################################################################### SAVING an ARchR project :

# save the current Arrow files to the designated outputDirectory so that they are exclusively associated with the new ArchRProject object.

saveArchRProject(ArchRProj = proj, outputDirectory = "Save-Proj", load = FALSE)


# This process does NOT automatically update the ArchRProject object that is active in your current R session.
# Specifically, the object named projHeme1 in the current R session will still point to the original location of the Arrow files,
# not the copied Arrow files that reside in the specified outputDirectory

###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
###################################################################################

# Filtering Doublets from an ArchRProject

## COMMAND : addDoubletScores()

## COMMAND : filterDoublets()  # in the

## filterRatio which is the maximum ratio of predicted doublets to filter based on the number of pass-filter cells.

# If you wanted to filter more cells from the ArchR Project, you would use a higher filterRatio.

## if we wanna to increase the STRINGENCY of FILTERING :
## projHemeTmp <- filterDoublets(proj, filterRatio = 1.5)

proj <- filterDoublets(proj, filterRatio = 1)

###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
##########################################; ; #########################################
####  Dimensionality Reduction with ArchR : ArchR’s LSI Implementation

#### To perform iterative LSI in ArchR, we use the addIterativeLSI() function.

#### addIterativeLSI()

### The most common parameters to tweak are iterations, varFeatures, and resolution.
### It is important to note that LSI is not deterministic.
### This means that even if you run LSI in exactly the same way with exactly the same parameters,
### you will not get exactly the same results.

projHeme2 = proj

projHeme2 <- addIterativeLSI(
    ArchRProj = projHeme2,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 2,
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2),
        sampleCells = 10000,
        n.start = 10
    ),
    varFeatures = 25000,
    dimsToUse = 1:30,
    force = TRUE
)

projHeme2@geneAnnotation
projHeme2@peakAnnotation
projHeme2@embeddings
projHeme2@genomeAnnotation
projHeme2@imputeWeights

# Estimated LSI is accessed in ArchR via the addIterativeLSI()
# function by setting the sampleCellsFinal and projectCellsPre parameters

# BATCH EFFECT CORRECTION : SHALL we have more than 2 SAMPLES :

projHeme2 <- addHarmony(
    ArchRProj = projHeme2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
################################################################################### CLUSTERING :

# Clustering using Seurat’s FindClusters() function

projHeme2 <- addClusters(
    input = projHeme2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force=TRUE
)

table(proj$Clusters)

# A confusion matrix across EACH SAMPLE :

cM <- confusionMatrix(paste0(proj$Clusters), paste0(proj$Sample))
cM

# > cM
#16 x 1 sparse Matrix of class "dgCMatrix"
#    atac_v1_pbmc_10k
#C12              123
#C4              2888
#C15              814
#C10              824
#C11              909
#C3               235

# A HEATMAP that represents each SAMPLE and CLUSTERS : it must have more than 2 CLUSTERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM),
    color = paletteContinuous("whiteBlue"),
    border_color = "black"
)
p

################################################################################### Single-cell EMBEDDINGS :
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
################################################################################### EMBEDDINGS

# We call these “embeddings” because they are strictly used to visualize the clusters and are not used to identify clusters which is done in
# an LSI sub-space as mentioned in previous chapters. The primary difference between UMAP and t-SNE is the interpretatino of the distance between cells or clusters.
# t-SNE is designed to preserve the local structure in the data while UMAP is designed to preserve both the local and most of the global structure in the data.

# UMAP :

# addUMAP
# plotEmbedding
# ggAlignPlots
# ggAlignPlots(p1, p2, type = "h")
# projHeme2@geneAnnotation
# projHeme2@peakAnnotation
# projHeme2@embeddings
# projHeme2@genomeAnnotation
# projHeme2@imputeWeights

projHeme2 <- addUMAP(
    ArchRProj = projHeme2,
    reducedDims = "IterativeLSI",
    name = "UMAP",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine", force = TRUE
)

p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

ggAlignPlots(p1, p2, type = "h")
 
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
 
# TSNE :
# addTSNE()
# plotEmbedding()
# ggAlignPlots
 
projHeme2 <- addTSNE(
    ArchRProj = projHeme2,
    reducedDims = "IterativeLSI",
    name = "TSNE",
    perplexity = 30
)

p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2, name = "Plot-TSNE-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

###################################################################################
###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
###################################################################################
################################################################################### HARMONY ::

# addHarmony() function, creating a reducedDims object named “Harmony”
# Dimensionality Reduction After Harmony

# addHarmony() function, creating a reducedDims object named “Harmony”.

projHeme2 <- addUMAP(
    ArchRProj = projHeme2,
    reducedDims = "Harmony",
    name = "UMAPHarmony",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine"
)

p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")

ggAlignPlots(p3, p4, type = "h")

plotPDF(p1,p2,p3,p4, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

projHeme2 <- addTSNE(
    ArchRProj = projHeme2,
    reducedDims = "Harmony",
    name = "TSNEHarmony",
    perplexity = 30
)

p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNEHarmony")
p4 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "TSNEHarmony")

ggAlignPlots(p3, p4, type = "h")

plotPDF(p1,p2,p3,p4, name = "Plot-TSNE2Harmony-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
###################################################################################  Gene Scores and Marker Genes with ArchR

# GENE SCORES

# Accessibility within the entire gene body contributes to the gene score.

# An exponential weighting function that accounts for the activity of putative distal regulatory elements
# in a distance-dependent fashion.

# Imposed gene boundaries that minimizes the contribution of unrelated regulatory elements to the gene score.

# addGeneScoreMat is set to TRUE - this is the default behavior.
# Alternatively, gene scores can be added to Arrow files at any time by using the addGeneScoreMatrix().

# MARKER FEATURES

# These features features can be anything -
# peaks,
# genes (based on gene scores),
# or transcription factor motifs (based on chromVAR deviations).

# getMarkerFeatures()

# useMatrix = "TileMatrix" would identify genomic regions that are highly specific to a certain cell group
# useMatrix = "PeakMatrix" would identify peaks that are highly specific to a certain cell group.

library(presto)
library(harmony)

################################################# MARKER GENES :

markersGS <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",                     ### to STRATIFY CELL GROUPS
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

# THE STRUCTURE of the object "markerGS"

# A SummarizedExperiment object contains one or more assays, each represented by a matrix-like object of numeric data,
# and metadata that applies to the rows or columns of the assays matrices.

markersGS@colData

# str(markersGS@colData)
# Formal class 'DFrame' [package "S4Vectors"] with 6 slots
#  ..@ rownames       : chr [1:16] "C1" "C2" "C3" "C4" ...

# markersGS@assays
# markersGS@elementMetadata
# markersGS@metadata

#markersGS@metadata$Params$groupBy
#[1] "Clusters"

#markersGS@metadata$Params$useGroups
#NULL

#markersGS@metadata$Params$bgdGroups
#NULL

#markersGS@metadata$Params$useMatrix
#[1] "GeneScoreMatrix"

#markersGS@metadata$Params$bias
#[1] "TSSEnrichment" "log10(nFrags)"

#markersGS@metadata$Params$normBy
#NULL

#markersGS@metadata$Params$testMethod
#[1] "wilcoxon"

#markersGS@metadata$Params$maxCells
#[1] 500

#markersGS@metadata$Params$scaleTo
#[1] 10000

# the CUTOFF by default is : cutOff = "FDR <= 0.1 & Log2FC >= 0.5

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

markerList$C1$name
markerList$C2$name
markerList$C3$name
markerList$C4$name
markerList$C5$name
markerList$C6$name
markerList$C7$name
markerList$C8$name
markerList$C9$name
markerList$C10$name
markerList$C11$name
markerList$C12$name
markerList$C13$name
markerList$C14$name

markerGenes  <- c(
    "CD34",  # Early Progenitor
    "GATA1", # Erythroid
    "PAX5", "MS4A1", "EBF1", "MME", # B-Cell Trajectory
    "CD14", "CEBPB", "MPO", # Monocytes
    "IRF8",
    "CD3D", "CD8A", "TBX21", "IL7R" # TCells
  )

library("magick")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE)


######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
################################################################################### VISUALIZE MARKER GENES on the EMBEDDINGS
################################################################################### TO CHANGE the SET of MARKER GENES that we do have 

# markerGenes  <- c(
#    "CD34",  #Early Progenitor
#    "GATA1", #Erythroid
#    "PAX5", "MS4A1", "MME", #B-Cell Trajectory
#    "CD14", "MPO", #Monocytes
#    "CD3D", "CD8A"#TCells
#  )

### SLECTING a SET of GENES that we may extract from the DATASET
 
markerGenes = c(markerList$C14$name[1], markerList$C13$name[1],
                markerList$C15$name[1], markerList$C16$name[1])  
               
# markerGenes
# [1] "RORC"   "TARP"   "NFATC2" "CCR7"

p <- plotEmbedding(
    ArchRProj = projHeme2,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)

p$CD14
p$TARP
p$NFATC2
p$CCR
p$RORC

p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

plotPDF(plotList = p,
    name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf",
    ArchRProj = projHeme2,
    addDOC = FALSE, width = 5, height = 5)


################################################################################### 
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
################################################################################### Marker Genes Imputation with MAGIC

# The sparsity of scATAC-seq data.
# We can use MAGIC to impute gene scores by smoothing signal across nearby cells.
# It greatly improves the visual interpretation of gene scores.

projHeme2 <- addImputeWeights(projHeme2)

# markerGenes  <- c(
#    "CD34",  #Early Progenitor
#    "GATA1", #Erythroid
#    "PAX5", "MS4A1", "MME", #B-Cell Trajectory
#    "CD14", "MPO", #Monocytes
#    "CD3D", "CD8A"#TCells
#  )

markerGenes = c(markerList$C14$name[1], markerList$C13$name[1],
                markerList$C15$name[1], markerList$C16$name[1])  
               

p <- plotEmbedding(
    ArchRProj = projHeme2,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme2)
)

p$CD14
p$TARP
p$NFATC2
p$CCR
p$RORC

### TO PLOT ALL MARKER GENES : #### Rearrange for grid PLOTTING
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

plotPDF(plotList = p,
    name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf",
    ArchRProj = projHeme2,
    addDOC = FALSE, width = 5, height = 5)

################################################################################### there is an ARCHBROWSER :
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
################################################################################### PLOTTING with ARCHBROWSER : plotBrowserTrack()
#### if needed, we do RELOAD the PROJECT :
# loadArchRProject(path = "./", force = FALSE, showLogo = TRUE)
# readRDS("/media/bogdan/easystore/ATAC_seq_scATAC_seq/DATASET_PBMC/HemeTutorial/Save-ArchR-Project.rds", refhook = NULL)


markerGenes  <- c(
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", #B-Cell Trajectory
    "CD14", #Monocytes
    "CD3D", "CD8A", "TBX21", "IL7R" #TCells
)
 
markerGenes = c(markerList$C14$name[1], markerList$C13$name[1],
                markerList$C15$name[1], markerList$C16$name[1])  
                 
p <- plotBrowserTrack(
    ArchRProj = projHeme2,
    groupBy = "Clusters",
    geneSymbol = markerGenes,
    upstream = 50000,
    downstream = 50000
)

### In the configuration of plotBrowserTrack : we can setup the UPSTREAM and DOWNSTREAM
# plotBrowserTrack(
#       ArchRProj = NULL,
#       region = NULL,
#       groupBy = "Clusters",
#       useGroups = NULL,
#       plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
#       sizes = c(10, 1.5, 3, 4),
#       features = getPeakSet(ArchRProj),
#       loops = getCoAccessibility(ArchRProj),
#       geneSymbol = NULL,
#       useMatrix = NULL,
#       log2Norm = TRUE,
#       upstream = 50000,
#       downstream = 50000,
#       tileSize = 250,
#       minCells = 25)
       

grid::grid.newpage()
grid::grid.draw(p$markerList$C14$name[1])
grid::grid.draw(p$markerList$C13$name[1])
grid::grid.draw(p$markerList$C15$name[1])
grid::grid.draw(p$markerList$C16$name[1])

plotPDF(plotList = p,
    name = "Plot-Tracks-Marker-Genes.pdf",
    ArchRProj = projHeme2,
    addDOC = FALSE, width = 5, height = 5)

### in order to startArchR Browser :

ArchRBrowser()

###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
################################################################################### Defining Cluster Identity with scRNA-seq

# The integration works is by directly aligning cells from scATAC-seq with cells from scRNA-seq
# by comparing the scATAC-seq gene score matrix with the scRNA-seq gene expression matrix

########################################################################## TO COMPARE : comparing

# the scATAC-seq gene score matrix with
# the scRNA-seq gene expression matrix

#  FindTransferAnchors()

################################################################################### CROSS-PLATFORM LINKAGE of scATAC-seq cells with scRNA-seq cells
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
################################################################################### to LOAD the existing saved PROJECT :

library(ArchR)
library(parallel)
library(hexbin)

set.seed(1)
addArchRThreads(threads = 1)
addArchRGenome("hg19")  ### it is a masked version of HG19 !!!

proj = loadArchRProject(path = "./", force = FALSE, showLogo = TRUE)

# readRDS("/media/bogdan/easystore/ATAC_seq_scATAC_seq/DATASET_PBMC/HemeTutorial/Save-ArchR-Project.rds", refhook = NULL)

readRDS("/media/bogdan/easystore/ATAC_seq_scATAC_seq_part4/Save-ProjHeme3/Save-ArchR-Project.rds", refhook = NULL)

seRNA <- readRDS("/media/bogdan/easystore/ATAC_seq_scATAC_seq/DATASET_PBMC/HemeTutorial/scRNA-Hematopoiesis-Granja-2019.rds")
seRNA

# class: RangedSummarizedExperiment
# dim: 20287 35582
# metadata(0):
# assays(1): counts
# rownames(20287): FAM138A OR4F5 ... S100B PRMT2
# rowData names(3): gene_name gene_id exonLength
# colnames(35582): CD34_32_R5:AAACCTGAGTATCGAA-1
# str(seRNA)
# str(seRNA@colData)
# str(seRNA@assays)
# str(seRNA@rowRanges)

# > colnames(colData(seRNA))
# [1] "Group"             "nUMI_pre"          "nUMI"            
# [4] "nGene"             "initialClusters"   "UMAP1"            
# [7] "UMAP2"             "Clusters"          "BioClassification"
# [10] "Barcode"

#> table(colData(seRNA)$BioClassification)
#
#        01_HSC 02_Early.Eryth  03_Late.Eryth  04_Early.Baso    05_CMP.LMPP
#          1425           1653            446            111           2260
#      06_CLP.1         07_GMP    08_GMP.Neut         09_pDC         10_cDC
#           903           2097           1050            544            325
#11_CD14.Mono.1 12_CD14.Mono.2   13_CD16.Mono         14_Unk       15_CLP.2
#          1800           4222            292            520            377
#      16_Pre.B           17_B      18_Plasma       19_CD8.N      20_CD4.N1
#           710           1711             62           1521           2470
#     21_CD4.N2       22_CD4.M      23_CD8.EM      24_CD8.CM          25_NK
#          2364           3539            796           2080           2143
#        26_Unk
#           161

colnames(colData(seRNA))
table(colData(seRNA)$BioClassification)

###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
################################################################################### from here we will describe the CONCEPTS and
################################################################################### will leave to add the R code at a later time

##################################  <> UNCONTRAINED :

# It would take all of the cells in your scATAC-seq experiment and attempt to align them
# to any of the cells in the scRNA-seq experiment

# We use the addGeneIntegrationMatrix() function. As mentioned above,
# this function accepts either a Seurat object or a RangedSummarizedExperiment object via the seRNA parameters.

projHeme2 = proj

projHeme2 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme2,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "BioClassification",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

################################# apparently this coresponds to the following DATASETS : BMMC :

# > seRNA
# class: RangedSummarizedExperiment
# dim: 20287 35582
# metadata(0):
# assays(1): counts
# rownames(20287): FAM138A OR4F5 ... S100B PRMT2
# rowData names(3): gene_name gene_id exonLength
# colnames(35582): CD34_32_R5:AAACCTGAGTATCGAA-1
#  CD34_32_R5:AAACCTGAGTCGTTTG-1 ...
#  BMMC_10x_GREENLEAF_REP2:TTTGTTGCATGTGTCA-1
#  BMMC_10x_GREENLEAF_REP2:TTTGTTGCATTGAAAG-1

##################################   <> CONSTRAINED :

# We use prior knowledge of the cell types to limit the search space of the alignment.

###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
################################################################################### adding PSEUDO ScRNA-profiles

projHeme3 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme2,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = TRUE,
    force= TRUE,
    # groupList = groupList,
    groupRNA = "BioClassification",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)

getAvailableMatrices(projHeme3)

projHeme3 <- addImputeWeights(projHeme3)

markerGenes  <- c(
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", #B-Cell Trajectory
    "CD14", #Monocytes
    "CD3D", "CD8A", "TBX21", "IL7R" #TCells
  )

p1 <- plotEmbedding(
    ArchRProj = projHeme3,
    colorBy = "GeneIntegrationMatrix",
    name = markerGenes,
    continuousSet = "horizonExtra",
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme3)
)

p2 <- plotEmbedding(
    ArchRProj = projHeme3,
    colorBy = "GeneScoreMatrix",
    continuousSet = "horizonExtra",
    name = markerGenes,
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme3)
)

p1c <- lapply(p1, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})

p2c <- lapply(p2, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})

do.call(cowplot::plot_grid, c(list(ncol = 3), p1c))

do.call(cowplot::plot_grid, c(list(ncol = 3), p2c))

plotPDF(plotList = p1,
    name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation.P1.pdf",
    ArchRProj = projHeme3,
    addDOC = FALSE, width = 5, height = 5)
   
plotPDF(plotList = p2,
    name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation.P2.pdf",
    ArchRProj = projHeme3,
    addDOC = FALSE, width = 5, height = 5)
   
plotPDF(plotList = p1c,
    name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation.P1c.pdf",
    ArchRProj = projHeme3,
    addDOC = FALSE, width = 5, height = 5)
   
plotPDF(plotList = p2c,
    name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation.P2c.pdf",
    ArchRProj = projHeme3,
    addDOC = FALSE, width = 5, height = 5)    
       
################################################################
#################### MAKING the READ COVERAGES :
################################################################ in order to call the peaks
################################################################ peaks can be called with MACS2 or with a TILE MATRIX
 
projHeme4 <- addGroupCoverages(ArchRProj = projHeme3, groupBy = "Clusters2")

###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
################################################################################### CALLING the PEAKS


pathToMacs2 = "/home/bogdan/Desktop/SOFTWARE/MACS2-2.2.7.1/bin/macs2" ### it does not work with MACS2 !!!

projHeme4 <- addReproduciblePeakSet(
    ArchRProj = projHeme4,
    groupBy = "Cluster",
    pathToMacs2 = pathToMacs2
)

getPeakSet(projHeme4)

###################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
################################################################################### CALLING the PEAKS with TILEMATRIX

projHemeTmp <- addReproduciblePeakSet(
    ArchRProj = projHeme4,
    groupBy = "Clusters",
    peakMethod = "Tiles",
    method = "p")

getPeakSet(projHemeTmp)

######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
