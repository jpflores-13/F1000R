## Install Packages
BiocManager::install(c("mariner",
                       "marinerData",
                       "InteractionSet",
                       "data.table",
                       "plyranges",
                       "apeglm",
                       "DESeq2",
                       "plotgardener",
                       "RColorBrewer"))

## Load packages
library(mariner)
library(marinerData)
library(InteractionSet)
library(data.table)
library(plyranges)
library(DESeq2)
library(plotgardener)
library(RColorBrewer)

## Access WT Hi-C data from GEO
wt_hicFiles <-
  c("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4259nnn/GSM4259896/suppl/GSM4259896_HEK_HiC_NUP_IDR_WT_A9_1_1_inter_30.hic",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4259nnn/GSM4259897/suppl/GSM4259897_HEK_HiC_NUP_IDR_WT_A9_1_2_inter_30.hic",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4259nnn/GSM4259898/suppl/GSM4259898_HEK_HiC_NUP_IDR_WT_A9_2_1_inter_30.hic",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4259nnn/GSM4259899/suppl/GSM4259899_HEK_HiC_NUP_IDR_WT_A9_2_2_inter_30.hic")

## Download all WT .hic files of interest using a `wget` command
for (i in seq_along(wt_hicFiles)){
  system(paste0("wget -c ", wt_hicFiles[i], " -P data/"))
}

## Access FS Hi-C data from GEO
fs_hicFiles <- 
  c("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4259nnn/GSM4259900/suppl/GSM4259900_HEK_HiC_NUP_IDR_FS_A9_1_1_inter_30.hic",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4259nnn/GSM4259901/suppl/GSM4259901_HEK_HiC_NUP_IDR_FS_A9_1_2_inter_30.hic",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4259nnn/GSM4259902/suppl/GSM4259902_HEK_HiC_NUP_IDR_FS_A9_2_1_inter_30.hic",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4259nnn/GSM4259903/suppl/GSM4259903_HEK_HiC_NUP_IDR_FS_A9_2_2_inter_30.hic")

## Download all FS .hic files of interest using a `wget` command
for (i in seq_along(fs_hicFiles)){
  system(paste0("wget -c ", fs_hicFiles[i], " -P data/"))
}

## Create a variable for .hic file file paths
hicFiles <- list.files("data",
                       pattern = "GSM4259*",
                       full.names = T)

## replace hicFile names with shorter easier to read names
names(hicFiles) <- c("WT_rep1","WT_rep2","WT_rep3","WT_rep4", "FS_rep1","FS_rep2","FS_rep3","FS_rep4")

## Access megaMap Hi-C data from GEO
megaMap_hicFiles <- 
  c("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143465/suppl/GSE143465%5FHEK%5FHiC%5FNUP%5FIDR%5FFS%5FA9%5FmegaMap%5Finter%5F30.hic",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143465/suppl/GSE143465%5FHEK%5FHiC%5FNUP%5FIDR%5FWT%5FA9%5FmegaMap%5Finter%5F30.hic")

## Download all megaMap .hic files of interest using a `wget` command
for (i in seq_along(megaMap_hicFiles)){
  system(paste0("wget -c ", megaMap_hicFiles[i], " -P data/"))
}

## Create a variable for megaMap .hic file file paths
megaMap_files <- list.files("data",
                            pattern = "*megaMap*",
                            full.names = T)

## Replace megaMap_files names with shorter easier to read names
names(megaMap_files) <- c("megaMap_FS",
                          "megaMap_WT")

## Convert loops into GInteractions objects & expand to 10kb resolution
wtLoops <- fread(marinerData::WT_5kbLoops.txt())
wtLoopsGI <- wtLoops |> 
  as_ginteractions() |> 
  snapToBins(10e3)

fsLoops <- fread(marinerData::FS_5kbLoops.txt())
fsLoopsGI <- fsLoops |> 
  as_ginteractions() |> 
  snapToBins(10e3)

## Show summaries of GInteractions
summary(wtLoopsGI)
summary(fsLoopsGI)

## change seqlevelsStyle such that "chr" is removed from seqnames() columns
seqlevelsStyle(wtLoopsGI) <- "ENSEMBL"
seqlevelsStyle(fsLoopsGI) <- "ENSEMBL"

## Create list of loops
loopList <- list("WT" = wtLoopsGI, "FS" = fsLoopsGI)
loopList

## Remove redundant loops (i.e. loops that are present in both WT & FS)
mergedLoops <-
  mergePairs(loopList,
             radius = 10e3,
             method = "manhattan",
             pos = "center")

## View number of loops before merging
summary(loopList)

## View number of loops after merging
summary(mergedLoops)

## Extract pixels
pixels <- pullHicPixels(
  x = mergedLoops,
  files = hicFiles,
  binSize = 10e3)

## Show attributes of `pixels`
pixels

## Show count matrix housed within `pixels`
counts(pixels)

## Filter out loops with low counts (at least 10 counts in at least 4 samples)
keep <- rowSums(counts(pixels) >= 10) >= 4
pixels <- pixels[keep,]

## Construct colData
colData <- data.frame(condition = factor(rep(c("WT", "FS"), each = 4)),
                      replicate = factor(rep(1:4, 2)))

## Add rownames to colData for DESeq object
rownames(colData) <- colnames(counts(pixels))

## Ensure the colnames of the count matrix is equal to the rownames of the colData
all(colnames(counts(pixels)) == rownames(colData))

## Build a DESeq Dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts(pixels),
  colData = colData,
  design = ~ replicate + condition)

## Perform differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
dds <- DESeq(dds)
dds

## Get shrunken results
res <- lfcShrink(dds,
                 coef = "condition_WT_vs_FS",
                 type = "apeglm")

summary(res, alpha = 0.05)

## Inspect the results with a PCA plot
varianceStabilizingTransformation(dds) |>
  plotPCA(intgroup="condition") +
  ggplot2::theme(aspect.ratio = 1)

## Inspect the results with an MA plot
plotMA(res,
       alpha = 0.05,
       ylim = c(-4,4))

## Add results to rowData of our `pixels` InteractionMatrix
rowData(pixels) <- res

## Separate WT/FS-specific loops 
fsLoops <- pixels[rowData(pixels)$padj <= 0.05 & rowData(pixels)$log2FoldChange > 0]

wtLoops <- pixels[rowData(pixels)$padj <= 0.05 & rowData(pixels)$log2FoldChange < 0]

## Initiate plotgardener page
pageCreate(width = 4, height = 3,
           showGuides = F)

## Define shared parameters
p <- pgParams(assembly = "hg19",
              resolution = 10e3,
              chrom = "9",
              chromstart = 128399598 - 100e3,
              chromend = 128839674 + 100e3,
              zrange = c(0,150),
              norm = "SCALE",
              x = 0.25,
              width = 3.5,
              length = 3.5,
              height = 1)


## Plot WT Hi-C Mega Map
wt_hic <- plotHicRectangle(data = megaMap_files[["megaMap_WT"]],                             
                           params = p,
                           x = 0.25,
                           y = 0.25)


## Plot FS Hi-C Map
fs_hic <- plotHicRectangle(data = megaMap_files[["megaMap_FS"]],                             
                           params = p,
                           x = 0.25,
                           y = 1.3)

## Add legend for WT Hi-C map
annoHeatmapLegend(plot = wt_hic,
                  x = 3.85,
                  y = 0.25,
                  width = 0.1,
                  height = 0.75,
                  fontcolor = 'black')

## Add legend for FS Hi-C map
annoHeatmapLegend(plot = fs_hic,
                  x = 3.85,
                  y = 1.3,
                  width = 0.1,
                  height = 0.75,
                  fontcolor = 'black')

## Annotate loops
annoPixels(plot = wt_hic,
           data = interactions(fsLoops),
           type = "arrow",
           col = 'blue')

annoPixels(plot = fs_hic,
           data = interactions(fsLoops),
           type = "arrow",
           col = 'blue')


## Add text labels
plotText(label = "WT",
         x = 0.25,
         y = 0.25,
         just = c('top', 'left'))

plotText(label = "FS",
         x = 0.25,
         y = 1.3,
         just = c('top', 'left'))

## Add Genes + Gene labels
plotGenes(width = 3.5,
          height = 0.5, 
          chrom = paste0("chr", p$chrom),
          params = p,
          x = 0.25,
          y = 2.325)

plotGenomeLabel(params = p,
                chrom = paste0("chr", p$chrom),
                x = 0.25,
                y = 2.825)


# Create Survey Plot ------------------------------------------------------

## Take Top 50 FS loops
fsLoops_50 <- head(fsLoops[order(rowData(fsLoops)$padj, decreasing = F)], 50)

## Convert to GRanges Object
fsLoops_gr <-
  GRanges(seqnames = as.character(seqnames(anchors(x = fsLoops, "first"))),
          ranges = IRanges(start = start(anchors(fsLoops, "first")),
                           end = end(anchors(fsLoops, "second"))),
          mcols = mcols(fsLoops))

## Add buffer 
buffer <- 200e3
fsLoops_gr_buffer <- fsLoops_gr + buffer

## Make pdf
pdf(file = "plots/Flores_F1000R.pdf",
    width = 5.75,
    height = 6)

## Loop through each region
for(i in seq_along(fsLoops_gr_buffer)){ 
  
  ## Initiate plotgardener page
  pageCreate(width = 5.75, height = 6,
             showGuides = F)
  
  ## Define shared parameters
  p <- pgParams(assembly = "hg19",
                resolution = 10e3,
                chrom = as.character(seqnames(fsLoops_gr_buffer))[i],
                chromstart = start(fsLoops_gr_buffer)[i],
                chromend = end(fsLoops_gr_buffer)[i],
                zrange = c(0,150),
                norm = "SCALE",
                x = 0.25,
                width = 5,
                length = 5,
                height = 2)
  
  ## Plot WT Hi-C Mega Map
  wt_hic <- plotHicRectangle(data = megaMap_files[["megaMap_WT"]],                             
                             params = p,
                             x = 0.25,
                             y = 0.25)
  
  ## Plot FS Hi-C Map
  fs_hic <- plotHicRectangle(data = megaMap_files[["megaMap_FS"]],                             
                             params = p,
                             x = 0.25,
                             y = 2.3)
  
  ## Add legend for WT Hi-C map
  annoHeatmapLegend(plot = wt_hic,
                    x = 5.35,
                    y = 0.25,
                    width = 0.1,
                    height = 1.5,
                    fontcolor = 'black')
  
  ## Add legend for FS Hi-C map
  annoHeatmapLegend(plot = fs_hic,
                    x = 5.35,
                    y = 2.5,
                    width = 0.1,
                    height = 1.5,
                    fontcolor = 'black')
  
  ## Annotate FS loops
  annoPixels(plot = wt_hic,
             data = interactions(fsLoops_50),
             type = "arrow",
             col = '#DC3220')
  
  annoPixels(plot = fs_hic,
             data = interactions(fsLoops_50),
             type = "arrow",
             col = '#DC3220')
  
  ## Add text labels
  plotText(label = "WT",
           x = 0.25,
           y = 0.25,
           just = c('top', 'left'))
  
  plotText(label = "FS",
           x = 0.25,
           y = 2.35,
           just = c('top', 'left'))
  
  ## Add Genes + Gene labels
  plotGenes(width = 5,
            chrom = paste0("chr", p$chrom),
            params = p,
            x = 0.25,
            y = 4.25,
            height = 1)
  
  plotGenomeLabel(params = p,
                  chrom = paste0("chr", p$chrom),
                  x = 0.25,
                  y = 5.25)
}
dev.off()

## Get session information
sessionInfo()
