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

## Create a variable for .hic file filepaths
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

## Create a variable for megaMap .hic file filepaths
megaMap_files <- list.files("data",
                            pattern = "*megaMap*",
                            full.names = T)

## Replace megaMap_files names with shorter easier to read names
names(megaMap_files) <- c("megaMap_FS",
                          "megaMap_WT")

## Read in WT & FS loop files from marinerData package
wtLoops <- read.table(WT_5kbLoops.txt(),
                      header = T)

fsLoops <- read.table(FS_5kbLoops.txt(),
                      header = T)

## Convert loops into GInteractions objects & expand to 10kb resolution
wtLoopsGI <- wtLoops |> 
  as_ginteractions() |> 
  snapToBins(10e3)

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
summary(loopList)

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

## View information about the pixels before merging
sets(mergedLoops)

## Extract pixels
pixels <- pullHicPixels(
  x = mergedLoops,
  files = hicFiles,
  binSize = 10e3 ## loops were originally called at 5kb resolution
)

## Show attributes of pixels
pixels

## Show count matrix housed within `pixels`
counts(pixels)

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
                 coef="condition_WT_vs_FS",
                 type= "apeglm")

summary(res, alpha = 0.05)

## Inspect the results with a PCA plot
plotMA(res, alpha = 0.05)

## Inspect the results with a PCA plot
varianceStabilizingTransformation(dds) |>
  plotPCA(intgroup="condition") +
  ggplot2::theme(aspect.ratio = 1)

## We can then add these results from DESeq2 to our InteractionMatrix object
rowData(pixels) <- res

## Separate WT/FS-specific loops 
fsLoops <- pixels[which(res$padj <= 0.05 & res$log2FoldChange > 0)]
wtLoops <- pixels[which(res$padj <= 0.05 & res$log2FoldChange < 0)]

## Initiate plotgardener page
pageCreate(width = 4, height = 3,
           showGuides = T)

## Define shared parameters
p <- pgParams(assembly = "hg19",
              resolution = 10e3,
              chrom = "9",
              chromstart = 128399598 - 100e3,
              chromend = 128839674 + 100e3,
              zrange = c(0,100),
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

## Annotate WT loops
annoPixels(plot = wt_hic,
           data = interactions(wtLoops),
           type = "arrow",
           col = '#005AB5')

annoPixels(plot = wt_hic,
           data = interactions(fsLoops),
           type = "arrow",
           col = '#DC3220')

## Annotate FS loops
annoPixels(plot = fs_hic,
           data = interactions(wtLoops),
           type = "arrow",
           col = '#005AB5')

annoPixels(plot = fs_hic,
           data = interactions(fsLoops),
           type = "arrow",
           col = '#DC3220')

## Plot CTCF in WT ChIP-seq data
#https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144076/suppl/GSE144076%5FHEK293%5FCTCF%5FCHIP%5FN%5FIDR%5FWT%5FA9%5Fcell.bw
plotSignal(data = ,)

## Plot CTCF in FS ChIP-seq data
#https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144076/suppl/GSE144076%5FHEK293%5FCTCF%5FCHIP%5FN%5FIDR%5FFS%5FA9%5Fcell.bw
plotSignal(data = ,)

## Plot NHA9 in WT ChIP-seq data
plotSignal(data = ,)

## Plot NHA9 in FS ChIP-seq data
plotSignal(data = ,)

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
                x = 0.25,
                y = 2.825,
                scale = "Mb")
