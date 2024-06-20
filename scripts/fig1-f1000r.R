## Figure 1 for F1000R Workflow

## Load packages
library(mariner)
library(marinerData)
library(InteractionSet)
library(data.table)
library(plyranges)
library(glue)

## Create a variable for .hic file filepaths
hicFile <- list.files("data",
                      pattern = "*megaMap*",
                      full.names = T)

## Read in WT & FS loop files from marinerData package
loops <- list(WT_5kbLoops.txt(),
              FS_5kbLoops.txt()) |> 
  lapply(read.table, header = T)

## Combine loop calls and convert into GInteractions objects
loopsGI <- do.call(rbind, loops) |> 
  as_ginteractions() |> 
  swapAnchors()

## Remove redundant loops
mergedLoops <-
  mergePairs(loopsGI,
             radius = 10e3) |> 
  swapAnchors()

## change seqlevelsStyle such that "chr" is removed from seqnames() columns
seqlevelsStyle(loopsGI) <- "ENSEMBL"
seqlevelsStyle(mergedLoops) <- "ENSEMBL"

## Find merged loops
unmerged_gr <-
  GRanges(seqnames = as.character(seqnames(anchors(x = loopsGI, "first"))),
          ranges = IRanges(start = start(anchors(loopsGI, "first")),
                           end = end(anchors(loopsGI, "second"))),
          mcols = mcols(loopsGI))

merged_gr <-
  GRanges(seqnames = as.character(seqnames(anchors(x = mergedLoops, "first"))),
          ranges = IRanges(start = start(anchors(mergedLoops, "first")),
                           end = end(anchors(mergedLoops, "second"))),
          mcols = mcols(mergedLoops))

consensusCalls <- setdiff(unmerged_gr, merged_gr)

## set buffer for plotting consensusCalls
buffer <- 300e3
consensusCalls_buff <- consensusCalls + buffer

## Download visualization packages
library(plotgardener)
library(RColorBrewer)
library(stringr)

pdf(width = 4.25,
    height = 3.5,
    file = "plots/consensus-survey-plots.pdf")

## Initiate plotgardener page
pageCreate(width = 4.25, height = 3.5,
           showGuides = F)

## Define shared parameters
p <- pgParams(assembly = "hg19",
              resolution = 10e3,
              chrom = as.character(seqnames(consensusCalls_buff))[18],
              chromstart = start(consensusCalls_buff)[18] + 125e3,
              chromend = end(consensusCalls_buff)[18] + 110e3,
              zrange = c(0,100),
              norm = "SCALE",
              x = 0.25,
              width = 3.5,
              length = 3.5,
              height = 1.25)

## Plot unmerged Hi-C Mega Map
unmerged_hic <- plotHicRectangle(data = hicFile[1],                             
                                 params = p,
                                 chrom = str_remove(p$chrom, "chr"),
                                 chromstart = p$chromstart,
                                 chromend = p$chromend,
                                 x = 0.25,
                                 y = 0.25)

## Add legend for unmerged Hi-C map
annoHeatmapLegend(plot = unmerged_hic,
                  x = 3.85,
                  y = 0.25,
                  width = 0.1,
                  height = 0.75,
                  fontcolor = 'black')

## Plot merged Hi-C Mega Map
merged_hic <- plotHicRectangle(data = hicFile[2],                             
                               params = p,
                               chrom = str_remove(p$chrom, "chr"),
                               chromstart = p$chromstart,
                               chromend = p$chromend,
                               x = 0.25,
                               y = 1.55)

## Add legend for merged Hi-C map
annoHeatmapLegend(plot = merged_hic,
                  x = 3.85,
                  y = 1.55,
                  width = 0.1,
                  height = 0.75,
                  fontcolor = 'black')

## Annotate loops before merging
annoPixels(plot = unmerged_hic,
           data = loopsGI,
           type = "arrow",
           col = '#DC3220')

## Annotate loops before merging
annoPixels(plot = unmerged_hic,
           data = mergedLoops,
           type = "arrow",
           col = '#005AB5')

## Annotate loops after merging
annoPixels(plot = merged_hic,
           data = loopsGI,
           type = "arrow",
           col = '#DC3220')

## Annotate loops after merging
annoPixels(plot = merged_hic,
           data = mergedLoops,
           type = "arrow",
           col = '#005AB5')

## Add text labels
plotText(label = "WT",
         x = 0.25,
         y = 0.25,
         just = c('top', 'left'))

plotText(label = "FS",
         x = 0.25,
         y = 1.55,
         just = c('top', 'left'))

## Add Genes + Gene labels
plotGenes(width = 3.5,
          height = 0.5, 
          chrom = paste0("chr", p$chrom),
          params = p,
          x = 0.25,
          y = 2.8)

plotGenomeLabel(params = p,
                x = 0.25,
                y = 3.35)

dev.off()