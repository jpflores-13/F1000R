## Download visualization packages
library(plotgardener)
library(RColorBrewer)
library(InteractionSet)

## Concatenate WT & FS loop GRanges
fsLoops_50 <- head(fsLoops[order(interactions(fsLoops)$padj, decreasing = F)], 50)

# top 50 FS loops
fsLoops_gr <-
  GRanges(seqnames = as.character(seqnames(anchors(x = fsLoops, "first"))),
          ranges = IRanges(start = start(anchors(fsLoops, "first")),
                           end = end(anchors(fsLoops, "second"))),
          mcols = mcols(fsLoops))

## Add buffer 
buffer <- 200e3
fsLoops_gr_buffer <- fsLoops_gr + buffer

##make pdf
pdf(file = "plots/Flores_F1000R.pdf",
    width = 5.75,
    height = 6)

## Loop through each region
for(i in seq_along(fsLoops_gr_buffer)){ 
  
  ## Initiate plotgardener page
  pageCreate(width = 5.75, height = 6,
             showGuides = F)
  
  ## Define shared parameters
  p <- pgParams(assembly = "hg38",
                resolution = 10e3,
                chrom = as.character(seqnames(fsLoops_gr_buffer))[i],
                chromstart = start(fsLoops_gr_buffer)[i],
                chromend = end(fsLoops_gr_buffer)[i],
                zrange = c(0,100),
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
                  y = 5.25,
                  scale = "Mb")
  
}
dev.off()