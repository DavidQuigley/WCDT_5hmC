###########################################################################################################
###########################################################################################################
# Figure 2C 
# 2022-01-05
# Martin Sjoestroem
###########################################################################################################
###########################################################################################################
# Directories and filenames
createFn <- function(fn, dirBase = dir_base) {paste0(dirBase, fn)}

dir_base <- "/home/martin/mnt/data1/projects/WCDT_5hmc_2019/publication_reproduce"

fn_load_5hmc_environment <- createFn("/data/hmc_211221.RData")
###########################################################################################################
# Load environment
###########################################################################################################
load(fn_load_5hmc_environment)
lapply(libraries_to_load,require,character.only = T)

###########################################################################################################
FCfilter <- 0
FDRfilter <- 0.00001

dataToPlot <- loc_mcrpc_res
dataToPlot <- dataToPlot[dataToPlot$FDR < FDRfilter & abs(dataToPlot$Fold) >= FCfilter, ]

###########################################################################################################
###########################################################################################################
# Idiogram
###########################################################################################################
###########################################################################################################

###########################################################################################################
## Prepare chromosome and gene density data
###########################################################################################################
# Load centromere and telomere info
centromere.coords <-as.data.frame(centromeres)
colnames(centromere.coords) <- c("chrom","start","end","width","strand")
telomere.coords <- as.data.frame(telomeres)
colnames(telomere.coords) <- c("chrom","start","end","width","strand")
telomere.coords <- telomere.coords[c(1:8,10,9,11:nrow(telomere.coords)),]
telomere.coords.5prime <- telomere.coords[seq(from=1,to=nrow(telomere.coords), by=2),]
telomere.coords.3prime <- telomere.coords[seq(from=2,to=nrow(telomere.coords), by=2),]

# Load protein coding gene info
genes.gencode.28 <- ensembl2sym[,c("chr","start","end","strand")]
all.genes <- makeGRangesFromDataFrame(genes.gencode.28)
pc.genes <- all.genes[ensembl2sym$type=="protein_coding",]

# Generate histogram bins tallying number of genes in each 1MB segment of the genome
chrlist <- list()
density_binsz <- 1000000
for(i in 1:nrow(hg38)) {
  currow <- hg38[i,]
  chrlen <- currow$length
  chrbreaks <- seq(from=1,to=chrlen, by=density_binsz)
  binstart <- chrbreaks
  binend <- chrbreaks + density_binsz - 1
  binend <- sapply(binend, min, chrlen)
  
  # Reformat chromosome names
  curchr <- currow$chrom
  if(curchr == 23) curchr <- 'X'
  if(curchr == 24) curchr <- 'Y'
  curchr <- paste('chr', curchr, sep='')
  
  chrtab <- cbind.data.frame(chr=curchr,start=binstart,end=binend)
  chrlist[[i]] <- chrtab
}

chrtab_full <- rbindlist(chrlist)
chrbins_full_gr <- makeGRangesFromDataFrame(chrtab_full)
gene_bins <- countOverlaps(chrbins_full_gr, pc.genes)

stopifnot(length(gene_bins) == nrow(chrtab_full))
chrtab_full$gene.counts <- gene_bins
chrtab_full$gened <- chrtab_full$gene.counts / (chrtab_full$end - chrtab_full$start + 1) 

# Scale 
chrtab_full$gened.n <- chrtab_full$gened / max(chrtab_full$gened)
chrtab_full$value <- 1

###########################################################################################################
## Prep the input dfata
###########################################################################################################
#
dataToPlot <- dataToPlot[,c("seqnames","start","end","Fold")]
colnames(dataToPlot) <- c("chr","start","end","value")
dataToPlot$intensity <- dataToPlot$value

mid <- dataToPlot

# Colors
scale_range <- c(-1.5,1.5)

mid$intensity[mid$intensity < scale_range[1]] <- scale_range[1]
mid$intensity[mid$intensity > scale_range[2]] <- scale_range[2]

pal <- colorNumeric("RdBu", domain = scale_range, reverse = T)
mid$color <- pal(mid$intensity)

###########################################################################################################
## Plot settings
##########################################################################################################
barwidth=10000
genes.highlight=NULL 
flank.width=0.25 
left.margin = 0.1
ymax <- 250000000/barwidth
numcol <- 1
featureWidth <- (1 - flank.width - left.margin) / numcol
chr_width <- 1
chrbar.width = chr_width - flank.width - left.margin

###########################################################################################################
## Plotting
##########################################################################################################
pdf(file = paste0(dir_out_figures, "/Figure_2C_idiogram.pdf"), height = 7, width = 12)

layout(matrix(1,1,1))
par(mar=c(2,4,1,1))

plot( -100, -100, xlim=c(1, 25.8), ylim=c(0, ymax), axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")

axis(2, at=seq(from=0, to=ymax, by=ymax/25), labels= seq(from=0, to=ymax, by=ymax/25) / 100, las=1)
axis(1, at=(1:24)+left.margin+chrbar.width/2, labels = c(1:22, 'X', 'Y'),las=1)
  
chrtab_full$gened.n <-  (chrtab_full$gened/max(chrtab_full$gened)) * flank.width * 1.5
chrtab_full$gened.n[chrtab_full$gened.n >= flank.width * 1.5] <- flank.width * 1.5 # Truncate outliers

for(i in 1:24) {
  chri <- paste0("chr",i)
  if(i==23) {
    chri <- 'chrX'
  }else if(i==24) {
    chri <- 'chrY'
  }
  print(paste('chromosome',chri))

  # Centromeric region
  centro_start <- centromere.coords[centromere.coords$chrom == paste(chri, sep=''), 'start']
  centro_end <- centromere.coords[centromere.coords$chrom == paste(chri, sep=''), 'end']
    
  # Ensure that centromere is at least some visible width
  if(centro_end - centro_start <= 200 * barwidth) {
    centro_start <- centro_end - 200 * barwidth
  }
    
  bins_centromere <- floor( round( centro_start / barwidth, 0 ) ) :
  ceiling( round( centro_end / barwidth, 0 ) )
  top_cent <- max(bins_centromere)
  bottom_cent <- min(bins_centromere)
  
  rowsmid <- mid$chr==chri
  print(paste('Number of mid bins =', sum(rowsmid)))
  
  midi <- mid[rowsmid,]
  
  gdi <- chrtab_full[chrtab_full$chr == chri,]
  
  #Don't draw gene density in centromeric regions (not reliable)
  gdi[gdi$start < centro_start & gdi$end > centro_start,'end'] <- centro_start 
  gdi[gdi$start >= centro_start & gdi$end <= centro_end,'gened.n'] <- 0
  gdi[gdi$start <= centro_end & gdi$end > centro_end,'start'] <- centro_end
    

  for(j in seq_len(nrow(midi))){
    rowindex_start <- (midi[j,'start'] / barwidth)
    rowindex_end <- max((midi[j,'end'] / barwidth))
    
    # ensure not to draw anything in telomeric regions
    rowindex_start <- max(rowindex_start, (telomere.coords.5prime[telomere.coords.5prime$chrom == paste('chr', chri, sep=''), 'chromEnd']/barwidth))
    rowindex_end <- min(rowindex_end, (telomere.coords.3prime[telomere.coords.3prime$chrom == paste('chr', chri, sep=''), 'chromStart']/barwidth))
    
    if(rowindex_end > rowindex_start) {
      rect( i + left.margin, 
            rowindex_start, 
            i+ left.margin + featureWidth, 
            rowindex_end, 
            col=midi[j,'color'], border=NA)
    }
  }

  # Plot mutation density
  for(j in seq_len(nrow(gdi))) {
    rowindex_start <- ceiling(gdi[j,'start']/barwidth)
    rowindex_end <- ceiling(gdi[j,'end']/barwidth)
    
    # ensure not to draw anything in telomeric regions
    rowindex_start <- max(rowindex_start, ceiling(telomere.coords.5prime[telomere.coords.5prime$chrom == paste('chr', chri, sep=''), 'chromEnd']/barwidth))
    rowindex_end <- min(rowindex_end, floor(telomere.coords.3prime[telomere.coords.3prime$chrom == paste('chr', chri, sep=''), 'chromStart']/barwidth))
    
    if(rowindex_end > rowindex_start) {
      rect( i+ left.margin + numcol*featureWidth, rowindex_start, i +left.margin + numcol*featureWidth + gdi[j, 'gened.n'], rowindex_end, col='lavenderblush', lwd=0.5)
    }
  }

  # Write out chromosome box, moved here to be above the lines.
  rect( i + left.margin, 0, i+chr_width - flank.width, 1 + (hg38[i,'length'] / barwidth), col=NA )
  
  # White out centromeric region
  rect(i+left.margin-0.05,bottom_cent+15,i+left.margin+chrbar.width+0.1, top_cent-15, col="white", border=NA)
  
  #polygon( c( i, i+1, mean(c(i, i+1))),
  polygon( c( i + left.margin + 0.1, i+left.margin + chrbar.width - 0.1, i + left.margin + chrbar.width/2),       # Make centromere triangle  
           c( top_cent, top_cent, mean(c(top_cent, bottom_cent)) ), col="gold")
  #polygon( c( i, i+1, mean(c(i, i+1))),
  polygon( c( i + left.margin + 0.1, i+left.margin + chrbar.width - 0.1, i + left.margin + chrbar.width/2),
           c( bottom_cent, bottom_cent, mean(c(top_cent, bottom_cent)) ), col="gold")
  lines( c(i+left.margin,i+1-flank.width), c(top_cent, top_cent), col="black")
  lines( c(i+left.margin,i+1-flank.width), c(bottom_cent, bottom_cent), col="black")
}

dev.off()
        
###########################################################################################################
## Output color key
##########################################################################################################
# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='', curwidth=10, sideText = "Differential methylation (%)") {
  scale = (length(lut)-1)/(max-min)
  
  #dev.new(width=1.75, height=5)
  plot(c(0,1), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', main=title, ylab = '')
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,curwidth,y+1/scale, col=lut[i], border=lut[i])
  }
  mtext(text = sideText,
        side = 2, #side 2 = left
        line = 2.5)
}

pal <- colorNumeric("RdBu", domain = scale_range, reverse = T)
min <- scale_range[1]
max <- scale_range[2]
pdf(file = paste0(dir_out_figures, "/Figure_2C_color_bar.pdf"), height = 4, width = 1.5)
color.bar(pal(seq(min,max, length.out = 200)), min, max, nticks=7, curwidth=1, title = "", sideText = "Differential 5hmC (fold change)")
dev.off()

###########################################################################################################
## Plot GREAT analysis
nTop = 10
qThresh <- 0.05
foldThresh <- 2

theDF <- subset(great_mcrpc_loc, BinomFdrQ < qThresh & HyperFdrQ < qThresh & RegionFoldEnrich > foldThresh)
theDF$BinomFdrQ[theDF$BinomFdrQ == 0] <-  1e-323 # set for appropriate rankning
theDF <- theDF[order(theDF$HyperFdrQ)[1:nTop],]

p <- ggplot(theDF, aes(x = reorder(Desc,-log10(HyperFdrQ)), y = -log10(HyperFdrQ))) + 
  geom_bar(stat = "identity", fill = "royalblue3") +
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 26),
        plot.title = element_text(size = 28, hjust = 1),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent")) + # get rid of legend panel bg) 
  ylab("-log10(Hypergeometric FDR)") +
  ggtitle("GO biological process using GREAT for top up regulated peaks in mCRPC")
p

pdf(paste0(dir_out_figures,"/Figure_2C_great.pdf"), height = 8, width = 15)
print(p)
dev.off()