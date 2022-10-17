###########################################################################################################
###########################################################################################################
# Figure 4,S7, tracks plotting
# 2022-05-04
# Martin Sjoestroem/George Zhao
###########################################################################################################
###########################################################################################################
# Directories and filenames
createFn <- function(fn, dirBase = dir_base) {paste0(dirBase, fn)}

dir_base <- "/home/martin/mnt/data1/projects/WCDT_5hmc_2019/publication_reproduce/cancer_research_revision_1/"

fn_load_5hmc_environment <- createFn("/data/hmc_220428.RData")
###########################################################################################################
# Load environment
load(fn_load_5hmc_environment)
lapply(libraries_to_load, library, character.only = T)

###########################################################################################################
###########################################################################################################
#y = start coord
#z = end coord
inrange <- function(x, y, z) {
  stopifnot(y<=z)
  return(x>=y & x<=z)
}

#y = start coord
#z = end coord
overlap <- function(x1, x2, y, z) {
  stopifnot(length(x1)==length(x2))
  if(length(x1)==1 && !is.na(x1) && !is.na(x2) && x1>x2) {
    x3 <- x2
    x2 <- x1
    x1 <- x3
  } else if(length(x1)>1) {
    rowswrong <- x1>x2
    rowswrong[is.na(rowswrong)] <- F
    x3 <- x2[rowswrong]
    x2[rowswrong] <- x1[rowswrong]
    x1[rowswrong] <- x3
  }
  stopifnot(y<=z)
  x1in <- x1>=y & x1<=z
  x2in <- x2>=y & x2<=z
  xencompasses <- x1<y & x2>z
  return(x1in | x2in | xencompasses)
}

###################################
peakflank <- 0
windowsize <- 1000

# data assumes chr1, pos1, chr2, pos2, svtype, sample_id
svwindow <- function(data, chr, start, end, window=windowsize, extra=peakflank, endsonly=T) {
  stopifnot((end-start)%%window==0)
  stopifnot(sum(is.na(data$chr1))==0)
  stopifnot(sum(colnames(data)!=c('chr1','pos1','chr2','pos2','svtype','sample_id'))==0)
  steps <- (end-start)/window
  output <- data.frame()
  for(i in 1:steps) {
    if(i %% 500 == 0) {
      print(paste(i, '/', steps))
    }
    left <- start+((i-1)*window)
    right <- left+window-1
    output[i,'chrom'] <- chr
    output[i,'start'] <- left
    output[i,'end'] <- right
    onechr <- is.na(data$chr2) | is.na(data$pos2)
    samechr <- data$chr1==data$chr2 & !onechr
    diffchr <- data$chr1!=data$chr2 & !onechr
    if(endsonly) {
      diffchr <- diffchr | samechr
      samechr <- rep(F,length(samechr))
    }
    left <- left-extra
    right <- right+extra
    rowsone <- onechr & inrange(data$pos1,left,right) & data$chr1==chr
    rowssame <- samechr & overlap(data$pos1,data$pos2,left,right) & data$chr1==chr
    rowsdiff <- diffchr &
      ((inrange(data$pos1,left,right) & data$chr1==chr) |
         (inrange(data$pos2,left,right) & data$chr2==chr))
    
    datarows <- rowsone | rowssame | rowsdiff
    stopifnot(sum(is.na(datarows))==0)
    output[i,'value'] <- length(unique(data[datarows,'sample_id']))
  }
  return(output)
}


# data assumes chrom, start, end, sample_id
slidingwindow <- function(data, chr, start, end, bysample=T, window=windowsize, extra=peakflank) {
  stopifnot((end-start)%%window==0)
  if(bysample) {
    stopifnot(sum(colnames(data)!=c('chrom','start','end','sample_id'))==0)
  } else {
    stopifnot(sum(colnames(data)!=c('chrom','start','end','value'))==0)
  }
  
  steps <- (end-start)/window
  output <- data.frame()
  for(i in 1:steps) {
    if(i %% 500 == 0) {
      print(paste(i, '/', steps))
    }
    left <- start+((i-1)*window)
    right <- left+window-1
    output[i,'chrom'] <- chr
    output[i,'start'] <- left
    output[i,'end'] <- right
    
    left <- left-extra
    right <- right+extra
    stopifnot((right-left+1) == (window+extra+extra))
    rows <- overlap(data$start,data$end,left,right) & data$chrom==chr
    stopifnot(sum(is.na(rows))==0)
    if(sum(rows)>0) {
      dataslice <- data[rows,]
      if(bysample) {
        output[i,'value'] <- length(unique(dataslice$sample_id))
      } else {
        #dataslice$start <- pmax(dataslice$start,start)
        #dataslice$end <- pmin(dataslice$start,end)
        #meansize <- mean(dataslice$end-dataslice$start+1)
        #dataslice$weight <- (dataslice$end-dataslice$start+1)/meansize
        #output[i,'value'] <- sum(dataslice$value*dataslice$weight,na.rm=T)
        output[i,'value'] <- sum(dataslice$value,na.rm=T)
      }
    } else {
      output[i,'value'] <- 0
    }
  }
  return(output)
}


#####
# Added MS
get_peak_enrich <- function(CHR, START, END, SAMPLES, PEAK_OBJ) {
  theRange <- GRanges(CHR,IRanges(START,END))
  thePeaks <- PEAK_OBJ[SAMPLES]
  subPeaks <- list()
  for(i in 1:length(thePeaks)){
    rows2keep <- countOverlaps(thePeaks[[i]],theRange)>0
    subPeaks[names(thePeaks[i])] <- thePeaks[[i]][rows2keep]
  }
  averageFC <- lapply(subPeaks, FUN = function(x){ max(x$fold_enrichment) })
  return(unlist(averageFC))
}

################################################################################
# Load/organize data and samples as in the Nat Gen 2020 paper.
rhmr <- tracks[["HMRsegs"]]
genomic_calls <- metadata_wcdt_genomics

upitt_benign_prostate <- c('Bis49AT','Bis159AT','Bis165AT','Bis171AT')
upitt_localized_tumor <- c('Bis49T','Bis158T','Bis159T','Bis165T','Bis171T')
upitt_benign_prostate <- toupper(upitt_benign_prostate)
upitt_localized_tumor <- toupper(upitt_localized_tumor)

########################################################################
######################### LOAD IN VARIOUS DATA #########################
########################################################################
ar_normal <- data.frame(pomerantz[['AR_Normal']]$consensusRanges)[,-4:-5]
colnames(ar_normal) <- c('chrom','start','end','value')
ar_normal$value <- ar_normal$value/max(ar_normal$value)
ar_primary <- data.frame(pomerantz[['AR_Tumour']]$consensusRanges)[,-4:-5]
colnames(ar_primary) <- c('chrom','start','end','value')
ar_primary$value <- ar_primary$value/max(ar_primary$value)
ar_met <- data.frame(pomerantz[['AR_PDX']]$consensusRanges)[,-4:-5]
colnames(ar_met) <- c('chrom','start','end','value')
ar_met$value <- ar_met$value/max(ar_met$value)

foxa1_normal <- data.frame(pomerantz[['FOXA_Normal']]$consensusRanges)[,-4:-5]
colnames(foxa1_normal) <- c('chrom','start','end','value')
foxa1_normal$value <- foxa1_normal$value/max(foxa1_normal$value)
foxa1_primary <- data.frame(pomerantz[['FOXA_Tumour']]$consensusRanges)[,-4:-5]
colnames(foxa1_primary) <- c('chrom','start','end','value')
foxa1_primary$value <- foxa1_primary$value/max(foxa1_primary$value)
foxa1_met <- data.frame(pomerantz[['FOXA_PDX']]$consensusRanges)[,-4:-5]
colnames(foxa1_met) <- c('chrom','start','end','value')
foxa1_met$value <- foxa1_met$value/max(foxa1_met$value)

hoxb13_normal <- data.frame(pomerantz[['HOXB13_Normal']]$consensusRanges)[,-4:-5]
colnames(hoxb13_normal) <- c('chrom','start','end','value')
hoxb13_normal$value <- hoxb13_normal$value/max(hoxb13_normal$value)
hoxb13_primary <- data.frame(pomerantz[['HOXB13_Tumour']]$consensusRanges)[,-4:-5]
colnames(hoxb13_primary) <- c('chrom','start','end','value')
hoxb13_primary$value <- hoxb13_primary$value/max(hoxb13_primary$value)
hoxb13_met <- data.frame(pomerantz[['HOXB13_PDX']]$consensusRanges)[,-4:-5]
colnames(hoxb13_met) <- c('chrom','start','end','value')
hoxb13_met$value <- hoxb13_met$value/max(hoxb13_met$value)

h3k27ac_normal <- data.frame(pomerantz[['H3k27ac_Normal']]$consensusRanges)[,-4:-5]
colnames(h3k27ac_normal) <- c('chrom','start','end','value')
h3k27ac_normal$value <- h3k27ac_normal$value/max(h3k27ac_normal$value)
h3k27ac_primary <- data.frame(pomerantz[['H3k27ac_Tumour']]$consensusRanges)[,-4:-5]
colnames(h3k27ac_primary) <- c('chrom','start','end','value')
h3k27ac_primary$value <- h3k27ac_primary$value/max(h3k27ac_primary$value)
h3k27ac_met <- data.frame(pomerantz[['H3k27ac_PDX']]$consensusRanges)[,-4:-5]
colnames(h3k27ac_met) <- c('chrom','start','end','value')
h3k27ac_met$value <- h3k27ac_met$value/max(h3k27ac_met$value)

h3k27me3_normal <- data.frame(pomerantz[['H3K27me3_Normal']]$consensusRanges)[,-4:-5]
colnames(h3k27me3_normal) <- c('chrom','start','end','value')
h3k27me3_normal$value <- h3k27me3_normal$value/max(h3k27me3_normal$value)
h3k27me3_primary <- data.frame(pomerantz[['H3K27me3_Tumour']]$consensusRanges)[,-4:-5]
colnames(h3k27me3_primary) <- c('chrom','start','end','value')
h3k27me3_primary$value <- h3k27me3_primary$value/max(h3k27me3_primary$value)

h3k4me2_normal <- data.frame(pomerantz[['H3K4me2_Normal']]$consensusRanges)[,-4:-5]
colnames(h3k4me2_normal) <- c('chrom','start','end','value')
h3k4me2_normal$value <- h3k4me2_normal$value/max(h3k4me2_normal$value)
h3k4me2_primary <- data.frame(pomerantz[['H3K4me2_Tumour']]$consensusRanges)[,-4:-5]
colnames(h3k4me2_primary) <- c('chrom','start','end','value')
h3k4me2_primary$value <- h3k4me2_primary$value/max(h3k4me2_primary$value)

dups = list_sv_m[list_sv_m$svtype=="TANDEM",]
dups = dups[, c('chrom_start','pos_start','chrom_end', 'pos_end', 'svtype', 'sample_id')] 
names(dups) = c('chr1','pos1','chr2','pos2','svtype','sample_id')

###########################################################################################################
###########################################################################################################
# Figure 4
###########################################################################################################
###########################################################################################################
rangeup <- 850000
rangedown <- 100000
gene2find <- 'AR'

# rangeup <- 500000
# rangedown <- 200000
# gene2find <- 'FOXA1'

numplots <- 11

rowensembl <- ensembl2sym$name==gene2find
ensemblid <- rownames(ensembl2sym)[rowensembl]
stopifnot(sum(rowensembl)==1)

chr <- ensembl2sym[rowensembl,'chr']
gene_start <- ensembl2sym[rowensembl,'start']
gene_end <- ensembl2sym[rowensembl,'end']

start <- round(gene_start-rangeup,-3)
end <- round(gene_end+rangedown,-3)
stopifnot(start<end)

grange <- GRanges(chr, IRanges(start, end))

name <- paste(gene2find,'loci',sep='_')
print(paste(name, ' - ', chr, ':', start, '-', end, sep=''))

group2compare <- rownames(genomic_calls)[genomic_calls$tSCNC]

bed_mantadup = svwindow(dups,chr,start,end,endsonly=FALSE, extra=0, window=1000)

## Load UMRs
all_hmrs <- list()
for(i in 1:length(hmr)) {
  samplename <- names(hmr)[i]
  rows2keep <- countOverlaps(hmr[[i]],grange)>0
  hmroverlap <- hmr[[i]][rows2keep]
  curtab <- data.frame(hmroverlap)
  rowindex <- dim(curtab)[1]+1
  curtab[rowindex,1] <- chr
  curtab[rowindex,2] <- end-10
  curtab[rowindex,3] <- end
  curtab$sample <- toupper(samplename)
  all_hmrs[[samplename]] <- curtab
}
hmrslice <- rbindlist(all_hmrs)[,c(1:3,11)]

## Load 5hmc
peaks_5hmc <- c(peaks_wcdt_tissue, peaks_localized)
all_5hmc <- list()
for(i in 1:length(peaks_5hmc)) {
  samplename <- names(peaks_5hmc)[i]
  rows2keep <- countOverlaps(peaks_5hmc[[i]],grange)>0
  peaksoverlap <- peaks_5hmc[[i]][rows2keep]
  curtab <- data.frame(peaksoverlap)
  rowindex <- dim(curtab)[1]+1
  curtab[rowindex,1] <- chr
  curtab[rowindex,2] <- end-10
  curtab[rowindex,3] <- end
  curtab$sample <- toupper(samplename)
  all_5hmc[[samplename]] <- curtab
}
fivehmcslice <- rbindlist(all_5hmc)[,c(1:3,13)]

## CpG islands
cgi <- as.data.frame(tracks[["cpg_island"]])[1:3]
colnames(cgi)[1] <- 'chr'
cgi$row <- 4

## correlate 5HMC peaks
rhmrlevel_granges <- makeGRangesFromDataFrame(rhmrlevel, keep.extra.columns = T)
rhmrtrack <- rhmrlevel_granges[countOverlaps(rhmrlevel_granges,grange)>0]
samples_5hmc_expr <- samples_mcrpc_uniq_all_data
geneexpr <- as.numeric(tpm[ensemblid,])
names(geneexpr) <- colnames(tpm)

###########################
around <- 1000
peaks <- data.frame(rhmrtrack)[,1:3]
for(i in 1:length(rhmrtrack)) {
  print(i)
  starti <- peaks[i,'start']
  endi <- peaks[i,'end']
  startiaround <- peaks[i,'start']-around
  endiaround <- peaks[i,'end']+around
  
  hmci <- get_peak_enrich(chr,startiaround,endiaround,samples_5hmc_expr,peaks_5hmc)
  
  corobj <- cor.test(hmci[samples_5hmc_expr],geneexpr[samples_5hmc_expr],method='spearman',na.rm=T)
  peaks[i,'hmc_p'] <- corobj$p.value
  peaks[i,'hmc_cor'] <- corobj$estimate

  wgbsi <- as.numeric(data.frame(rhmrtrack[i,samples_5hmc_expr])[-1:-5])
  corobj <- cor.test(	wgbsi,geneexpr[samples_5hmc_expr],method='spearman',na.rm=T)
  peaks[i,'wgbs_p'] <- corobj$p.value
  peaks[i,'wgbs_cor'] <- corobj$estimate
}

peaks$expr_fdr <- p.adjust(peaks$hmc_p,method='fdr')
peaks$wgbs_fdr <- p.adjust(peaks$wgbs_p,method='fdr')
colnames(peaks)[1] <- 'chr'

#######################################################################
############################# PLOT TRACKS #############################
#######################################################################
pdf(file=paste0(dir_out_figures,"/Figure_4.pdf"), onefile=T, width=10, height=numplots/1.2)

matrows <- c(1,1,2,2,2,2,3,3,3,3,4,5,5,6:numplots)
layout(matrix(matrows, length(matrows), 1, byrow=T))

## Plot genes
par(mai=c(0.4,1,0.05,0.5))
exonrows <- bedgene$chr==chr & overlap(bedgene$start, bedgene$end, start, end)
if(sum(exonrows)>0) {
  plotGenes(bedgene[exonrows,],chr,start,end)
  labelgenome(chr,start,end,n=20,scale="Mb")
} else {
  plot.new()
}

colsglobal <- c(met='firebrick3',metgroup='orange2',loc='goldenrod3',nt='springgreen3',benign='steelblue3')

par(mai=c(0.05,1,0.05,0.5))

## Plot 5hmc
fivehmcslice$start <- pmax(fivehmcslice$start,start+1)
fivehmcslice$end <- pmin(fivehmcslice$end,end-1)
colnames(fivehmcslice) <- c('chrom','start','end','sample')
fivehmcslice$score <- 1
fivehmcslice$strand <- 1

fivehmcslice$color <- colsglobal['met']
fivehmcslice[fivehmcslice$sample %in% group2compare,'color'] <- colsglobal['metgroup']
fivehmcslice[fivehmcslice$sample %in% samples_localized_tissue,'color'] <- colsglobal['loc']
fivehmcslice[fivehmcslice$sample %in% samples_wcdt_tissue_nt,'color'] <- colsglobal['nt']

rows2keep <- fivehmcslice$sample %in% c(samples_localized_tissue,samples_mcrpc_uniq_all_data,samples_wcdt_tissue_nt)
fivehmcslice <- fivehmcslice[rows2keep,]

roworder <- as.numeric(factor(fivehmcslice$sample,levels=
                                c(samples_wcdt_tissue_nt,samples_localized_tissue,group2compare,setdiff(samples_mcrpc_uniq_all_data,group2compare))))
plotBed(fivehmcslice,chr,start,end,row='supplied',rownumber=roworder,color=fivehmcslice$color)
axis(side=2,las=2,tcl=.2)
mtext('5hmc',side=2,line=4,cex=0.5)

## Plot HMR
hmrslice$start <- pmax(hmrslice$start,start+1)
hmrslice$end <- pmin(hmrslice$end,end-1)
colnames(hmrslice) <- c('chrom','start','end','sample')
hmrslice$score <- 1
hmrslice$strand <- 1

hmrslice$color <- colsglobal['met']
hmrslice[hmrslice$sample %in% group2compare,'color'] <- colsglobal['metgroup']
hmrslice[hmrslice$sample %in% upitt_localized_tumor,'color'] <- colsglobal['loc']
hmrslice[hmrslice$sample %in% upitt_benign_prostate,'color'] <- colsglobal['benign']
hmrslice[hmrslice$sample %in% samples_normal,'color'] <- colsglobal['nt']

rows2keep <- hmrslice$sample %in% c(samples_mcrpc_uniq_all_data,upitt_localized_tumor,upitt_benign_prostate,samples_normal)
hmrslice <- hmrslice[rows2keep,]

roworder <- as.numeric(factor(hmrslice$sample,levels=
                                c(samples_normal,upitt_benign_prostate,upitt_localized_tumor,
                                  group2compare,setdiff(samples_mcrpc_uniq_all_data,group2compare))))
plotBed(hmrslice,chr,start,end,row='supplied',rownumber=roworder,color=hmrslice$color)
axis(side=2,las=2,tcl=.2)
mtext('HMR',side=2,line=4,cex=0.5)

plotBed(cgi,chr,start,end,type="region",row='supplied',rownumber=1)
axis(side=2,las=2,tcl=.2)
mtext('CGI',side=2,line=4,cex=0.5)
abline(h=0, col='black')
legend("topleft",inset=0.01,legend=c('met','sub','loc','NT','benign'),
       fill=opaque(colsglobal),border=colsglobal,cex=0.5)

cols <- c(colsglobal['benign'],colsglobal['loc'],colsglobal['met'])
legend("topright",inset=0,legend=c('prostate','PCa','PDX'),
       fill=opaque(cols),border=cols,cex=0.5)

## Plot correlation with expression
cor_range <- c(-1,1)
colscor <- c('darkred','dodgerblue3')
plotBedgraph(peaks[,c('chr','start','end','hmc_cor')],chr,start,end,color=colscor[2],range=cor_range,transparency=0.6)
plotBedgraph(peaks[,c('chr','start','end','wgbs_cor')],chr,start,end,color=colscor[1],range=cor_range,transparency=0.6,overlay=T)
legend("topright",inset=0,legend=c('WGBS','5hMC'),
       fill=opaque(colscor),border=colscor,cex=0.5)
axis(side=2,las=2,tcl=.2)
mtext('Rho',side=2,line=4,cex=0.5)
abline(h=0, col='black')

## Plot chipseq
plotBedgraph(h3k27ac_normal,chr,start,end,color=cols[1],range=c(0,1))
plotBedgraph(h3k27ac_primary,chr,start,end,color=cols[2],overlay=T)
plotBedgraph(h3k27ac_met,chr,start,end,color=cols[3],overlay=T)
axis(side=2,las=2,tcl=.2)
mtext("H3K27ac",side=2,line=4,cex=0.5)
abline(h=0, col='black')

plotBedgraph(ar_normal,chr,start,end,color=cols[1],range=c(0,1))
plotBedgraph(ar_primary,chr,start,end,color=cols[2],overlay=T)
plotBedgraph(ar_met,chr,start,end,color=cols[3],overlay=T)
axis(side=2,las=2,tcl=.2)
mtext("AR",side=2,line=4,cex=0.5)
abline(h=0, col='black')

plotBedgraph(foxa1_normal,chr,start,end,color=cols[1],range=c(0,1))
plotBedgraph(foxa1_primary,chr,start,end,color=cols[2],overlay=T)
plotBedgraph(foxa1_met,chr,start,end,color=cols[3],overlay=T)
axis(side=2,las=2,tcl=.2)
mtext("FOXA1",side=2,line=4,cex=0.5)
abline(h=0, col='black')

plotBedgraph(hoxb13_normal,chr,start,end,color=cols[1],range=c(0,1))
plotBedgraph(hoxb13_primary,chr,start,end,color=cols[2],overlay=T)
plotBedgraph(hoxb13_met,chr,start,end,color=cols[3],overlay=T)
axis(side=2,las=2,tcl=.2)
mtext("HOXB13",side=2,line=4,cex=0.5)
abline(h=0, col='black')

cols <- c('darkmagenta',colsglobal['met'])
plotBedgraph(ergvcap,chr,start,end,color=cols[1])
plotBedgraph(erg,chr,start,end,color=cols[2],overlay=T)
axis(side=2,las=2,tcl=.2)
mtext('ERG',side=2,line=4,cex=0.5)
abline(h=0, col='black')
legend("topright",inset=0,legend=c('VCaP','Met'),
       fill=opaque(cols),border=cols,cex=0.5)

plotBedgraph(bed_mantadup,chr,start,end,color='green')
axis(side=2,las=2,tcl=.2)
mtext("DUP",side=2,line=4,cex=0.5)

dev.off()

###########################################################################################################
###########################################################################################################
# Figure S7
###########################################################################################################
###########################################################################################################
rangeup <- 500000
rangedown <- 200000
gene2find <- 'FOXA1'

numplots <- 11

rowensembl <- ensembl2sym$name==gene2find
ensemblid <- rownames(ensembl2sym)[rowensembl]
stopifnot(sum(rowensembl)==1)

chr <- ensembl2sym[rowensembl,'chr']
gene_start <- ensembl2sym[rowensembl,'start']
gene_end <- ensembl2sym[rowensembl,'end']

start <- round(gene_start-rangeup,-3)
end <- round(gene_end+rangedown,-3)
stopifnot(start<end)

grange <- GRanges(chr, IRanges(start, end))

name <- paste(gene2find,'loci',sep='_')
print(paste(name, ' - ', chr, ':', start, '-', end, sep=''))

group2compare <- rownames(genomic_calls)[genomic_calls$tSCNC]

bed_mantadup = svwindow(dups,chr,start,end,endsonly=FALSE, extra=0, window=1000)

## Load UMRs
all_hmrs <- list()
for(i in 1:length(hmr)) {
  samplename <- names(hmr)[i]
  rows2keep <- countOverlaps(hmr[[i]],grange)>0
  hmroverlap <- hmr[[i]][rows2keep]
  curtab <- data.frame(hmroverlap)
  rowindex <- dim(curtab)[1]+1
  curtab[rowindex,1] <- chr
  curtab[rowindex,2] <- end-10
  curtab[rowindex,3] <- end
  curtab$sample <- toupper(samplename)
  all_hmrs[[samplename]] <- curtab
}
hmrslice <- rbindlist(all_hmrs)[,c(1:3,11)]

## Load 5hmc
peaks_5hmc <- c(peaks_wcdt_tissue, peaks_localized)
all_5hmc <- list()
for(i in 1:length(peaks_5hmc)) {
  samplename <- names(peaks_5hmc)[i]
  rows2keep <- countOverlaps(peaks_5hmc[[i]],grange)>0
  peaksoverlap <- peaks_5hmc[[i]][rows2keep]
  curtab <- data.frame(peaksoverlap)
  rowindex <- dim(curtab)[1]+1
  curtab[rowindex,1] <- chr
  curtab[rowindex,2] <- end-10
  curtab[rowindex,3] <- end
  curtab$sample <- toupper(samplename)
  all_5hmc[[samplename]] <- curtab
}
fivehmcslice <- rbindlist(all_5hmc)[,c(1:3,13)]

## CpG islands
cgi <- as.data.frame(tracks[["cpg_island"]])[1:3]
colnames(cgi)[1] <- 'chr'
cgi$row <- 4

## correlate 5HMC peaks
rhmrlevel_granges <- makeGRangesFromDataFrame(rhmrlevel, keep.extra.columns = T)
rhmrtrack <- rhmrlevel_granges[countOverlaps(rhmrlevel_granges,grange)>0]
samples_5hmc_expr <- samples_mcrpc_uniq_all_data
geneexpr <- as.numeric(tpm[ensemblid,])
names(geneexpr) <- colnames(tpm)

###########################
around <- 1000
peaks <- data.frame(rhmrtrack)[,1:3]
for(i in 1:length(rhmrtrack)) {
  print(i)
  starti <- peaks[i,'start']
  endi <- peaks[i,'end']
  startiaround <- peaks[i,'start']-around
  endiaround <- peaks[i,'end']+around
  
  hmci <- get_peak_enrich(chr,startiaround,endiaround,samples_5hmc_expr,peaks_5hmc)
  
  corobj <- cor.test(hmci[samples_5hmc_expr],geneexpr[samples_5hmc_expr],method='spearman',na.rm=T)
  peaks[i,'hmc_p'] <- corobj$p.value
  peaks[i,'hmc_cor'] <- corobj$estimate
  
  wgbsi <- as.numeric(data.frame(rhmrtrack[i,samples_5hmc_expr])[-1:-5])
  corobj <- cor.test(	wgbsi,geneexpr[samples_5hmc_expr],method='spearman',na.rm=T)
  peaks[i,'wgbs_p'] <- corobj$p.value
  peaks[i,'wgbs_cor'] <- corobj$estimate
}

peaks$expr_fdr <- p.adjust(peaks$hmc_p,method='fdr')
peaks$wgbs_fdr <- p.adjust(peaks$wgbs_p,method='fdr')
colnames(peaks)[1] <- 'chr'

#######################################################################
############################# PLOT TRACKS #############################
#######################################################################
pdf(file=paste0(dir_out_figures,"/Figure_S7.pdf"), onefile=T, width=10, height=numplots/1.2)

matrows <- c(1,1,2,2,2,2,3,3,3,3,4,5,5,6:numplots)
layout(matrix(matrows, length(matrows), 1, byrow=T))

## Plot genes
par(mai=c(0.4,1,0.05,0.5))
exonrows <- bedgene$chr==chr & overlap(bedgene$start, bedgene$end, start, end)
if(sum(exonrows)>0) {
  plotGenes(bedgene[exonrows,],chr,start,end)
  labelgenome(chr,start,end,n=20,scale="Mb")
} else {
  plot.new()
}

colsglobal <- c(met='firebrick3',metgroup='orange2',loc='goldenrod3',nt='springgreen3',benign='steelblue3')

par(mai=c(0.05,1,0.05,0.5))

## Plot 5hmc
fivehmcslice$start <- pmax(fivehmcslice$start,start+1)
fivehmcslice$end <- pmin(fivehmcslice$end,end-1)
colnames(fivehmcslice) <- c('chrom','start','end','sample')
fivehmcslice$score <- 1
fivehmcslice$strand <- 1

fivehmcslice$color <- colsglobal['met']
fivehmcslice[fivehmcslice$sample %in% group2compare,'color'] <- colsglobal['metgroup']
fivehmcslice[fivehmcslice$sample %in% samples_localized_tissue,'color'] <- colsglobal['loc']
fivehmcslice[fivehmcslice$sample %in% samples_wcdt_tissue_nt,'color'] <- colsglobal['nt']

rows2keep <- fivehmcslice$sample %in% c(samples_localized_tissue,samples_mcrpc_uniq_all_data,samples_wcdt_tissue_nt)
fivehmcslice <- fivehmcslice[rows2keep,]

roworder <- as.numeric(factor(fivehmcslice$sample,levels=
                                c(samples_wcdt_tissue_nt,samples_localized_tissue,group2compare,setdiff(samples_mcrpc_uniq_all_data,group2compare))))
plotBed(fivehmcslice,chr,start,end,row='supplied',rownumber=roworder,color=fivehmcslice$color)
axis(side=2,las=2,tcl=.2)
mtext('5hmc',side=2,line=4,cex=0.5)

## Plot HMR
hmrslice$start <- pmax(hmrslice$start,start+1)
hmrslice$end <- pmin(hmrslice$end,end-1)
colnames(hmrslice) <- c('chrom','start','end','sample')
hmrslice$score <- 1
hmrslice$strand <- 1

hmrslice$color <- colsglobal['met']
hmrslice[hmrslice$sample %in% group2compare,'color'] <- colsglobal['metgroup']
hmrslice[hmrslice$sample %in% upitt_localized_tumor,'color'] <- colsglobal['loc']
hmrslice[hmrslice$sample %in% upitt_benign_prostate,'color'] <- colsglobal['benign']
hmrslice[hmrslice$sample %in% samples_normal,'color'] <- colsglobal['nt']

rows2keep <- hmrslice$sample %in% c(samples_mcrpc_uniq_all_data,upitt_localized_tumor,upitt_benign_prostate,samples_normal)
hmrslice <- hmrslice[rows2keep,]

roworder <- as.numeric(factor(hmrslice$sample,levels=
                                c(samples_normal,upitt_benign_prostate,upitt_localized_tumor,
                                  group2compare,setdiff(samples_mcrpc_uniq_all_data,group2compare))))
plotBed(hmrslice,chr,start,end,row='supplied',rownumber=roworder,color=hmrslice$color)
axis(side=2,las=2,tcl=.2)
mtext('HMR',side=2,line=4,cex=0.5)

plotBed(cgi,chr,start,end,type="region",row='supplied',rownumber=1)
axis(side=2,las=2,tcl=.2)
mtext('CGI',side=2,line=4,cex=0.5)
abline(h=0, col='black')
legend("topleft",inset=0.01,legend=c('met','sub','loc','NT','benign'),
       fill=opaque(colsglobal),border=colsglobal,cex=0.5)

cols <- c(colsglobal['benign'],colsglobal['loc'],colsglobal['met'])
legend("topright",inset=0,legend=c('prostate','PCa','PDX'),
       fill=opaque(cols),border=cols,cex=0.5)

## Plot correlation with expression
cor_range <- c(-1,1)
colscor <- c('darkred','dodgerblue3')
plotBedgraph(peaks[,c('chr','start','end','hmc_cor')],chr,start,end,color=colscor[2],range=cor_range,transparency=0.6)
plotBedgraph(peaks[,c('chr','start','end','wgbs_cor')],chr,start,end,color=colscor[1],range=cor_range,transparency=0.6,overlay=T)
legend("topright",inset=0,legend=c('WGBS','5hMC'),
       fill=opaque(colscor),border=colscor,cex=0.5)
axis(side=2,las=2,tcl=.2)
mtext('Rho',side=2,line=4,cex=0.5)
abline(h=0, col='black')

## Plot chipseq
plotBedgraph(h3k27ac_normal,chr,start,end,color=cols[1],range=c(0,1))
plotBedgraph(h3k27ac_primary,chr,start,end,color=cols[2],overlay=T)
plotBedgraph(h3k27ac_met,chr,start,end,color=cols[3],overlay=T)
axis(side=2,las=2,tcl=.2)
mtext("H3K27ac",side=2,line=4,cex=0.5)
abline(h=0, col='black')

plotBedgraph(ar_normal,chr,start,end,color=cols[1],range=c(0,1))
plotBedgraph(ar_primary,chr,start,end,color=cols[2],overlay=T)
plotBedgraph(ar_met,chr,start,end,color=cols[3],overlay=T)
axis(side=2,las=2,tcl=.2)
mtext("AR",side=2,line=4,cex=0.5)
abline(h=0, col='black')

plotBedgraph(foxa1_normal,chr,start,end,color=cols[1],range=c(0,1))
plotBedgraph(foxa1_primary,chr,start,end,color=cols[2],overlay=T)
plotBedgraph(foxa1_met,chr,start,end,color=cols[3],overlay=T)
axis(side=2,las=2,tcl=.2)
mtext("FOXA1",side=2,line=4,cex=0.5)
abline(h=0, col='black')

plotBedgraph(hoxb13_normal,chr,start,end,color=cols[1],range=c(0,1))
plotBedgraph(hoxb13_primary,chr,start,end,color=cols[2],overlay=T)
plotBedgraph(hoxb13_met,chr,start,end,color=cols[3],overlay=T)
axis(side=2,las=2,tcl=.2)
mtext("HOXB13",side=2,line=4,cex=0.5)
abline(h=0, col='black')

cols <- c('darkmagenta',colsglobal['met'])
plotBedgraph(ergvcap,chr,start,end,color=cols[1])
plotBedgraph(erg,chr,start,end,color=cols[2],overlay=T)
axis(side=2,las=2,tcl=.2)
mtext('ERG',side=2,line=4,cex=0.5)
abline(h=0, col='black')
legend("topright",inset=0,legend=c('VCaP','Met'),
       fill=opaque(cols),border=cols,cex=0.5)

plotBedgraph(bed_mantadup,chr,start,end,color='green')
axis(side=2,las=2,tcl=.2)
mtext("DUP",side=2,line=4,cex=0.5)

dev.off()