library(stringr)
library(Sushi)
library(rtracklayer)
library(data.table)
library(ggplot2)
library(ggpubr)
source('helper_peaks.R')

outputdir <- 'C:/Users/Admin/Desktop/'

rangeup <- 850000
rangedown <- 100000
gene2find <- 'AR'

#rangeup <- 200000
#rangedown <- 200000
#gene2find <- 'MYC'

#rangeup <- -3100000
#rangedown <- 2930000
#rangeup <- 100000
#rangedown <- 100000
#gene2find <- 'ERG'

#rangeup <- 500000
#rangedown <- 200000
#gene2find <- 'FOXA1'

#rangeup <- 100000
#rangedown <- 100000
#gene2find <- 'TMPRSS2'
#gene2find <- 'BRCA2'
#gene2find <- 'RB1'
#gene2find <- 'TP53'

#gene2find <- 'SCHLAP1'
#gene2find <- 'TACSTD2'

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
filename <- paste(outputdir, name,'.pdf',sep='')


if(gene2find=='ERG' || gene2find=='TMPRSS2') {
	group2compare <- rownames(genomic_calls)[genomic_calls$ERG_act_sv]
} else if(gene2find=='MYC') {
	group2compare <- rownames(genomic_calls)[genomic_calls$MYC_amp]
} else if(gene2find=='TP53') {
	group2compare <- rownames(genomic_calls)[genomic_calls$TP53_2loss]
} else if(gene2find=='BRCA2') {
	group2compare <- rownames(genomic_calls)[genomic_calls$BRCA2_2loss]
} else if(gene2find=='RB1') {
	group2compare <- rownames(genomic_calls)[genomic_calls$RB1_2loss]
} else { 
	group2compare <- rownames(genomic_calls)[genomic_calls$tSCNC]
}



########################################################################
######################### LOAD IN VARIOUS DATA #########################
########################################################################
bedgene <- exons[,c(5,6,7,3,1,4)]
bedgene$strand <- as.numeric(bedgene$strand=='+')
bedgene[bedgene$strand==0,'strand'] <- -1

## CHIPseq tracks
erg <- data.frame(tracks[['chip_ERG']])[,-4:-5]
colnames(erg) <- c('chrom','start','end','value')
ergvcap <- data.frame(tracks[['chip_ERG_VCaP']])[,-4:-5]
colnames(ergvcap) <- c('chrom','start','end','value')

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

## MANTA tracks
print('Sliding windows for SVs')
fusionrows <- (data_manta_slim$chr1==chr & inrange(data_manta_slim$pos1,start,end)) |
	(data_manta_slim$chr2==chr & inrange(data_manta_slim$pos2,start,end)) |
	(data_manta_slim$chr1==data_manta_slim$chr2 & data_manta_slim$chr1==chr &
		overlap(data_manta_slim$pos1,data_manta_slim$pos2,start,end)) &
	data_manta_slim$sample_id %in% samples_wgbs
fusionrows <- fusionrows & rows_manta_pass
fusionrows[is.na(fusionrows)] <- F

cols2keep <- c('chr1','pos1','chr2','pos2','svtype','sample_id')
plot_sv <- data_manta_slim[fusionrows,cols2keep]
rowsdel <- plot_sv$svtype=='MantaDEL'
rowsinv <- plot_sv$svtype=='MantaINV'
rowsins <- plot_sv$svtype=='MantaINS'
rowsbnd <- plot_sv$svtype=='MantaBND'
rowsdup <- plot_sv$svtype=='MantaDUP'
bed_mantadel <- svwindow(plot_sv[rowsdel,],chr,start,end,endsonly=F)
bed_mantainv <- svwindow(plot_sv[rowsinv,],chr,start,end)
bed_mantabnd <- svwindow(plot_sv[rowsbnd,],chr,start,end)
bed_mantadup <- svwindow(plot_sv[rowsdup,],chr,start,end,endsonly=F)
bed_mantains <- svwindow(plot_sv[rowsins,],chr,start,end)

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
cgi <- as.data.frame(tracks[['CGI']])[1:3]
colnames(cgi)[1] <- 'chr'
cgi$row <- 4

## correlate 5HMC peaks
rhmrtrack <- tracks[['rHMR']][countOverlaps(tracks[['rHMR']],grange)>0]
samples_5hmc_expr <- intersect(colnames(tpm),names(peaks_5hmc))
geneexpr <- as.numeric(tpm[ensemblid,])
names(geneexpr) <- colnames(tpm)

around <- 1000
peaks <- data.frame(rhmrtrack)[,1:3]
for(i in 1:length(rhmrtrack)) {
	print(i)
	starti <- peaks[i,'start']
	endi <- peaks[i,'end']
	startiaround <- peaks[i,'start']-around
	endiaround <- peaks[i,'end']+around
	
	hmci <- get_peak_enrich(chr,startiaround,endiaround,samples_5hmc_expr,peaks_5hmc,'fold_enrichment')
	corobj <- cor.test(hmci,geneexpr[samples_5hmc_expr],method='spearman',na.rm=T)
	peaks[i,'hmc_p'] <- corobj$p.value
	peaks[i,'hmc_cor'] <- corobj$estimate
	
	wgbsi <- as.numeric(data.frame(rhmrtrack[i,samples_wgbs])[-1:-5])
	corobj <- cor.test(	wgbsi,geneexpr[samples_wgbs],method='spearman',na.rm=T)
	peaks[i,'wgbs_p'] <- corobj$p.value
	peaks[i,'wgbs_cor'] <- corobj$estimate
}
peaks$expr_fdr <- p.adjust(peaks$hmc_p,method='fdr')
peaks$wgbs_fdr <- p.adjust(peaks$wgbs_p,method='fdr')
colnames(peaks)[1] <- 'chr'


#######################################################################
############################# PLOT TRACKS #############################
#######################################################################
pdf(file=filename, onefile=T, width=10, height=numplots/1.2)
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
fivehmcslice[fivehmcslice$sample %in% samples_5hmc_loc,'color'] <- colsglobal['loc']
fivehmcslice[fivehmcslice$sample %in% samples_5hmc_nt,'color'] <- colsglobal['nt']

rows2keep <- fivehmcslice$sample %in% c(samples_5hmc_loc,samples_5hmc_mcrpc,samples_5hmc_nt)
fivehmcslice <- fivehmcslice[rows2keep,]

roworder <- as.numeric(factor(fivehmcslice$sample,levels=
	c(samples_5hmc_nt,samples_5hmc_loc,group2compare,setdiff(samples_5hmc_mcrpc,group2compare))))
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

rows2keep <- hmrslice$sample %in% c(samples_wgbs,upitt_localized_tumor,upitt_benign_prostate,samples_normal)
hmrslice <- hmrslice[rows2keep,]

roworder <- as.numeric(factor(hmrslice$sample,levels=
	c(samples_normal,upitt_benign_prostate,upitt_localized_tumor,
	group2compare,setdiff(samples_wgbs,group2compare))))
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
