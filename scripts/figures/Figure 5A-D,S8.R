###########################################################################################################
###########################################################################################################
# Figure 5,S8
# 2022-06-02
# Martin Sjoestroem
###########################################################################################################
###########################################################################################################
# Directories and filenames
createFn <- function(fn, dirBase = dir_base) {paste0(dirBase, fn)}

dir_base <- "/home/martin/mnt/data1/projects/WCDT_5hmc_2019/publication_reproduce/cancer_research_revision_1/"

fn_load_5hmc_environment <- createFn("/data/hmc_220428.RData")
###########################################################################################################
# Load environment
###########################################################################################################
load(fn_load_5hmc_environment)
lapply(libraries_to_load,require,character.only = T)
###########################################################################################################
matrix_hmc <- counts_wcdt_tissue_tpm
matrix_hmc <- as.matrix(matrix_hmc[,samples_mcrpc_uniq_all_data])
matrix_gex <- as.matrix(tpm[,samples_mcrpc_uniq_all_data])

nameToEnsemblV28 <- function(gene){
  ensembl <- rownames(ensembl2sym)[which(ensembl2sym$name == gene)]
  if(length(ensembl)>1){
    print("Warning! More than one match!")
  }
  return(ensembl)
}

getGeneRegion <- function(gene, upstream = 10000, downstream = 10000){
  idx <- which(ensembl2sym$name==gene)
  stopifnot(length(idx)==1)
  chrom_name <- ensembl2sym$chr[idx] 
  chrom_start <- ensembl2sym$start[idx] - upstream
  chrom_end <- ensembl2sym$end[idx] + downstream
  return(c(chrom_name,chrom_start,chrom_end))
}

makeDataFrame <- function(gr){
  new.df <- as.data.frame(gr)
  new.df <- new.df[,!colnames(new.df) %in% c("width","strand")]
  return(new.df)
}

getChromosomePeakCoverage <- function(peaks = peaks_wcdt_tissue, theSamples, chromName){
  if( ! class(peaks)=="CompressedGRangesList" ){
    peaks <- GRangesList(peaks)
  }
  cov <- coverage(peaks[theSamples]) 
  gr <- GRanges(seqnames = chromName, ranges = IRanges(start = start(cov[[chromName]]), end = end(cov[[chromName]])), strand = "*")
  gr$value <- runValue(cov[[chrom_name]])
  df <- makeDataFrame(gr)
  colnames(df) <- c("chrom","start","end","count")
  df$count <- (df$count/length(theSamples))*100
  return(df)
}

boxplot_jitter <- function(data_vector, group_vector, Title = "Title", xLab = NULL, yLab = NULL, ftLevels = NULL, plotPval = F, pXpos = 1.5, pYpos = 10, Colors = c("red","green"), 
                           savePlot = F, file_name = NULL,fileW = 10, fileH = 7){
  
  data <- data.frame(value = data_vector, 
                     name = group_vector)
  
  if(!is.null(ftLevels)){
    data$name[data$name == "FALSE"] <- ftLevels[1]
    data$name[data$name == "TRUE"] <- ftLevels[2]
  }
  
  # Add pval
  pval <- signif(wilcox.test(data$value~data$name)$p.value,2)
  
  p <- ggplot(data, aes(x=name, y=value, fill=name)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(color="black", size=1, alpha=0.9, width = 0.2) +
    ggtitle(Title) +
    xlab(xLab) + 
    ylab(yLab) +
    scale_fill_manual(values = Colors) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(color="Black", size=14, face="bold", hjust = 0.5))
  if(plotPval){
    p <- p + 
      geom_text(aes(x = pXpos, y = pYpos, label = paste0("P = ",pval)))
  }
  print(p)
  
  if(savePlot){
    for(fileType in c("png","pdf")){
      if(fileType == "pdf"){
        pdf(file = paste0(file_name,".",fileType), width = fileW, height = fileH)
      }else if(fileType == "png"){
        png(file = paste0(file_name,".",fileType), width = fileW, height = fileH, units = "in", res = 300)
      }
      print(p)
      dev.off()
    }
  }
}

plot_gene_track <- function( chrom_name, chrom_start, chrom_end){
  hg_tmp = hg38_models[hg38_models$chrom==chrom_name & 
                         hg38_models$txStart>chrom_start-500000 & 
                         hg38_models$txEnd< chrom_end + 500000,]
  if( dim(hg_tmp)[1] > 0 ){
    tmp_starts = c()
    tmp_ends = c()
    tmp_names = c()
    tmp_strands = c()
    for(i in 1:dim(hg_tmp)[1]){
      new_starts = as.numeric( strsplit(hg_tmp$exonStarts[i], ",", fixed=TRUE)[[1]] )
      new_ends = as.numeric( strsplit(hg_tmp$exonEnds[i], ",", fixed=TRUE)[[1]] )
      tmp_strands = c(tmp_strands, rep( hg_tmp$strand[i], length(new_starts)))
      tmp_starts = c(tmp_starts, new_starts)
      tmp_ends = c(tmp_ends, new_ends)
      tmp_names = c(tmp_names, rep(hg_tmp$name2[i], length(new_starts)))
    }
    tmp_chroms = rep(chrom_name, length(tmp_names))
    hg38_exons_tmp = data.frame(
      chrom=tmp_chroms,
      start=tmp_starts,
      end=tmp_ends,
      name=tmp_names,
      score = rep('.', length(tmp_chroms)),
      strand=tmp_strands,
      stringsAsFactors = FALSE)
    
    
    pg=plotGenes(geneinfo=hg38_exons_tmp, chrom_name, chrom_start, chrom_end,
                 plotgenetype="box", type=rep("exon", dim(hg38_exons_tmp)[1]),
                 bheight=0.15, 
                 labeltext=TRUE, 
                 packrow=FALSE,
                 fontsize=1, labeloffset = 0.5,
                 fonttype=3)
  }else{
    plot(-1,-1, xlim=c(0,10), ylim=c(0,10), xlab="", ylab="", axes=FALSE)   
  }
}

plot_peak_track <- function( peaktrack, chrom_name, chrom_start, chrom_end, label_text="",range = NULL, textyPos=NULL, Col =  "darkblue", plotLines = FALSE){
  data = data.frame( chrom=peaktrack$chrom,
                     start=peaktrack$start,
                     end=peaktrack$end,
                     value=peaktrack$count,
                     stringsAsFactors=FALSE)
  
  plotBedgraph(data,chrom_name,chrom_start,chrom_end, range = range,
               transparency=.50,
               color=Col, flip=FALSE)
  axis(side=2,las=2,tcl=.2)
  if(plotLines){
    for( x in seq(from=10,to=100,by=10)){
      abline( x,0, col="#0000ff33")   
    }
  }
  
  abline( 0,0, col="black")   
  
  if( length(label_text)>0 ){
    maxy = 0.9
    if( dim(data)[1]>0 ){
      dvals = data$value[data$start>=chrom_start & data$end<=chrom_end & data$chrom==chrom_name]
      if( length(dvals)>0 ){
        maxy = max(dvals, na.rm=TRUE) * 0.9
      }
    }
    
    if(is.null(textyPos)){
      textyPos = maxy
    }
    
    text( chrom_start + (chrom_end-chrom_start)/50, textyPos, 
          adj=0, label_text, font=2, col=Col)
  }
  data
}

################################################################################################################################
# Plot boxplots
################################################################################################################################
# ERG
boxplot_jitter(log2(matrix_gex[nameToEnsemblV28("ERG"),samples_mcrpc_uniq_all_data]+1),metadata_wcdt_genomics[samples_mcrpc_uniq_all_data,"TMPRSS2_ERG"], plotPval = T, pYpos = 6 ,Title = "ERG GEX", ftLevels = c("T2E-neg","T2E-pos"), yLab = "log2(TPM+1)", Colors = col_neg_pos,savePlot = F, fileW = 2, fileH = 5,  file_name = paste0(dir_out_figures,"/Figure_5D_gex"))

boxplot_jitter(log2(matrix_hmc[nameToEnsemblV28("ERG"),samples_mcrpc_uniq_all_data]+1), metadata_wcdt_genomics[samples_mcrpc_uniq_all_data,"TMPRSS2_ERG"], plotPval = T, pYpos = 3 ,Title = "ERG 5hmC", ftLevels = c("T2E-neg","T2E-pos"), yLab = "log2(TPM+1)", Colors = col_neg_pos, savePlot = F, fileW = 2, fileH = 5, file_name = paste0(dir_out_figures,"/Figure_5D_hmc"))

# KCNS3
boxplot_jitter(log2(matrix_gex[nameToEnsemblV28("KCNS3"),samples_mcrpc_uniq_all_data]+1),metadata_wcdt_genomics[samples_mcrpc_uniq_all_data,"TMPRSS2_ERG"], plotPval = T, pYpos = 5.2 ,Title = "KCNS3 GEX", ftLevels = c("T2E-neg","T2E-pos"), yLab = "log2(TPM+1)", Colors = col_neg_pos, savePlot = T, fileW = 2, fileH = 5,  file_name = paste0(dir_out_figures,"/Figure_S8B_gex"))
boxplot_jitter(log2(matrix_hmc[nameToEnsemblV28("KCNS3"),samples_mcrpc_uniq_all_data]+1), metadata_wcdt_genomics[samples_mcrpc_uniq_all_data,"TMPRSS2_ERG"], plotPval = T, pYpos = 3.6 ,Title = "KCNS3 5hmC", ftLevels = c("T2E-neg","T2E-pos"), yLab = "log2(TPM+1)", Colors = col_neg_pos, savePlot = T, fileW = 2, fileH = 5, file_name = paste0(dir_out_figures,"/Figure_S8B_hmc"))

################################################################################################################################
# Plot tracksplot
################################################################################################################################
# ERG
################################################################################################################################
theGene <- "ERG"
gene_start <- ensembl2sym[ensembl2sym$name==theGene,"start"]
gene_end <- ensembl2sym[ensembl2sym$name==theGene,"end"]

gene_loc <- getGeneRegion(theGene, upstream = 100000, downstream = 100000)

chrom_name <- gene_loc[1]
plot_region_start <- as.numeric(gene_loc[2])
plot_region_end <- as.numeric(gene_loc[3])
  
data_start <- plot_region_start-100000
data_end <- plot_region_end+100000

group_neg <- as.character(metadata_wcdt_genomics$sample[metadata_wcdt_genomics$sample %in% samples_mcrpc_uniq_all_data & metadata_wcdt_genomics$TMPRSS2_ERG == FALSE])
group_pos <-  as.character(metadata_wcdt_genomics$sample[metadata_wcdt_genomics$sample %in% samples_mcrpc_uniq_all_data & metadata_wcdt_genomics$TMPRSS2_ERG == TRUE])

group_neg_cov <- getChromosomePeakCoverage(peaks = peaks_wcdt_tissue, theSamples = group_neg, chromName = chrom_name)
group_pos_cov <- getChromosomePeakCoverage(peaks = peaks_wcdt_tissue, theSamples = group_pos, chromName = chrom_name)

group_neg_cov_hmr <- getChromosomePeakCoverage(peaks = hmr, theSamples = group_neg, chromName = chrom_name)
group_pos_cov_hmr <- getChromosomePeakCoverage(peaks = hmr, theSamples = group_pos, chromName = chrom_name)

breakends <- curated_sv[curated_sv$threeprime=="ERG" & 
                          curated_sv$sample %in% metadata_wcdt_genomics$sample[metadata_wcdt_genomics$TMPRSS2_ERG==T] &
                          curated_sv$sample %in% samples_mcrpc_uniq_all_data,"pos3"]

png(filename = paste0(dir_out_figures,"/Figure_5A.png"),width = 7, height = 5, units = "in", res = 300)
nrow = 5
heights = NULL
if(is.null(heights)){
  heights = rep(1, nrow)
  heights[1] = 2
  heights[ length(heights) ] = 2
}

layout(matrix(1:nrow,nrow,1), heights=heights)
par(mar=c(0,3,1,1))
plot_gene_track(chrom_name = chrom_name, chrom_start = plot_region_start, chrom_end = plot_region_end)
points(breakends,rep(2,length(breakends)), col = "red", pch = 5, cex = 1)
par(mar=c(0,3,1,1))
d=plot_peak_track(group_pos_cov, chrom_name = chrom_name, chrom_start = plot_region_start, chrom_end = plot_region_end, 
                  label_text="5hmC peaks, T2E-pos", range = c(0,100),textyPos = 90)
par(mar=c(0,3,1,1))
d=plot_peak_track(group_neg_cov, chrom_name = chrom_name, chrom_start = plot_region_start, chrom_end = plot_region_end, 
                  label_text="5hmC peaks, T2E-neg", range = c(0,100),textyPos = 90)
par(mar=c(0,3,1,1))
d=plot_peak_track(group_pos_cov_hmr, chrom_name = chrom_name, chrom_start = plot_region_start, chrom_end = plot_region_end, 
                  label_text="HMRs, T2E-pos", range = c(0,100),textyPos = 90, Col = "chartreuse3")
par(mar=c(5,3,1,1))
d=plot_peak_track(group_neg_cov_hmr, chrom_name = chrom_name, chrom_start = plot_region_start, chrom_end = plot_region_end, 
                  label_text="HMRs, T2E-neg", range = c(0,100),textyPos = 90, Col = "chartreuse3")

labelgenome(chrom = chrom_name, chromstart = plot_region_start, chromend = plot_region_end,n=10,scale="Mb")
dev.off()

# Plot density around breakend
adjustedPeaks <- GRangesList() 
for(i in 1:length(group_pos)){
  theSample <- group_pos[i]
  thisData <- peaks_wcdt_tissue[theSample][[1]]
  thisData <- thisData[seqnames(thisData) == chrom_name & start(thisData) > data_start & end(thisData) < data_end]
  breakend = curated_sv[curated_sv$threeprime=="ERG" & curated_sv$sample == theSample,"pos3"]
  # adjust peaks to position relative to breakend, add 1Mb to avoid Granges problem of not handling negataive values. 1Mb will be the breakpoint
  newData <- IRanges::shift(thisData, (-breakend+1000000))
  adjustedPeaks[[theSample]] <- newData
}

adj_cov <- getChromosomePeakCoverage(peaks = adjustedPeaks, theSamples = names(adjustedPeaks), chromName = chrom_name)

nrow = 1
heights = NULL
if(is.null(heights)){
  heights = rep(1, nrow)
  heights[1] = 2
  heights[ length(heights) ] = 2
}

png(filename = paste0(dir_out_figures,"/Figure_5C.png"),width = 7, height = 5, units = "in", res = 300)
layout(matrix(1:nrow,nrow,1), heights=heights)
par(mar=c(5,4,3,3))
d=plot_peak_track(adj_cov, chrom_name = chrom_name, chrom_start = 900000, chrom_end = 1100000, label_text="", range = c(0,35),textyPos = 70)
abline(v = 1000000, lty = 2, col = "red", lwd =2)
mtext("0",1,1,cex = 1)
mtext("-100kb",1,1,cex = 1, at = 900000)
mtext("+100kb",1,1,cex = 1, at = 1100000)
mtext("Position relative to breakend",1,3,cex = 1)
mtext("% with peak",2,2,cex = 1)
title("Location of 5hmC peaks in ERG relative to the fusion breakend")
dev.off()

################################################################################################################################
# TMPRSS2
################################################################################################################################
theGene <- "TMPRSS2"
gene_start <- 41464305 #from ucsc, as our ensemble annotation has long transcript
gene_end <- 41508158

chrom_name <- "chr21"
plot_region_start <- gene_start -100000
plot_region_end <- gene_end + 100000

data_start <- plot_region_start-100000
data_end <- plot_region_end+100000

group_neg_cov <- getChromosomePeakCoverage(peaks = peaks_wcdt_tissue, theSamples = group_neg, chromName = chrom_name)
group_pos_cov <- getChromosomePeakCoverage(peaks = peaks_wcdt_tissue, theSamples = group_pos, chromName = chrom_name)

group_neg_cov_hmr <- getChromosomePeakCoverage(peaks = hmr, theSamples = group_neg, chromName = chrom_name)
group_pos_cov_hmr <- getChromosomePeakCoverage(peaks = hmr, theSamples = group_pos, chromName = chrom_name)

breakends <- curated_sv[curated_sv$fiveprime=="TMPRSS2" & 
                          curated_sv$sample %in% metadata_wcdt_genomics$sample[metadata_wcdt_genomics$TMPRSS2_ERG==T] &
                          curated_sv$sample %in% samples_mcrpc_uniq_all_data,"pos5"]

png(filename = paste0(dir_out_figures,"/Figure_5B.png"),width = 7, height = 5, units = "in", res = 300)
nrow = 5
heights = NULL
if(is.null(heights)){
  heights = rep(1, nrow)
  heights[1] = 2
  heights[ length(heights) ] = 2
}

layout(matrix(1:nrow,nrow,1), heights=heights)
par(mar=c(0,3,1,1))
plot_gene_track(chrom_name = chrom_name, chrom_start = plot_region_start, chrom_end = plot_region_end)
points(breakends,rep(2,length(breakends)), col = "red", pch = 5, cex = 1)
par(mar=c(0,3,1,1))
d=plot_peak_track(group_pos_cov, chrom_name = chrom_name, chrom_start = plot_region_start, chrom_end = plot_region_end, 
                  label_text="5hmC peaks, T2E-pos", range = c(0,100),textyPos = 90)
par(mar=c(0,3,1,1))
d=plot_peak_track(group_neg_cov, chrom_name = chrom_name, chrom_start = plot_region_start, chrom_end = plot_region_end, 
                  label_text="5hmC peaks, T2E-neg", range = c(0,100),textyPos = 90)
par(mar=c(0,3,1,1))
d=plot_peak_track(group_pos_cov_hmr, chrom_name = chrom_name, chrom_start = plot_region_start, chrom_end = plot_region_end, 
                  label_text="HMRs, T2E-pos", range = c(0,100),textyPos = 90, Col = "chartreuse3")
par(mar=c(5,3,1,1))
d=plot_peak_track(group_neg_cov_hmr, chrom_name = chrom_name, chrom_start = plot_region_start, chrom_end = plot_region_end, 
                  label_text="HMRs, T2E-neg", range = c(0,100),textyPos = 90, Col = "chartreuse3")

labelgenome(chrom = chrom_name, chromstart = plot_region_start, chromend = plot_region_end,n=10,scale="Mb")
dev.off()

###########################################################################################################
# KCNS3
###########################################################################################################

differential_consensus_peaks <- read.delim(fn_diffbind_t2e, sep = "\t", header = T, stringsAsFactors = F)
differential_consensus_peaks_for_plot <- differential_consensus_peaks[, c("seqnames", "start", "end", "p.value")]
differential_consensus_peaks_for_plot$count <- (-1)*log10(differential_consensus_peaks_for_plot$p.value)
differential_consensus_peaks_for_plot$p.value <- NULL
colnames(differential_consensus_peaks_for_plot) <- c("chrom", "start", "end", "count")
differential_consensus_peaks_for_plot <- differential_consensus_peaks_for_plot[with(differential_consensus_peaks_for_plot, order(chrom, start)),]
row.names(differential_consensus_peaks_for_plot) <- NULL

differential_consensus_peaks <- differential_consensus_peaks[, !str_detect(colnames(differential_consensus_peaks), "DTB")]
differential_consensus_peaks$loci <- paste0(differential_consensus_peaks$seqnames, ":", differential_consensus_peaks$start, "-", differential_consensus_peaks$end)

erg_vcap <- read.delim(fn_chip_erg_vcap, header = F)
colnames(erg_vcap) <- c("chrom","start","end","count")
erg_vcap_gr <- makeGRangesFromDataFrame(erg_vcap)

motif_annotation <- read.delim(fn_ERG_motif_annotation)

motif_annotation$PeakID..cmd.annotatePeaks.pl..home.meng.Desktop.marlowe_data1.projects.WCDT_5hmc_2019.analysis_by_meng.differential_5hmc_diffbind.differential_5hmc_consensus_peaks_no_control.ERG_fusion.significant_up_peaks_ERG_fusion.bed.hg38..m..home.meng.Desktop.marlowe_data1.projects.WCDT_5hmc_2019.analysis_by_meng.motif_search.consensus_peaks_no_control.ERG_fusion.up.knownResults.known1.motif..size.given..mbed..home.meng.Desktop.marlowe_data1.projects.WCDT_5hmc_2019.analysis_by_meng.differential_5hmc_diffbind.differential_5hmc_consensus_peaks_no_control.ERG_fusion.significant_up_peaks_ERG_fusion_hommer_ERG_motif_annotation.bed. <- NULL
motif_annotation$motif <- motif_annotation[, colnames(motif_annotation)[length(colnames(motif_annotation))]]
motif_annotation$loci <- paste0(motif_annotation$Chr, ":", (motif_annotation$Start-1), "-", motif_annotation$End)

motif_bed <- read.delim(fn_ERG_motif_bed, skip = 1, header = F)
colnames(motif_bed) <- c("chrom", "start", "end", "motif", "score", "strand")
motif_bed$count <- 1
motif_bed <- motif_bed[, c("chrom", "start", "end", "count", "motif", "strand", "score")]
motif_bed <- motif_bed[with(motif_bed, order(start)),]
motif_bed_plot <- motif_bed[,c("chrom","start","end","score")]
colnames(motif_bed_plot) <- c("chrom","start","end","count")

differential_consensus_peaks$erg_motif <- motif_annotation$motif[match(differential_consensus_peaks$loci, motif_annotation$loci)]

# Plot
theGene <- "KCNS3"
gene_start <- ensembl2sym$start[ensembl2sym$name==theGene]
gene_end <- ensembl2sym$end[ensembl2sym$name==theGene] # a longer transcript here

gene_loc <- getGeneRegion(theGene, upstream = 50000, downstream = 10000)

chrom_name <- gene_loc[1]
plot_region_start <- as.numeric(gene_loc[2])
plot_region_end <- as.numeric(gene_loc[3])-300000

group_neg <- as.character(metadata_wcdt_genomics$sample[metadata_wcdt_genomics$sample %in% samples_mcrpc_uniq_all_data & metadata_wcdt_genomics$TMPRSS2_ERG == FALSE])
group_pos <-  as.character(metadata_wcdt_genomics$sample[metadata_wcdt_genomics$sample %in% samples_mcrpc_uniq_all_data & metadata_wcdt_genomics$TMPRSS2_ERG == TRUE])

group_neg_cov <- getChromosomePeakCoverage(peaks = peaks_wcdt_tissue, theSamples = group_neg, chromName = chrom_name)
group_pos_cov <- getChromosomePeakCoverage(peaks = peaks_wcdt_tissue, theSamples = group_pos, chromName = chrom_name)

png(filename = paste0(dir_out_figures,"/Figure_S8A.png"),width = 7, height = 6, units = "in", res = 300)
nrow = 7
heights = NULL
if(is.null(heights)){
  heights = rep(1, nrow)
  heights[1] = 2
  heights[ length(heights) ] = 2
}
layout(matrix(1:nrow,nrow,1), heights=heights)
par(mar=c(0,3,1,1))
plot_gene_track(chrom_name = chrom_name, chrom_start = plot_region_start, chrom_end = plot_region_end)

par(mar=c(0,3,1,1))
d=plot_peak_track(group_pos_cov, chrom_name = chrom_name, chrom_start = plot_region_start, chrom_end = plot_region_end, 
                  label_text="5hmC peaks, T2E-pos", range = c(0,100),textyPos = 90)

par(mar=c(0,3,1,1))
d=plot_peak_track(group_neg_cov, chrom_name = chrom_name, chrom_start = plot_region_start, chrom_end = plot_region_end, 
                  label_text="5hmC peaks, T2E-neg", range = c(0,100),textyPos = 90)

par(mar=c(0,3,1,1))
d=plot_peak_track(differential_consensus_peaks_for_plot, chrom_name = chrom_name, chrom_start = plot_region_start, chrom_end = plot_region_end, 
                  label_text="-log10(p-val) diff. 5hmC", range = c(0,11),textyPos = (10))

par(mar=c(0,3,1,1))
d=plot_peak_track(motif_bed_plot, chrom_name = chrom_name, chrom_start = plot_region_start, chrom_end = plot_region_end, 
                  label_text="ERG motif score in 5hmC peaks", range = c(0,15),textyPos = (12))

par(mar=c(0,3,1,1))
d=plot_peak_track(erg_vcap, chrom_name = chrom_name, chrom_start = plot_region_start, chrom_end = plot_region_end,
                  label_text="ERG ChIP VCaP", range = c(0,35),textyPos = (30))
labelgenome(chrom = chrom_name, chromstart = plot_region_start, chromend = plot_region_end,n=10,scale="Mb")
dev.off()