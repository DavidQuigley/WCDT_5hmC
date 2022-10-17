###########################################################################################################
# Prepare environment for analysis
# 2022-05-02
# Martin Sjoestroem
###########################################################################################################
###########################################################################################################
# Libraries
libraries_to_load <- c("stringr",
                       "stringi",
                       "dplyr",
                       "plyr",
                       "tidyverse",
                       "factoextra",
                       "readxl",
                       "data.table",
                       "ggplot2",
                       "ggrepel",
                       "RColorBrewer",
                       "viridis",
                       "ggthemes",
                       "ggformula",
                       "reshape2",
                       "ggsci",
                       "ggpubr",
                       "cowplot",
                       "gridExtra",
                       "ggfortify",
                       "gridExtra",
                       "circlize",
                       "rtracklayer",
                       "Sushi",
                       "rCGH",
                       "plotrix",
                       "leaflet",
                       "rstatix",
                       "matrixStats",
                       "ComplexHeatmap",
                       "TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "GenomicFeatures",
                       "abind",
                       "foreach",
                       "scales",
                       "DESeq2",
                       "fgsea",
                       "msigdbr",
                       "MASS",
                       "qdapTools",
                       "survival",
                       "rms",
                       "survminer",
                       "ChIPseeker",
                       "ChIPpeakAnno")

lapply(libraries_to_load,require,character.only = T)

###########################################################################################################
###########################################################################################################
# Directories and filenames
createFn <- function(fn, dirBase = dir_base) {paste0(dirBase, fn)}

dir_base <- "/home/martin/mnt/data1/projects/WCDT_5hmc_2019/publication_reproduce/cancer_research_revision_1/" 

dir_out_figures <- createFn("/figures")
fn_environment_out <- createFn("/data/hmc_220428.RData")

## WGS data
fn_wgs_environment <- createFn("/data/WGS_rdata_NO_UPLOAD/2018_04_25_build_image.RData")
fn_cn_gencodev28 <- createFn("/data/WGS_additional_data/copy_number_gencode_v28_93.txt") 


## WGBS data
fn_wgbs_environment <- createFn("/data/WGBS_rdata_NO_UPLOAD/2019_05_15_prepare_environment.RData")
fn_hmrseg_upd <- createFn("/data/WGBS_additional_data/hmrseg_upd.txt") 
fn_hmrlevels <- createFn("/data/WGBS_additional_data/rhmr.tsv")

fn_promoter_methylation <- createFn("/data/WGBS_additional_data/promoter_methylation_wcdt_tissue.txt")

## ChIP-seq
fn_chip_erg <- createFn("/data/hmc_secondary_data/processed_chipseq/erg_hg38.bed")
fn_chip_erg_vcap <- createFn("/data/hmc_secondary_data/processed_chipseq/erg_vcap_hg38.bed")

fn_pomerantz <- createFn("/data/hmc_secondary_data/processed_chipseq/pomerantz_2020/allsamples_consensusPeaks_withCounts.RDS")

## 5hmC secondary data
# Peaks
fn_peaks_wcdt_tissue <- createFn("/data/hmc_secondary_data/peaks/peaks_wcdt_tissue_incl.rds")
fn_peaks_wcdt_cfdna <- createFn("/data/hmc_secondary_data/peaks/peaks_wcdt_cfdna_incl.rds")
fn_peaks_localized <- createFn("/data/hmc_secondary_data/peaks/peaks_localized_incl.rds")
fn_peaks_ubc_cfdna <- createFn("/data/hmc_secondary_data/peaks/peaks_ubc_cfdna_incl.rds")

# Counts
fn_counts_wcdt_tissue_raw <- createFn("/data/hmc_secondary_data/counts/counts_wcdt_tissue_103_raw.txt")
fn_counts_wcdt_cfdna_raw <- createFn("/data/hmc_secondary_data/counts/counts_wcdt_cfdna_15_raw.txt")
fn_counts_localized_raw <- createFn("/data/hmc_secondary_data/counts/counts_localized_52_raw.txt")
fn_counts_ubc_cfdna_raw <- createFn("/data/hmc_secondary_data/counts/counts_ubc_cfdna_64_raw.txt")
fn_counts_benign_prostate_raw <- createFn("/data/hmc_secondary_data/counts/counts_benign_prostate_raw.txt")

fn_counts_wcdt_tissue_tpm <- createFn("/data/hmc_secondary_data/counts/counts_wcdt_tissue_103_tpm.txt")
fn_counts_wcdt_cfdna_tpm <- createFn("/data/hmc_secondary_data/counts/counts_wcdt_cfdna_15_tpm.txt")
fn_counts_localized_tpm <- createFn("/data/hmc_secondary_data/counts/counts_localized_52_tpm.txt")
fn_counts_ubc_cfdna_tpm <- createFn("/data/hmc_secondary_data/counts/counts_ubc_cfdna_64_tpm.txt")
fn_counts_benign_prostate_tpm <- createFn("/data/hmc_secondary_data/counts/counts_benign_prostate_tpm.txt")

# Tissue Map
fn_tissue_map_wcdt_tissue <- createFn("/data/hmc_secondary_data/tissue_map/tissue_map_wcdt_tissue_103.txt")
fn_tissue_map_wcdt_cfdna <- createFn("/data/hmc_secondary_data/tissue_map/tissue_map_wcdt_cfdna_15.txt")
fn_tissue_map_localized <- createFn("/data/hmc_secondary_data/tissue_map/tissue_map_localized_52.txt")
fn_tissue_map_ubc_cfdna <- createFn("/data/hmc_secondary_data/tissue_map/tissue_map_ubc_cfdna_64.txt")

# NGSplot enrichemnt data
fn_ngs_enrichment <- createFn("/data/hmc_secondary_data/ngs_plot/ngs_plot_enrichment_data.RData")
fn_ngs_chip_enrichment <- createFn("/data/hmc_secondary_data/ngs_plot/NGS.peaks.RData")

# Diffbind output
fn_diffbind_mcrpc_loc <- createFn("/data/hmc_secondary_data/diffbind/differential_peaks_mCRPC_localized.txt")
fn_diffbind_t2e <- createFn("/data/hmc_secondary_data/diffbind/differential_peaks_ERG_fusion.txt")
fn_peaks_up_t2ePos <- createFn("/data/hmc_secondary_data/diffbind/significant_up_peaks_ERG_fusion.bed")

# Homer
fn_ERG_motif_annotation <- createFn("/data/hmc_secondary_data/homer/all_consensus_peaks_homer_ERG_motif_annotation.txt")
fn_ERG_motif_bed <- createFn("/data/hmc_secondary_data/homer/all_consensus_peaks_homer_ERG_motif_annotation.bed")
fn_homer_t2e_up <- createFn("/data/hmc_secondary_data/homer/homer_t2e_up_knownResults.txt")

# GREAT output
fn_great_mcrpc_loc <- createFn("/data/hmc_secondary_data/great/greatExportAll_mcrpc_vs_loc_wo_subtraction_20210626.tsv")

## Annotations
# Gene annotation Gencode v.28
fn_gencodev28_gtf <- createFn("/data/annotation_data/gene_models/gencode.v28.annotation.gtf.gz")

fn_counts_anno <- createFn("/data/annotation_data/gene_models/counts_gencode28_annot.txt")
fn_genes_mapping_gencode_v28 <- createFn("/data/annotation_data/gene_models/gencode.v28.mapping.txt")

# Gene model for tracks plotting
fn_gene_model_hg38 <- createFn("/data/annotation_data/gene_models/hg38_models.txt")
fn_exons_gencode28 <- createFn("/data/annotation_data/gene_models/gencode28_exons.txt")

# Gene lists
fn_genelist_pca_specific_overexpr <- createFn("/data/annotation_data/gene_lists/genelist_pca_specific_overexpr.txt")
fn_genes_zhang <- createFn("/data/annotation_data/gene_lists/zhang_nat_com_2016_fig1i.txt")
fn_genes_beltran <- createFn("/data/annotation_data/gene_lists/41591_2016_BFnm4045_MOESM29_ESM.xlsx")
fn_genes_pc_gi <- createFn("/data/annotation_data/gene_lists/1-s2.0-S1535610817304610-mmc2.xlsx")

## Metadata
# WCDT
fn_metadata_wcdt_genomics <- createFn("/data/sample_metadata/metadata_wcdt_genomics.txt")
fn_pred_ctfrac_wcdt_cfdna <- createFn("/data/sample_metadata/predicted_ctfrac_cfdna_wcdt.txt")

# Localized
fn_metadata_localized <- createFn("/data/sample_metadata/metadata_cpcg.txt")

# UBC
fn_metadat_ubc <- createFn("/data/sample_metadata_NO_UPLOAD/clinical_ubc.txt")

###########################################################################################################
###########################################################################################################
# Samples
samples_mcrpc_uniq_all_data <- c("DTB-003-BL","DTB-005-BL","DTB-008-BL","DTB-011-BL","DTB-018-BL","DTB-019-PRO","DTB-020-BL","DTB-021-BL","DTB-022-BL","DTB-023-BL","DTB-024-PRO","DTB-034-BL", "DTB-036-BL","DTB-037-BL","DTB-040-BL","DTB-042-BL","DTB-053-BL","DTB-055-PRO","DTB-059-BL","DTB-060-BL","DTB-061-BL","DTB-063-BL","DTB-064-BL","DTB-067-PRO","DTB-069-BL","DTB-071-BL","DTB-074-BL","DTB-077-PRO","DTB-080-BL","DTB-083-BL","DTB-085-BL","DTB-089-BL","DTB-090-PRO","DTB-091-BL","DTB-092-BL","DTB-094-BL","DTB-097-PRO","DTB-098-PRO2","DTB-100-BL","DTB-101-BL","DTB-104-BL","DTB-111-PRO","DTB-119-PRO","DTB-121-BL","DTB-124-BL","DTB-126-BL","DTB-127-PRO","DTB-128-BL","DTB-129-BL","DTB-132-BL","DTB-137-PRO","DTB-138-BL","DTB-140-BL","DTB-141-BL","DTB-143-BL","DTB-146-BL","DTB-149-BL","DTB-151-BL","DTB-156-BL","DTB-159-BL","DTB-165-PRO","DTB-167-PRO","DTB-170-BL","DTB-172-BL","DTB-173-BL","DTB-175-BL","DTB-176-BL","DTB-186-BL","DTB-187-BL","DTB-188-BL","DTB-190-BL","DTB-194-PRO","DTB-201-PRO","DTB-202-BL","DTB-204-BL","DTB-205-BL","DTB-206-BL","DTB-210-BL","DTB-213-BL","DTB-214-BL","DTB-216-PRO","DTB-222-BL","DTB-223-BL","DTB-232-PRO","DTB-234-BL","DTB-251-BL","DTB-252-BL","DTB-255-BL","DTB-258-BL","DTB-260-BL","DTB-261-BL","DTB-265-PRO","DTB-266-BL")

samples_wcdt_tissue_nt <- c("NT-DTB-003-BL","NT-DTB-021-BL","NT-DTB-037-BL","NT-DTB-074-BL","NT-DTB-077-PRO","NT-DTB-188-BL","NT-DTB-194-PRO")

samples_wcdt_cfdna <-  c("DTB-119-PRO_cfdna","DTB-127-PRO_cfdna","DTB-140-BL_cfdna","DTB-149-BL_cfdna","DTB-149-PRO_cfdna","DTB-165-PRO_cfdna","DTB-202-BL_cfdna","DTB-214-BL_cfdna","DTB-216-BL_cfdna","DTB-216-PRO_cfdna","DTB-234-BL_cfdna","DTB-252-BL_cfdna","DTB-258-BL_cfdna","DTB-261-BL_cfdna","DTB-265-BL_cfdna")

samples_wcdt_tissue_mcrpc <- unique(c(samples_mcrpc_uniq_all_data,gsub("_cfdna","",samples_wcdt_cfdna))) 

samples_localized_tissue <- c("CPCG0196","CPCG0208","CPCG0210","CPCG0235","CPCG0236","CPCG0238","CPCG0241","CPCG0246","CPCG0248","CPCG0249","CPCG0250","CPCG0256","CPCG0258","CPCG0260" ,"CPCG0263","CPCG0266","CPCG0324","CPCG0331","CPCG0336","CPCG0341","CPCG0342","CPCG0344","CPCG0345","CPCG0346","CPCG0348","CPCG0349","CPCG0350","CPCG0352","CPCG0353","CPCG0354","CPCG0355","CPCG0356","CPCG0357","CPCG0358","CPCG0360","CPCG0361","CPCG0362","CPCG0368","CPCG0369","CPCG0371","CPCG0372","CPCG0373","CPCG0374","CPCG0375","CPCG0377","CPCG0379","CPCG0380","CPCG0382","CPCG0391","CPCG0394","CPCG0458","CPCG0545")

samples_ubc_cfdna <- c("ubc1","ubc10","ubc11","ubc12","ubc13","ubc14","ubc15","ubc17","ubc18","ubc19","ubc2","ubc20","ubc21","ubc22","ubc23","ubc25","ubc26","ubc27","ubc28","ubc29","ubc3","ubc31","ubc32","ubc33","ubc34","ubc35","ubc36","ubc37","ubc38","ubc39","ubc4","ubc40","ubc41","ubc42","ubc43","ubc44","ubc45","ubc46","ubc47","ubc48","ubc5","ubc51","ubc52","ubc53","ubc54","ubc55","ubc56","ubc57","ubc58","ubc59","ubc6","ubc60","ubc61","ubc62","ubc63","ubc64","ubc65","ubc66","ubc67","ubc68","ubc7","ubc70","ubc8","ubc9")


###########################################################################################################
###########################################################################################################
# Load WGS data
###########################################################################################################
###########################################################################################################
# Load data used in Quigley et al. Cell 2018. Copied from /data1/projects/WCDT_WGS_2018/WCDT/2018_04_25_build_image.RData on 2021-12-21.
# Data and scripts available at https://github.com/DavidQuigley/WCDT
load(fn_wgs_environment)

matrix_CNA_gencode28 <- as.matrix(read.table(fn_cn_gencodev28, header = T, stringsAsFactors = F))
colnames(matrix_CNA_gencode28) <- gsub("\\.","-",colnames(matrix_CNA_gencode28))

###########################################################################################################
###########################################################################################################
# Load WGBS data
###########################################################################################################
###########################################################################################################
# Load data used in the Zhao et al. Nat Gen 2020. Created with /data1/projects/WCDT_WGBS_2019/WCDT_WGBS/scripts/2019_05_15_prepare_environment.r on 2021-12-21
# Data and script available on https://github.com/DavidQuigley/WCDT_WGBS
load(fn_wgbs_environment)

# Updated rHMR file
hmrseg <- read.table(fn_hmrseg_upd, stringsAsFactors = F, header = T)
tracks[["HMRsegs"]] <- makeGRangesFromDataFrame(hmrseg)

# Methylation levels per sample in rHMRs
rhmrlevel <- read.delim(fn_hmrlevels,header=T,stringsAsFactors=F,sep='\t',strip.white=T,check.names=F)

# Promoter methylation
promoter_methylation <- read.table(fn_promoter_methylation, header = T, stringsAsFactors = F, check.names = F)

###########################################################################################################
###########################################################################################################
# Load additional external data 
erg = read.delim(fn_chip_erg, header=FALSE, sep='\t', stringsAsFactors=FALSE, strip.white=T)
colnames(erg) = c('chrom','start','end','value')

ergvcap = read.delim(fn_chip_erg_vcap, header=FALSE, sep='\t', stringsAsFactors=FALSE, strip.white=T)
colnames(ergvcap) = c('chrom','start','end','value')

pomerantz <- readRDS(fn_pomerantz)

###########################################################################################################
###########################################################################################################
# Load 5hmC data
###########################################################################################################
###########################################################################################################
# Peaks
peaks_wcdt_tissue <- readRDS(fn_peaks_wcdt_tissue)
peaks_wcdt_cfdna <- readRDS(fn_peaks_wcdt_cfdna)
peaks_localized <- readRDS(fn_peaks_localized)
peaks_ubc_cfdna <- readRDS(fn_peaks_ubc_cfdna)

# Counts
counts_wcdt_tissue_raw <- read.table(fn_counts_wcdt_tissue_raw, header = T, stringsAsFactors = F, sep = "\t", check.names = F)
counts_wcdt_cfdna_raw <- read.table(fn_counts_wcdt_cfdna_raw, header = T, stringsAsFactors = F, sep = "\t", check.names = F)
counts_localized_raw <- read.table(fn_counts_localized_raw, header = T, stringsAsFactors = F, sep = "\t", check.names = F)
counts_ubc_cfdna_raw <- read.table(fn_counts_ubc_cfdna_raw, header = T, stringsAsFactors = F, sep = "\t", check.names = F)
counts_benign_prostate_raw <- read.table(fn_counts_benign_prostate_raw, header = T, stringsAsFactors = F, sep = "\t", check.names = F)

counts_wcdt_tissue_tpm <- read.table(fn_counts_wcdt_tissue_tpm, header = T, stringsAsFactors = F, sep = "\t", check.names = F)
counts_wcdt_cfdna_tpm <- read.table(fn_counts_wcdt_cfdna_tpm, header = T, stringsAsFactors = F, sep = "\t", check.names = F)
counts_localized_tpm <- read.table(fn_counts_localized_tpm, header = T, stringsAsFactors = F, sep = "\t", check.names = F)
counts_ubc_cfdna_tpm <- read.table(fn_counts_ubc_cfdna_tpm, header = T, stringsAsFactors = F, sep = "\t", check.names = F)
counts_benign_prostate_tpm <- read.table(fn_counts_benign_prostate_tpm, header = T, stringsAsFactors = F, sep = "\t", check.names = F)

# Tissue Map
tissue_map_wcdt_tissue <- read.table(fn_tissue_map_wcdt_tissue, header = T, stringsAsFactors = F)
tissue_map_wcdt_cfdna <- read.table(fn_tissue_map_wcdt_cfdna, header = T, stringsAsFactors = F)
tissue_map_localized <- read.table(fn_tissue_map_localized, header = T, stringsAsFactors = F)
tissue_map_ubc <- read.table(fn_tissue_map_ubc_cfdna, stringsAsFactors = F) 

gi_sig <- c("gastric","colon_tran","colon_sig","liver","pancreas") 

# NGSplot enrichment data
load(fn_ngs_enrichment)
load(fn_ngs_chip_enrichment)

# Diffbind consensuspeaks mcrpc vs localized
loc_mcrpc_res <- read.table(fn_diffbind_mcrpc_loc, header = T, stringsAsFactors = F)

# Homer
homer_res_erg <- read.delim(fn_homer_t2e_up, stringsAsFactors = F, header = T, sep = "\t")

# Great analysis 
great_mcrpc_loc <- read.delim(fn_great_mcrpc_loc, stringsAsFactors =  F, sep = "\t", comment.char = "", skip = 3)

###########################################################################################################
########################################################################################################### 
# Load gene model and annotation data
gencode.v28.txdb <- makeTxDbFromGFF(file = fn_gencodev28_gtf) # This makes the annotation transcript based

gencode.v28 <- read.table(fn_genes_mapping_gencode_v28, header = T, stringsAsFactors = F)
gencode.v28.pc <- gencode.v28[gencode.v28$gene_type == "protein_coding",]

hg38_models <- read.table(fn_gene_model_hg38, sep = "\t", header = T, stringsAsFactors = F)

bedgenes = ensembl2sym[,c('chr','start','end','name','type','strand')]
bedgenes$strand = as.numeric(bedgenes$strand=='+')
bedgenes[bedgenes$strand==0,'strand'] = -1
colnames(bedgenes) = c('chrom','start','stop','gene','score','strand')

exons = read.delim(file=fn_exons_gencode28, header=TRUE, stringsAsFactors = F) 
colnames(exons)[1] = 'exon_ID'
bedgene = exons[,c(5,6,7,3,1,4)]
bedgene$strand = as.numeric(bedgene$strand=='+')
bedgene[bedgene$strand==0,'strand'] = -1

# Genelists
genelist_pca_specific_overexpr <- read.table(fn_genelist_pca_specific_overexpr,stringsAsFactors = F, header = T, sep = "\t")

# Msigdb gene sets
msigdb_df <- msigdbr(species = "Homo sapiens")
msigdb_h <- msigdb_df[msigdb_df$gs_cat=="H",]
hallmark_pathways <- msigdb_h %>% split(x = .$gene_symbol, f = .$gs_name)

msigdb_c2 <- msigdb_df[msigdb_df$gs_cat=="C2",]
wp_pathways <- msigdb_c2[msigdb_c2$gs_subcat == "CP:WIKIPATHWAYS",]
wp_pathways <- wp_pathways %>% split(x = .$gene_symbol, f = .$gs_name)

genes_coding <- gencode.v28$gene_name[gencode.v28$gene_type=="protein_coding"]

genes_zhang <- read.table(fn_genes_zhang, header = T, stringsAsFactors = F)
genes_zhang_luminal <- as.vector(genes_zhang$Luminal)
genes_zhang_basal <- as.vector(genes_zhang$Basal)

genes_beltran <- as.data.frame(read_excel(fn_genes_beltran, sheet = 10, skip = 1))
genes_beltran_nepc_up <- genes_beltran$`HGNC ID`[genes_beltran$`RNA (CRPC-NE vs CRPC-Adeno)` == "Over-expressed"]
genes_beltran_nepc_dn <- genes_beltran$`HGNC ID`[genes_beltran$`RNA (CRPC-NE vs CRPC-Adeno)` == "Under-expressed"]

genes_pc_gi <- as.data.frame(read_excel(fn_genes_pc_gi, sheet = 1, skip = 1))
genes_pc_gi <- genes_pc_gi[-nrow(genes_pc_gi),]
genes_pc_gi <- unique(as.character(genes_pc_gi[,"GENE_SYMBOL"])) 

all_pathways <- c(hallmark_pathways,wp_pathways)

all_pathways[["Beltran_NEPC_up"]] <- genes_beltran_nepc_up
all_pathways[["Beltran_NEPC_dn"]] <- genes_beltran_nepc_dn
all_pathways[["Zhang_Luminal"]] <- genes_zhang_luminal
all_pathways[["Zhang_Basal"]] <- genes_zhang_basal
all_pathways[["PCa_GI"]] <- genes_pc_gi

for(i in 1:length(all_pathways)){
  theDF <- data.frame(gs_name = names(all_pathways[i]),
                      gene_name = all_pathways[[i]])
  if(i == 1){
    m_t2g <- theDF
  } else {
    m_t2g <- rbind(m_t2g, theDF)
  } 
}

###########################################################################################################
###########################################################################################################
# Colors
cols_disease_stage <- c("Benign Prostate" = "steelblue3",
                        "Adjacent tissue" = "springgreen3",
                        "Localized PCa" = "goldenrod3",
                        "mCRPC" = "firebrick3")

cols_histology_type <- c("t-SCNC" = "orange2",
                         "Adeno" = "black")

cols_met_site <- c("Bone" = "bisque3",
                   "Liver" = "sienna3",
                   "Lymph_node" = "olivedrab4",
                   "Other" = "grey50")

cols_data_types <- c("darkgreen","darkred","dodgerblue3")
names(cols_data_types) <-c("CNA","PM","HMC")

cols_biopsy_type <- c("cfDNA" = "red1",
                      "tissue" = "tan4")

cols_wo_5hmc <- c("black","dodgerblue3")
cols_false_true <- c("gray","darkred")

cols_true_false <- c("TRUE" = "black",
                     "FALSE" = "white")

cols_5hmc_clusters <- c("1" = "orange2",
                        "2" = "skyblue2",
                        "3" = "darkred")

cols_select_mutations <- c("FALSE" = "white",
                           "BRAF" = "orange",
                           "IDH1" = "green",
                           "TET2" = "purple",
                           "DNMT3B" = "deeppink",
                           "CDH1" = "gold",
                           "SPOP" = "red")

col_neg_pos <- c("darkblue","orangered3")
col_3 <- c("darkblue","goldenrod2","orangered3")

# Color functions for heatmaps etc
cols_blue_red = colorRampPalette(c('blue','white','red'))(n = 1000)
col_fun_RdYlBl = colorRamp2(c(seq(0,1,0.1)), rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(11))) # This keeps 0 and 1 constant, otherwise scaled to min/max value.
col_fun_wb = colorRamp2(c(0,100),c("white","black"))
col_fun_purity = colorRamp2(c(0,100),c("white","red"))
col_fun_wHMC = colorRamp2(c(0,1),c("white",cols_data_types["HMC"]))

##
cols2 <- c("blue","red")
cols3 <- c("black","blue","red")
cols4 <- c("green","blue","orange","red")

###########################################################################################################
###########################################################################################################
# Metadata
###########################################################################################################
###########################################################################################################
# WCDT
###########################################################################################################
metadata_wcdt_genomics <- read.table(fn_metadata_wcdt_genomics, stringsAsFactors = F, header = T, sep = "\t")
pred_ctfrac_wcdt_cfdna <- read.table(fn_pred_ctfrac_wcdt_cfdna)

###########################################################################################################
# Localized
metadata_localized <- read.table(fn_metadata_localized)


###########################################################################################################
# UBC cfDNA
clin_ubc <- read.table(fn_metadat_ubc, stringsAsFactors = F, header = T)

###########################################################################################################
###########################################################################################################
# Output RData object
save.image(file = fn_environment_out)
