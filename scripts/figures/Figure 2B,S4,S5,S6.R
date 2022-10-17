###########################################################################################################
###########################################################################################################
# Figure 2B,S4,S5,S6, differential 5hmC analysis
# 2022-05-03
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
# Functions
plotGSEA <- function(theDF, colStat = "stat", colName = "gene_name", pathways = all_pathways, theSeed = 12, topUp = 10, topDn = 10,
                     Title = "", fdrCut = 0.05, savePlot = FALSE, fileName = NULL, plotH = 7, plotW = 7, nProc = 6){
        
        set.seed(12)
        
        # Get ranked genes
        genes_rank <- theDF[,colStat]
        names(genes_rank) <-theDF[,colName]  
        
        # Remove dups and NAs
        genes_rank <- genes_rank[!is.na(genes_rank)]
        genes_rank <- genes_rank[!duplicated(genes_rank)]
        
        # Perform GSEA
        fgseaRes <- fgseaMultilevel(pathways = all_pathways, stats = genes_rank, minSize = 5, maxSize = 500, nproc = nProc)
        
        topPathwaysUp <- fgseaRes[ES > 0][head(order(NES, decreasing = T), n=topUp), pathway]
        topPathwaysDown <- fgseaRes[ES < 0][head(order(NES), n=topDn), pathway]
        topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
        
        # Subset and order for pathways
        plotDF <- as.data.frame(fgseaRes)
        rownames(plotDF) <- plotDF$pathway
        plotDF <- plotDF[topPathways,]
        
        # Plot results
        thePlot <- ggplot(plotDF, aes(reorder(pathway, NES), NES)) +
                geom_col(aes(fill = padj < fdrCut)) +
                coord_flip() +
                labs(x = "Pathway", y = "Normalized Enrichment Score",
                     title = Title,fill= paste0("FDR < ",fdrCut)) +
                theme_few() + 
                scale_fill_manual(values = c("slategray3","forestgreen"),
                                  breaks = c("FALSE","TRUE"))
        
        print(thePlot)
        
        if(savePlot){
                pdf(file = fileName, width = plotW, height = plotH)
                print(thePlot)
                dev.off()
        }
}


plotVolcano <- function(theDF, padjCutUp = 40, padjCutDn = 40, l2fcCutUp = 2, l2fcCutDn = 2, Use = "or",
                        Title = "", savePlot = FALSE, fileName = NULL, plotH = 7, plotW = 7){
        
        if(Use == "or"){
                toLabel <- (-log10(theDF$padj) > padjCutUp & theDF$log2FoldChange > 0  & !is.na(theDF$padj)) |
                        (-log10(theDF$padj) > padjCutDn & theDF$log2FoldChange < 0  & !is.na(theDF$padj)) |
                        (theDF$log2FoldChange > l2fcCutUp & !is.na(theDF$log2FoldChange)) |
                        (theDF$log2FoldChange < -l2fcCutDn & !is.na(theDF$log2FoldChange))
        } else if(Use == "and"){
                toLabel <- ((-log10(theDF$padj) > padjCutUp & theDF$log2FoldChange > 0  & !is.na(theDF$padj)) & 
                                    (theDF$log2FoldChange > l2fcCutUp & !is.na(theDF$log2FoldChange)))    |
                        ((-log10(theDF$padj) > padjCutDn & theDF$log2FoldChange < 0  & !is.na(theDF$padj)) &
                                 (theDF$log2FoldChange < -l2fcCutDn & !is.na(theDF$log2FoldChange)))
        }else{
                stop("Use must be and/or!")
        }
        
        
        
        theDF$theLabel <- NA
        theDF$theLabel[toLabel] <- theDF$gene_name[toLabel]
        
        p <- ggplot(theDF, aes(x = log2FoldChange, y = -log10(padj), label = theLabel)) + 
                geom_point(color = ifelse(toLabel, "red","black"), size = 3, 
                           alpha = ifelse(toLabel, 1, 0.2)) +
                geom_label_repel(box.padding = 0.5, max.overlaps = 100) + 
                labs(x = "log2FC", y = "-log10(padj)", title = Title) +
                theme_classic() +
                theme(axis.text = element_text(size = 14),
                      axis.title = element_text(size = 16))
        
        print(p)
        
        if(savePlot){
                pdf(file = fileName, width = plotW, height = plotH)
                print(p)
                dev.off()
        } 
}


###########################################################################################################
# Set up data
###########################################################################################################
# 5hmC data
wcdt_hmc_counts_raw <- counts_wcdt_tissue_raw[,samples_mcrpc_uniq_all_data]
wcdt_hmc_counts <- wcdt_hmc_counts_raw[gencode.v28.pc$gene_ensembl,]

all_hmc_counts_raw <- cbind(counts_wcdt_tissue_raw, counts_localized_raw, counts_benign_prostate_raw)
colnames(all_hmc_counts_raw) <- gsub("_5hmc","",colnames(all_hmc_counts_raw))
all_hmc_counts_raw <- all_hmc_counts_raw[,c(samples_mcrpc_uniq_all_data,samples_localized_tissue,"BSGD54","BSGD55","BSGD56","BSGD57","BSGD58")]
all_hmc_counts <- all_hmc_counts_raw[gencode.v28.pc$gene_ensembl,]

# Annotation data
wcdt_annotation_data <- as.data.frame(metadata_wcdt_genomics[samples_mcrpc_uniq_all_data,])
wcdt_annotation_data$tSCNC <- as.factor(wcdt_annotation_data$tSCNC)
wcdt_annotation_data$site <- as.factor(wcdt_annotation_data$site)
stopifnot(all(rownames(wcdt_annotation_data)==colnames(wcdt_hmc_counts)))


all_annotation_data <- data.frame(sample = colnames(all_hmc_counts),
                                  stage = c(rep("mCRPC",length(samples_mcrpc_uniq_all_data)),rep("localized_PCa",length(samples_localized_tissue)),
                                            rep("benign",5)), 
                                  stringsAsFactors = TRUE)

rownames(all_annotation_data) <- all_annotation_data$sample
stopifnot(all(rownames(all_annotation_data)==colnames(all_hmc_counts)))

###########################################################################################################
# Perform diff analyses
###########################################################################################################
## tSCNC vs. adeno
cnts <- wcdt_hmc_counts
colData <- wcdt_annotation_data

dds_tscnc <- DESeqDataSetFromMatrix(countData = cnts, colData = colData, design = ~ site + tSCNC)
dds_tscnc <- DESeq(dds_tscnc)
res_tscnc <- results(dds_tscnc)

resDF_tscnc <- as.data.frame(res_tscnc)
resDF_tscnc$ens <- rownames(resDF_tscnc)
resDF_tscnc <- merge(resDF_tscnc, gencode.v28.pc, by.x = "ens", by.y = "gene_ensembl")

## mCRPC vs. loc and loc. vs benign
cnts <- all_hmc_counts
colData <- all_annotation_data

dds_all <- DESeqDataSetFromMatrix(countData = cnts, colData = colData, design = ~ stage)
dds_all <- DESeq(dds_all)

## mCRPC vs. loc
res_mcrpc_loc <- results(dds_all, contrast = c("stage","mCRPC","localized_PCa"))
resDF_mcrpc_loc <- as.data.frame(res_mcrpc_loc)
resDF_mcrpc_loc$gene_ensembl <- rownames(resDF_mcrpc_loc)
resDF_mcrpc_loc <- merge(resDF_mcrpc_loc, gencode.v28.pc, by = "gene_ensembl")

## loc. vs benign
res_loc_ben <- results(dds_all, contrast = c("stage","localized_PCa","benign"))
resDF_loc_ben <- as.data.frame(res_loc_ben)
resDF_loc_ben$gene_ensembl <- rownames(resDF_loc_ben)
resDF_loc_ben <- merge(resDF_loc_ben, gencode.v28.pc, by = "gene_ensembl")


#### Tissue site in mCRPC vs loc
### Bone
tissue <- "Bone"
not_tissue <- rownames(wcdt_annotation_data)[!wcdt_annotation_data$site %in% tissue]

cnts <- all_hmc_counts[,!colnames(all_hmc_counts) %in% not_tissue]
colData <- all_annotation_data[!rownames(all_annotation_data) %in% not_tissue,]

dds_bone <- DESeqDataSetFromMatrix(countData = cnts, colData = colData, design = ~ stage)
dds_bone <- DESeq(dds_bone)

res_bone <- results(dds_bone, contrast = c("stage","mCRPC","localized_PCa"))
resDF_bone <- as.data.frame(res_bone)
resDF_bone$ens <- rownames(resDF_bone)
resDF_bone <- merge(resDF_bone, gencode.v28.pc, by.x = "ens", by.y = "gene_ensembl")

### Lymph node
tissue <- "Lymph_node"
not_tissue <- rownames(wcdt_annotation_data)[!wcdt_annotation_data$site %in% tissue]

cnts <- all_hmc_counts[,!colnames(all_hmc_counts) %in% not_tissue]
colData <- all_annotation_data[!rownames(all_annotation_data) %in% not_tissue,]

dds_ln <- DESeqDataSetFromMatrix(countData = cnts, colData = colData, design = ~ stage)
dds_ln <- DESeq(dds_ln)

res_ln <- results(dds_ln, contrast = c("stage","mCRPC","localized_PCa"))
resDF_ln <- as.data.frame(res_ln)
resDF_ln$ens <- rownames(resDF_ln)
resDF_ln <- merge(resDF_ln, gencode.v28.pc, by.x = "ens", by.y = "gene_ensembl")

### Liver
tissue <- "Liver"
not_tissue <- rownames(wcdt_annotation_data)[!wcdt_annotation_data$site %in% tissue]

cnts <- all_hmc_counts[,!colnames(all_hmc_counts) %in% not_tissue]
colData <- all_annotation_data[!rownames(all_annotation_data) %in% not_tissue,]

dds_liver <- DESeqDataSetFromMatrix(countData = cnts, colData = colData, design = ~ stage)
dds_liver <- DESeq(dds_liver)

res_liver <- results(dds_liver, contrast = c("stage","mCRPC","localized_PCa"))
resDF_liver <- as.data.frame(res_liver)
resDF_liver$ens <- rownames(resDF_liver)
resDF_liver <- merge(resDF_liver, gencode.v28.pc, by.x = "ens", by.y = "gene_ensembl")

### Other
tissue <- "Other"
not_tissue <- rownames(wcdt_annotation_data)[!wcdt_annotation_data$site %in% tissue]

cnts <- all_hmc_counts[,!colnames(all_hmc_counts) %in% not_tissue]
colData <- all_annotation_data[!rownames(all_annotation_data) %in% not_tissue,]

dds_other <- DESeqDataSetFromMatrix(countData = cnts, colData = colData, design = ~ stage)
dds_other <- DESeq(dds_other)

res_other <- results(dds_other, contrast = c("stage","mCRPC","localized_PCa"))
resDF_other <- as.data.frame(res_other)
resDF_other$ens <- rownames(resDF_other)
resDF_other <- merge(resDF_other, gencode.v28.pc, by.x = "ens", by.y = "gene_ensembl")

###########################################################################################################
# Visualize results
###########################################################################################################
# Plot combined per stage, Figure 2B
###########################################################################################################
set.seed(12)
genes_rank <- resDF_loc_ben$stat
names(genes_rank) <-resDF_loc_ben$gene_name  
genes_rank <- genes_rank[!is.na(genes_rank)]
genes_rank <- genes_rank[!duplicated(genes_rank)]

fgseaRes <- fgseaMultilevel(pathways = all_pathways, stats = genes_rank, maxSize = 500, nproc = 6)
theRes_loc_ben <- fgseaRes

# mCRPC Localized 
set.seed(12)
genes_rank <- resDF_mcrpc_loc$stat
names(genes_rank) <-resDF_mcrpc_loc$gene_name  
genes_rank <- genes_rank[!is.na(genes_rank)]
genes_rank <- genes_rank[!duplicated(genes_rank)]

fgseaRes <- fgseaMultilevel(pathways = all_pathways, stats = genes_rank, maxSize = 500, nproc = 6)
theRes_mcrpc_loc <-fgseaRes

theRes_mcrpc_loc_df <- as.data.frame(theRes_mcrpc_loc)
theRes_loc_ben_df <- as.data.frame(theRes_loc_ben)

rownames(theRes_mcrpc_loc_df) <- theRes_mcrpc_loc_df$pathway
rownames(theRes_loc_ben_df) <- theRes_loc_ben_df$pathway

# Select pathways to plot
all(theRes_loc_ben$pathway == theRes_mcrpc_loc$pathway)

significant_pathways <- theRes_mcrpc_loc$pathway[theRes_mcrpc_loc$padj < 0.0001 | theRes_loc_ben$padj < 0.0001]

theDF <- data.frame(row.names = significant_pathways,
                    mCRPC = theRes_mcrpc_loc_df[significant_pathways,"NES"],
                    loc = theRes_loc_ben_df[significant_pathways,"NES"])

thePs <- data.frame(row.names = significant_pathways,
                    mCRPC = theRes_mcrpc_loc_df[significant_pathways,"padj"],
                    loc = theRes_loc_ben_df[significant_pathways,"padj"])

stopifnot(all(rownames(theDF)==rownames(thePs)))

rownames(theDF)[rownames(theDF)=="WP_SUDDEN_INFANT_DEATH_SYNDROME_SIDS_SUSCEPTIBILITY_PATHWAYS"] <- "WP_SIDS_SUSCEPTIBILITY_PATHWAYS" # cannot be plotted otherwise

theMat = as.matrix(theDF)

col_fun <- colorRamp2(c(-2.5,-1,0,1,2.5),c("blue","white","white","white","red"))

hm <- Heatmap(theMat,
              col = col_fun,
              rect_gp = gpar(col = "black", lwd = 1),
              row_order = order(theMat[,"mCRPC"],decreasing = T),
              column_order = c(1,2),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              column_split = factor(c("mCRPC vs. Loc","Loc vs. benign"), levels = c("mCRPC vs. Loc","Loc vs. benign")),
              cluster_column_slices = FALSE,
              show_column_names = FALSE,
              heatmap_legend_param = list(at = c(-3, 0, 3),
                                          title = "NES"),
              heatmap_width = unit(6, "in"),
              heatmap_height = unit(6,"in"),
              cell_fun = function(j,i,x,y,width,height, fill) {
                      grid.text(signif(thePs[i,j],2), x, y, gp = gpar(fontsize = 14))
                      }
)

pdf(paste0(dir_out_figures,"/Figure_2B.pdf"), height = 6.5, width = 10)
draw(hm, heatmap_legend_side = "left")
dev.off()

###########################################################################################################
# mCRPC vs loc
## All 
# Volcano
plotVolcano(theDF = resDF_mcrpc_loc, padjCutUp = 25, padjCutDn = 25,  l2fcCutUp = 0.75, l2fcCutDn = 0.75, Use = "and",
            Title = "mCRPC vs. Localized PCa", savePlot = TRUE,  
            fileName = paste0(dir_out_figures,"/Figure_S4A.pdf"), 
            plotH = 7, plotW = 7)
# GSEA
plotGSEA(theDF = resDF_mcrpc_loc, Title = "mCRPC vs. Localized PCa", topUp = 15, topDn = 5,
         savePlot = TRUE,  fileName = paste0(dir_out_figures,"/Figure_S4B.pdf"), plotH = 5, plotW = 10)


# Loc vs. benign
# Volcano
plotVolcano(theDF = resDF_loc_ben, padjCutUp = 1, padjCutDn = 4,  l2fcCutUp = 1.5, l2fcCutDn = 1.7, Use = "or",
            Title = "Localized PCa vs. benign prostate", savePlot = TRUE,  
            fileName = paste0(dir_out_figures,"/Figure_S4C.pdf"), 
            plotH = 7, plotW = 7)
# GSEA
plotGSEA(theDF = resDF_loc_ben, Title = "Localized PCa vs. benign prostate", topUp = 15, topDn = 5,
         savePlot = TRUE,  fileName = paste0(dir_out_figures,"/Figure_S4D.pdf"), plotH = 5, plotW = 10)


## Bone
# Volcano
plotVolcano(theDF = resDF_bone, padjCutUp = 25, padjCutDn = 25,  l2fcCutUp = 0.75, l2fcCutDn = 0.75, Use = "and",
            Title = "mCRPC vs. Localized PCa, Bone only", savePlot = TRUE,  
            fileName = paste0(dir_out_figures,"/Figure_S5A.pdf"), 
            plotH = 7, plotW = 7)
# GSEA
plotGSEA(theDF = resDF_bone, Title = "mCRPC vs. Localized PCa, Bone only", topUp = 15, topDn = 5,
         savePlot = TRUE,  fileName = paste0(dir_out_figures,"/Figure_S5C.pdf"), plotH = 5, plotW = 12)

## Lymph node
# Volcano
plotVolcano(theDF = resDF_ln, padjCutUp = 25, padjCutDn = 25,  l2fcCutUp = 0.75, l2fcCutDn = 0.75, Use = "and",
            Title = "mCRPC vs. Localized PCa, LN only", savePlot = TRUE,  
            fileName = paste0(dir_out_figures,"/Figure_S5B.pdf"), 
            plotH = 7, plotW = 7)
# GSEA
plotGSEA(theDF = resDF_ln, Title = "mCRPC vs. Localized PCa, LN only", topUp = 15, topDn = 5,
         savePlot = TRUE,  fileName = paste0(dir_out_figures,"/Figure_S5D.pdf"), plotH = 5, plotW = 12)

## Liver
# Volcano
plotVolcano(theDF = resDF_liver, padjCutUp = 20, padjCutDn = 20,  l2fcCutUp = 1.5, l2fcCutDn = 1.5, Use = "or",
            Title = "mCRPC vs. Localized PCa, Liver only", savePlot = TRUE,  
            fileName = paste0(dir_out_figures,"/Figure_S5E.pdf"), 
            plotH = 7, plotW = 7)
# GSEA
plotGSEA(theDF = resDF_liver, Title = "mCRPC vs. Localized PCa, Liver only", topUp = 15, topDn = 5,
         savePlot = TRUE,  fileName = paste0(dir_out_figures,"/Figure_S5G.pdf"), plotH = 5, plotW = 12)

## Other
# Volcano
plotVolcano(theDF = resDF_other, padjCutUp = 25, padjCutDn = 25,  l2fcCutUp = 1.5, l2fcCutDn = 2, Use = "or",
            Title = "mCRPC vs. Localized PCa, Other only", savePlot = TRUE,  
            fileName = paste0(dir_out_figures,"/Figure_S5F.pdf"), 
            plotH = 7, plotW = 7)
# GSEA
plotGSEA(theDF = resDF_other, Title = "mCRPC vs. Localized PCa, Other only", topUp = 15, topDn = 5,
         savePlot = TRUE,  fileName = paste0(dir_out_figures,"/Figure_S5H.pdf"), plotH = 5, plotW = 12)

# tSCNC
# GSEA
plotGSEA(theDF = resDF_tscnc, Title = "t-SCNC vs. adenocarcinoma", 
         savePlot = FALSE,  fileName = paste0(dir_out_figures,"/Figure_S6.pdf"), plotH = 5, plotW = 10)

set.seed(12)
genes_rank <- resDF_tscnc$stat
names(genes_rank) <-resDF_tscnc$gene_name  
genes_rank <- genes_rank[!is.na(genes_rank)]
genes_rank <- genes_rank[!duplicated(genes_rank)]

theRes_tscnc <- fgseaMultilevel(pathways = all_pathways, stats = genes_rank, maxSize = 500, nproc = 6)
theRes_tscnc[pathway=="HALLMARK_ANDROGEN_RESPONSE",] # padj 1.462359e-07
theRes_tscnc[pathway=="Beltran_NEPC_up",] # padj 4.378389e-06

