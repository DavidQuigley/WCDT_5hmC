###########################################################################################################
###########################################################################################################
# Figure 3,6A,7A
# 2022-05-02
# Martin Sjoestroem
###########################################################################################################
###########################################################################################################
# Directories and filenames
createFn <- function(fn, dirBase = dir_base) {paste0(dirBase, fn)}

dir_base <- "/home/martin/mnt/data1/projects/WCDT_5hmc_2019/publication_reproduce/cancer_research_revision_1"

fn_load_5hmc_environment <- createFn("/data/hmc_220428.RData")
###########################################################################################################
# Load environment
###########################################################################################################
load(fn_load_5hmc_environment)
lapply(libraries_to_load,require,character.only = T)
###########################################################################################################

###########################################################################################################
# Heatmap settings
###########################################################################################################
tissuemap_hm_colors <- col_fun_RdYlBl

row_order_tm_tissue <- c("prostate","colon_tran","colon_sig","gastric","pancreas","liver","lung","skin","breast","bladder","brain","heart","hypothalamus","kidney","ovary","placenta","testis","lymph.node","bone.marrow")
row_order_tm_cfdna <- c("bone.marrow","prostate","colon_tran","colon_sig","gastric","pancreas","liver","lung","skin","breast","bladder","brain","heart","hypothalamus","kidney","ovary","placenta","testis","lymph.node")

FONTSIZE <- 6

###########################################################################################################
###########################################################################################################
# Heatmap plots
###########################################################################################################
###########################################################################################################
# WCDT tissue
# site
site.mcrpc <- metadata_wcdt_genomics[,"site",drop=F]
site.nt <- site.mcrpc
rownames(site.nt) <- paste0("NT-",rownames(site.nt))
anno_df <- rbind(site.mcrpc,site.nt)[c(samples_mcrpc_uniq_all_data, samples_wcdt_tissue_nt),,drop = F]
colnames(anno_df) <- "Site_of_metastasis"
anno_df$Disease_stage <- ifelse(grepl("NT",rownames(anno_df)),"Adjacent tissue","mCRPC")
anno_df$Tumor_purity <- metadata_wcdt_genomics[rownames(anno_df),"purity"]
anno_df$Tumor_subtype <- ifelse(rownames(anno_df) %in% samples_tscnc, "t-SCNC","Adeno")

anno_df <- anno_df[,c("Disease_stage","Site_of_metastasis","Tumor_purity","Tumor_subtype")]

column_ha <- HeatmapAnnotation(df = anno_df,
                               show_annotation_name = T,
                               gap = unit(0.05, "in"),
                               border  = TRUE,
                               gp = gpar(col = "black",lwd = 0.5),
                               col = list(Site_of_metastasis = cols_met_site,
                                          Disease_stage = cols_disease_stage,
                                          Tumor_purity = col_fun_purity,
                                          Tumor_subtype = cols_histology_type),
                               annotation_legend_param = list(Site_of_metastasis = list(title = "Site of metastasis"),
                                                              Disease_stage = list(title = "Disease_stage"),
                                                              Tumor_purity = list(title = "Tumor purity"),
                                                              Tumor_subtype = list(title = "Tumor subtype")))

mat <- t(as.matrix(tissue_map_wcdt_tissue[,!colnames(tissue_map_wcdt_tissue) %in% c("SeqID","SampleID")]))
mat <- mat[,c(samples_mcrpc_uniq_all_data,samples_wcdt_tissue_nt)]

hm1 <- Heatmap(mat, name = "TissueMap score", col = tissuemap_hm_colors, rect_gp = gpar(col = "darkgray",lwd = 0.05),
               column_title = "mCRPC",
               cluster_rows = F, cluster_columns = F,
               top_annotation = column_ha,
               show_row_dend = F,
               row_order = row_order_tm_tissue,
               column_order = order(anno_df$Disease_stage,anno_df$Site_of_metastasis, decreasing = T),
               show_column_names = F,
               column_names_gp = grid::gpar(fontsize = FONTSIZE))

###########################################################################################################
# Localized
mat <- t(as.matrix(tissue_map_localized[,!colnames(tissue_map_localized) %in% c("SeqID","SampleID")]))
mat <- mat[,samples_localized_tissue]

anno_df <- data.frame(sample_id = colnames(mat))
anno_df <- merge(anno_df,metadata_localized)
rownames(anno_df) <- anno_df$sample_id
colnames(anno_df) <- c("sample_id","Tumor_purity")
anno_df <- anno_df[,c("Tumor_purity"),drop = F]
anno_df$Tumor_purity <- anno_df$Tumor_purity*100
anno_df$Disease_stage <- "Localized PCa"
anno_df$Site_of_metastasis <- NA # note that plotting only the localized with only NA in one variable will throw an error. Works when combining in a list with other heatmaps though.
anno_df$Tumor_subtype <- NA

anno_df <- anno_df[,c("Disease_stage","Site_of_metastasis","Tumor_purity","Tumor_subtype")]

column_ha <- HeatmapAnnotation(df = anno_df,
                               show_annotation_name = F,
                               gap = unit(0.05, "in"),
                               border  = TRUE,
                               gp = gpar(col = "black",lwd = 0.5),
                               col = list(Site_of_metastasis = cols_met_site,
                                          Disease_stage = cols_disease_stage,
                                          Tumor_purity = col_fun_purity,
                                          Tumor_subtype = cols_histology_type),
                               annotation_legend_param = list(Site_of_metastasis = list(title = "Metastatic site"),
                                                              Disease_stage = list(title = "Disease stage"),
                                                              Tumor_purity = list(title = "Tumor purity"),
                                                              Tumor_subtype = list(title = "Tumor subtype")))

hm2 <- Heatmap(mat, name = "TissueMap score", col = tissuemap_hm_colors, rect_gp = gpar(col = "darkgray",lwd = 0.05),
               column_title = "Localized prostate cancer",
               top_annotation = column_ha,
               cluster_rows = F, 
               cluster_columns = F,
               show_row_dend = F,
               row_order = row_order_tm_tissue,
               column_order = order(mat["prostate",],decreasing = T),
               show_column_names = F,
               column_names_gp = grid::gpar(fontsize = FONTSIZE))

###########################################################################################################
# Figure 3A
###########################################################################################################
ht_list <- hm2 + hm1
#pdf(file = paste0(dir_out_figures,"/Figure_3A.pdf"),height = 6, width = 13)
draw(ht_list, ht_gap = unit(0.3,"in"))
#dev.off()

###########################################################################################################
# Figure 6A
mat <- t(as.matrix(tissue_map_wcdt_cfdna[,!colnames(tissue_map_wcdt_cfdna) %in% c("SeqID","SampleID","predTissue")]))

ctfrac <- pred_ctfrac_wcdt_cfdna[colnames(mat),,drop=F]
colnames(ctfrac) <- "ct-fraction (5hmC)"

column_ha <- HeatmapAnnotation(df = ctfrac,
                               gap = unit(0.05, "in"),
                               border  = TRUE,
                               gp = gpar(col = "black",lwd = 0.5),
                               show_annotation_name = T,
                               col = list("ct-fraction (5hmC)" = col_fun_purity),
                               annotation_legend_param = list("ct-fraction (5hmC)" = list(title = "ct-fraction (5hmC)")))


hm <- Heatmap(mat, name = "TissueMap score", col = tissuemap_hm_colors, rect_gp = gpar(col = "darkgray",lwd = 0.05), 
              cluster_rows = F, 
              cluster_columns = F,
              show_row_dend = F,
              show_column_names = F,
              row_order = row_order_tm_cfdna,
              top_annotation = column_ha,
              column_order = order(mat["prostate",],decreasing = T))

#pdf(file = paste0(dir_out_figures,"/Figure_6A.pdf"),height = 5, width = 6)
print(hm)
#dev.off()

cor.test(ctfrac[,"ct-fraction (5hmC)"],mat["prostate",], method = "spearman")
cor.test(ctfrac[,"ct-fraction (5hmC)"],mat["bone.marrow",], method = "spearman")

####################### 
# Figure 7A
toPlot <- rownames(clin_ubc) 
mat <- t(as.matrix(tissue_map_ubc[toPlot,!colnames(tissue_map_ubc) %in% 
                                    c("SeqID","SampleID","predTissue","cohortType","seq_id","source_id","Outlier","sample_id")]))

theCTfrac <- clin_ubc[colnames(mat),c("est_ctfrac_targeted","est_ctfrac_5hmc"), drop = F]
theCTfrac$est_ctfrac_targeted <- theCTfrac$est_ctfrac_targeted*100
colnames(theCTfrac) <- c("ct-fraction (DNA-seq)","ct-fraction (5hmC)")

column_ha <- HeatmapAnnotation(df = theCTfrac,
                               show_annotation_name = T,
                               gap = unit(0.05, "in"),
                               border  = TRUE,
                               gp = gpar(col = "black",lwd = 0.5),
                               col = list("ct-fraction (DNA-seq)" = col_fun_purity,
                                          "ct-fraction (5hmC)" = col_fun_purity),
                               annotation_legend_param = list("ct-fraction (DNA-seq)" = list(title = "ct-fraction (DNA-seq)"),
                                                              "ct-fraction (5hmC)" = list(title = "ct-fraction (5hmC)")))

hm <- Heatmap(mat, name = "TissueMap score", col = tissuemap_hm_colors, rect_gp = gpar(col = "darkgray",lwd = 0.05), 
              top_annotation = column_ha,
              cluster_rows = F, 
              cluster_columns = F,
              column_names_gp = grid::gpar(fontsize = 8),
              show_row_dend = F,
              row_order = row_order_tm_cfdna,
              show_column_names = F,
              column_order = order(mat["prostate",],decreasing = T))

#pdf(file = paste0(dir_out_figures,"/Figure_7A.pdf"),height = 5, width = 9)
print(hm)
#dev.off()

rownames(theCTfrac)==colnames(mat)
cor.test(theCTfrac[,"ct-fraction (5hmC)"],theCTfrac[,"ct-fraction (DNA-seq)"], method = "spearman")
cor.test(theCTfrac[,"ct-fraction (DNA-seq)"],mat["prostate",], method = "spearman")
cor.test(theCTfrac[,"ct-fraction (DNA-seq)"],mat["bone.marrow",], method = "spearman")

###########################################################################################################
# Figure 3B-C
###########################################################################################################
theScore <- "prostate"
theJitter <- 0.1
theColors <- c("goldenrod3","firebrick3","orange2")
theColorsDark <- c("goldenrod4","red4","orange4")
theNames <- c("Localized PC","mCRPC, adeno","mCRPC, t-SCNC")

theData <- list(Loc = tissue_map_localized[samples_localized_tissue,theScore],
                adeno = tissue_map_wcdt_tissue[samples_mcrpc_uniq_all_data[!samples_mcrpc_uniq_all_data %in% samples_tscnc],theScore],
                nepc = tissue_map_wcdt_tissue[samples_mcrpc_uniq_all_data[samples_mcrpc_uniq_all_data %in% samples_tscnc],theScore])

p_loc_adeno <- signif(wilcox.test(theData[["Loc"]],theData[["adeno"]])$p.value,2)
p_loc_nepc <- signif(wilcox.test(theData[["Loc"]],theData[["nepc"]])$p.value,2)
p_loc_mcrpc <- signif(wilcox.test(theData[["Loc"]],c(theData[["adeno"]],theData[["nepc"]]))$p.value,2)
p_adeno_nepc <- signif(wilcox.test(theData[["adeno"]],theData[["nepc"]])$p.value,2)

png(filename = paste0(dir_out_figures,"/Figure_3B.png"),height = 7, width = 5, units = "in", res = 300)
par(mar = c(7,4,2,2) + 0.2)
boxplot(theData[[1]],
        theData[[2]],
        theData[[3]],
        names = theNames,
        ylab = "5hmC Tissue Map prostate score",
        xlab = NULL,
        xaxt = "n",
        col = theColors,
        outline = F,
        ylim = c(), # To kepp outlier points in the plot.
        main = "",
        #cex.main = 2,
        cex.lab = 1.5,
        las = 2)
text(x = c(1,2,3)-0.3, y = -0.18, labels = theNames, xpd = T, srt = 45, cex = 1.3)
text(x = c(3), y = 0.6, labels = paste0("Loc - mCRPC\n p = ",p_loc_mcrpc,"\nAdeno - t-SCNC\n p = ",p_adeno_nepc))

points(x = jitter(rep(1, length(theData[[1]])), amount = theJitter), y = theData[[1]], pch = 23, col = "black", bg = theColorsDark[1])
points(x = jitter(rep(2, length(theData[[2]])), amount = theJitter), y = theData[[2]], pch = 23, col = "black", bg = theColorsDark[2])
points(x = jitter(rep(3, length(theData[[3]])), amount = theJitter), y = theData[[3]], pch = 23, col = "black", bg = theColorsDark[3])
dev.off()

##################
gi_tissue <- gi_sig
theScore <- gi_tissue
theJitter <- 0.1

theColors <- c("goldenrod3","firebrick3")
theColorsDark <- c("goldenrod4","red4")
theNames <- c("Localized PC","mCRPC")

theData <- list(Loc = rowSums(tissue_map_localized[samples_localized_tissue,theScore]),
                mcrpc = rowSums(tissue_map_wcdt_tissue[samples_mcrpc_uniq_all_data,theScore]))

p_loc_mcrpc <- signif(wilcox.test(theData[["Loc"]],theData[["mcrpc"]])$p.value,2)

png(filename = paste0(dir_out_figures,"/Figure_3C.png"),height = 7, width = 4, units = "in", res = 300)
par(mar = c(7,6,2,2) + 0.2)
boxplot(theData[[1]],
        theData[[2]],
        names = theNames,
        ylab = "5hmC Tissue Map GI score",
        xaxt = "n",
        col = theColors,
        outline = F,
        ylim = c(0,1), # To kepp outlier points in the plot.
        main = "",
        #cex.main = 2,
        cex.lab = 1.5,
        las = 2)
abline(h = 0.25, lty = 2, lwd = 3)
text(x = c(1,2)-0.3, y = -0.18, labels = theNames, xpd = T, srt = 45, cex = 1.3)
text(x = c(1), y = 0.6, labels = paste0("p = ",p_loc_mcrpc))

points(x = jitter(rep(1, length(theData[[1]])), amount = theJitter), y = theData[[1]], pch = 23, col = "black", bg = theColorsDark[1])
points(x = jitter(rep(2, length(theData[[2]])), amount = theJitter), y = theData[[2]], pch = 23, col = "black", bg = theColorsDark[2])
dev.off()

table(theData$mcrpc > 0.25)
32/93 # 34% 
table(theData$Loc > 0.25)
2/52 # 4%


###########################################################################################################
# Figure 3D
###########################################################################################################
tm <- t(as.matrix(tissue_map_wcdt_tissue[,!colnames(tissue_map_wcdt_tissue) %in% c("SeqID","SampleID")]))
tm <- tm[,c(samples_mcrpc_uniq_all_data)]
gex <- log2(as.matrix(tpm[,samples_mcrpc_uniq_all_data])+1)

stopifnot(all(colnames(tm)==colnames(gex)))
stopifnot(all(rownames(gex)==gencode.v28$gene_ensembl))

gex <- gex[gencode.v28$gene_type == "protein_coding",]
stopifnot(all(rownames(gex) == gencode.v28.pc$gene_ensembl))

theDF <- data.frame(matrix(NA,nrow = nrow(gex), ncol = nrow(tm)))
colnames(theDF) <- rownames(tm)
theDF$gene <- gencode.v28.pc$gene_name
theDF$ens <- gencode.v28.pc$gene_ensembl
rownames(theDF) <- theDF$ens

for(i in 1:nrow(theDF)){
  if(i %% 100 == 0){
    print(i)
  }
  
  theEns <- theDF[i,"ens"]
  
  for(j in 1:nrow(tm)){
    theTissue <- rownames(tm)[j]
    theC <- cor(tm[theTissue,],gex[theEns,], use = "na.or.complete", method = "spearman")
    theDF[theEns,theTissue] <- theC
  }
}

# Add GI score
theSigs <- c("pancreas","gastric","liver","colon_sig","colon_tran")
gi_score <- colSums(tm[theSigs,])

theDF$gi_score <- NA
all(names(gi_score)==colnames(gex))
for(i in 1:nrow(theDF)){
  if(i %% 100 == 0){
    print(i)
  }
  
  theEns <- rownames(theDF)[i]
  theC <- cor(gi_score,gex[theEns,],method = "spearman",use = "na.or.complete")
  theDF[theEns,"gi_score"] <- theC
}

genes_rank <- theDF[,"gi_score"]
names(genes_rank) <-theDF[,"gene"]  
genes_rank <- genes_rank[!is.na(genes_rank)]
genes_rank <- genes_rank[!duplicated(genes_rank)]

set.seed(12)
fgseaRes <- fgseaMultilevel(pathways = all_pathways, stats = genes_rank, minSize = 5, maxSize = 500, nproc = 6)

pval <- signif(as.numeric(fgseaRes[fgseaRes$pathway == "PCa_GI","padj"]),2)

p <- plotEnrichment(pathway = all_pathways[["PCa_GI"]], stats = genes_rank) +
  labs(y = "Enrichment Score", x = "Gene Rank", title = "PCa_GI") + 
  annotate("text", x = 13000, y = 0.4, label  = paste0("P (adj.) = ",pval, "\nNES = ", round(fgseaRes[fgseaRes$pathway == "PCa_GI","NES"],2)), size = 6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

p

pdf(paste0(dir_out_figures, "/Figure_3D.pdf"), width = 6, height = 4)
print(p)
dev.off()