###########################################################################################################
# Figure 6B-F
# 2022-01-06
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
plotScatter <- function(var1,var2,xLab = "xlab", yLab = "ylab", Main = "main", pXpos = 1, pYpos = 1, CEX = 1.5, plotLine = T){
  C <- signif(cor(var1,var2, use = "complete.obs", method = "spearman"),2)
  P <- signif(cor.test(var1,var2,use = "complete.obs", method = "spearman")$p.value,2)
  
  fit <-  lm(var2~var1)
  A <- fit$coefficients[1]
  B <- fit$coefficients[2]
  
  plot(var1,var2,xlab = xLab, ylab = yLab, main = Main,
       pch = 16, cex = CEX, col = "darkblue")
  text(x = pXpos, y = pYpos, labels = paste0("Rho = ",C,"\nP = ",P))
  if(plotLine){
    abline(a=A,b=B,lty = 2, lwd = 2, col = 2)
  }
}

###########################################################################################################
#Calculate correlation tisuse-cfdna
tpm_cfdna <- as.matrix(counts_wcdt_cfdna_tpm[,samples_wcdt_cfdna])
tpm_cfdna_full_name <- tpm_cfdna
colnames(tpm_cfdna) <- gsub("_cfdna","",colnames(tpm_cfdna))
tpm_tissue <- as.matrix(counts_wcdt_tissue_tpm[,gsub("_cfdna","",samples_wcdt_cfdna)])

all(rownames(tpm_tissue)==rownames(tpm_cfdna))
all(colnames(tpm_tissue)==colnames(tpm_cfdna))

samples_use <- colnames(tpm_tissue)

correlation_mat <- matrix(NA, nrow = nrow(tpm_tissue), ncol = 2)
rownames(correlation_mat) <- rownames(tpm_tissue)
colnames(correlation_mat) <- c("C","P")

for(i in 1:nrow(correlation_mat)){

  if(i %% 1000 == 0){
    print(i)
  }
  gene <- rownames(correlation_mat)[i]

  C <- cor(tpm_tissue[gene,samples_use],tpm_cfdna[gene,samples_use],method = "spearman", use = "complete.obs")
  P <- cor.test(tpm_tissue[gene,samples_use],tpm_cfdna[gene,samples_use],method = "spearman", use = "complete.obs")$p.value

  correlation_mat[gene,"C"] <- C
  correlation_mat[gene,"P"] <- P

}

correlation_mat <- as.data.frame(correlation_mat)
correlation_mat$gene_ensembl <- rownames(correlation_mat)
all(correlation_mat$gene_ensembl == gencode.v28$gene_ensembl)
correlation_mat_full <- merge(correlation_mat, gencode.v28, by = "gene_ensembl")

###########################################################################################################
# Figure 6B
###########################################################################################################
theSig <- gi_sig
samps <- samples_use

tissuemap_gi_score_wcdt_tissue <- rowSums(tissue_map_wcdt_tissue[samps,theSig,drop=F])   
tissuemap_gi_score_wcdt_cfdna <- rowSums(tissue_map_wcdt_cfdna[samps,theSig,drop=F])   

pdf(paste0(dir_out_figures,"/Figure_6B.pdf"), height = 5, width = 5)
plotScatter(var1 = tissuemap_gi_score_wcdt_tissue[samps], var2 = tissuemap_gi_score_wcdt_cfdna[samps], 
            xLab = "TissueMap GI-score in tissue", yLab = "TissueMap GI-score in cfDNA", Main = "",pXpos = 0.5, pYpos = 0.1, plotLine = F)
dev.off()

###########################################################################################################
# Figure 6C
###########################################################################################################
theDF <- as.data.frame(correlation_mat_full[correlation_mat_full$gene_type == "protein_coding",])
theCol <- cols_data_types["HMC"]

p <- ggplot(theDF, aes(x = C)) +
  geom_density(alpha = 0.4, fill = theCol) +
  geom_vline(xintercept = 0, color = "black", alpha = 0.5,linetype = "dashed", size = 2) +
  geom_vline(xintercept = median(theDF$C,na.rm=T), color = theCol, alpha = 0.5, linetype = "dashed", size = 2) +
  labs(title = "Correlation 5hmC in tissue and cfDNA for protein coding genes",
       y = "Density", x = "Spearman's Rho")+ 
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 16),
        axis.line.x = element_line(),
        axis.line.y = element_line()) +
  theme(legend.position = "none")

pdf(paste0(dir_out_figures,"/Figure_6C.pdf"), width = 8, height = 5)
print(p)
dev.off()


###########################################################################################################
# Format data
theData <- cbind(counts_wcdt_tissue_tpm,counts_wcdt_cfdna_tpm)

theSampleData <- data.frame(original_sample = colnames(theData),stringsAsFactors = F)
rownames(theSampleData) <- theSampleData$original_sample
theSampleData$base_sample <- theSampleData$original_sample
theSampleData$base_sample <- gsub("_cfdna$","",theSampleData$base_sample)
theSampleData$patient <- gsub("-.{2,4}$", "",theSampleData$base_sample)
theSampleData$is.cfdna <- grepl("cfdna",theSampleData$original_sample)
theSampleData$matched <- theSampleData$base_sample %in% gsub("_cfdna","",samples_wcdt_cfdna)
sampleData <- theSampleData[theSampleData$matched==T,]

theData_matched <- theData[,sampleData$original_sample]
ensembl_protein_coding <- gencode.v28$gene_ensembl[gencode.v28$gene_type == "protein_coding"]
theDF_matched <- theData_matched[ensembl_protein_coding,]
theDF_matched <- log2(theDF_matched+1)

###########################################################################################################
# Select top varying genes, scale per cohort
###########################################################################################################
cutOff <- 0.9

theVar <- apply(theDF_matched,1,var)
topVar <- theVar > quantile(theVar,cutOff)
theDF_matched <- t(theDF_matched[topVar,])

all(rownames(theDF_matched)==rownames(sampleData))

theDF_scale <- rbind(scale(theDF_matched[sampleData$is.cfdna == T,]),
                     scale(theDF_matched[sampleData$is.cfdna == F,]))

sampleData <- sampleData[rownames(theDF_scale),]
all(rownames(theDF_scale)==rownames(sampleData))

# Remove genes with NAs 
theDF_scale <- theDF_scale[,!apply(theDF_scale,2,function(x) sum(!is.finite(x))) > 0]

###########################################################################################################
# Figure 6D
###########################################################################################################
mat <- as.matrix(t(theDF_scale)) 
res <- cor(mat,method = "spearman") 

basenames <- unique(sampleData$base_sample)

corrs <- matrix(NA,ncol = 2, nrow = length(basenames))
rownames(corrs) <- basenames
colnames(corrs) <- c("C_self","C_others")

for(i in 1:length(basenames)){
  
  theBasename <- basenames[i]
  tissueName <- paste0(theBasename)
  cfdnaName <- paste0(theBasename,"_cfdna")
  
  thisCorr <- res[tissueName,cfdnaName]
  
  cfdna_names <- colnames(res)[grepl("cfdna",colnames(res))]
  other_cfdna_names <- cfdna_names[!cfdna_names==cfdnaName]
  
  avOtherCorr <- mean(res[tissueName,other_cfdna_names])
  
  corrs[theBasename,"C_self"] <- thisCorr
  corrs[theBasename,"C_others"] <- avOtherCorr
}

theDF <- data.frame(sample = rep(rownames(corrs),2),
                    group = c(rep("Matched",15),rep("Others",15)),
                    corr = c(corrs[,1],corrs[,2]),
                    stringsAsFactors = T)

theDF$group <- factor(theDF$group,levels = c("Others","Matched"))

p <- ggplot(theDF, aes(x = group, y = corr, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(width = 0.1, alpha = 0.6, size = 3) +
  labs(x="", y = "Spearman's Rho", title = "") +
  scale_fill_manual(values = c("red3","steelblue3")) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 10, hjust = 0),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent")) + 
  theme(legend.position = "none")

p

pdf(paste0(dir_out_figures,"/Figure_6D.pdf"),height = 5,width = 3)
print(p)
dev.off()

###########################################################################################################
# Figure 6E
###########################################################################################################
theDF_sub <- theDF[theDF$group=="Matched",] 
rownames(theDF_sub) <- theDF_sub$sample
samps2use <- gsub("_cfdna","",samples_wcdt_cfdna)

pdf(paste0(dir_out_figures,"/Figure_6E.pdf"),height = 5,width = 5)
plotScatter(var1 = pred_ctfrac_wcdt_cfdna[samps2use,"ctfrac_by_5hmc_pred"], var2 = theDF_sub[samps2use,"corr"],
            xLab = "ct-fraction (5hmC)", yLab = "Sample similarity (Spearman's correlation)", Main = "", 
            pXpos = 40, pYpos = 0.1)
dev.off()

###########################################################################################################
# Figure 6F
###########################################################################################################
method.distance <- "euclidean" 
method.clust <- "ward.D2"

mat <- as.matrix(t(theDF_scale))
all(colnames(mat)==rownames(sampleData))

theAnno <- data.frame(Sample_type = ifelse(sampleData$is.cfdna,"cfDNA","tissue"),
                      patient = sampleData$patient)

rownames(theAnno) <- sampleData$original_sample
theAnno$patient <- gsub("DTB-","",theAnno$patient)

theAnno$'ct-fraction/purity'[theAnno$Sample_type == "cfDNA"] <- pred_ctfrac_wcdt_cfdna[gsub("_cfdna","",rownames(theAnno)[theAnno$Sample_type == "cfDNA"]),"ctfrac_by_5hmc_pred"]
sampsWpurity <- intersect(rownames(theAnno), rownames(metadata_wcdt_genomics)) 
theAnno[sampsWpurity,"ct-fraction/purity"] <- metadata_wcdt_genomics[sampsWpurity,"purity"]

sampleColors <- colorRampPalette(brewer.pal(9, "Set3"))(length(unique(theAnno$patient)))
names(sampleColors) <- unique(theAnno$patient) 

column_ha <- HeatmapAnnotation(df = theAnno[,c("Sample_type","patient","ct-fraction/purity")],
                               Patient = anno_text(theAnno$patient, rot = 0, location = 3, just = "center",  gp = gpar(fontsize = 10)),
                               gap = unit(c(0,0,-0.1), "in"),
                               border  = TRUE,
                               gp = gpar(col = "black",lwd = 0.5),
                               show_annotation_name = T,
                               col = list(Sample_type = cols_biopsy_type,
                                          patient = sampleColors,
                                          'ct-fraction/purity' = col_fun_purity),
                               show_legend = c(T,F,T),
                               annotation_legend_param = list("Sample_type" = list(title = "Sample type")))

hm <- Heatmap(mat, 
              name = "5hmC levels",
              top_annotation = column_ha,
              clustering_distance_columns = method.distance,
              clustering_distance_rows = method.distance, 
              clustering_method_columns = method.clust,
              clustering_method_rows = method.clust,
              show_row_names = F,
              show_column_names = F,
              show_row_dend = FALSE)

pdf(paste0(dir_out_figures,"/Figure_6F.pdf"),width = 11, height = 5)
print(hm)
dev.off()
