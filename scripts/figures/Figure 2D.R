###########################################################################################################
###########################################################################################################
# Figure 2D 
# 2022-01-05
# Arian Lundberg/Martin Sjoestroem
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

###########################################################################################################library(tidyverse)
levelNames = c("Normal adjacent tissue", "Localized prostate cancer","mCRPC")
startPos = 41 
endPos = 61
xLab = "Genomic position" 
xBreaks = c(0,21,41,61,81,101) 
xLabels = c("-5kb","-2.5kb","Start","End","+2.5kb","+5kb")
legendXpos = 0.1 
legendYpos = 0.8
yLab = "5hmC enrichment" 
theColors = c("springgreen3","goldenrod3","firebrick3")
titles <- c("AR mCRPC", "AR Localized", "AR Normal",
            "FOXA1 mCRPC", "FOXA1 Localized", "FOXA1 Normal",
            "H3K27ac mCRPC", "H3K27ac Localized", "H3K27ac Normal",
            "HOXB13 mCRPC", "HOXB13 Localized", "HOXB13 Normal",
            "H3K4me2 Localized", "H3K4me2 Normal",
            "H3K4me3 Localized", "H3K4me3 Normal",
            "H3K27me3 Localized", "H3K27me3 Normal")
myTitle <- foreach(i=titles) %do% {
  paste0(i,", average 5hmC enrichment")};names(myTitle) <- titles


samples.mcrpc <- lapply(all.samples.list,'[',c('MCRPC'))
samples.local <- lapply(all.samples.list,'[',c('LOCAL'))
samples.normal <- lapply(all.samples.list,'[',c('NORMAL'))


MCRPC_AR.data <- data.frame(samples='mCRPC',gene='AR',samples.mcrpc[c("AR_NORMAL_All","AR_LOCAL_All","AR_MCRPC_All")]);colnames(MCRPC_AR.data)[-c(1:2)] <- c('Normal','Primary','mCRPC')
Local_AR.data <- data.frame(samples='Localized',gene='AR',samples.local[c("AR_NORMAL_All","AR_LOCAL_All","AR_MCRPC_All")]);colnames(Local_AR.data)[-c(1:2)] <- c('Normal','Primary','mCRPC')
Normal_AR.data <- data.frame(samples='Normal',gene='AR',samples.normal[c("AR_NORMAL_All","AR_LOCAL_All","AR_MCRPC_All")]);colnames(Normal_AR.data)[-c(1:2)] <- c('Normal','Primary','mCRPC')

MCRPC_FOXA1.data <- data.frame(samples='mCRPC',gene='FOXA1',samples.mcrpc[c("FOXA1_NORMAL_All","FOXA1_LOCAL_All","FOXA1_MCRPC_All")]);colnames(MCRPC_FOXA1.data)[-c(1:2)] <- c('Normal','Primary','mCRPC')
Local_FOXA1.data <- data.frame(samples='Localized',gene='FOXA1',samples.local[c("FOXA1_NORMAL_All","FOXA1_LOCAL_All","FOXA1_MCRPC_All")]);colnames(Local_FOXA1.data)[-c(1:2)] <- c('Normal','Primary','mCRPC')
Normal_FOXA1.data <- data.frame(samples='Normal',gene='FOXA1',samples.normal[c("FOXA1_NORMAL_All","FOXA1_LOCAL_All","FOXA1_MCRPC_All")]);colnames(Normal_FOXA1.data)[-c(1:2)] <- c('Normal','Primary','mCRPC')

MCRPC_H3K27AC.data <- data.frame(samples='mCRPC',gene='H3K27ac',samples.mcrpc[c("H3K27AC_NORMAL_All","H3K27AC_LOCAL_All","H3K27AC_MCRPC_All")]);colnames(MCRPC_H3K27AC.data)[-c(1:2)] <- c('Normal','Primary','mCRPC')
Local_H3K27AC.data <- data.frame(samples='Localized',gene='H3K27ac',samples.local[c("H3K27AC_NORMAL_All","H3K27AC_LOCAL_All","H3K27AC_MCRPC_All")]);colnames(Local_H3K27AC.data)[-c(1:2)] <- c('Normal','Primary','mCRPC')
Normal_H3K27AC.data <- data.frame(samples='Normal',gene='H3K27ac',samples.normal[c("H3K27AC_NORMAL_All","H3K27AC_LOCAL_All","H3K27AC_MCRPC_All")]);colnames(Normal_H3K27AC.data)[-c(1:2)] <- c('Normal','Primary','mCRPC')

MCRPC_HOXB13.data <- data.frame(samples='mCRPC',gene='HOXB13',samples.mcrpc[c("HOXB13_NORMAL_All","HOXB13_LOCAL_All","HOXB13_MCRPC_All")]);colnames(MCRPC_HOXB13.data)[-c(1:2)] <- c('Normal','Primary','mCRPC')
Local_HOXB13.data <- data.frame(samples='Localized',gene='HOXB13',samples.local[c("HOXB13_NORMAL_All","HOXB13_LOCAL_All","HOXB13_MCRPC_All")]);colnames(Local_HOXB13.data)[-c(1:2)] <- c('Normal','Primary','mCRPC')
Normal_HOXB13.data <- data.frame(samples='Normal',gene='HOXB13',samples.normal[c("HOXB13_NORMAL_All","HOXB13_LOCAL_All","HOXB13_MCRPC_All")]);colnames(Normal_HOXB13.data)[-c(1:2)] <- c('Normal','Primary','mCRPC')

MCRPC_H3K4ME2.data <- data.frame(samples='mCRPC',gene='H3K4me2',samples.mcrpc[c("H3K4ME2_NORMAL_All","H3K4ME2_LOCAL_All")]);colnames(MCRPC_H3K4ME2.data)[-c(1:2)] <- c('Normal','Primary')
Local_H3K4ME2.data <- data.frame(samples='Localized',gene='H3K4me2',samples.local[c("H3K4ME2_NORMAL_All","H3K4ME2_LOCAL_All")]);colnames(Local_H3K4ME2.data)[-c(1:2)] <- c('Normal','Primary')
Normal_H3K4ME2.data <- data.frame(samples='Normal',gene='H3K4me2',samples.normal[c("H3K4ME2_NORMAL_All","H3K4ME2_LOCAL_All")]);colnames(Normal_H3K4ME2.data)[-c(1:2)] <- c('Normal','Primary')

MCRPC_H3K4ME3.data <- data.frame(samples='mCRPC',gene='H3K4me3',samples.mcrpc[c("H3K4ME3_NORMAL_All","H3K4ME3_LOCAL_All")]);colnames(MCRPC_H3K4ME3.data)[-c(1:2)] <- c('Normal','Primary')
Local_H3K4ME3.data <- data.frame(samples='Localized',gene='H3K4me3',samples.local[c("H3K4ME3_NORMAL_All","H3K4ME3_LOCAL_All")]);colnames(Local_H3K4ME3.data)[-c(1:2)] <- c('Normal','Primary')
Normal_H3K4ME3.data <- data.frame(samples='Normal',gene='H3K4me3',samples.normal[c("H3K4ME3_NORMAL_All","H3K4ME3_LOCAL_All")]);colnames(Normal_H3K4ME3.data)[-c(1:2)] <- c('Normal','Primary')

MCRPC_H3K27ME3.data <- data.frame(samples='mCRPC',gene='H3K27me3',samples.mcrpc[c("H3K27ME3_NORMAL_All","H3K27ME3_LOCAL_All")]);colnames(MCRPC_H3K27ME3.data)[-c(1:2)] <- c('Normal','Primary')
Local_H3K27ME3.data <- data.frame(samples='Localized',gene='H3K27me3',samples.local[c("H3K27ME3_NORMAL_All","H3K27ME3_LOCAL_All")]);colnames(Local_H3K27ME3.data)[-c(1:2)] <- c('Normal','Primary')
Normal_H3K27ME3.data <- data.frame(samples='Normal',gene='H3K27me3',samples.normal[c("H3K27ME3_NORMAL_All","H3K27ME3_LOCAL_All")]);colnames(Normal_H3K27ME3.data)[-c(1:2)] <- c('Normal','Primary')


Fig2D.data <- cbind(x=as.factor(rep(1:101,3)),rbind(MCRPC_AR.data,Local_AR.data,Normal_AR.data,
                                                           MCRPC_FOXA1.data,Local_FOXA1.data,Normal_FOXA1.data,
                                                           MCRPC_H3K27AC.data,Local_H3K27AC.data,Normal_H3K27AC.data,
                                                           MCRPC_HOXB13.data,Local_HOXB13.data,Normal_HOXB13.data))[,c("x","samples","gene","mCRPC")]
Fig2D.data <- melt(Fig2D.data)

Fig2D.data$x <- as.numeric(Fig2D.data$x)
Fig2D.data$samples <- factor(Fig2D.data$samples,levels=c('Normal','Localized','mCRPC'))
Fig2D.data$gene <- as.factor(Fig2D.data$gene)

name_vector <- c("mCRPC" = "")

Fig2D.plot <- ggplot(Fig2D.data, aes(x =x,color=samples)) +
  geom_vline(xintercept = startPos, linetype = "dotted", color = "red", size = 0.5, alpha = 0.5) +
  geom_vline(xintercept = endPos, linetype = "dotted", color = "red", size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype='solid', color = "black", size = 0.5, alpha = 0.2) +
  labs(title = '', y = yLab, x = xLab) +
  theme_few() +
  geom_spline(aes(y = value), size = 1) +
  scale_color_manual(values=c('springgreen3','goldenrod3','firebrick3'),name = "Sample type",guide = "legend",labels = levelNames)+
  scale_x_continuous(breaks = c(0,21,41,61,81,101),labels = c("-5kb","-2.5kb","Start","End","+2.5kb","+5kb")) +
  facet_wrap(~variable~gene,nrow=1, labeller = labeller(variable  = name_vector))+
  theme(legend.position = 'left',
        plot.title = element_text(size=18),axis.text = element_text(size=12,angle=45,hjust = c(1,1)),axis.title = element_text(size=15),legend.title = element_text(size=17), legend.text = element_text(size=15),
        legend.background = element_rect(),strip.text.x = element_text(size = 15),
        legend.key = element_rect(),strip.text.y = element_text(size = 12, color = "black", face = "bold"),
        plot.tag=element_text(angle=-90,face = 'bold'),
        plot.tag.position=c(1.02, 0.5),plot.margin = unit(c(1, 1, 1, 1), "cm"))

pdf(file=paste0(dir_out_figures,"/Figure_2D.pdf"), paper = 'special',onefile = T, width = 15, height = 6,pointsize = 10, colormodel = 'cmyk')
Fig2D.plot
dev.off()