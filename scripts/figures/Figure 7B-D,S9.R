###########################################################################################################
# Figure 5E
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
# Settings and setup
endpoint <- "dead"
endpoint_time <- "OS.months"

all(rownames(clin_ubc)==colnames(counts_ubc_cfdna_tpm))
pats2use <- colnames(counts_ubc_cfdna_tpm)

clin_ubc <- clin_ubc[pats2use,]
hmc_ubc <- as.matrix(counts_ubc_cfdna_tpm[,pats2use])

theSig <- gi_sig
tissuemap_gi_score_ubc <- rowSums(tissue_map_ubc[pats2use,theSig,drop=F])   
all(names(tissuemap_gi_score_ubc)==rownames(clin_ubc))
clin_ubc$tissuemap_gi_score <- tissuemap_gi_score_ubc

clin_covars <- c("age.at.CRPC","PSA.1L","type.1L","docetaxel.for.mHSPC","time.to.CRPC.from.ADT","Hb.1L","visceral.met")
clin_covars_ctfrac_targeted <- c("age.at.CRPC","PSA.1L","type.1L","docetaxel.for.mHSPC","time.to.CRPC.from.ADT","Hb.1L","visceral.met","est_ctfrac_targeted")

###########################################################################################################
#Functions
plotGGsurv <- function(theGGdata, colEvent, colTime, colGroup, clinCovars = NULL, coxGroup = NULL,
                       groupLevels = c(FALSE,TRUE),
                       groupLabels = c("FALSE","TRUE"),
                       baseName = NULL,
                       Title = "Title",
                       xLab = "Months",
                       yLab = "Overall survival",
                       xLim = c(0,42),
                       yLim = c(0,1),
                       Colors = cols2,
                       textXpos = 25,
                       textYpos = 0.8,
                       legendXpos = 0.13,
                       legendYpos = 0.3,
                       adjustHR = FALSE,
                       uni_P  = c("wald","logrank"),
                       add_default_pvalue = TRUE,
                       savePlot = FALSE,
                       plotH = 8,
                       plotW = 9)
{
  
  # Data checks
  # Check factors
  if(!all(groupLevels %in% theGGdata[,colGroup])){
    stop("Not all group leveles are found in variable!")
  }
  
  if(!(all(unique(theGGdata[,colGroup]) %in% groupLevels))){
    stop("Not all group levels in variable are supplied!")
  }
  
  if(sum(is.na(theGGdata[,colGroup])) > 1){
    stop("groupLevels contains NAs!")
  }
  
  if(length(unique(groupLevels)) > 2){
    TREND = TRUE
  } else{
    TREND = FALSE
  }
  
  theGGdata[,colGroup] <- factor(theGGdata[,colGroup],levels = groupLevels, labels = groupLabels)
  theGGdata[,colGroup] <- as.numeric(theGGdata[,colGroup])
  
  theSurv <- Surv(time = theGGdata[,colTime], event = theGGdata[,colEvent])

  theFit <- surv_fit(theSurv ~ theGGdata[,colGroup], data = theGGdata) 
  
  if(is.null(coxGroup)){
    theCoxSummary <- summary(coxph(theSurv ~ ., data = theGGdata[,colGroup, drop = F]))
    
    print(theCoxSummary)
    
    P_wald_uni <- signif(theCoxSummary$coefficients[colGroup,"Pr(>|z|)"],2)
    P_wald_uni
    HR <- signif(theCoxSummary$coefficients[colGroup,"exp(coef)"],2) 
    HRci <- paste0("[",
                   signif(theCoxSummary$conf.int[colGroup,"lower .95"],2),
                   "-",
                   signif(theCoxSummary$conf.int[colGroup,"upper .95"],2),
                   "]")
    
    if(length(groupLevels) == 2){
      survdiff_obj <- survdiff(theSurv~theGGdata[,colGroup])
      pval_logrank <- 1 - pchisq(survdiff_obj$chisq, length(survdiff_obj$n) - 1)
      pval_logrank <- signif(pval_logrank,2)
    } else {
      pval_logrank <- "NA, not currently available for more than 2 groups!"
    }
    
    if(uni_P == "wald"){
      theText = paste0("HR = ",HR," ",HRci,"\nP (Wald) = ",P_wald_uni) 
    }else if(uni_P=="logrank"){
      theText = paste0("HR = ",HR," ",HRci,"\nP (logrank) = ",pval_logrank) 
    }else{
      stop("Univariate p-value must be specified as logrank or wald")
    }
    
    
    if(adjustHR){
      theAdjustedCoxSummary <- summary(coxph(theSurv~., data = theGGdata[,c(colGroup,clinCovars)]))
      
      print(theAdjustedCoxSummary)
      
      Padj <- signif(theAdjustedCoxSummary$coefficients[colGroup,"Pr(>|z|)"],2)
      HRadj <- signif(theAdjustedCoxSummary$coefficients[colGroup,"exp(coef)"],2) 
      HRciadj <- paste0("[",
                        signif(theAdjustedCoxSummary$conf.int[colGroup,"lower .95"],2),
                        "-",
                        signif(theAdjustedCoxSummary$conf.int[colGroup,"upper .95"],2),
                        "]")
      theText = paste0(theText,"\n", "Adj. HR = ",HRadj," ",HRciadj,"\nAdj. P (Wald) = ",Padj)
    }
    
  }else{
    
    theCoxSummary <- summary(coxph(theSurv ~ ., data = theGGdata[,coxGroup, drop = F]))
    P_wald_uni <- signif(theCoxSummary$coefficients[coxGroup,"Pr(>|z|)"],2)
    P_wald_uni
    HR <- signif(theCoxSummary$coefficients[coxGroup,"exp(coef)"],2) 
    HRci <- paste0("[",
                   signif(theCoxSummary$conf.int[coxGroup,"lower .95"],2),
                   "-",
                   signif(theCoxSummary$conf.int[coxGroup,"upper .95"],2),
                   "]")
    
    theText = paste0("HR = ",HR," ",HRci,"\nP (Wald) = ",P_wald_uni) 

    
    if(adjustHR){
      theAdjustedCoxSummary <- summary(coxph(theSurv~., data = theGGdata[,c(coxGroup,clinCovars)]))

      print(theAdjustedCoxSummary)
      
      Padj <- signif(theAdjustedCoxSummary$coefficients[coxGroup,"Pr(>|z|)"],2)
      HRadj <- signif(theAdjustedCoxSummary$coefficients[coxGroup,"exp(coef)"],2) 
      HRciadj <- paste0("[",
                        signif(theAdjustedCoxSummary$conf.int[coxGroup,"lower .95"],2),
                        "-",
                        signif(theAdjustedCoxSummary$conf.int[coxGroup,"upper .95"],2),
                        "]")
      theText = paste0(theText,"\n", "Adj. HR = ",HRadj," ",HRciadj,"\nAdj. P (Wald) = ",Padj)
    }
  }

  
  #print(theText)
  #survplot(theFit)
  #plot(theFit)
  
  theGGsurv <- ggsurvplot(theFit,
                          #data = theGGdata, # does not matter when it is inside the function?
                          pval = add_default_pvalue, 
                          pval.coord = c(2,0.05),
                          test.for.trend = TREND,
                          risk.table = TRUE,
                          risk.table.pos = "out",
                          palette = Colors,
                          size = 1,
                          risk.table.height = 0.2,
                          xlab = xLab,
                          ylab = yLab,
                          risk.table.y.text.col = TRUE,
                          risk.table.y.text = FALSE,
                          tables.theme = theme_cleantable(plot.title = element_text(size = 14)),
                          legend.labs = groupLabels,
                          conf.int = FALSE,
                          legend.title = "",
                          title = Title,
                          risk.table.title = "Number at risk",
                          legend = c(legendXpos,legendYpos),
                          censor.size = 10,
                          font.title = c(18,"bold"),
                          font.x = c(18),
                          font.y = c(18),
                          font.tickslab = c(14),
                          risk.table.font.size = 16,
                          axes.offset = FALSE,
                          xlim = xLim,
                          ylim = yLim
  )
  
  theGGsurv$plot <- theGGsurv$plot + ggplot2::annotate("text", x = textXpos, y = textYpos, label = theText, hjust = 0, cex = 5)
  theGGsurv$plot <- theGGsurv$plot + theme(plot.title = element_text(hjust = 0.5),
                                           legend.text = element_text(size = 14))
  
  
  print(theGGsurv)
  
  ### SAVE as PDF and PNG
  if(savePlot){
    pdf(file = paste0(baseName,".pdf"), width = plotW, height = plotH, onefile = FALSE)
    print(theGGsurv)
    dev.off()
    
    png(filename = paste0(baseName,".png"), width = plotW, height = plotH, units = "in", res = 300)
    print(theGGsurv)
    dev.off()
    
  }
}

plotScatter <- function(var1,var2,xLab = "xlab", yLab = "ylab", Main = "main", pXpos = 1, pYpos = 1,...){
  C <- signif(cor(var1,var2, use = "complete.obs", method = "spearman"),2)
  P <- signif(cor.test(var1,var2,use = "complete.obs", method = "spearman")$p.value,2)
  if(P==0){
    P <- "< 2.2e-16"
  }
  fit <-  lm(var2~var1)
  A <- fit$coefficients[1]
  B <- fit$coefficients[2]
  
  plot(var1,var2,xlab = xLab, ylab = yLab, main = Main,
       pch = 16, cex = 1.5, col = "darkblue")
  text(x = pXpos, y = pYpos, labels = paste0("Rho = ",C,"\nP = ",P))
  abline(a=A,b=B,lty = 2, lwd = 2, col = 2)
}


getTertiles <- function(theScore){
  tertiles <- rep(NA, length(theScore))
  tertiles[theScore < quantile(theScore,0.33, na.rm=T)] <- 1
  tertiles[theScore >= quantile(theScore,0.33, na.rm=T)] <- 2
  tertiles[theScore >= quantile(theScore,0.67, na.rm=T)] <- 3
  return(tertiles)
}

###########################################################################################################
# Figure S9A
###########################################################################################################
mat <- as.matrix(log2(hmc_ubc+1))
all(rownames(mat)==rownames(ensembl2sym))
mat <- mat[ensembl2sym$type=="protein_coding",]

theVars <- apply(mat, 1, var)
highVar <- theVars > quantile(theVars, 0.8)
mat <- mat[highVar,]

res.pca <- prcomp(t(mat), scale = T)
table(rownames(res.pca$x)==rownames(clin_ubc))
X <- as.data.frame(res.pca$x)

myPalette <- colorRampPalette(viridis(12,option = "E"))
colorVar <- clin_ubc$est_ctfrac_targeted

p <- ggplot(X, aes(x = PC1, y = PC2, label = rownames(X))) +
  geom_point(aes(colour = colorVar), size = 4) +
  #geom_text_repel(max.overlaps = 10) +
  scale_colour_gradientn(colours = myPalette(100), limits=c(min(colorVar),max(colorVar)),name = "ct-fraction (DNA-seq)") +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        axis.line.x = element_line(),
        axis.line.y = element_line())
p

pdf(paste0(dir_out_figures,"/Figure_S9A.pdf"), height = 7, width = 7)
p
dev.off()


###########################################################################################################
# Figure S9B
###########################################################################################################
pdf(paste0(dir_out_figures,"/Figure_S9B.pdf"), height = 7, width = 7)
plotScatter(clin_ubc$est_ctfrac_targeted*100,clin_ubc$est_ctfrac_5hmc, yLab = "ct-fraction predicted by 5hmC (%)", 
            xLab = "ct-fraction by targeted DNA-seq (%)", Main = "",pXpos = 75, pYpos = 25)
dev.off()

cor.test(clin_ubc$est_ctfrac_targeted*100,clin_ubc$prostate, method = "spearman")
cor.test(clin_ubc$est_ctfrac_targeted*100,clin_ubc$bone.marrow, method = "spearman")

plotScatter(clin_ubc$est_ctfrac_targeted*100,clin_ubc$prostate, yLab = "Tissue Map prostate score", 
            xLab = "ct-fraction by targeted DNA-seq (%)", Main = "5hmC prostate score vs. ct-fraction",pXpos = 50, pYpos = 0.3)
dev.off()

plotScatter(clin_ubc$est_ctfrac_targeted*100,clin_ubc$bone.marrow, yLab = "Tissue Map bone marrow", 
            xLab = "ct-fraction by targeted DNA-seq (%)", Main = "5hmC bone marrow score vs. ct-fraction",pXpos = 70, pYpos = 0.8)
dev.off()

###########################################################################################################
# Figure S9C
###########################################################################################################
clin_ubc$est_ctfrac_targeted_tertiles <- getTertiles(clin_ubc$est_ctfrac_targeted)

plotGGsurv(theGGdata = clin_ubc, colEvent = "dead", colTime = "os", colGroup = "est_ctfrac_targeted_tertiles", clinCovars = clin_covars, adjustHR = TRUE,
           groupLevels = c(1,2,3), groupLabels = c("Low","Medium","High"), Colors = cols3, 
           baseName = paste0(dir_out_figures,"/Figure_S9C"), savePlot = FALSE,
           Title = "ct-fraction (targeted DNA-seq)",
           add_default_pvalue = FALSE)

###########################################################################################################
# Figure 7B
###########################################################################################################
clin_ubc$est_ctfrac_5hmc_tertiles <- getTertiles(clin_ubc$est_ctfrac_5hmc)

plotGGsurv(theGGdata = clin_ubc, colEvent = "dead", colTime = "os", colGroup = "est_ctfrac_5hmc_tertiles", clinCovars = clin_covars, adjustHR = TRUE,
           groupLevels = c(1,2,3), groupLabels = c("Low","Medium","High"), Colors = cols3, 
           baseName = paste0(dir_out_figures,"/Figure_7B"), savePlot = FALSE,
           Title = "ct-fraction, by 5hmC",
           add_default_pvalue = FALSE)

###########################################################################################################
theGenes <- c("AR","MYC","FOXA1","PIK3CA","TP53","RB1","BRCA2","PTEN","NKX3-1","ERG","BRAF","NCOA2","TOP2A","EZH2","MKI67","AURKA", "SPOP","CLU","ATM","CHD1","SPINK1", "CTNNB1","AKT1","MYCN","ZBTB16","NCOR1","NCOR2","SPOP")
geneTypes <- c("gain","gain","gain","gain","loss","loss","loss","loss","loss", "gain","gain","gain","gain","gain","gain","gain","gain","loss","loss", "loss", "gain","gain","gain","gain","loss","gain","gain","loss")
#View(cbind(theGenes, geneTypes))

stopifnot(all(rownames(clin_ubc)==colnames(hmc_ubc)))

for(i in 1:length(theGenes)){
  
  theGeneSymbol <- theGenes[i]
  geneType <- geneTypes[i]
  theGene <- rownames(ensembl2sym)[ensembl2sym$name==theGeneSymbol]
  stopifnot(length(theGene)==1)
  
  cutoff <- NA
  if(geneType == "gain"){
    cutoff <- 0.75
  }else{
    cutoff <- 0.25
  }
  
  hmc_score <- hmc_ubc[theGene,]
  hmc_high <- hmc_score > quantile(hmc_score,cutoff)
  
  # Handle dashes (e.g. NKX3-1)
  hmc_name <- gsub("-","_",paste0(theGeneSymbol,"_hmc_high"))
  wgs_name <- gsub("-","_",paste0(theGeneSymbol,"_wgs_high"))
  
  clin_ubc[,hmc_name] <- hmc_high
}

###########################################################################################################
# Figure 7D
###########################################################################################################
clin_ubc$ezh2_top2a_comb <- clin_ubc$EZH2_hmc_high + clin_ubc$TOP2A_hmc_high

plotGGsurv(theGGdata = clin_ubc, colEvent = "dead", colTime = "os", colGroup = "ezh2_top2a_comb", clinCovars = clin_covars_ctfrac_targeted, adjustHR = TRUE,
           groupLevels = c(0,1,2), groupLabels = c("None","Either","Both"), Colors = cols3, 
           baseName = paste0(dir_out_figures,"/Figure_7D"), savePlot = FALSE,
           Title = "EZH2 and TOP2A",
           add_default_pvalue = FALSE, legendXpos = 0.10, legendYpos = 0.15)

###########################################################################################################
# Prep data
###########################################################################################################
tsgs <- c("TP53", "RB1","PTEN","BRCA1","BRCA2","APC","ATM","CDK12","CHD1","NKX3-1","SPOP","ZFHX3","CLU","CHD1")
oncg <- c("AR","MYC","PIK3CA","FOXA1","NCOA2","ERG","BRAF","MED12","AKT1")

cutoff_tsg <- 0.25
cutoff_oncogene <- 0.75

stopifnot(all(rownames(clin_ubc)==colnames(hmc_ubc)))

for(i in 1:length(tsgs)){
  theGeneSymbol <- tsgs[i]
  theGene <- rownames(ensembl2sym)[ensembl2sym$name==theGeneSymbol]
  stopifnot(length(theGene)==1)
  
  hmc_score <- hmc_ubc[theGene,]
  hmc_low <- hmc_score < quantile(hmc_score,cutoff_tsg)
  hmc_name <- gsub("-","_",paste0(theGeneSymbol,"_hmc_loss"))
  
  clin_ubc[,hmc_name] <- hmc_low
}

for(i in 1:length(oncg)){
  theGeneSymbol <- oncg[i]
  theGene <- rownames(ensembl2sym)[ensembl2sym$name==theGeneSymbol]
  stopifnot(length(theGene)==1)
  
  hmc_score <- hmc_ubc[theGene,]
  hmc_high <- hmc_score > quantile(hmc_score,cutoff_oncogene)
  hmc_name <- gsub("-","_",paste0(theGeneSymbol,"_hmc_gain"))

  clin_ubc[,hmc_name] <- hmc_high
}

#Top oncogenes: AR, MYC, NCOA2, 
#Top TSGs: TP53, PTEN, NKX3-1, RB1, BRCA2

clin_ubc$nb_alt_hmc <- rowSums(clin_ubc[,c("TP53_hmc_loss","RB1_hmc_loss","AR_hmc_gain", "MYC_hmc_gain", "BRCA2_hmc_loss", "NCOA2_hmc_gain", "PTEN_hmc_loss", "NKX3_1_hmc_loss")])
clin_ubc$nb_alt_targ <- rowSums(clin_ubc[,c("tp53_loss_mut","rb1_loss_mut","brca2_1loss_mut","AR_ampl","myc_gain","nkx31_loss","ncoa2_gain","pten_loss_mut")])

get3levels <- function(theScore){
  comb_score <- theScore
  comb_score[comb_score < 2] <- 0
  comb_score[comb_score %in% c(2:3)] <- 1
  comb_score[comb_score > 3 ] <- 2
  return(comb_score)
}

clin_ubc$nb_alt_hmc_3 <- get3levels(clin_ubc$nb_alt_hmc)
clin_ubc$nb_alt_targ_3 <- get3levels(clin_ubc$nb_alt_targ)

###########################################################################################################
# Figure 7C
###########################################################################################################
# NB this analysis used full variable with number of alterations to keep info. I.e. the HR is not for the plotted curves, which are grouped to avoid overplotting.

plotGGsurv(theGGdata = clin_ubc, colEvent = "dead", colTime = "os", colGroup = "nb_alt_hmc_3", coxGroup = "nb_alt_hmc", clinCovars = clin_covars_ctfrac_targeted, adjustHR = TRUE,
           groupLevels = c(0,1,2), groupLabels = c("0-1","2-3","4+"), Colors = cols3, 
           baseName = paste0(dir_out_figures,"/Figure_7C"), savePlot = FALSE,
           Title = "Number of alterations by 5hmC",
           add_default_pvalue = FALSE)


###########################################################################################################
# Figure S9D
###########################################################################################################
# NB this analysis used full variable with number of alterations to keep info. I.e. the HR is not for the plotted curves, which are grouped to avoid overplotting.
plotGGsurv(theGGdata = clin_ubc, colEvent = "dead", colTime = "os", colGroup = "nb_alt_targ_3", coxGroup = "nb_alt_targ", clinCovars = clin_covars_ctfrac_targeted, adjustHR = TRUE,
           groupLevels = c(0,1,2), groupLabels = c("0-1","2-3","4+"), Colors = cols3, 
           baseName = paste0(dir_out_figures,"/Figure_S9D"), savePlot = FALSE,
           Title = "Number of alterations by targeted DNA-seq",
           add_default_pvalue = FALSE)