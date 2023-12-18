
#   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #  
######### Packages ####################
#   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #
options(rlib_downstream_check = FALSE)
library(minfi)
library(wateRmelon)
library(dplyr)
library(stringr)
library(corrplot)
library(cowplot)
library(rgl)
library(lm.beta)
library(WGCNA)

setwd("/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/")
#   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #  
######### Load relevant files #######
#   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #

load("EWAS/WGCNA/regrValidation/SN_CBlockwise_max_10000_SignedNetwork_cord_Rerun_Power8.Rdata")
load("EWAS/WGCNA/regrValidation/CP_CBlockwise_max_10000_SignedNetwork_cord_Rerun_Power9.Rdata")
load("EWAS/WGCNA/regrValidation/FC_CBlockwise_max_10000_SignedNetwork_cord_Rerun_Power12.Rdata")

load("EWAS/WGCNA/regrValidation/SN_OutRem_WGCNA.RData")
load("EWAS/WGCNA/regrValidation/FC_OutRem_WGCNA.RData")
load("EWAS/WGCNA/regrValidation/CP_OutRem_WGCNA.RData")


load("phenosForWGCNA.Rdata")
load("EWAS/WGCNA/phenoBatch.Rdata")
## Create phenotype matrix (bPheno) with binarised/continuous variables for testing with a joining variable consistent with region specific pheno files
bPheno <- read.csv("Pheno/Phenotypic files/PhenoCSV.csv", header = T, stringsAsFactors = F)
colnames(bPheno)[1] <- "PATNO"
for(x in colnames(bPheno)[5:18]){
  y <- str_sub(bPheno[,x], start = 1L, end = 1L)
  y <- as.factor(y)
  bPheno[,x] <- y
}
bPheno$AAO <- as.numeric(str_sub(bPheno$AAO,-2L,-1L))
bPheno <- bPheno[,-which(colnames(bPheno) %in% c("Sex","L.Dopa_responsive","Family.hx","Paranoia","Delusion","Memory.loss"))]

######### Substantia nigra testing
# Create association pheno files
phenoSN$batch <- pheno[rownames(phenoSN),"Batch"]
phenoSN$batch <- as.numeric(phenoSN$batch)
pheno.sn <- phenoSN[rownames(dat3),c("PATNO","Age","PMI","predSex3","NeuN_pos","batch")] 
pheno.sn <- left_join(pheno.sn,bPheno,by = "PATNO")
pheno.sn$YrsDisease <- pheno.sn$AgeAtDeath -as.numeric(pheno.sn$AAO)
pheno.sn <- pheno.sn[,-which(colnames(pheno.sn) %in% c("Age","AAO"))]
pheno.sn[which(is.na(pheno.sn$PMI)),"PMI"] <- mean(pheno.sn$PMI,na.rm = TRUE) 
pheno.sn[which(is.na(pheno.sn$YrsDisease)),"YrsDisease"] <- mean(pheno.sn$YrsDisease,na.rm = TRUE)

# Source the reference function for testing correlation/association
source("/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/Scripts/traitModCor.R")

#Define primary association labels "labs" (phenotypes of interest) and define whether they are confounder/covariates "cov"
labs <- colnames(pheno.sn)[c(2:6,10:16)]
cov <- c(rep(TRUE,5),rep(FALSE,7))

#Run correlation associaiton test
traitModCor(tVec = labs, tFile = pheno.sn, mFile = bwnet.sn, rDat = dat3)

# Create P value matrix
moduleTraitP <- as.data.frame(moduleTraitPvalue)

#define modules to remove with evidence of confounding association
rem <- which(moduleTraitP$predSex3 < 0.05 | moduleTraitP$batch < 0.05 | moduleTraitP$NeuN_pos < 0.05| moduleTraitP$AgeAtDeath < 0.05
      | moduleTraitP$PMI < 0.05 | rownames(moduleTraitP) == "MEgrey")

#Remove confounded mocules
moduleTraitPvalue <- moduleTraitPvalue[-rem,]
moduleTraitCor<- moduleTraitCor[-rem,]
labs[2] <- "Sex"

# create Heatmap
# Will display correlations and their p-values
pdf("EWAS/WGCNA/regrValidation/SNcorr_ConfRemoved.pdf",height = 10)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5,5,3,1));
# Display the correlation values within a heatmap plot

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = labs,
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               cex.lab = 0.4,
               zlim = c(-1,1),
               main = paste("SN Module-trait relationships"))          
dev.off()

# Remove non-primary phenotypes
moduleTraitPvalue <- moduleTraitPvalue[,-c(1:5)]
moduleTraitCor <- moduleTraitCor[,-c(1:5)]

# Create final correlation coefficient and p-value matrices
corSN <- moduleTraitCor
pSN <- moduleTraitPvalue
#Defined FDR significant associations
FDRpSN <- pSN
for(x in colnames(pSN)){
  FDRpSN[,x] <- p.adjust(pSN[,x], method = "fdr")
}

#Join into larger plotting frame corAll and pAll
corAll <- melt(corSN)
pAll <- melt(pSN)
pFDR <- melt(FDRpSN)
pAll$Cor <- abs(corAll$value)
pAll$FDR <- pFDR$value < 0.05
hist(pAll$FDR)

#Plot sideways association plot
snAssociation <- ggplot(pAll,aes(x = as.factor(Var2),y = -log10(value), fill = as.factor(Var2)))+
  geom_hline(yintercept = -log10(c(0.05,(.05/nrow(pSN)))),col = c("grey","black"),linetype = "dashed")+
  geom_jitter(aes(size = Cor,color = FDR,stroke = FDR),shape = 21)+
  coord_flip()+
  theme_cowplot()+
  scale_color_manual(values = c("white","black"))+ guides(fill = "none")
  

###### Do the same as above but for Frontal cortex and Caudate Nucleus

phenoFC$batch <- pheno[rownames(phenoFC),"Batch"]
phenoFC$batch <- as.numeric(phenoFC$batch)
pheno.fcx <- phenoFC[rownames(dat),c("PATNO","Age","PMI","predSex3","NeuN_pos","batch")] 
pheno.fcx <- left_join(pheno.fcx,bPheno,by = "PATNO")
pheno.fcx$YrsDisease <- pheno.fcx$AgeAtDeath -as.numeric(pheno.fcx$AAO)

pheno.fcx <-pheno.fcx[,-which(colnames(pheno.fcx) %in% c("Age","AAO"))]
pheno.fcx[which(is.na(pheno.fcx$PMI)),"PMI"] <- mean(pheno.fcx$PMI,na.rm = TRUE) 
pheno.fcx[which(is.na(pheno.fcx$YrsDisease)),"YrsDisease"] <- mean(pheno.fcx$YrsDisease,na.rm = TRUE)
names(bwnet.fcx)
moduleColors = labels2colors(bwnet.fcx$colors)



source("/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/Scripts/traitModCor.R")

labs <- colnames(pheno.fcx)[c(2:6,10:16)]
cov <- c(rep(TRUE,5),rep(FALSE,7))

traitModCor(tVec = labs, tFile = pheno.fcx, mFile = bwnet.fcx, rDat = dat)



moduleTraitP <- as.data.frame(moduleTraitPvalue)
rem <- which(moduleTraitP$predSex3 < 0.05 | moduleTraitP$batch < 0.05 | moduleTraitP$NeuN_pos < 0.05| moduleTraitP$AgeAtDeath < 0.05
             | moduleTraitP$PMI < 0.05 | rownames(moduleTraitP) == "MEgrey")


moduleTraitPvalue <- moduleTraitPvalue[-rem,]
moduleTraitCor<- moduleTraitCor[-rem,]

labs[2] <- "Sex"

pdf("EWAS/WGCNA/regrValidation/FCcorr_ConfRemoved.pdf",height = 10)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5,5,3,1));
# Display the correlation values within a heatmap plot

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = labs,
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               cex.lab = 0.4,
               zlim = c(-1,1),
               main = paste("FC Module-trait relationships"))          
dev.off()

moduleTraitPvalue <- moduleTraitPvalue[,-c(1:5)]
moduleTraitCor <- moduleTraitCor[,-c(1:5)]

dim(moduleTraitCor)

min(p.adjust(moduleTraitPvalue,method = "holm"))

corFC <- moduleTraitCor
pFC <- moduleTraitPvalue

FDRpFC <- pFC

for(x in colnames(pFC)){
  FDRpFC[,x] <- p.adjust(pFC[,x], method = "fdr")
}


colnames(FDRpFC) <- colnames(pFC)
rownames(FDRpFC)<- rownames(pFC)
corAll <- melt(corFC)
pAll <- melt(pFC)
pFDR <- melt(FDRpFC)
pAll$Cor <- abs(corAll$value)
pAll$FDR <- pFDR$value < 0.05
hist(pAll$FDR)
fcAssociation <- ggplot(pAll,aes(x = as.factor(Var2),y = -log10(value), fill = as.factor(Var2)))+
  geom_hline(yintercept = -log10(c(0.05,(0.05/nrow(pFC)))),col = c("grey","black"),linetype = "dashed")+
  geom_jitter(aes(size = Cor,color = FDR,stroke = FDR),shape = 21)+
  coord_flip()+
  theme_cowplot()+
  scale_color_manual(values = c("white","black"))+ guides(fill = "none")

snAssociation+ylim(0,3.5)
fcAssociation+ylim(0,3.5)


phenoCP$batch <- pheno[rownames(phenoCP),"Batch"]
phenoCP$batch <- as.numeric(phenoCP$batch)
pheno.cp <- phenoCP[rownames(dat2),c("PATNO","Age","PMI","predSex3","NeuN_pos","batch")] 
pheno.cp <- left_join(pheno.cp,bPheno,by = "PATNO")
pheno.cp$YrsDisease <- pheno.cp$AgeAtDeath -as.numeric(pheno.cp$AAO)

pheno.cp <-pheno.cp[,-which(colnames(pheno.cp) %in% c("Age","AAO"))]
pheno.cp[which(is.na(pheno.cp$PMI)),"PMI"] <- mean(pheno.cp$PMI,na.rm = TRUE) 
pheno.cp[which(is.na(pheno.cp$YrsDisease)),"YrsDisease"] <- mean(pheno.cp$YrsDisease,na.rm = TRUE)
names(bwnet.cp)
moduleColors = labels2colors(bwnet.cp$colors)



source("/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/Scripts/traitModCor.R")

labs <- colnames(pheno.cp)[c(2:6,10:16)]
cov <- c(rep(TRUE,5),rep(FALSE,7))

traitModCor(tVec = labs, tFile = pheno.cp, mFile = bwnet.cp, rDat = dat2)



moduleTraitP <- as.data.frame(moduleTraitPvalue)
rem <- which(moduleTraitP$predSex3 < 0.05 | moduleTraitP$batch < 0.05 | moduleTraitP$NeuN_pos < 0.05| moduleTraitP$AgeAtDeath < 0.05
             | moduleTraitP$PMI < 0.05 | rownames(moduleTraitP) == "MEgrey")


moduleTraitPvalue <- moduleTraitPvalue[-rem,]
moduleTraitCor<- moduleTraitCor[-rem,]

labs[2] <- "Sex"
pdf("EWAS/WGCNA/regrValidation/CNcorr_ConfRemoved.pdf",height = 10)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5,5,3,1));
# Display the correlation values within a heatmap plot

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = labs,
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               cex.lab = 0.4,
               zlim = c(-1,1),
               main = paste("CN Module-trait relationships"))          
dev.off()


moduleTraitPvalue <- moduleTraitPvalue[,-c(1:5)]
moduleTraitCor <- moduleTraitCor[,-c(1:5)]

dim(moduleTraitCor)

min(p.adjust(moduleTraitPvalue,method = "holm"))

corCP <- moduleTraitCor
pCP <- moduleTraitPvalue
FDRpCP <- pCP

for(x in colnames(pCP)){
  FDRpCP[,x] <- p.adjust(pCP[,x], method = "fdr")
}


colnames(FDRpCP) <- colnames(pCP)
rownames(FDRpCP)<- rownames(pCP)
corAll <- melt(corCP)
pAll <- melt(pCP)function(tVec,tFile,mFile,rDat){
  message(paste("You've started the function, at least"))
  # Define your variables
  nGenes = ncol(rDat)
pFDR <- melt(FDRpCP)
pAll$Cor <- abs(corAll$value)
pAll$FDR <- pFDR$value < 0.05  

labs <- c("Dementia","Hallucinations","Depression","Anxiety","Aggression","Sleep disturbance","Years of Disease")
cpAssociation <- ggplot(pAll,aes(x = as.factor(Var2),y = -log10(value), fill = as.factor(Var2)))+
  geom_hline(yintercept = -log10(c(0.05,(.05/nrow(pCP)))),col = c("grey","black"),linetype = "dashed")+
  geom_jitter(aes(size = Cor,color = FDR,stroke = FDR),shape = 21)+
  coord_flip()+
  theme_cowplot()+
  scale_color_manual(values = c("white","black"))+ guides(fill = "none")


#Plot all together using cowplot::plot_grid()  
legendPlot <- get_legend(snAssociation+ylim(0,3.5)+xlab(NULL)+ggtitle("Substantia Nigra")+
  scale_size_binned(breaks = c(0.1,0.2,0.3),
                    limits = c(0, .4),range = c(0,6))+labs(size="Correlation Coefficient", color = "FDR p-value < 0.05"))
snPlot <- snAssociation+ylim(0,3.5)+theme(legend.position="none")+xlab(NULL)+ggtitle("Substantia Nigra")+ 
  scale_size_binned(breaks = c(0.1,0.2,0.3),
                    limits = c(0, .4))+ scale_x_discrete(labels=labs)+ylab("-log10(p-value)")
fcPlot <- fcAssociation+ylim(0,3.5)+theme(legend.position="none")+xlab(NULL)+ggtitle("Prefrontal Cortex")+ 
  scale_size_binned(breaks = c(0.1,0.2,0.3),
                    limits = c(0, .4))+ scale_x_discrete(labels=labs)+ylab("-log10(p-value)")
cpPlot <- cpAssociation+ylim(0,3.5)+theme(legend.position="none")+xlab(NULL)+ggtitle("Caudate Nucleus")+ 
  scale_size_binned(breaks = c(0.1,0.2,0.3),
                    limits = c(0, .4))+ scale_x_discrete(labels=labs)+ylab("-log10(p-value)")

pdf("EWAS/WGCNA/SummaryAssociationPlots.pdf",height = 7, width = 12)
plot_grid(snPlot,fcPlot,cpPlot,legendPlot)
dev.off()
