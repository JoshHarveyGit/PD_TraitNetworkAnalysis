#   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #  
######### Introduction ################
#   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   # 

# Aim: Run WGCNA for module generation
# Run By: Joshua Harvey
# Date: 16/11/2021
#     

#   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #  
######### Packages ####################
#   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #

library(minfi)
library(wateRmelon)
library(dplyr)
library(stringr)
library(corrplot)
library(cowplot)
library(rgl)
library(WGCNA,options(rlib_downstream_check = FALSE))
#   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #  
######### Load relevant files #######
#   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #

#Load in QCed mSet / RGSet / pheno
load("/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/FinalData/CardiffDataQCed.R")
CellCounts <- read.csv("/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/FinalData/Cellcounts.csv",header = T, stringsAsFactors = F)
pheno <- cbind(pheno,CellCounts[,c(2,3)])

#   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #  
######### Filter samples and probes #######
#   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #

length(unique(pheno$Region_ID))
pheno <- pheno[-which(duplicated(pheno$Region_ID)),]

head(pheno)
betas <- betas(mSet1) # subset betas
betas <- betas[,which(colnames(betas) %in% pheno$Sentrix_ID_Position)] # subset retained samples
unique(colnames(betas) == pheno$Sentrix_ID_Position) # sanity check pheno and beta matrix match

#Remove crosshhyb
crosshyb <- read.csv("/mnt/data1/450K_reference/CrossHybridisingProbesPriceORWeksberg.csv", row.names = 1)
probes<-read.csv("/mnt/data1/450K_reference/SNPsinProbesAnno.csv", row.names = 1)
## remove cross hybridising probes
remove<-match(crosshyb[,1], rownames(betas))
remove<-remove[which(is.na(remove) != TRUE)]
betas<-betas[-remove,]
## remove SNP probes
probes<-probes[row.names(betas),]
betas<-betas[which(probes$Weksburg_CommonSNP_EuAf_within10bpSBE == ""),]
betas<-betas[-grep("rs", rownames(betas)),] 

# remove sex probes
ANNOT450<-read.csv("/mnt/data1/450K_reference/AdditionalAnnotation450K_Price_Great_SNPsinProbeSequence.csv")
rownames(ANNOT450)<-ANNOT450[,1]
sex_probes <- rownames(ANNOT450[which(ANNOT450$CHR %in% c("X","Y")),])
betas <- betas[-which(rownames(betas) %in% sex_probes),]

pcAll <- prcomp(t(betas))

phenoPC <- cbind(pheno,pcAll$x[,c(1:8)])
phenoPC <- phenoPC[,-c(12,13)]

plot(phenoPC$PC3,phenoPC$PC5)

ggplot(data = phenoPC, aes(x = as.factor(Sentrix_ID), y = PC3, color = Batch))+
  geom_point()+
  geom_boxplot()+
  geom_abline(intercept = -4.78,slope = 0)+
  coord_flip()


phenoPC$Batch <- phenoPC$PC3 > -4.78


pheno <- phenoPC 
#   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #  
######### Format pheno, add relevant variables #######
#   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #
# Subset regions and annotate as a numeric
pheno$BulkRegion <- as.factor(pheno$BulkRegion)
#Levels: CP:1 FC:2 SN:3
pheno$BulkRegion <- as.numeric(pheno$BulkRegion)
PATNO <- as.numeric(as.factor(pheno$PATNO)) # reformat PATNO for joining and downstream work

pheno <- pheno[which(pheno$Disease == "PD"),]
bPheno <- read.csv("/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/Pheno/Phenotypic files/PhenoCSV.csv", header = T, stringsAsFactors = F)
# Remove those with a long course of PD, these appear to be non-normal / asyn driven PD cases
bPheno$YrsDisease <- as.numeric(bPheno$AgeAtDeath) - as.numeric(bPheno$AAO)
bPheno <- bPheno[-which(bPheno$YrsDisease > 40),] # Two samples removed for long disease course: "PD071" "PD274"
colnames(bPheno)[1] <- "PATNO"
pheno <- pheno[which(pheno$PATNO %in% bPheno$PATNO),] # this removes controls and samples missing clinical annoation
pheno$PMI[which(is.na(pheno$PMI))] <- mean(pheno$PMI,na.rm = TRUE)
# Subset regions per beta matrix
phenoFC <- pheno[which(pheno$BulkRegion == 2),]
phenoCP <- pheno[which(pheno$BulkRegion == 1),]
phenoSN <- pheno[which(pheno$BulkRegion == 3),]
betasFC <- betas[,which(colnames(betas) %in% phenoFC$Sentrix_ID_Position)]
betasCP <- betas[,which(colnames(betas) %in% phenoCP$Sentrix_ID_Position)]
betasSN <- betas[,which(colnames(betas) %in% phenoSN$Sentrix_ID_Position)]


#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### Regress out covariates ########
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #
## Regress out covariates of age, sex and NeuN+ cell proportions
# Function to extract residuals


resid <- function(x,age,sex,cet,batch,PMI){
  fit <-try(lm(as.numeric(x) ~ age + sex + cet + batch + PMI), silent=T)
  if(inherits(fit,'try-error')) return(rep(NA,length(sex)))
  return(fit$residuals + fit$coefficients[1])
} 
# Frontal Cortex
sex    	<- 	as.factor(phenoFC$predSex3)
age    	<- 	phenoFC$AgeAtDeath 
cet		<- as.numeric(phenoFC$NeuN_pos)
batch <- as.numeric(phenoFC$Batch)
PMI <- as.numeric(phenoFC$PMI)

cl<-makeCluster(8)
normMat <- t(parApply(cl,betasFC,1,resid,age,sex,cet,batch,PMI))
stopCluster(cl)
colnames(normMat) <- colnames(betasFC)
rownames(normMat) <- rownames(betasFC) 
normFC <- normMat

# Caudate Nucleus
sex    	<- 	as.factor(phenoCP$predSex3)
age    	<- 	phenoCP$AgeAtDeath 
cet		<- as.numeric(phenoCP$NeuN_pos)
batch <- as.numeric(phenoCP$Batch)
PMI <- as.numeric(phenoCP$PMI)

cl<-makeCluster(8)
normMat <- t(parApply(cl,betasCP,1,resid,age,sex,cet,batch,PMI))
stopCluster(cl)
colnames(normMat) <- colnames(betasCP)
rownames(normMat) <- rownames(betasCP)
normCP <- normMat

# Substantia Nigra
sex    	<- 	as.factor(phenoSN$predSex3)
age    	<- 	phenoSN$AgeAtDeath 
cet		<- as.numeric(phenoSN$NeuN_pos)
batch <- as.numeric(phenoSN$Batch)
PMI <- as.numeric(phenoSN$PMI)

cl<-makeCluster(8)
normMat <- t(parApply(cl,betasSN,1,resid,age,sex,cet,batch,PMI))
stopCluster(cl)
colnames(normMat) <- colnames(betasSN)
rownames(normMat) <- rownames(betasSN)
normSN <- normMat


#remove the next three once you are done fiddling
saveFC <- normFC
saveCP <- normCP
saveSN <- normSN


#   #   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #   #  
######### Filter probes #################
#   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #

# Filter to most variable residualised probes, as calculated by median absolute deviation over a median(MAD)
varTest <- function(x){
  mad <- mad(x)
  med <- median(x)
  return(c(mad,med))
}

cl<- makeCluster(16)
resVarFC <-t(parApply(cl, normFC, 1, varTest))
stopCluster(cl)

cl<- makeCluster(16)
resVarCP <-t(parApply(cl, normCP, 1, varTest))
stopCluster(cl)

cl<- makeCluster(16)
resVarSN <-t(parApply(cl, normSN, 1, varTest))
stopCluster(cl)

medFC <- rownames(normFC[which(resVarFC[,1] > median(resVarFC[,1])),])
medCP <- rownames(normCP[which(resVarCP[,1] > median(resVarCP[,1])),])
medSN <- rownames(normSN[which(resVarSN[,1] > median(resVarSN[,1])),])

allProbes <- unique(c(medFC,medCP,medSN))

#remove non nariable probes (probes with variance <=0.0005092303
t(normFC[allProbes,])->dat
dim(dat)
#    89 243783
t(normCP[allProbes,])->dat2
dim(dat2)
#  82 243783
t(normSN[allProbes,])->dat3
dim(dat3)
# 88 243783
# 


# save(dat, file = "/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation//FC_dat_Reg_prepro.Rdata")
load("/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/FC_dat_Reg_prepro.Rdata")

# save(dat2, file = "/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/CP_dat_Reg_prepro.Rdata")
load("/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/CP_dat_Reg_prepro.Rdata")

# save(dat3, file = "/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/SN_dat_Reg_prepro.Rdata")
load(file = "/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/SN_dat_Reg_prepro.Rdata")



#   #   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #   #  
######### Check for extreme outliers ########
#   #   #   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #   #   #

#Generate hierarachical clustered sample tree to visually assess for outliers (as in the WGCNA documentation)
# This part is for Frontal cortex

# load(file = "/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/Promoters/FC_dat_Reg_prepro.Rdata")
sampleTree = hclust(dist(dat), method = "average");
# # Clustering dendrogram:
pdf(file = "/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/Promoters/SampleClusteringFC.pdf")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTreeFC, main = "Sample clustering to detect outliers", sub="", xlab="",labels=F, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2,cex=0.4)
dev.off()
save(sampleTreeFC, file = "/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/Promoters/FCsampleTree.Rdata" )
# 
# 
# # This part is for Caudate
load(file = "/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/Promoters/CP_dat_Reg_prepro.Rdata")
sampleTreeCP = hclust(dist(dat2), method = "average");
# # Clustering dendrogram:
pdf(file = "/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/Promoters/SampleClusteringCP.pdf")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTreeCP, main = "Sample clustering to detect outliers", sub="", xlab="",labels=F, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2,cex=0.4)
dev.off()
save(sampleTreeCP, file = "/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/Promoters/CPsampleTree.Rdata" )
# 
# 
# This part is for Substantia
load(file = "/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/SN_dat_Reg_prepro.Rdata")
sampleTreeSN = hclust(dist(dat3), method = "average");
# Clustering dendrogram:
pdf(file = "/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/SampleClusteringSN.pdf")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTreeSN, main = "Sample clustering to detect outliers", sub="", xlab="",labels=F, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2,cex=0.4)
dev.off()
save(sampleTreeSN, file = "/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/SNsampleTree.Rdata" )




load("/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/FCsampleTree.Rdata")
load("/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/CPsampleTree.Rdata")
load("/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/SNsampleTree.Rdata")


#Generate principal components for each dataset
pcFC <- prcomp(dat)
pcCP <- prcomp(dat2)
pcSN <- prcomp(dat3)


#Hierarchacal clustering based method, this code generates 
cutNkeep <- function(sampleTree, cutnum, region, data){
  clust <- cutreeStatic(sampleTree, cutHeight = cutnum, minSize = 10)
  message(paste(table(clust)))
  keepSamples = (clust==1)
  keepFrame <- data[keepSamples,]
  message(paste(region,"has dimensions", dim(keepFrame)))
  assign(paste("dat",region,sep = ""), keepFrame, envir = parent.frame())
  assign(paste("keepSamples",region,sep = ""), keepSamples, envir = parent.frame() )
}

cutNkeep(sampleTreeFC,30 ,"FC",dat)
cutNkeep(sampleTreeCP, 26,"CP",dat2)
cutNkeep(sampleTreeSN, 27,"SN",dat3)

# Visually assess removed samples via the sample tree based method as compered to the PC based outliers. 

pdf(file = "/gpfs/mrc0/projects/Research_Project-T112069/Meth/EWAS/WGCNA/SampleClustering.pdf",height = 10,width = 15)
par(cex = 0.6)
par(mar = c(0,4.5,2,0))
plot(sampleTreeSN, main = "Sample clustering to detect outliers: SN", sub="", xlab="",labels=F, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2,cex=0.4)
abline(h = 27, col = "red")
dev.off()

par(mfrow=c(2,2),mar = c(4,4,4,4))
plot(pcSN$x[,1],pcSN$x[,2],col = as.factor(rownames(dat3) %in% rownames(datSN)),pch = 4,xlab = "PC1", ylab = "PC2")
legend("topleft", c("Exclude", "Keep"), col = c("black", "red"), pch = 1)
plot(pcSN$x[,2],pcSN$x[,3],col = as.factor(rownames(dat3) %in% rownames(datSN)),pch = 4,xlab = "PC2", ylab = "PC3")
plot(pcSN$x[,1],pcSN$x[,4],col = as.factor(rownames(dat3) %in% rownames(datSN)),pch = 4,xlab = "PC1", ylab = "PC4")
plot(pcSN$x[,3],pcSN$x[,1],col = as.factor(rownames(dat3) %in% rownames(datSN)),pch = 4,xlab = "PC3", ylab = "PC1")
dev.off


#Keep final filtered data, with sample outliers removed
dat <- datFC
dat2 <- datCP
dat3 <- datSN
save(dat, file="/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/FC_OutRem_WGCNA.RData")
save(dat2, file="/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/CP_OutRem_WGCNA.RData")
save(dat3, file="/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/SN_OutRem_WGCNA.RData")
