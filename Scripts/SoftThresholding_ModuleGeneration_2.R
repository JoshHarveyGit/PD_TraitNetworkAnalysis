# Rscript to be run with nohup for Picking Soft Threshold

library(WGCNA)
#   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #  
######### Soft Thresholding #######
#   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #

# This part is for Substantia Nigra

load(file="/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/SN_OutRem_WGCNA.RData")
enableWGCNAThreads(8)
powers = 1:20
##Call the network topology analysis function, depending on whether you want a signed or unsigned network this can be changed (bold). I would recommend unsigned as this allows for correlations which are negative and positive
sft = pickSoftThreshold(dat3, powerVector = powers, verbose = 5, networkType = "unsigned")
# Removed , blockSize = 10000
save(sft,file="/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/SN_SFT_WGCNA.RData")

##plot results
pdf('/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/StepByStepConstructionFigs_Unsigned_SN.pdf')
## Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, Unsigned R^2", type = "n", main = paste("Scale independence"),ylim=c(0.2,1))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels= powers, cex = 0.8, col = "red")
##this line corresponds to using an R^2 cut off of h
abline(h = 0.9, col = "red")
abline(h = 0.8, col = "orange")
##Mean connectivity as
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex=0.8, col = "red")
dev.off()

# This part is for Caudate 

load(file="/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/CP_OutRem_WGCNA.RData")
enableWGCNAThreads(8)
powers = 1:20
##Call the network topology analysis function, depending on whether you want a signed or unsigned network this can be changed (bold). I would recommend unsigned as this allows for correlations which are negative and positive
sft = pickSoftThreshold(dat2, powerVector = powers, verbose = 5, networkType = "unsigned")
# Removed , blockSize = 10000
save(sft,file="/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/CP_SFT_WGCNA.RData")

##plot results
pdf('/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/StepByStepConstructionFigs_Unsigned_CP.pdf')
## Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, Unsigned R^2", type = "n", main = paste("Scale independence"),ylim=c(0.2,1))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels= powers, cex = 0.8, col = "red")
##this line corresponds to using an R^2 cut off of h
abline(h = 0.9, col = "red")
abline(h = 0.8, col = "orange")
##Mean connectivity as
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex=0.8, col = "red")
dev.off()


# This part is for Frontal cortex

##Set your powers to investigate
load(file="/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/FC_OutRem_WGCNA.RData")
enableWGCNAThreads(8)
powers = 1:20
##Call the network topology analysis function, depending on whether you want a signed or unsigned network this can be changed (bold). I would recommend unsigned as this allows for correlations which are negative and positive
sft = pickSoftThreshold(dat, powerVector = powers, verbose = 5, networkType = "unsigned")
# Removed , blockSize = 10000
save(sft,file="/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/FC_SFT_WGCNA.RData")

##plot results
pdf('/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/StepByStepConstructionFigs_Unsigned_FCX.pdf')
## Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, Unsigned R^2", type = "n", main = paste("Scale independence"),ylim=c(0.2,1))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels= powers, cex = 0.8, col = "red")
##this line corresponds to using an R^2 cut off of h
abline(h = 0.9, col = "red")
abline(h = 0.8, col = "orange")
##Mean connectivity as
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex=0.8, col = "red")
dev.off()
#   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #  
######### Module Construction #######
#   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #

load("/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/FC_OutRem_WGCNA.RData")

###After looking at your plot, choose the lowest number power that is closest to the 0.9 line
softPower = 12;
##use this function for large dataset. NOTE maxBlockSize depends on computer memory. Also depending on what your minimum module size required is. This will be a long part
enableWGCNAThreads(16)
bwnet.fcx = blockwiseModules(dat, maxBlockSize = 10000, power = softPower, TOMType = "unsigned", minModuleSize = 100, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, saveTOMs= FALSE, verbose = 3)
save(bwnet.fcx, file = "/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/FC_CBlockwise_max_10000_SignedNetwork_cord_Rerun_Power12.Rdata")

load("/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/SN_OutRem_WGCNA.RData")
###After looking at your plot, choose the lowest number power that is closest to the 0.9 line
softPower = 8;
enableWGCNAThreads(16)
bwnet.sn = blockwiseModules(dat3, maxBlockSize = 10000, power = softPower, TOMType = "unsigned", minModuleSize = 100, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, saveTOMs= FALSE, verbose = 3)
save(bwnet.sn, file = "/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/SN_CBlockwise_max_10000_SignedNetwork_cord_Rerun_Power8.Rdata")


load("/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/CP_OutRem_WGCNA.RData")
###After looking at your plot, choose the lowest number power that is closest to the 0.9 line
softPower = 9;
##use this function for large dataset. NOTE maxBlockSize depends on computer memory. Also depending on what your minimum module size required is. This will be a long part
enableWGCNAThreads(16)
bwnet.cp = blockwiseModules(dat2, maxBlockSize = 10000, power = softPower, TOMType = "unsigned", minModuleSize = 100, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, saveTOMs= FALSE, verbose = 3)
save(bwnet.cp, file = "/mnt/data1/Josh/Cardiff_PD/PDUK.IDAT/EWAS/WGCNA/regrValidation/CP_CBlockwise_max_10000_SignedNetwork_cord_Rerun_Power9.Rdata")


