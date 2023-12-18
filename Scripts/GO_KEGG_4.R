#Scripts to run GO and KEGG enrichment analysis on prioritised modules using missMethyl

#   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #  
######### Packages ####################
#   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #
.libPaths("/gpfs/mrc0/projects/Research_Project-T112069/packages")
setwd("/gpfs/mrc0/projects/Research_Project-T112069/Cardiff/scRNA/")
library(DropletUtils)
library(EWCE)
library(ewceData)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(cowplot)
library(missMethyl)

wgcnaLab <- read.table("/gpfs/mrc0/projects/Research_Project-T112069/Cardiff/scRNA/SNSummary.txt", header = T)

targetLabs <- rownames(wgcnaLab[which(wgcnaLab$color == "magenta"),])
allLabs <- rownames(wgcnaLab)

go_PP <- gometh(sig.cpg=targetLabs, all.cpg=allLabs, collection="GO",
              plot.bias=F,prior.prob = T)

kegg_PP <- gometh(sig.cpg=targetLabs, all.cpg=allLabs, collection="KEGG",
              plot.bias=F,prior.prob = T)
			  
topGO_PP		<- topGSA(go_PP,n = 200)

topKEGG_PP		<- topGSA(kegg_PP,n = 200)

save(topGO_PP, topKEGG_PP, file = "/gpfs/mrc0/projects/Research_Project-T112069/Cardiff/Methylation/Ontology/MM_SNdep.Rdata")
load(file = "/gpfs/mrc0/projects/Research_Project-T112069/Cardiff/Methylation/Ontology/MM_SNdep.Rdata")


cpLabs <- read.table("/gpfs/mrc0/projects/Research_Project-T112069/Cardiff/Methylation/CPSummary.Rdata", header = T)
cpLabs <- rownames(cpLabs[which(cpLabs$color == "magenta"),])

allLabs <- rownames(cpLabs)

go_PPcp <- gometh(sig.cpg=cpLabs, all.cpg=allLabs, collection="GO",
              plot.bias=F,prior.prob = T)

kegg_PPcp <- gometh(sig.cpg=cpLabs, all.cpg=allLabs, collection="KEGG",
              plot.bias=F,prior.prob = T)
			  
topGO_PPcp		<- topGSA(go_PPcp,n = 200)
topKEGG_PPcp		<- topGSA(kegg_PPcp,n = 200)

save(topGO_PPcp,topKEGG_PPcp,  file = "/gpfs/mrc0/projects/Research_Project-T112069/Cardiff/Methylation/Ontology/MM_CPagg.Rdata")

