#Scripts to run EWCE using the Kamath dataset

#   #   #   #   #   #   #   #   #   #  
#   #   #   #   #   #   #   #   #   #  
######### Packages ####################
#   #   #   #   #   #   #   #   #   # 
#   #   #   #   #   #   #   #   #   #
.libPaths("lustre/projects/Research_Project-T112069/packages")
setwd("/lustre/projects/Research_Project-T112069/Cardiff")
library(DropletUtils)
library(EWCE)
library(ewceData)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(cowplot)

data_dir <- paste("/gpfs/mrc0/projects/Research_Project-T112069/Cardiff/scRNA/SCP/")

print(data_dir)

sce <- Seurat::Read10X(data_dir)

metaSCE <- read.table("/gpfs/mrc0/projects/Research_Project-T112069/Cardiff/scRNA/SCP/METADATA_PD.tsv.gz")

umapMeta <- read.table("/gpfs/mrc0/projects/Research_Project-T112069/Cardiff/scRNA/UMAPS/UMAP_meta.tsv",header = T)

rownames(umapMeta) <- umapMeta$NAME


sceSub <- sce[,umapMeta$NAME]

CellData <- CreateSeuratObject(counts = sceSub,meta.data = umapMeta)
SE <-SummarizedExperiment(assays=list(counts=sceSub), colData=umapMeta)

gene <- "GFAP"
cellExpDist <- data.frame(Expression=assay(SE)[gene,],
                          Celltype=SE$Cell_Cat)
						  
png("/gpfs/mrc0/projects/Research_Project-T112069/Cardiff/scRNA/TestGFAP_Expression.png",res = 150, height = 1000, width = 1600)						  
try({
  ggplot(cellExpDist) + 
  geom_boxplot(aes(x=Celltype, y=Expression)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
})
dev.off()
SE <-SummarizedExperiment(assays=list(counts=sceSub), colData=umapMeta)


# This bit of code below is for plotting raw gene expression across cell types
gene <- "TMEM200A"
cellExpDist <- data.frame(Expression=assay(SE)[gene,which(SE$Cell_Cat == "DA_neurons")],
                          Celltype=SE$Cell_Type[which(SE$Cell_Cat == "DA_neurons")])
						  
png("/gpfs/mrc0/projects/Research_Project-T112069/Cardiff/scRNA/TestTMEM200A_Expression.png",res = 150, height = 1000, width = 1600)						  
try({
  ggplot(cellExpDist) + 
  geom_boxplot(aes(x=Celltype, y=Expression)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
})
dev.off()

# Fix any potential misannotated gene names
SNsce = fix_bad_hgnc_symbols(SE)

# Generate cell type data dropping uninformative genes with low expression, minimum human orthologues and no DEG between cell types
ExpSNdrop <- drop_uninformative_genes(exp=SE,
                                                   drop_nonhuman_genes = T,
                                                   input_species = "human",
                                                   level2annot=SE$Cell_Type,no_cores = 20)

#Check 41625Check 387483
#20 core(s) assigned as workers (12reserved).
#Converting to sparse matrix.
#Checking for non-expressed genes.
#-4236 / 37389 non-expressed genes dropped
#Checking for cells with no expressed genes.
#DGE:: Limma...
#3,532 / 37,389 genes dropped @ DGE adj_pval_thresh < 1e-05
#Time difference of 1.468923 hours
#Warning message:
#In asMethod(object) :
 # sparse->dense coercion: allocating vector of size 107.9 GiB




annotLevels = list(level1class=SE$Cell_Cat,
                   level2class=SE$Cell_Type)

fNames_SN <- generate_celltype_data(exp=ExpSNdrop,
                                           annotLevels=annotLevels,
                                           groupName="scSN",
                                           savePath="/gpfs/mrc0/projects/Research_Project-T112069/Cardiff/scRNA/", no_cores = 20) 
										   
# 20 core(s) assigned as workers (12reserved).
#Converting to sparse matrix.
#+ Calculating normalized mean expression.
#Converting to sparse matrix.
#Converting to sparse matrix.
#+ Calculating normalized specificity.
#Converting to sparse matrix.
#Converting to sparse matrix.
#Converting to sparse matrix.
#Converting to sparse matrix.
#Loading required namespace: ggdendro

ctd <- EWCE::load_rdata(fNames_SN) 

plt <- EWCE::plot_ctd(ctd = ctd,
                      level = 1,
                      genes = c("RBFOX3","GFAP","TH","CSF1R","OLIG2"),
                      metric = "mean_exp")

					 
plt2 <- EWCE::plot_ctd(ctd = ctd,
                      level = 1,
                      genes = c("RBFOX3","GFAP","TH","CSF1R","OLIG2"),
                      metric = "specificity")
					  
					  
pdf("/gpfs/mrc0/projects/Research_Project-T112069/Cardiff/scRNA/CellMarkerTestsEWCE.pdf", height = 6, width = 10)					  
plot_grid(plt,plt2)
dev.off()


hits <- c("TMEM200A","SNCA","GPNMB","SULT1C2","GEM","AGTR1")

unconditional_results <- EWCE::bootstrap_enrichment_test(
  sct_data = ctd,
  hits = hits,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = 100,
  annotLevel = 2, no_cores = 20)
  
  
EWCE::plot_ctd(ctd = ctd,
                      level = 1,
                      genes = c("RBFOX3","GFAP","TH","CSF1R","OLIG2"),
                      metric = "mean_exp")  
genes <- c("PRKAR1B",
"SEPSECS",
"SPTBN2",
"LYSMD4",
"C19orf21",
"KIF5C",
"CCDC165",
"KIAA0889",
"GRAMD1B",
"LOC285577")

plot1 <- EWCE::plot_ctd(ctd = ctd,
                      level = 1,
                      genes = genes,
                      metric = "mean_exp")
					  
plot2 <- EWCE::plot_ctd(ctd = ctd,
                      level = 1,
                      genes = genes,
                      metric = "specificity")
					  
pdf("/gpfs/mrc0/projects/Research_Project-T112069/Cardiff/scRNA/magentaHubs.pdf", height = 6, width = 10)					  
plot_grid(plot1,plot2)
dev.off()

pdf("NeuronSpecific_subcellLocalisation.pdf", width = 15)
EWCE::plot_ctd(ctd = ctd,
                      level = 2,
                      genes = c("PRKAR1B","CADPS2"),
                      metric = "mean_exp")
EWCE::plot_ctd(ctd = ctd,
                      level = 2,
                      genes = c("PRKAR1B","CADPS2"),
                      metric = "specificity")
dev.off()					  

magentaNames <- read.table("magentaNames.txt")

wgcnaLab <- read.table("/gpfs/mrc0/projects/Research_Project-T112069/Cardiff/Methylation/wgcnaAnnotations.txt", header = T)


unconditional_results <- EWCE::bootstrap_enrichment_test(
  sct_data = ctd,
  hits = magentaNames$x,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = 100,
  annotLevel = 1, no_cores = 20)


# 20 core(s) assigned as workers (12reserved).
# Retrieving all genes using: homologene.
# Retrieving all organisms available in homologene.
# Mapping species name: human
# Common name mapping found for human
# 1 organism identified from search: 9606
# Gene table with 19,129 rows retrieved.
# Returning all 19,129 genes from human.
# Standardising CellTypeDataset
# Converting to sparse matrix.
# Converting to sparse matrix.
# Checking gene list inputs.
# Retrieving all genes using: homologene.
# Retrieving all organisms available in homologene.
# Mapping species name: human
# Common name mapping found for human
# 1 organism identified from search: 9606
# Gene table with 19,129 rows retrieved.
# Returning all 19,129 genes from human.
# Standardising sct_data.
# Converting gene list input to standardised human genes.
# Running without gene size control.
# 673 hit genes remain after filtering.
# Computing summed proportions.
# Testing for enrichment in 7 cell types...
# Sorting results by p-value.
# Computing BH-corrected q-values.
# 2 significant cell type enrichment results @ q<0.05 :
#   
#   CellType annotLevel p fold_change sd_from_mean q
# 1 NonDA_neurons          1 0    1.164850     5.057850 0
# 2    DA_neurons          1 0    1.101451     3.208384 0

unconditional_results2 <- EWCE::bootstrap_enrichment_test(
  sct_data = ctd,
  hits = magentaNames$x,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = 100,
  annotLevel = 2, no_cores = 20)
  
  
  # 20 core(s) assigned as workers (12reserved).
# Retrieving all genes using: homologene.
# Retrieving all organisms available in homologene.
# Mapping species name: human
# Common name mapping found for human
# 1 organism identified from search: 9606
# Gene table with 19,129 rows retrieved.
# Returning all 19,129 genes from human.
# Standardising CellTypeDataset
# Converting to sparse matrix.
# Converting to sparse matrix.
# Checking gene list inputs.
# Retrieving all genes using: homologene.
# Retrieving all organisms available in homologene.
# Mapping species name: human
# Common name mapping found for human
# 1 organism identified from search: 9606
# Gene table with 19,129 rows retrieved.
# Returning all 19,129 genes from human.
# Standardising sct_data.
# Converting gene list input to standardised human genes.
# Running without gene size control.
# 673 hit genes remain after filtering.
# Computing summed proportions.
# Testing for enrichment in 68 cell types...
# Sorting results by p-value.
# Computing BH-corrected q-values.
# 10 significant cell type enrichment results @ q<0.05 :
#   
#   CellType annotLevel p fold_change sd_from_mean q
# 1               Ex_POSTN          2 0    1.302864     7.427478 0
# 2               Ex_OPRD1          2 0    1.314015     6.279050 0
# 3         Ex_LAMP5_NTNG2          2 0    1.234171     5.648140 0
# 4               Ex_SATB2          2 0    1.223266     5.047599 0
# 5              CALB1_GEM          2 0    1.196823     4.603890 0
# 6               Ex_MYO5B          2 0    1.138018     4.122592 0
# 7  Inh_PRLR_RP11_384J4_2          2 0    1.126418     3.830443 0
# 8               Inh_SIX3          2 0    1.129655     3.474751 0
# 9         Inh_PAX5_CCBE1          2 0    1.101558     3.207944 0
# 10            Inh_IGFBP5          2 0    1.130331     3.000949 0


pdf("magentaAnalysis/EWCE_HigRes.pdf",width = 15)
EWCE::ewce_plot(unconditional_results2$results, 
                             mtc_method = "BH") 
dev.off()

pdf("magentaAnalysis/EWCE_LowRes.pdf",width = 10)
EWCE::ewce_plot(unconditional_results$results, 
                             mtc_method = "BH") 
dev.off()
