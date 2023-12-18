# Function to perform module trait associations
# tVec: Vector of labels for trait assocation, must be in form binary or quantitative
# tFile: Data frame with tVec labels as colnames (pheno file)
# mFile: module file containing module eigengenes
# rDat: "Raw" base gene expression / methylation matrix used to construct module



traitModCor <- function(tVec,tFile,mFile,rDat){
  message(paste("You've started the function, at least"))
  # Define your variables
  nGenes = ncol(rDat)
  nSamples= nrow(rDat)
  moduleColors = labels2colors(mFile$colors)
  MEs0= moduleEigengenes(rDat, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  
  message(paste("You've gotten to the pre-check"))
  # 
  if(!identical(rownames(rDat),rownames(MEs))){
    
    stop("Raw data matrix and module eigengenes are discordant")
    
  } 
  
  if(length(unique(tVec %in% colnames(tFile))) > 1 | unique(tVec %in% colnames(tFile)) == FALSE){
    
    stop("Label vector and trait file are discordant")
    
  }
  
  message(paste("You've gotten to the results file creation"))
  
  moduleTraitCor <- matrix(data = NA, ncol = length(tVec), nrow = ncol(MEs))
  colnames(moduleTraitCor)<- tVec
  rownames(moduleTraitCor) <- colnames(MEs)
  
  moduleTraitPvalue <- matrix(data = NA, ncol = length(tVec), nrow = ncol(MEs))
  colnames(moduleTraitPvalue)<-  tVec
  rownames(moduleTraitPvalue) <- colnames(MEs)

  message(paste("You've gotten to the correlation loop"))
  
  for(t in tVec){
    if(length(unique(tFile[,t])) > 2){ # if Binary run spearman cor
      for(m in colnames(MEs)){
        
        try(res <-cor.test(MEs[,m] , as.numeric(tFile[,t]),method = "pearson"), silent = TRUE)
        if(class(res) != "try-error")
          moduleTraitCor[m,t]<-as.numeric(res$estimate)
          moduleTraitPvalue[m,t]<-res$p.value
      }
      
    } else { # do the same as above but if quantitative then run different pearson cor
      for(m in colnames(MEs)){
      try(res <-cor.test(MEs[,m] , as.numeric(tFile[,t]),method = "spearman"), silent = TRUE)
      if(class(res) != "try-error")
        moduleTraitCor[m,t]<-as.numeric(res$estimate)
        moduleTraitPvalue[m,t]<-res$p.value
      } 
    }
    
  message(paste("Testing",t))
  }
  assign("moduleTraitCor", moduleTraitCor, envir = parent.frame() )
  assign("moduleTraitPvalue", moduleTraitPvalue, envir = parent.frame() )
}




