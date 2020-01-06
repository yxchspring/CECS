
# This code is used to conduct differentially expression analysis (normal VS each stage)

setwd("Your current path")
rm(list = ls())

#################### Step1: data preparation
### file path for counts stage data
Counts_Stage_path <- "~/expression/workspace_one/express_COUNT_Stage_new"

#################### Step 2: Conduct statistic analysis using edgeR
### 1) edgeR function
edgeR.de.test <- function(df1, df2) {
  #' Use edgeR to perform gene differential expression analysis.
  #'
  #' Require edgeR and dplyr to work.
  #' 
  #' @param df1 First data frame / matrix. Must be counts value!!
  #' @param df2 Second. The result comes as df2:df1.
  #browser()
  library(edgeR)
  df1 <- df1[order(rownames(df1)), ]
  df2 <- df2[order(rownames(df2)), ]
  group <- c(rep("control", ncol(df1)), rep("experiment", ncol(df2)))
  design <- model.matrix(~0+group)
  x <- cbind(df1, df2)
  y <- DGEList(counts=x, group=group)
  y <- y[rowSums(cpm(y) > 1) >= 2, , keep.lib.sizes=F]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design=design)
  ## logFC logCPM PValue
  et <- exactTest(y)$table
  colnames(et) <- c("log2_fold_change", "log2_CPM", "p_value")
  et$symbol <- rownames(et)
  rownames(et) <- NULL
  return(et)
}

### 2) Obtain the projects names
file_name <- list.files(Counts_Stage_path,pattern = "RData$")
projects  <- sub("(.*)\\..*$", "\\1", file_name)
for (project in projects[2:9]) {
  ### a) load the count stage data
  load(paste0(Counts_Stage_path,'/',project,'.RData'))
  ### b) edgeR starts!
  ### b1) Normal VS Stage I
  edgeR_Stage_I <- edgeR.de.test(datan,Data_FPKM_i)
  cat('The Stage I edgeR analysis has been completed for ',project,'\n')
  
  ### b1) Normal VS Stage II
  edgeR_Stage_II <- edgeR.de.test(datan,Data_FPKM_ii)
  cat('The Stage II edgeR analysis has been completed for ',project,'\n')
  
  ### b1) Normal VS Stage III
  edgeR_Stage_III <- edgeR.de.test(datan,Data_FPKM_iii)
  cat('The Stage III edgeR analysis has been completed for ',project,'\n')
  
  ### b1) Normal VS Stage IV
  edgeR_Stage_IV <- edgeR.de.test(datan,Data_FPKM_iv)
  cat('The Stage IV edgeR analysis has been completed for ',project,'\n')
  
  save(edgeR_Stage_I,edgeR_Stage_II,edgeR_Stage_III,edgeR_Stage_IV,
       file = paste0("./DEGs/",project,"_edgeR_Stages.Rdata"))
  cat("The edgeR analysis has been completed for ",project,"\n")
}