
# This code is used to Otain the specific location PCC matrix for up- and down-upregulated, and normal groups, respectively.

setwd("Your current path")
# rm(list = ls())

PCC_cutoff=c(0.60, 0.70, 0.80)
FPKM_Stage_path <- "~/expression/workspace_one/express_FPKM_Stage_new"

### Some functions which are needed to be declared in advance
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

#################### Step 2: Extract the fpkm data
file_name <- list.files(FPKM_Stage_path,pattern = "RData$")
projects  <- sub("(.*)\\..*$", "\\1", file_name)

for(pcc_cutoff_item in PCC_cutoff){
  
  cat("The PCC_cutoff for",pcc_cutoff_item,"starts!\n")
  
  for (project in projects[2:9]) {
    cat("The work for ",project," has started!\n")
    load(paste0("./CEGs_Part1/",project,"_CEGs_Stages_Matrix_Part1.Rdata"))
    load(paste0("./DEGs_Ori/",project,"_DEGs_Stages_df.Rdata"))
    #################### Step 3: Find the significantly differential co-expression edges
    ####################
    #################### Step 2.3.1: To determine the significantly differential  co-expression edges for stage i
    ########################## Here, we get the result for stage i
    # The DEGs for stage i: DEGs_Stage_I
    # The PCC for normal tissues: Ref_PCC_Stage_N_I
    # The PCC for stage i: Stage_I_PCC_Net
    # And the corresponding the SDCN: sig_loc_stage_i
    
    ### compute the co-expression location (gene pairs) for normal tissues
    CEG_Loc_N_I <- which(Ref_PCC_Stage_N_I >= pcc_cutoff_item)
    ### compute the co-expression location (gene pairs) for stage i
    CEG_Loc_Stage_I_up <- which(PCC_Stage_I_up >= pcc_cutoff_item)
    CEG_Loc_Stage_I_down <- which(PCC_Stage_I_down >= pcc_cutoff_item)
    
    #################### Step b: construct the matrix form for normal and stage i
    ### for normal
    nlen_DEGs_Stage_I <- length(DEGs_Stage_I_All)
    CEG_Matrix_N_I<- matrix(0,nrow = nlen_DEGs_Stage_I,ncol = nlen_DEGs_Stage_I)
    colnames(CEG_Matrix_N_I) <- DEGs_Stage_I_All
    rownames(CEG_Matrix_N_I) <- DEGs_Stage_I_All
    CEG_Matrix_N_I[CEG_Loc_N_I] <- 100
    CEG_Matrix_N_I <- get_upper_tri(CEG_Matrix_N_I)
    
    ### for stage i
    # for up regulated
    nlen_DEGs_Stage_I_up <- length(DEGs_Stage_I_upregulated)
    CEG_Matrix_Stage_I_up <- matrix(0,nrow = nlen_DEGs_Stage_I_up,ncol = nlen_DEGs_Stage_I_up)
    colnames(CEG_Matrix_Stage_I_up) <- DEGs_Stage_I_upregulated
    rownames(CEG_Matrix_Stage_I_up) <- DEGs_Stage_I_upregulated
    CEG_Matrix_Stage_I_up[CEG_Loc_Stage_I_up] <- 100
    CEG_Matrix_Stage_I_up <- get_upper_tri(CEG_Matrix_Stage_I_up)
    
    # for down regulated
    nlen_DEGs_Stage_I_down <- length(DEGs_Stage_I_downregulated)
    CEG_Matrix_Stage_I_down <- matrix(0,nrow = nlen_DEGs_Stage_I_down,ncol = nlen_DEGs_Stage_I_down)
    colnames(CEG_Matrix_Stage_I_down) <- DEGs_Stage_I_downregulated
    rownames(CEG_Matrix_Stage_I_down) <- DEGs_Stage_I_downregulated
    CEG_Matrix_Stage_I_down[CEG_Loc_Stage_I_down] <- 100
    CEG_Matrix_Stage_I_down <- get_upper_tri(CEG_Matrix_Stage_I_down)
    # up
    ### save this matrix: CEG_Matrix_N_I and CEG_Matrix_Stage_I
    
    ####################
    #################### Step 2.3.2: To determine the significantly differential  co-expression edges for stage ii
    ########################## Here, we get the result for stage ii
    # The DEGs for stage ii: DEGs_Stage_II
    # The PCC for normal tissues: Ref_PCC_Stage_N_I
    # The PCC for stage ii: Stage_II_PCC_Net
    # And the corresponding the SDCN: sig_loc_stage_ii
    
    ### compute the co-expression location (gene pairs) for normal tissues
    CEG_Loc_N_II <- which(Ref_PCC_Stage_N_II >= pcc_cutoff_item)
    ### compute the co-expression location (gene pairs) for stage ii
    CEG_Loc_Stage_II_up <- which(PCC_Stage_II_up >= pcc_cutoff_item)
    CEG_Loc_Stage_II_down <- which(PCC_Stage_II_down >= pcc_cutoff_item)
    #################### Step b: construct the matrix form for normal and stage ii
    ### for normal
    nlen_DEGs_Stage_II <- length(DEGs_Stage_II_All)
    CEG_Matrix_N_II<- matrix(0,nrow = nlen_DEGs_Stage_II,ncol = nlen_DEGs_Stage_II)
    colnames(CEG_Matrix_N_II) <- DEGs_Stage_II_All
    rownames(CEG_Matrix_N_II) <- DEGs_Stage_II_All
    CEG_Matrix_N_II[CEG_Loc_N_II] <- 100
    CEG_Matrix_N_II <- get_upper_tri(CEG_Matrix_N_II)
    
    ### for stage ii
    # for up-regulated
    nlen_DEGs_Stage_II_up <- length(DEGs_Stage_II_upregulated)
    CEG_Matrix_Stage_II_up <- matrix(0,nrow = nlen_DEGs_Stage_II_up,ncol = nlen_DEGs_Stage_II_up)
    colnames(CEG_Matrix_Stage_II_up) <- DEGs_Stage_II_upregulated
    rownames(CEG_Matrix_Stage_II_up) <- DEGs_Stage_II_upregulated
    CEG_Matrix_Stage_II_up[CEG_Loc_Stage_II_up] <- 100
    CEG_Matrix_Stage_II_up <- get_upper_tri(CEG_Matrix_Stage_II_up)
    
    # for down-regulated
    nlen_DEGs_Stage_II_down <- length(DEGs_Stage_II_downregulated)
    CEG_Matrix_Stage_II_down <- matrix(0,nrow = nlen_DEGs_Stage_II_down,ncol = nlen_DEGs_Stage_II_down)
    colnames(CEG_Matrix_Stage_II_down) <- DEGs_Stage_II_downregulated
    rownames(CEG_Matrix_Stage_II_down) <- DEGs_Stage_II_downregulated
    CEG_Matrix_Stage_II_down[CEG_Loc_Stage_II_down] <- 100
    CEG_Matrix_Stage_II_down <- get_upper_tri(CEG_Matrix_Stage_II_down)
    # up
    ### save this matrix: CEG_Matrix_N_II and CEG_Matrix_Stage_II
    
    ####################
    #################### Step 2.3.3: To determine the significantly differential  co-expression edges for stage iii
    ########################## Here, we get the result for stage iii
    # The DEGs for stage iii: DEGs_Stage_III
    # The PCC for normal tissues: Ref_PCC_Stage_N_I
    # The PCC for stage iii: Stage_III_PCC_Net
    # And the corresponding the SDCN: sig_loc_stage_iii
    
    ### compute the co-expression location (gene pairs) for normal tissues
    CEG_Loc_N_III <- which(Ref_PCC_Stage_N_III >= pcc_cutoff_item)
    ### compute the co-expression location (gene pairs) for stage III
    CEG_Loc_Stage_III_up <- which(PCC_Stage_III_up >= pcc_cutoff_item)
    CEG_Loc_Stage_III_down <- which(PCC_Stage_III_down >= pcc_cutoff_item)
    #################### Step b: construct the matrix form for normal and stage III
    ### for normal
    nlen_DEGs_Stage_III <- length(DEGs_Stage_III_All)
    CEG_Matrix_N_III<- matrix(0,nrow = nlen_DEGs_Stage_III,ncol = nlen_DEGs_Stage_III)
    colnames(CEG_Matrix_N_III) <- DEGs_Stage_III_All
    rownames(CEG_Matrix_N_III) <- DEGs_Stage_III_All
    CEG_Matrix_N_III[CEG_Loc_N_III] <- 100
    CEG_Matrix_N_III <- get_upper_tri(CEG_Matrix_N_III)
    
    ### for stage III
    # for up-regulated
    nlen_DEGs_Stage_III_up <- length(DEGs_Stage_III_upregulated)
    CEG_Matrix_Stage_III_up <- matrix(0,nrow = nlen_DEGs_Stage_III_up,ncol = nlen_DEGs_Stage_III_up)
    colnames(CEG_Matrix_Stage_III_up) <- DEGs_Stage_III_upregulated
    rownames(CEG_Matrix_Stage_III_up) <- DEGs_Stage_III_upregulated
    CEG_Matrix_Stage_III_up[CEG_Loc_Stage_III_up] <- 100
    CEG_Matrix_Stage_III_up <- get_upper_tri(CEG_Matrix_Stage_III_up)
    
    # for down-regulated
    nlen_DEGs_Stage_III_down <- length(DEGs_Stage_III_downregulated)
    CEG_Matrix_Stage_III_down <- matrix(0,nrow = nlen_DEGs_Stage_III_down,ncol = nlen_DEGs_Stage_III_down)
    colnames(CEG_Matrix_Stage_III_down) <- DEGs_Stage_III_downregulated
    rownames(CEG_Matrix_Stage_III_down) <- DEGs_Stage_III_downregulated
    CEG_Matrix_Stage_III_down[CEG_Loc_Stage_III_down] <- 100
    CEG_Matrix_Stage_III_down <- get_upper_tri(CEG_Matrix_Stage_III_down)
    # up
    ### save this matrix: CEG_Matrix_N_III and CEG_Matrix_Stage_III  
    # III
    
    
    ####################
    #################### Step 2.3.4: To determine the significantly differential  co-expression edges for stage iv
    ########################## Here, we get the result for stage iv
    # The DEGs for stage iv: DEGs_Stage_IV
    # The PCC for normal tissues: Ref_PCC_Stage_N_I
    # The PCC for stage iv: Stage_IV_PCC_Net
    # And the corresponding the SDCN: sig_loc_stage_iv
    
    
    ### compute the co-expression location (gene pairs) for normal tissues
    CEG_Loc_N_IV <- which(Ref_PCC_Stage_N_IV >= pcc_cutoff_item)
    ### compute the co-expression location (gene pairs) for stage IV
    CEG_Loc_Stage_IV_up <- which(PCC_Stage_IV_up >= pcc_cutoff_item)
    CEG_Loc_Stage_IV_down <- which(PCC_Stage_IV_down >= pcc_cutoff_item)
    #################### Step b: construct the matrix form for normal and stage IV
    ### for normal
    nlen_DEGs_Stage_IV <- length(DEGs_Stage_IV_All)
    CEG_Matrix_N_IV<- matrix(0,nrow = nlen_DEGs_Stage_IV,ncol = nlen_DEGs_Stage_IV)
    colnames(CEG_Matrix_N_IV) <- DEGs_Stage_IV_All
    rownames(CEG_Matrix_N_IV) <- DEGs_Stage_IV_All
    CEG_Matrix_N_IV[CEG_Loc_N_IV] <- 100
    CEG_Matrix_N_IV <- get_upper_tri(CEG_Matrix_N_IV)
    
    ### for stage IV
    # for up-regulated
    nlen_DEGs_Stage_IV_up <- length(DEGs_Stage_IV_upregulated)
    CEG_Matrix_Stage_IV_up <- matrix(0,nrow = nlen_DEGs_Stage_IV_up,ncol = nlen_DEGs_Stage_IV_up)
    colnames(CEG_Matrix_Stage_IV_up) <- DEGs_Stage_IV_upregulated
    rownames(CEG_Matrix_Stage_IV_up) <- DEGs_Stage_IV_upregulated
    CEG_Matrix_Stage_IV_up[CEG_Loc_Stage_IV_up] <- 100
    CEG_Matrix_Stage_IV_up <- get_upper_tri(CEG_Matrix_Stage_IV_up)
    
    # for down-regulated
    nlen_DEGs_Stage_IV_down <- length(DEGs_Stage_IV_downregulated)
    CEG_Matrix_Stage_IV_down <- matrix(0,nrow = nlen_DEGs_Stage_IV_down,ncol = nlen_DEGs_Stage_IV_down)
    colnames(CEG_Matrix_Stage_IV_down) <- DEGs_Stage_IV_downregulated
    rownames(CEG_Matrix_Stage_IV_down) <- DEGs_Stage_IV_downregulated
    CEG_Matrix_Stage_IV_down[CEG_Loc_Stage_IV_down] <- 100
    CEG_Matrix_Stage_IV_down <- get_upper_tri(CEG_Matrix_Stage_IV_down)
    # up
    ### save this matrix: CEG_Matrix_N_IV and CEG_Matrix_Stage_IV  
    # IV
    
    #################### Step 4: Save the resulsts
    ####################
    save(CEG_Matrix_N_I, CEG_Matrix_N_II, CEG_Matrix_N_III, CEG_Matrix_N_IV,
         CEG_Matrix_Stage_I_up, CEG_Matrix_Stage_II_up, CEG_Matrix_Stage_III_up, CEG_Matrix_Stage_IV_up,
         CEG_Matrix_Stage_I_down, CEG_Matrix_Stage_II_down, CEG_Matrix_Stage_III_down, CEG_Matrix_Stage_IV_down,
         file = paste0("./CEGs_Part2/",pcc_cutoff_item,"/",project,"_CEGs_Stages_Matrix.Rdata"))
    cat("The significantly differenatilal location for",project,"has completed!\n")
  }
}

cat("The whole task has completed!\n")









