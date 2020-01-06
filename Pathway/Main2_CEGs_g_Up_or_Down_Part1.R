
# This code is used to obtain the PCC for control (normal) groups (i.e. ctrl) up- and down-upregulated, and normal groups, respectively.

setwd("Your current path")
rm(list = ls())

# load the required packages
# library(CancerSubtypes)
library(caret)
library(Hmisc) # for cor.test
library(BSDA) # for z.test
library(DMwR)
library(ROSE)
#################### Step 1: data preparation
### file path for counts stage data
# Counts_Stage_path <- "~/expression/workspace_one/express_COUNT_Stage_new"
### file path for FPKM stage data
FPKM_Stage_path <- "~/expression/workspace_one/express_FPKM_Stage_new"
# dis_method <- "cos" # Euclidean, cos
dis_method <- "Euclidean" # Euclidean, cos
PCC_cutoff=0.50
pvalue_cutoff=0.05
delta_PCC_cutoff <- 0.05
# split_rate_ref <- 0.3 # 30% as the reference set
split_rate_train <- 0.7 # 70% as the traning set

### some variables which are needed to be declared in advance
# the specific sampling methods
# smp_method <- c("downSample","upSample","smote")
smp_method <- c("smote")
classifier_method <- c("nb","treebag","C5.0","rf","rFerns","wsrf")


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

# Function calculate_PCC_network: Conduct the PCC for project using traning set
calculate_PCC_network <- function(Data_DEGs_topk_stage,train_sample_names_stage,if_abs=TRUE){
  # inuput:
  # output
  # note: 
  ### output variables include:
  ### 1) PCC Network: PCC_cor_stage_upper
  ### 2) significant location: location_intersect_stage
  ### 3) representive vector: PCC_stage_vector
  ####################
  ### e) construct PCC network 
  
  # rcorr(x, type="pearson")
  rcorr_stage <- rcorr(t(Data_DEGs_topk_stage[,train_sample_names_stage]),type = "pearson")
  PCC_cor_stage <- rcorr_stage$r
  diag(PCC_cor_stage) <- 0 
  
  # the correlation is set to be positive
  if(if_abs==TRUE){
    PCC_cor_stage <- abs(PCC_cor_stage) 
  }
  # the elements who are NA is set to be 0
  PCC_cor_stage[is.na(PCC_cor_stage)] <- 0 
  PCC_cor_stage_upper <- get_upper_tri(PCC_cor_stage)
  # diag(PCC_cor_stage_upper) <- NA # comment it
  # PCC_stage_vector <- PCC_cor_stage_upper[which(is.na(PCC_cor_stage_upper)!=TRUE)] #comment this line of code
  PCC_stage_vector <- as.vector(PCC_cor_stage_upper)  # add this line of code
  # PCC_stage_vector <- as.vector(PCC_cor_stage_upper[which(is.na(PCC_cor_stage_upper)!=TRUE)])  # add this line of code and modify it.
  
  # PCC_stage_vector # The return value
  
  # get the observation number
  obs_cor_stage <- rcorr_stage$n
  obs_cor_stage <- get_upper_tri(obs_cor_stage)
  
  # get the pvalues
  pvalue_cor_stage <- rcorr_stage$P
  pvalue_cor_stage[is.na(pvalue_cor_stage)] <- 1 # assign a quite large value for the na location
  pvalue_cor_stage_upper <- get_upper_tri(pvalue_cor_stage)
  
  # set the return value
  PCC_stage_vector # The return value
}

# Function getSampleFea: get the feature vector for the newsample
getSampleFea <- function(newsample_df,Final_PCC_Stage,Data_DEGs_topk_stage,vec_sigma_ref,if_abs=TRUE){
  # make sure the newsample has the similiar genes with the Data_DEGs_topk_stage
  newsample <- newsample_df[rownames(Data_DEGs_topk_stage),1]
  Data_DEGs_topk_stage_U_s <- cbind(Data_DEGs_topk_stage,newsample)
  
  # the merged sample set: R U {s}
  rcorr_stage <- rcorr(t(Data_DEGs_topk_stage_U_s),type = "pearson")
  PCC_cor_stage <- rcorr_stage$r
  diag(PCC_cor_stage) <- 0 
  # the correlation is set to be positive
  if(if_abs==TRUE){
    PCC_cor_stage <- abs(PCC_cor_stage) 
  }
  
  # the elements who are NA is set to be 0
  PCC_cor_stage[is.na(PCC_cor_stage)] <- 0 
  PCC_cor_stage_upper <- get_upper_tri(PCC_cor_stage)
  # diag(PCC_cor_stage_upper) <- NA #comment it
  # PCC_stage_vector <- PCC_cor_stage_upper[which(is.na(PCC_cor_stage_upper)!=TRUE)] #comment this line of code
  PCC_stage_vector <- as.vector(PCC_cor_stage_upper)  # add this line of code
  # PCC_stage_vector <- as.vector(PCC_cor_stage_upper[which(is.na(PCC_cor_stage_upper)!=TRUE)])  # add this line of code and modify it.
  
  
  # get the observation number
  obs_cor_stage <- rcorr_stage$n
  obs_cor_stage <- get_upper_tri(obs_cor_stage)
  
  # get the pvalues
  pvalue_cor_stage <- rcorr_stage$P
  pvalue_cor_stage[is.na(pvalue_cor_stage)] <- 1 # assign a quite large values
  pvalue_cor_stage_upper <- get_upper_tri(pvalue_cor_stage)
  
  # get the final feature for newsample
  # PCC_stage_vector <- (PCC_cor_stage_upper/vec_sigma_ref)[Final_PCC_Stage]
  PCC_stage_vector <- (PCC_stage_vector/vec_sigma_ref)[Final_PCC_Stage]
  PCC_stage_vector # returen value
}


FSbyVar <- function(Data, cut.type = "topk", value){
  vars = apply(Data, 1, var)
  feature_num = length(vars)
  hist(vars, breaks = feature_num * 0.1, col = "red", main = "Expression (Variance) distribution", 
       xlab = "The Variance of feature")
  if (cut.type == "topk") {
    index = sort(vars, decreasing = TRUE, index.return = TRUE)
    if (value > nrow(Data)) {
      value = nrow(Data)
      cat("Warning: the feature selection number is beyond the original feature numnber")
    }
    cutoff = index$x[value]
    abline(v = cutoff, col = "blue", lty = 5, lwd = 1.5)
    index = index$ix[1:value]
    selectData = Data[index, ]
  }
  if (cut.type == "cutoff") {
    abline(v = value, col = "blue", lty = 5, lwd = 1.5)
    index = which(vars > value)
    selectData = Data[index, ]
  }
  selectData
}

data.normalization <- function(Data, type = "feature_Median", log2 = FALSE){
  if (log2) {
    data = log2(Data + 1)
  }
  else {
    data = Data
  }
  if (type == "feature_Median") {
    result = sweep(data, 1, apply(data, 1, function(x) median(x, 
                                                              na.rm = TRUE)))
  }
  else if (type == "feature_Mean") {
    result = sweep(data, 1, apply(data, 1, function(x) mean(x, 
                                                            na.rm = TRUE)))
  }
  else if (type == "feature_zscore") {
    var_row = apply(data, 1, var)
    index = which(var_row < 1e-10)
    if (length(index) > 0) {
      data = data[-index, ]
      cat("The features with the zero variance have been removed.")
    }
    result = t(scale(t(data)))
  }
  else if (type == "sample_zscore") {
    var_col = apply(data, 2, var)
    index = which(var_col < 1e-10)
    if (length(index) > 0) {
      data = data[-index, ]
      cat("The samples with the zero variance have been removed.")
    }
    result = scale(data)
  }
  result
}

#################### Step 2: Extract the fpkm data
file_name <- list.files(FPKM_Stage_path,pattern = "RData$")
projects  <- sub("(.*)\\..*$", "\\1", file_name)
for (project in projects[2:9]) {
  cat("The work for ",project," has started!\n")
  #################### Step 2.1: data preparation
  ####################
  ### load the fpkm stage data for project
  load(paste0(FPKM_Stage_path,'/',project,'.RData'))
  
  ### load the DEGs for project
  load(paste0("./DEGs_Ori/",project,"_DEGs_Stages_df.Rdata"))

  Data_FPKM_all_labels <- c(rep("Stage I",dim(Data_FPKM_i)[2]),
                            rep("Stage II",dim(Data_FPKM_ii)[2]),
                            rep("Stage III",dim(Data_FPKM_iii)[2]),
                            rep("Stage IV",dim(Data_FPKM_iv)[2]))
  
  
  Data_FPKM_all_tissue_name <- c(colnames(Data_FPKM_i),
                                 colnames(Data_FPKM_ii),
                                 colnames(Data_FPKM_iii),
                                 colnames(Data_FPKM_iv)
  )
  
  ### 1) For Stage I VS Normal
  # Data_DEGs_topk_i <- log2(Data_FPKM_i[DEGs_Stage_I,] + 1)
  Data_DEGs_topk_i <- Data_FPKM_i[DEGs_Stage_I_All,]
  
  ### 2) For Stage II VS Normal
  # Data_DEGs_topk_ii <- log2(Data_FPKM_ii[DEGs_Stage_II,] + 1)
  Data_DEGs_topk_ii <- Data_FPKM_ii[DEGs_Stage_II_All,]
  
  ### 3) For Stage III VS Normal
  # Data_DEGs_topk_iii <- log2(Data_FPKM_iii[DEGs_Stage_III,] + 1)
  Data_DEGs_topk_iii <- Data_FPKM_iii[DEGs_Stage_III_All,]
  
  ### 4) For Stage IV VS Normal
  # Data_DEGs_topk_iv <- log2(Data_FPKM_iv[DEGs_Stage_IV,] + 1)
  Data_DEGs_topk_iv <- Data_FPKM_iv[DEGs_Stage_IV_All,]
  
  cat("The PFKM data processing has completed!\n")
  ####################
  #################### Step 2.2: Obtain the spliting Reference + training  + testing tissues
  ### split the sample for each stage of project into three parts: Reference + training  + testing
  ###    note: please save the sample names for each stage
  # set.seed(256)
  ### 30% reference + 40% traning + 30% testing
  # length_tissues <- dim(Data_DEGs_topk_n)[2]
  # index_tissues <- sample(length_tissues)
  
  ### Obtain the Reference_labels
  Reference_labels <- rep("Normal",dim(datan)[2])
  Reference_tissues <- colnames(datan)

  #################### Step 2.3: Construct the reference network for each stage
  ####################
  ### we should define a fucntion which can build the pcc network
  cat("The calculation for PCC reference network of Normal starts!\n")
  
  #################### Step 2.3.1: Construct the reference network for stage i
  Data_DEGs_topk_n_i <- datan[rownames(Data_DEGs_topk_i),]
  # normalization
  Data_DEGs_topk_n_i <- data.normalization(Data_DEGs_topk_n_i,type="feature_Median",log2=TRUE)
  # Data_DEGs_topk_n_i <- log2(Data_DEGs_topk_n_i + 1)
  Ref_PCC_Stage_N_I <- calculate_PCC_network(Data_DEGs_topk_n_i,Reference_tissues,TRUE)
  
  #################### Step 2.3.2: Construct the reference network for stage ii
  Data_DEGs_topk_n_ii <- datan[rownames(Data_DEGs_topk_ii),]
  # normalization
  Data_DEGs_topk_n_ii <- data.normalization(Data_DEGs_topk_n_ii,type="feature_Median",log2=TRUE)
  # Data_DEGs_topk_n_ii <- log2(Data_DEGs_topk_n_ii + 1)
  Ref_PCC_Stage_N_II <- calculate_PCC_network(Data_DEGs_topk_n_ii,Reference_tissues,TRUE)
  
  #################### Step 2.3.3: Construct the reference network for stage iii
  Data_DEGs_topk_n_iii <- datan[rownames(Data_DEGs_topk_iii),]
  # normalization
  Data_DEGs_topk_n_iii <- data.normalization(Data_DEGs_topk_n_iii,type="feature_Median",log2=TRUE)
  # Data_DEGs_topk_n_iii <- log2(Data_DEGs_topk_n_iii + 1)
  Ref_PCC_Stage_N_III <- calculate_PCC_network(Data_DEGs_topk_n_iii,Reference_tissues,TRUE)
  
  #################### Step 2.3.4: Construct the reference network for stage iv
  Data_DEGs_topk_n_iv <- datan[rownames(Data_DEGs_topk_iv),]
  # normalization
  Data_DEGs_topk_n_iv <- data.normalization(Data_DEGs_topk_n_iv,type="feature_Median",log2=TRUE)
  # Data_DEGs_topk_n_iv <- log2(Data_DEGs_topk_n_iv + 1)
  Ref_PCC_Stage_N_IV <- calculate_PCC_network(Data_DEGs_topk_n_iv,Reference_tissues,TRUE)
  
  #################### Step 2.2: Conduct the PCC for each Stage of project
  ####################
  #################### Step 2.2.1: Construct the co-expression network networks for stage i
  # normalization
  Data_DEGs_topk_i <- data.normalization(Data_DEGs_topk_i,type="feature_Median",log2=TRUE)
  # Data_DEGs_topk_i <- log2(Data_DEGs_topk_i + 1)
  PCC_Stage_I_up <- calculate_PCC_network(Data_DEGs_topk_i[DEGs_Stage_I_upregulated,],colnames(Data_DEGs_topk_i),TRUE)
  PCC_Stage_I_down <- calculate_PCC_network(Data_DEGs_topk_i[DEGs_Stage_I_downregulated,],colnames(Data_DEGs_topk_i),TRUE)
  
  
  #################### Step 2.2.2: Construct the co-expression network networks for stage ii
  # normalization
  Data_DEGs_topk_ii <- data.normalization(Data_DEGs_topk_ii,type="feature_Median",log2=TRUE)
  # Data_DEGs_topk_ii <- log2(Data_DEGs_topk_ii + 1)
  PCC_Stage_II_up <- calculate_PCC_network(Data_DEGs_topk_ii[DEGs_Stage_II_upregulated,],colnames(Data_DEGs_topk_ii),TRUE)
  PCC_Stage_II_down <- calculate_PCC_network(Data_DEGs_topk_ii[DEGs_Stage_II_downregulated,],colnames(Data_DEGs_topk_ii),TRUE)

  
  
  #################### Step 2.2.3: Construct the co-expression network networks for stage iii
  # normalization
  Data_DEGs_topk_iii <- data.normalization(Data_DEGs_topk_iii,type="feature_Median",log2=TRUE)
  # Data_DEGs_topk_iii <- log2(Data_DEGs_topk_iii + 1)
  PCC_Stage_III_up <- calculate_PCC_network(Data_DEGs_topk_iii[DEGs_Stage_III_upregulated,],colnames(Data_DEGs_topk_iii),TRUE)
  PCC_Stage_III_down <- calculate_PCC_network(Data_DEGs_topk_iii[DEGs_Stage_III_downregulated,],colnames(Data_DEGs_topk_iii),TRUE)

  
  
  #################### Step 2.2.4: Construct the co-expression network networks for stage iv
  # normalization
  Data_DEGs_topk_iv <- data.normalization(Data_DEGs_topk_iv,type="feature_Median",log2=TRUE)
  # Data_DEGs_topk_iv <- log2(Data_DEGs_topk_iv + 1)
  PCC_Stage_IV_up <- calculate_PCC_network(Data_DEGs_topk_iv[DEGs_Stage_IV_upregulated,],colnames(Data_DEGs_topk_iv),TRUE)
  PCC_Stage_IV_down <- calculate_PCC_network(Data_DEGs_topk_iv[DEGs_Stage_IV_downregulated,],colnames(Data_DEGs_topk_iv),TRUE)
  
  
  #################### Step 3 Save the variables
  ####################
  save(Ref_PCC_Stage_N_I, Ref_PCC_Stage_N_II, Ref_PCC_Stage_N_III, Ref_PCC_Stage_N_IV,
       PCC_Stage_I_up, PCC_Stage_II_up, PCC_Stage_III_up, PCC_Stage_IV_up,
       PCC_Stage_I_down, PCC_Stage_II_down, PCC_Stage_III_down, PCC_Stage_IV_down,
       file = paste0("./CEGs_Part1/",project,"_CEGs_Stages_Matrix_Part1.Rdata"))
  
  cat("The significantly differenatilal location for",project,"has completed!\n")
}





