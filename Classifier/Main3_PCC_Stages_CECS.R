
# This code is used to conduct the co-expression networks(i.e. PCC) and perturbed PCC for each Stage and complete the classification task

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
PCC_cutoff=0.7
pvalue_cutoff=0.05
delta_PCC_cutoff <- 0.05
split_rate_ref <- 0.3 # 30% as the reference set
split_rate_train <- 0.4 # 70% as the traning set

### some variables which are needed to be declared in advance
# the specific sampling methods
# smp_method <- c("downSample","upSample","smote")
smp_method <- c("smote")
# classifier_method <- c("nb","treebag","C5.0","rf","rFerns","wsrf")
classifier_method <- c("C5.0")

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
  diag(PCC_cor_stage_upper) <- NA
  PCC_stage_vector <- PCC_cor_stage_upper[which(is.na(PCC_cor_stage_upper)!=TRUE)] #comment this line of code
  # PCC_stage_vector <- as.vector(PCC_cor_stage_upper)  # add this line of code
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
  diag(PCC_cor_stage_upper) <- NA
  # PCC_stage_vector <- PCC_cor_stage_upper[which(is.na(PCC_cor_stage_upper)!=TRUE)] #comment this line of code
  PCC_stage_vector <- as.vector(PCC_cor_stage_upper)  # add this line of code
  
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
for (project in projects[c(2:3,5:7,9)]) {
  cat("The work for ",project," has started!\n")
  #################### Step 2.1: data preparation
  ####################
  ### load the fpkm stage data for project
  load(paste0(FPKM_Stage_path,'/',project,'.RData'))
  
  ### load the DEGs for project
  load(paste0("./DEGs/",project,"_DEGs_Stages_df.Rdata"))
  
  # Data_FPKM_all <- cbind(datan[order(rownames(datan)),],
  #                        Data_FPKM_i[rownames(Data_FPKM_i),],
  #                        Data_FPKM_ii[rownames(Data_FPKM_ii),],
  #                        Data_FPKM_iii[rownames(Data_FPKM_iii),],
  #                        Data_FPKM_iv[rownames(Data_FPKM_iv),]
  # )
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
  ### Obatin the top M DEGs
  gene_num_small <- min(length(DEGs_Stage_I),length(DEGs_Stage_II),
                        length(DEGs_Stage_III),length(DEGs_Stage_IV))
  # Data_DEGs_topk_i=FSbyVar(Data_FPKM_i[DEGs_Stage_I,], cut.type="topk",value=1000)
  Data_DEGs_topk_i=FSbyVar(Data_FPKM_i[DEGs_Stage_I,], cut.type="topk",value=gene_num_small)
  Data_DEGs_topk_i <- na.omit(Data_DEGs_topk_i) # delete the rows whose values are all NAs
  # normalization
  Data_DEGs_topk_i <- data.normalization(Data_DEGs_topk_i,type="feature_Median",log2=TRUE)
  # Data_DEGs_topk_i <- log(Data_DEGs_topk_i + 1)
  norm_para_top_i <- apply(Data_DEGs_topk_i, 1, function(x) median(x,na.rm = TRUE))
  
  ### 2) For Stage II VS Normal
  # Data_DEGs_topk_ii=FSbyVar(Data_FPKM_ii[DEGs_Stage_II,], cut.type="topk",value=1000)
  Data_DEGs_topk_ii=FSbyVar(Data_FPKM_ii[DEGs_Stage_II,], cut.type="topk",value=gene_num_small)
  Data_DEGs_topk_ii <- na.omit(Data_DEGs_topk_ii) # delete the rows whose values are all NAs
  Data_DEGs_topk_ii <- data.normalization(Data_DEGs_topk_ii,type="feature_Median",log2=TRUE)
  # Data_DEGs_topk_ii <- log(Data_DEGs_topk_ii + 1)
  norm_para_top_ii <- apply(Data_DEGs_topk_ii, 1, function(x) median(x,na.rm = TRUE))
  
  ### 3) For Stage III VS Normal
  # Data_DEGs_topk_iii=FSbyVar(Data_FPKM_iii[DEGs_Stage_III,], cut.type="topk",value=1000)
  Data_DEGs_topk_iii=FSbyVar(Data_FPKM_iii[DEGs_Stage_III,], cut.type="topk",value=gene_num_small)
  Data_DEGs_topk_iii <- na.omit(Data_DEGs_topk_iii) # delete the rows whose values are all NAs
  Data_DEGs_topk_iii <- data.normalization(Data_DEGs_topk_iii,type="feature_Median",log2=TRUE)
  # Data_DEGs_topk_iii <- log(Data_DEGs_topk_iii + 1)
  norm_para_top_iii <- apply(Data_DEGs_topk_iii, 1, function(x) median(x,na.rm = TRUE))
  
  ### 4) For Stage IV VS Normal
  # Data_DEGs_topk_iv=FSbyVar(Data_FPKM_iv[DEGs_Stage_IV,], cut.type="topk",value=1000)
  Data_DEGs_topk_iv=FSbyVar(Data_FPKM_iv[DEGs_Stage_IV,], cut.type="topk",value=gene_num_small)
  Data_DEGs_topk_iv <- na.omit(Data_DEGs_topk_iv) # delete the rows whose values are all NAs
  Data_DEGs_topk_iv <- data.normalization(Data_DEGs_topk_iv,type="feature_Median",log2=TRUE)
  # Data_DEGs_topk_iv <- log(Data_DEGs_topk_iv + 1)
  norm_para_top_iv <- apply(Data_DEGs_topk_iv, 1, function(x) median(x,na.rm = TRUE))
  
  
  ### summary again
  min_num_genes <- min(dim(Data_DEGs_topk_i)[1],dim(Data_DEGs_topk_ii)[1],
                       dim(Data_DEGs_topk_iii)[1],dim(Data_DEGs_topk_iv)[1])
  
  Data_DEGs_topk_i = Data_DEGs_topk_i[1:min_num_genes,]
  Data_DEGs_topk_ii = Data_DEGs_topk_ii[1:min_num_genes,]
  Data_DEGs_topk_iii = Data_DEGs_topk_iii[1:min_num_genes,]
  Data_DEGs_topk_iv = Data_DEGs_topk_iv[1:min_num_genes,]
  cat("The PFKM data processing has completed!\n")
  ####################
  #################### Step 2.2: Obtain the spliting Reference + training  + testing tissues
  ### split the sample for each stage of project into three parts: Reference + training  + testing
  ###    note: please save the sample names for each stage
  # set.seed(256)
  ### 30% reference + 40% traning + 30% testing
  # length_tissues <- dim(Data_DEGs_topk_n)[2]
  # index_tissues <- sample(length_tissues)

  
  ### Obtain the Ref_labels
  Ref_labels_split <- createDataPartition(Data_FPKM_all_labels, p=split_rate_ref, list=FALSE)
  Ref_labels <- Data_FPKM_all_labels[Ref_labels_split]
  Ref_tissues <- Data_FPKM_all_tissue_name[Ref_labels_split]
  ### Obtain the the remaining the labels and tissue_names except Ref_labels
  Data_FPKM_Train_Test_labels = Data_FPKM_all_labels[-Ref_labels_split]
  Data_FPKM_Train_Test_tissue_name = Data_FPKM_all_tissue_name[-Ref_labels_split]
  
  
  ### Obtain the Train_labels
  Train_labels_split <- createDataPartition(Data_FPKM_Train_Test_labels, p=split_rate_train/(1-split_rate_ref), list=FALSE)
  Train_labels <- Data_FPKM_Train_Test_labels[Train_labels_split]
  Train_tissues <- Data_FPKM_Train_Test_tissue_name[Train_labels_split]
  
  ### Obtain the Test_labels
  Test_labels <- Data_FPKM_Train_Test_labels[-Train_labels_split]
  Test_tissues <- Data_FPKM_Train_Test_tissue_name[-Train_labels_split]
  
  save(Ref_labels,Ref_tissues,Train_labels,Train_tissues, Test_labels,Test_tissues,
       file = paste0("./labels_tissues/",project,"_labels_tissues.Rdata"))
  cat("The save work for sampling method: for project: ",project," has been completed!\n")
  
  
  cat("The samples dividion has completed!\n")
  #################### Step 2.3: Construct the reference network for each stage
  ####################
  ### we should define a fucntion which can build the pcc network
  cat("The calculation for PCC reference network of Normal starts!\n")
  
  #################### Step 2.3.1: Construct the reference network for stage i
  Ref_PCC_Stage_N_I <- calculate_PCC_network(Data_DEGs_topk_i,Ref_tissues[which(Ref_labels=="Stage I")],TRUE)
  
  #################### Step 2.3.2: Construct the reference network for stage ii
  Ref_PCC_Stage_N_II <- calculate_PCC_network(Data_DEGs_topk_ii,Ref_tissues[which(Ref_labels=="Stage II")],TRUE)
  
  #################### Step 2.3.3: Construct the reference network for stage iii
  Ref_PCC_Stage_N_III <- calculate_PCC_network(Data_DEGs_topk_iii,Ref_tissues[which(Ref_labels=="Stage III")],TRUE)
  
  #################### Step 2.3.4: Construct the reference network for stage iv
  Ref_PCC_Stage_N_IV <- calculate_PCC_network(Data_DEGs_topk_iv,Ref_tissues[which(Ref_labels=="Stage IV")],TRUE)
  
  #################### Step 2.2: Conduct the PCC for each Stage of project using traning set
  ####################
  #################### Step 2.2.1: Construct the co-expression network (DCN) networks for stage i
  # Obatin the stage i data in tranining set
  train_tissue_stage_i <- Train_tissues[which(Train_labels=="Stage I")]
  row_num_stage_i = length(train_tissue_stage_i)
  col_num_stage_i = length(Ref_PCC_Stage_N_I)
  PCC_tr_stage_i <- matrix(NA,nrow = row_num_stage_i,ncol = col_num_stage_i*4)
  # Construct the PCC for reference network + one training sample
  cat("Start to calculate the PCC for stage i of traning set!\n")
  for (tr_stage_i in c(1:row_num_stage_i)) {
    # cat("Start to calculate the tissue: ",tr_stage_i, " of stage i of traning set!\n")
    # Data_DEGs_topk_stage <- cbind(Data_DEGs_topk_n_i,Data_DEGs_topk_i)
    PCC_tr_stage_i[tr_stage_i,1:col_num_stage_i] <- calculate_PCC_network(Data_DEGs_topk_i,
                                                                          c(Ref_tissues[which(Ref_labels=="Stage I")],train_tissue_stage_i[tr_stage_i]),TRUE)
    
    PCC_tr_stage_i[tr_stage_i,(col_num_stage_i+1):(col_num_stage_i*2)] <- calculate_PCC_network(cbind(Data_DEGs_topk_ii,Data_DEGs_topk_i),
                                                                                                c(Ref_tissues[which(Ref_labels=="Stage II")],train_tissue_stage_i[tr_stage_i]),TRUE)
    
    PCC_tr_stage_i[tr_stage_i,(col_num_stage_i*2+1):(col_num_stage_i*3)] <- calculate_PCC_network(cbind(Data_DEGs_topk_iii,Data_DEGs_topk_i), 
                                                                                                  c(Ref_tissues[which(Ref_labels=="Stage III")],train_tissue_stage_i[tr_stage_i]),TRUE)
    
    PCC_tr_stage_i[tr_stage_i,(col_num_stage_i*3+1):(col_num_stage_i*4)] <- calculate_PCC_network(cbind(Data_DEGs_topk_iv,Data_DEGs_topk_i),
                                                                                                  c(Ref_tissues[which(Ref_labels=="Stage IV")],train_tissue_stage_i[tr_stage_i]),TRUE)
    
  }
  
  #################### Step 2.2.2: Construct the co-expression network (DCN) networks for stage ii
  # Obatin the stage ii data in tranining set
  train_tissue_stage_ii <- Train_tissues[which(Train_labels=="Stage II")]
  row_num_stage_ii = length(train_tissue_stage_ii)
  col_num_stage_ii = length(Ref_PCC_Stage_N_II)
  PCC_tr_stage_ii <- matrix(NA,nrow = row_num_stage_ii,ncol = col_num_stage_ii*4)
  # Construct the PCC for reference network + one training sample
  cat("Start to calculate the PCC for stage ii of traning set!\n")
  
  for (tr_stage_ii in c(1:row_num_stage_ii)) {
    # cat("Start to calculate the tissue: ",tr_stage_ii, " of stage ii of traning set!\n")
    PCC_tr_stage_ii[tr_stage_ii,1:col_num_stage_ii] <- calculate_PCC_network(cbind(Data_DEGs_topk_i,Data_DEGs_topk_ii),
                                                                          c(Ref_tissues[which(Ref_labels=="Stage I")],train_tissue_stage_ii[tr_stage_ii]),TRUE)
    
    PCC_tr_stage_ii[tr_stage_ii,(col_num_stage_ii+1):(col_num_stage_ii*2)] <- calculate_PCC_network(cbind(Data_DEGs_topk_ii),
                                                                                            c(Ref_tissues[which(Ref_labels=="Stage II")],train_tissue_stage_ii[tr_stage_ii]),TRUE)
    
    PCC_tr_stage_ii[tr_stage_ii,(col_num_stage_ii*2+1):(col_num_stage_ii*3)] <- calculate_PCC_network(cbind(Data_DEGs_topk_iii,Data_DEGs_topk_ii),
                                                                                              c(Ref_tissues[which(Ref_labels=="Stage III")],train_tissue_stage_ii[tr_stage_ii]),TRUE)
    
    PCC_tr_stage_ii[tr_stage_ii,(col_num_stage_ii*3+1):(col_num_stage_ii*4)] <- calculate_PCC_network(cbind(Data_DEGs_topk_iv,Data_DEGs_topk_ii),
                                                                                              c(Ref_tissues[which(Ref_labels=="Stage IV")],train_tissue_stage_ii[tr_stage_ii]),TRUE)
    
  }

  #################### Step 2.2.3: Construct the co-expression network (DCN) networks for stage iii
  # Obatin the stage iii data in tranining set
  train_tissue_stage_iii <- Train_tissues[which(Train_labels=="Stage III")]
  row_num_stage_iii = length(train_tissue_stage_iii)
  col_num_stage_iii = length(Ref_PCC_Stage_N_III)
  PCC_tr_stage_iii <- matrix(NA,nrow = row_num_stage_iii,ncol = col_num_stage_iii*4)
  # Construct the PCC for reference network + one training sample
  cat("Start to calculate the PCC for stage iii of traning set!\n")
  
  for (tr_stage_iii in c(1:row_num_stage_iii)) {
    # cat("Start to calculate the tissue: ",tr_stage_iii, " of stage iii of traning set!\n")
    PCC_tr_stage_iii[tr_stage_iii,1:col_num_stage_iii] <- calculate_PCC_network(cbind(Data_DEGs_topk_i,Data_DEGs_topk_iii),
                                                                             c(Ref_tissues[which(Ref_labels=="Stage I")],train_tissue_stage_iii[tr_stage_iii]),TRUE)
    
    PCC_tr_stage_iii[tr_stage_iii,(col_num_stage_iii+1):(col_num_stage_iii*2)] <- calculate_PCC_network(cbind(Data_DEGs_topk_ii,Data_DEGs_topk_iii),
                                                                                                c(Ref_tissues[which(Ref_labels=="Stage II")],train_tissue_stage_iii[tr_stage_iii]),TRUE)
    
    PCC_tr_stage_iii[tr_stage_iii,(col_num_stage_iii*2+1):(col_num_stage_iii*3)] <- calculate_PCC_network(Data_DEGs_topk_iii,
                                                                                                  c(Ref_tissues[which(Ref_labels=="Stage III")],train_tissue_stage_iii[tr_stage_iii]),TRUE)
    
    PCC_tr_stage_iii[tr_stage_iii,(col_num_stage_iii*3+1):(col_num_stage_iii*4)] <- calculate_PCC_network(cbind(Data_DEGs_topk_iv,Data_DEGs_topk_iii),
                                                                                                  c(Ref_tissues[which(Ref_labels=="Stage IV")],train_tissue_stage_iii[tr_stage_iii]),TRUE)
    
  }
  
  
  #################### Step 2.2.4: Construct the co-expression network (DCN) networks for stage iv
  # Obatin the stage iv data in tranining set
  train_tissue_stage_iv <- Train_tissues[which(Train_labels=="Stage IV")]
  row_num_stage_iv = length(train_tissue_stage_iv)
  col_num_stage_iv = length(Ref_PCC_Stage_N_IV)
  PCC_tr_stage_iv <- matrix(NA,nrow = row_num_stage_iv,ncol = col_num_stage_iv*4)
  # Construct the PCC for reference network + one training sample
  cat("Start to calculate the PCC for stage iv of traning set!\n")
  
  for (tr_stage_iv in c(1:row_num_stage_iv)) {
    # cat("Start to calculate the tissue: ",tr_stage_iv, " of stage iii of traning set!\n")
    PCC_tr_stage_iv[tr_stage_iv,1:col_num_stage_iv] <- calculate_PCC_network(cbind(Data_DEGs_topk_i,Data_DEGs_topk_iv),
                                                                                c(Ref_tissues[which(Ref_labels=="Stage I")],train_tissue_stage_iv[tr_stage_iv]),TRUE)
    
    PCC_tr_stage_iv[tr_stage_iv,(col_num_stage_iv+1):(col_num_stage_iv*2)] <- calculate_PCC_network(cbind(Data_DEGs_topk_ii,Data_DEGs_topk_iv),
                                                                                                    c(Ref_tissues[which(Ref_labels=="Stage II")],train_tissue_stage_iv[tr_stage_iv]),TRUE)
    
    PCC_tr_stage_iv[tr_stage_iv,(col_num_stage_iv*2+1):(col_num_stage_iv*3)] <- calculate_PCC_network(cbind(Data_DEGs_topk_iii,Data_DEGs_topk_iv),
                                                                                                      c(Ref_tissues[which(Ref_labels=="Stage III")],train_tissue_stage_iv[tr_stage_iv]),TRUE)
    
    PCC_tr_stage_iv[tr_stage_iv,(col_num_stage_iv*3+1):(col_num_stage_iv*4)] <- calculate_PCC_network(Data_DEGs_topk_iv,
                                                                                                      c(Ref_tissues[which(Ref_labels=="Stage IV")],train_tissue_stage_iv[tr_stage_iv]),TRUE)
    
  }
  
  
  ####################
  #################### Step 3: Obtain the testing data
  #################### Step 3.1: Obtain the testing data for stage i
  #################### Step 3.1.1: Get the the PCC for stage i of testing data
  test_tissue_stage_i <- Test_tissues[which(Test_labels=="Stage I")]
  row_num_stage_i = length(test_tissue_stage_i)
  col_num_stage_i = length(Ref_PCC_Stage_N_I)
  PCC_te_stage_i <- matrix(NA,nrow = row_num_stage_i,ncol = col_num_stage_i*4)
  # Construct the PCC for reference network + one testing sample
  cat("Start to calculate the PCC for stage i of testing set!\n")
  
  for (te_stage_i in c(1:row_num_stage_i)) {
    # cat("Start to calculate the tissue: ",te_stage_i, " of stage i of traning set!\n")
    # Data_DEGs_topk_stage <- cbind(Data_DEGs_topk_n_i,Data_DEGs_topk_i)
    PCC_te_stage_i[te_stage_i,1:col_num_stage_i] <- calculate_PCC_network(Data_DEGs_topk_i,
                                                                          c(Ref_tissues[which(Ref_labels=="Stage I")],test_tissue_stage_i[te_stage_i]),TRUE)
    
    PCC_te_stage_i[te_stage_i,(col_num_stage_i+1):(col_num_stage_i*2)] <- calculate_PCC_network(cbind(Data_DEGs_topk_ii,Data_DEGs_topk_i),
                                                                                            c(Ref_tissues[which(Ref_labels=="Stage II")],test_tissue_stage_i[te_stage_i]),TRUE)
    
    PCC_te_stage_i[te_stage_i,(col_num_stage_i*2+1):(col_num_stage_i*3)] <- calculate_PCC_network(cbind(Data_DEGs_topk_iii,Data_DEGs_topk_i),
                                                                                              c(Ref_tissues[which(Ref_labels=="Stage III")],test_tissue_stage_i[te_stage_i]),TRUE)
    
    PCC_te_stage_i[te_stage_i,(col_num_stage_i*3+1):(col_num_stage_i*4)] <- calculate_PCC_network(cbind(Data_DEGs_topk_iv,Data_DEGs_topk_i),
                                                                                              c(Ref_tissues[which(Ref_labels=="Stage IV")],test_tissue_stage_i[te_stage_i]),TRUE)
    
  }
  
  
  #################### Step 3.2.2: Construct the co-expression network (DCN) networks for stage ii
  # Obatin the stage ii data in tranining set
  test_tissue_stage_ii <- Test_tissues[which(Test_labels=="Stage II")]
  row_num_stage_ii = length(test_tissue_stage_ii)
  col_num_stage_ii = length(Ref_PCC_Stage_N_II)
  PCC_te_stage_ii <- matrix(NA,nrow = row_num_stage_ii,ncol = col_num_stage_ii*4)
  # Construct the PCC for reference network + one testing sample
  cat("Start to calculate the PCC for stage ii of testing set!\n")
  
  for (te_stage_ii in c(1:row_num_stage_ii)) {
    # cat("Start to calculate the tissue: ",te_stage_ii, " of stage ii of testing set!\n")
    PCC_te_stage_ii[te_stage_ii,1:col_num_stage_ii] <- calculate_PCC_network(cbind(Data_DEGs_topk_i,Data_DEGs_topk_ii),
                                                                             c(Ref_tissues[which(Ref_labels=="Stage I")],test_tissue_stage_ii[te_stage_ii]),TRUE)
    
    PCC_te_stage_ii[te_stage_ii,(col_num_stage_ii+1):(col_num_stage_ii*2)] <- calculate_PCC_network(cbind(Data_DEGs_topk_ii),
                                                                                                c(Ref_tissues[which(Ref_labels=="Stage II")],test_tissue_stage_ii[te_stage_ii]),TRUE)
    
    PCC_te_stage_ii[te_stage_ii,(col_num_stage_ii*2+1):(col_num_stage_ii*3)] <- calculate_PCC_network(cbind(Data_DEGs_topk_iii,Data_DEGs_topk_ii),
                                                                                                  c(Ref_tissues[which(Ref_labels=="Stage III")],test_tissue_stage_ii[te_stage_ii]),TRUE)
    
    PCC_te_stage_ii[te_stage_ii,(col_num_stage_ii*3+1):(col_num_stage_ii*4)] <- calculate_PCC_network(cbind(Data_DEGs_topk_iv,Data_DEGs_topk_ii),
                                                                                                  c(Ref_tissues[which(Ref_labels=="Stage IV")],test_tissue_stage_ii[te_stage_ii]),TRUE)
    
  }
  
  #################### Step 3.2.3: Construct the co-expression network (DCN) networks for stage iii
  # Obatin the stage iii data in tranining set
  test_tissue_stage_iii <- Test_tissues[which(Test_labels=="Stage III")]
  row_num_stage_iii = length(test_tissue_stage_iii)
  col_num_stage_iii = length(Ref_PCC_Stage_N_III)
  PCC_te_stage_iii <- matrix(NA,nrow = row_num_stage_iii,ncol = col_num_stage_iii*4)
  # Construct the PCC for reference network + one testing sample
  cat("Start to calculate the PCC for stage iii of testing set!\n")
  
  for (te_stage_iii in c(1:row_num_stage_iii)) {
    # cat("Start to calculate the tissue: ",te_stage_iii, " of stage iii of testing set!\n")
    PCC_te_stage_iii[te_stage_iii,1:col_num_stage_iii] <- calculate_PCC_network(cbind(Data_DEGs_topk_i,Data_DEGs_topk_iii),
                                                                                c(Ref_tissues[which(Ref_labels=="Stage I")],test_tissue_stage_iii[te_stage_iii]),TRUE)
    
    PCC_te_stage_iii[te_stage_iii,(col_num_stage_iii+1):(col_num_stage_iii*2)] <- calculate_PCC_network(cbind(Data_DEGs_topk_ii,Data_DEGs_topk_iii),
                                                                                                    c(Ref_tissues[which(Ref_labels=="Stage II")],test_tissue_stage_iii[te_stage_iii]),TRUE)
    
    PCC_te_stage_iii[te_stage_iii,(col_num_stage_iii*2+1):(col_num_stage_iii*3)] <- calculate_PCC_network(Data_DEGs_topk_iii,
                                                                                                      c(Ref_tissues[which(Ref_labels=="Stage III")],test_tissue_stage_iii[te_stage_iii]),TRUE)
    
    PCC_te_stage_iii[te_stage_iii,(col_num_stage_iii*3+1):(col_num_stage_iii*4)] <- calculate_PCC_network(cbind(Data_DEGs_topk_iv,Data_DEGs_topk_iii),
                                                                                                      c(Ref_tissues[which(Ref_labels=="Stage IV")],test_tissue_stage_iii[te_stage_iii]),TRUE)
    
  }
  
  
  #################### Step 3.2.4: Construct the co-expression network (DCN) networks for stage iv
  # Obatin the stage iv data in tranining set
  test_tissue_stage_iv <- Test_tissues[which(Test_labels=="Stage IV")]
  row_num_stage_iv = length(test_tissue_stage_iv)
  col_num_stage_iv = length(Ref_PCC_Stage_N_IV)
  PCC_te_stage_iv <- matrix(NA,nrow = row_num_stage_iv,ncol = col_num_stage_iv*4)
  # Construct the PCC for reference network + one testing sample
  cat("Start to calculate the PCC for stage iv of testing set!\n")
  
  for (te_stage_iv in c(1:row_num_stage_iv)) {
    # cat("Start to calculate the tissue: ",te_stage_iv, " of stage iii of testing set!\n")
    PCC_te_stage_iv[te_stage_iv,1:col_num_stage_iv] <- calculate_PCC_network(cbind(Data_DEGs_topk_i,Data_DEGs_topk_iv),
                                                                             c(Ref_tissues[which(Ref_labels=="Stage I")],test_tissue_stage_iv[te_stage_iv]),TRUE)
    
    PCC_te_stage_iv[te_stage_iv,(col_num_stage_iv+1):(col_num_stage_iv*2)] <- calculate_PCC_network(cbind(Data_DEGs_topk_ii,Data_DEGs_topk_iv),
                                                                                                c(Ref_tissues[which(Ref_labels=="Stage II")],test_tissue_stage_iv[te_stage_iv]),TRUE)
    
    PCC_te_stage_iv[te_stage_iv,(col_num_stage_iv*2+1):(col_num_stage_iv*3)] <- calculate_PCC_network(cbind(Data_DEGs_topk_iii,Data_DEGs_topk_iv),
                                                                                                  c(Ref_tissues[which(Ref_labels=="Stage III")],test_tissue_stage_iv[te_stage_iv]),TRUE)
    
    PCC_te_stage_iv[te_stage_iv,(col_num_stage_iv*3+1):(col_num_stage_iv*4)] <- calculate_PCC_network(Data_DEGs_topk_iv,
                                                                                                  c(Ref_tissues[which(Ref_labels=="Stage IV")],test_tissue_stage_iv[te_stage_iv]),TRUE)
    
  }
  

  ####################
  #################### Step 4: Construct the Nearest Neighbor (NN) classifier
  #################### Step 4.1: Obtain the final data and lables for training and testing set
  ### training set 
  df_train_final <- rbind(PCC_tr_stage_i,PCC_tr_stage_ii,PCC_tr_stage_iii,PCC_tr_stage_iv)
  label_train_final <- Train_labels
  
  ### testing set 
  df_test_final <- rbind(PCC_te_stage_i,PCC_te_stage_ii,PCC_te_stage_iii,PCC_te_stage_iv)
  label_test_final <- Test_labels
  
  ### Conduct pca operation
  # row binding
  df_all_final <- rbind(df_train_final,df_test_final)
  
  # remove the columns whose values are all 0
  df_all_final_new <- as.data.frame(df_all_final[,apply(df_all_final,2,function(x) !all(x==0))])
  df_all_final_new <- as.data.frame(df_all_final_new[,apply(df_all_final_new, 2, function(x) all(is.finite(x)))])
  
  # df_final_pca = prcomp(df_all_final_new[,c(1:50000)], center = TRUE, scale. = TRUE)
  df_final_pca = prcomp(df_all_final_new, center = TRUE, scale. = TRUE)
  expl.var <- df_final_pca$sdev^2/sum(df_final_pca$sdev^2)*100
  expl.var_cum <- cumsum(expl.var)
  pc_idx <- which(expl.var_cum >= 95)[1]
  
  ### Obatin the training and testing data after conducting dimensionality reduction
  df_train_final_pca <- as.data.frame(df_final_pca$x[c(1:length(label_train_final)),c(1:pc_idx)])
  df_test_final_pca <- as.data.frame(df_final_pca$x[c( (length(label_train_final)+1):(length(label_train_final)+length(label_test_final))),
                                                    c(1:pc_idx)])
  
  #################### Step 4.2: Start to conduct sampling operation
  ## reference website: https://topepo.github.io/caret/subsampling-for-class-imbalances.html
  df_train_final_pca$Class <- as.factor(label_train_final)# class label for training set
  # df_train_final_pca$Class <- factor(df_train_final_pca$Class)
  df_test_final_pca$Class <- label_test_final # class label for testing set
  for(smp in smp_method){
    cat("The sampling method: ",smp," begins!\n")
    ### declare a variable which can save the classifing results
    results_classifier_list <- list()
    pred_results_list <- list()
    mtd_counter <- 0
    for(mtd in classifier_method){
      mtd_counter <- mtd_counter + 1
      # smp_method <- "smote"
      switch(smp,
             "downSample"={
               #################### down sampling
               df_train_final_pca_sampling <- downSample(x = df_train_final_pca[,-ncol(df_train_final_pca)],
                                                         y = df_train_final_pca$Class)
             },
             "upSample"= {
               df_train_final_pca_sampling <- upSample(x = df_train_final_pca[,-ncol(df_train_final_pca)],
                                                       y = df_train_final_pca$Class)
             },
             "smote"= {
               df_train_final_pca_sampling <- SMOTE(Class ~ ., data = df_train_final_pca)
             }, 
             stop("error"))
      
      
      
      #################### Step 4.2: Start to construct classifier using caret
      # classifier_method <- "nb"
      control <- trainControl(method = "repeatedcv",
                              number = 3,
                              repeats = 100)
      
      ### Reference information: https://topepo.github.io/caret/available-models.html
      switch(mtd,
             "nb"={
               model_train <- train(x = df_train_final_pca_sampling[,-ncol(df_train_final_pca_sampling)],
                                    y = df_train_final_pca_sampling$Class,
                                    trConstrol=control,
                                    preProcess = c("scale", "center"),
                                    method="nb", tuneLength = 15)
             },
             
             "treebag"={
               model_train <- train(x = df_train_final_pca_sampling[,-ncol(df_train_final_pca_sampling)],
                                    y = df_train_final_pca_sampling$Class,
                                    trConstrol=control,
                                    preProcess = c("scale", "center"),
                                    method="treebag", tuneLength = 15)
             },
             
             "C5.0"={
               model_train <- train(x = df_train_final_pca_sampling[,-ncol(df_train_final_pca_sampling)],
                                    y = df_train_final_pca_sampling$Class,
                                    trConstrol=control,
                                    preProcess = c("scale", "center"),
                                    method="C5.0", tuneLength = 15)
             },
             
             "rf"={
               model_train <- train(x = df_train_final_pca_sampling[,-ncol(df_train_final_pca_sampling)],
                                    y = df_train_final_pca_sampling$Class,
                                    trConstrol=control,
                                    preProcess = c("scale", "center"),
                                    method="rf", tuneLength = 15)
             },
             "rFerns"={
               model_train <- train(x = df_train_final_pca_sampling[,-ncol(df_train_final_pca_sampling)],
                                    y = df_train_final_pca_sampling$Class,
                                    trConstrol=control,
                                    preProcess = c("scale", "center"),
                                    method="rFerns", tuneLength = 15)
             },
             
             "wsrf"={
               model_train <- train(x = df_train_final_pca_sampling[,-ncol(df_train_final_pca_sampling)],
                                    y = df_train_final_pca_sampling$Class,
                                    trConstrol=control,
                                    preProcess = c("scale", "center"),
                                    method="wsrf", tuneLength = 15)
             },
             
             stop("error"))
      
      #################### Step 4.3: Obatin the testing results
      pred_results <- predict(model_train,df_test_final_pca[,-ncol(df_test_final_pca)])
      cfn_Mtx <- confusionMatrix(pred_results, as.factor(df_test_final_pca$Class))
      cat("The result of ",mtd," for ",project,"is: ",cfn_Mtx$table,"\n")
      
      acc_rate <-sum(diag(cfn_Mtx$table))/length(label_test_final)
      cat("The final accuracy of ",mtd," for ", project, " is: ",acc_rate,"\n")
      
      results_classifier_list[[mtd_counter]] <- cfn_Mtx
      pred_results_list[[mtd_counter]] <- pred_results
    }
    names(results_classifier_list) <- classifier_method
    names(pred_results_list) <- classifier_method
    save(results_classifier_list,pred_results_list, file = paste0("./Results_Co_Expres/",smp,"/",project,"_classifier_results.Rdata"))
    cat("The save work for sampling method: ",smp," for project: ",project," has been completed!\n")
  }
  
  
  ####################
  
  cat("The work for ",project," has ended!\n")
  
}
