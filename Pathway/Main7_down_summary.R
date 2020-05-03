
# This code is used to conduct the result analysis (i.e. summary) for down-regulated group.

setwd("Your current path")
# library(readxl)
library(openxlsx)

projects <- c("BLCA_Stage","BRCA_Stage","COAD_Stage","HNSC_Stage","KIRC_Stage","KIRP_Stage","LUAD_Stage","STAD_Stage","THCA_Stage")
stages <- c("i","ii","iii","iv")

PCC_cutoff=c(0.60, 0.70, 0.80)
for(pcc_cutoff_item in PCC_cutoff){
  path_down_regulated <- paste0("./085/",pcc_cutoff_item,"/","generalized_Down/")
  ########################## step 1: load the data
  
  file_name_down <- list.files(path_down_regulated,pattern = "xlsx$")
  
  down_regulated_list <- list()
  projects_count <- 0
  for (ii in projects[2:9]) {
    projects_count <- projects_count + 1
    stage_list_type <- list()
    
    file_count <- 0
    for (jj in file_name_down) {
      file_count <- file_count + 1
      stage_list_type[[file_count]] <- read.xlsx(paste(path_down_regulated,jj,sep = ""),sheet =ii)
      cat("compute for the file: ",file_count,"\n")
      
    }
    df_list_type <- NULL
    for (kk in c(1:4)) {
      temp_term <- stage_list_type[[kk]]
      dim_term <- dim(temp_term)
      col_name <- rep(stages[kk],dim_term[1])
      temp_term_new <- data.frame(temp_term,Stages=col_name)
      
      df_list_type <- rbind(df_list_type,temp_term_new)
    }
    down_regulated_list[[projects_count]] <- df_list_type
    cat("The project ",ii," is completed!","\n")
    
  }
  names(down_regulated_list) <- projects[2:9]
  # write.xlsx(down_regulated_list,file = "./generalized_Down_summary.xlsx")
  write.xlsx(down_regulated_list,file = paste0("./085/",pcc_cutoff_item,"/","generalized_Down_summary_",pcc_cutoff_item,".xlsx"))
}








