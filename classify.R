# this script is:
# do manual classifications to Maternal/Zygotic/Maternal-Zygotic: 
# first take the three criteria for classification
# second plot the distribution of these three criteria
# third classify genes

# the data needed:
# log_exp_tbl.csv file

# setwd(getwd())
# rm(list = ls())

# ---------------------------------------------------------------------------------
# 1.
# this part takes the MZT hour of this organism
# read the filtered normalized data from csv file
# cut the data to mature and precursor 

read_relev_data <- function(file_name) {
  
  normalized_exp <- read.csv(file = file_name, header = TRUE, sep = ",")[,-1]
  return(normalized_exp)
  
}

get_mzt <- function() {
  
  mzt <- readline(prompt = "Enter hour of the MZT: ")
  return(mzt)
}


return_mature_only <- function(normal_exp) {
  # take the mature transcripts
  mature <- cbind(normal_exp[,1:2], normal_exp[,grep("transcript", colnames(normal_exp))])
  return(mature)
}

return_precursor_only <- function(normal_exp) {
  # take the pre-mRNA
  precursor <- cbind(normal_exp[,1:2], normal_exp[,grep("precursor", colnames(normal_exp))])
  return(precursor)
}


get_info_saperate_mat_pre <- function() {
  
  file_name <- "normalized_exp.csv"
  normalized_exp <- read_relev_data("normalized_exp.csv")
  
  # for zebrafish:
  # MZT <- 3 
  MZT <- as.numeric(get_mzt())
  
  mature <- return_mature_only(normalized_exp)
  precursor <- return_precursor_only(normalized_exp)
  
  return(list(mature, precursor, MZT))
}


# ----------------------------------------------------------------------------------
# 2.
# this part is taking the three criteria 

# ----------------
# first criterion:
# take maximum hour for mature transcript before the mzt:
mature_bf_MZT <- function(mature_trans, hours_bf_mzt) {
  
  before_mzt_col_name <- mature_trans[,1:(length(hours_bf_mzt)+2)] 
  before_mzt_col_name$max_ <- unlist(apply(before_mzt_col_name[,-c(1,2)], 1, max))
  
  return(before_mzt_col_name)
}

# -----------------
# second criterion:
# take the ratio of maximum expression of pre-mRNA after MZT to before MZT
pre_ratio_aft_bf <- function(precursors, MZT) {
  
  pre_hours <- as.numeric(gsub("precursor_|_\\d+","",colnames(precursors)[-c(1,2)]))
  pre_hour_names_aft_mzt <- pre_hours[pre_hours > MZT]
  pre_hour_names_bf_mzt <- pre_hours[pre_hours < MZT]
  
  pre_hour_bf_mzt <- precursors[,1:(length(pre_hour_names_bf_mzt)+2)]
  pre_hour_aft_mzt <- precursors[,c(1,2,(length(pre_hour_names_bf_mzt)+3):(ncol(precursors)))]
  pre_hour_bf_mzt$max_ <- unlist(apply(pre_hour_bf_mzt[,-c(1,2)], 1, max))
  pre_hour_aft_mzt$max_ <- unlist(apply(pre_hour_aft_mzt[,-c(1,2)], 1, max))
  
  # becuase the maximum expression of precursor before the mzt is the deliminator, it shouldn't be zero 
  if (any(pre_hour_bf_mzt$max_ == 0))
    pre_hour_bf_mzt$max_[which(pre_hour_bf_mzt$max_ ==0)] <- 0.1^5
  
  ratio_pre <- cbind(pre_hour_bf_mzt[,1:2], ratio_precursor = pre_hour_aft_mzt$max_ - pre_hour_bf_mzt$max_)
  
  return(ratio_pre)
}

# ----------------
# third criterion:
# ratio of mature log max hour after mzt to mature log max hour before mzt
mat_ratio_aft_bf <- function(matures, hours_bf_mzt, before_mzt_col_name) {
  
  mature_hour_aft_mzt <- matures[,c(1,2,(length(hours_bf_mzt)+3):(ncol(matures)))]
  mature_hour_aft_mzt$max_ <- unlist(apply(mature_hour_aft_mzt[,-c(1,2)], 1, max))
  
  # becuase the maximum expression of mature before the mzt is the deliminator, it shouldn't be zero 
  if (any(before_mzt_col_name$max_ == 0))
    before_mzt_col_name$max_[which(before_mzt_col_name$max_ ==0)] <- 0.1^5
  
  ratio_mature <- cbind(mature_hour_aft_mzt[,1:2], ratio_mature = mature_hour_aft_mzt$max_ - before_mzt_col_name$max_)
  
}

set_criteria_thresholds <- function() {
  
  mat_pre <- get_info_saperate_mat_pre()
  mature <- mat_pre[[1]]
  precursor <- mat_pre[[2]]
  MZT <- mat_pre[[3]]
  
  gene_hours <- as.numeric(gsub("transcript_|_\\d+","",colnames(mature)[-c(1,2)]))
  hours_bf_mzt <- gene_hours[gene_hours < MZT]
  
  first_mat_bf <- mature_bf_MZT(mature, hours_bf_mzt)
  second_pre_bf_aft <- pre_ratio_aft_bf(precursor, MZT)
  third_mat_bf_aft <- mat_ratio_aft_bf(mature, hours_bf_mzt, first_mat_bf)
  
  # combine all criteria to one data frame
  three_df <- cbind(third_mat_bf_aft[,1:2],
                    mature_bf_mzt = first_mat_bf$max_ , 
                    ratio_precursor = second_pre_bf_aft$ratio_precursor, 
                    ratio_mature = third_mat_bf_aft$ratio_mature)
  
  return(three_df)
  
}

# ------------------------------------------------------------------------------------------
# 3.
# combine all criteria to one data frame
# determine threshold for the 3 criteria


determine_3_thresholds <- function(df_criteria) {
  
  threshold_mat_bf_mzt <- (mean(df_criteria$mature_bf_mzt) - 1*sd(df_criteria$mature_bf_mzt))
  threshold_ratio_precursor <- log2(1.25)
  threshold_ratio_mature <- log2(0.9)
  
  thresholds_3 <- c(threshold_mat_bf_mzt = threshold_mat_bf_mzt, 
                    threshold_ratio_precursor = threshold_ratio_precursor,
                    threshold_ratio_mature = threshold_ratio_mature)
  
  return(thresholds_3)
}



# -------------------------------------------------------------------------------------------
# 4.
# create histograms of each one from the three criteria in PDF file

pdf_thresholds <- function(three_df, thresholds_3) {
  
  pdf("three_criteria_graph.pdf")
  
  # maximum hour for mature transcript before the mzt from section 1.
  plot(ggplot(data = three_df, aes(x = mature_bf_mzt)) + 
         geom_histogram(binwidth = 0.1, color="plum4", fill="plum4")  + 
         coord_cartesian(ylim = c(0,700)) + 
         geom_vline(aes(xintercept = thresholds_3["threshold_mat_bf_mzt"]), col = "orange", show.legend = F) + 
         labs(title = "Maximum Exp for Mature Transcripts Before MZT", x = "log(FPKM)") + 
         theme(text = element_text(size = 16), axis.text = element_text(size = 18), legend.text = element_text(size = 18)))
  
  # ratio precursor after the mzt to before max exp from section 2.
  plot(ggplot(data = three_df, aes(x = ratio_precursor)) + 
         geom_histogram(binwidth = 0.03, color="plum4", fill="plum4")  + 
         coord_cartesian(ylim = c(0,700)) + 
         geom_vline(aes(xintercept = thresholds_3["threshold_ratio_precursor"]), col = "orange", show.legend = F) + 
         labs(title = "Ratio of Precursor After MZT to Before MZT", x = "log(FPKM) ratio") + 
         theme(text = element_text(size = 16), axis.text = element_text(size = 18), legend.text = element_text(size = 18)))
  
  # ratio of mature max exp after mzt to before mzt from section 3.
  plot(ggplot(data = three_df, aes(x = ratio_mature)) + 
         geom_histogram(binwidth = 0.1, color="plum4", fill="plum4")  + 
         coord_cartesian(xlim = c(- 10, 10), ylim = c(0,700)) + 
         geom_vline(aes(xintercept = thresholds_3["threshold_ratio_mature"]), col = "orange", show.legend = F) + 
         labs(title = "Ratio of Mature Transcript After MZT to Before MZT", x = "log(FPKM) ratio") + 
         theme(text = element_text(size = 16), axis.text = element_text(size = 18), legend.text = element_text(size = 18)))
  
  dev.off()
  
}


# ----------------------------------------------------------------------------------------
# 5.
# caracterize the genes by the three criteria and three thresholds
# classify by its characterization

above_below_threshold <- function(three_df, thresholds_3) {
  
  three_df$class_mat_bf <- "-"
  three_df$class_mat_bf[which(three_df$mature_bf_mzt >= thresholds_3["threshold_mat_bf_mzt"])] <- "+"
  
  three_df$class_ratio_pre <- "-"
  three_df$class_ratio_pre[which(three_df$ratio_precursor >= thresholds_3["threshold_ratio_precursor"])] <- "+"
  
  three_df$class_ratio_mat <- "-"
  three_df$class_ratio_mat[which(three_df$ratio_mature >= thresholds_3["threshold_ratio_mature"])] <- "+"
  
  return(three_df)
}

decision_tree <- function(three_df) {
  
  three_df$classification <- "NONE"
  # for maternal class is should:
  # have high expression before the MZT (+)
  # have low or no change of precursor expression (-)
  # have negative change of mature expression, because it suppose to degrade after the MZT (-)
  three_df$classification[(three_df$class_mat_bf == "+") & 
                            (three_df$class_ratio_pre == "-") & 
                            (three_df$class_ratio_mat == "-")] <- "M"
  
  
  # for zygotic class is should:
  # have low to no expression before the MZT (-)
  # have high change of precursor expression (+)
  # have high positive change of mature expression, because it suppose to express only after the MZT (+)
  three_df$classification[(three_df$class_mat_bf == "-") &
                            (three_df$class_ratio_pre == "+") &
                            (three_df$class_ratio_mat == "+")] <- "Z"
  
  
  # for maternal-zygotic class is should:
  # have high expression before the MZT (+) 
  # have high change of precursor expression (+)
  # have high positive (or not) change of mature expression, because it suppose to express only after the MZT (+/-)
  three_df$classification[(three_df$class_mat_bf == "+") &
                            (three_df$class_ratio_pre == "+") &
                            (three_df$class_ratio_mat == "+")] <- "MZ"
  
  three_df$classification[(three_df$class_mat_bf == "+") &
                            (three_df$class_ratio_pre == "+") &
                            (three_df$class_ratio_mat == "-")] <- "MZ"
  
  saveRDS(three_df, file = "class_criteria.RDS") 
  
  return(three_df)
}


determine_classifications <- function() {
  
  df_3_criteria <- set_criteria_thresholds()
  thresholds <- determine_3_thresholds(df_3_criteria)
  pdf_thresholds(df_3_criteria, thresholds)
  df_criteria <- above_below_threshold(df_3_criteria, thresholds)
  df_classes <- decision_tree(df_criteria)
  
  return(df_classes)
}

# -----------------------------------------------------------------------------------------
# 6.
# plot a graph of classification division

classification_count_graph <- function(three_df) {
  
  pdf("class_count_graph.pdf")
  
  plot(ggplot(three_df, aes(classification)) + 
         geom_bar() + 
         labs(title = "manual classification - counts of each category") + 
         theme(text = element_text(size = 18)))
  
  dev.off()
  
}

save_class_fun <- function(normalized_exp, three_df) {
  
  # save the expression data frame with the classification
  exp_classified <- normalized_exp %>%
    inner_join(three_df %>%
                 select(gene_id, gene_name, classification), by = c("gene_id", "gene_name"))
  
  saveRDS(exp_classified, file = "exp_df_classified.RDS")
  
  write_csv(exp_classified, "exp_df_classified")
  
}

write_csv <- function(tble, file_name) {
  
  # save the expression table in csv file
  write.csv(tble, file = paste0(file_name,".csv"))
  
}

create_exp_class_tbl <- function() {
  
  df_gene_classes <- determine_classifications()
  
  normalized_exp <- read_relev_data("normalized_exp.csv")
  save_class_fun(normalized_exp, df_gene_classes)
  classification_count_graph(df_gene_classes)
  
}

