# this script is:
# calculate linear regression model for the maternals only:
# first filter for the classified maternal genes
# second calculate mathematical models for the degradation
# third choose the better model with likelihood ratio test


#setwd(getwd())
#rm(list = ls())

# model calculating -----------------------------------------------------------------------
# 4.
# find the best onset for each hour

# already exist in the first script
is_include_repetitions <- function(hours) {
  
  occur_hours <- data.frame(table(hours))
  if (any(occur_hours$Freq > 1))
    return(TRUE)
  return(FALSE)
  
}

# this function is getting a gene with expression and hours
# it calculates linear models for each start degradation time 
# it computes r square for each model to determine the best t1
# it returns the gene with hours manipulated according to t1
calc_t1 <- function(gene, to_hour) {
  
  best_t1 <- choose_t1(gene, to_hour)
  gene$ind[gene$ind < best_t1] <- best_t1
  return(gene)
  
}

# calculate all linear models for changing t1 (start checking from 0)
# take a gene with expression and hours
# returns the t1 for which the r square was the maximum
choose_t1 <- function(gene, to_hour) {
  
  r_square <- 0
  t1 <- 0
  for (h in seq(0,to_hour, 0.25)) { # gene$ind[nrow(gene)]
    gene_h <- gene
    gene_h$ind[gene_h$ind < h] <- h
    lm_h <- calculate_lm(gene_h)
    if (is.null(lm_h)) 
      next
    if (r_square < summary(lm_h[[1]])$r.squared) {
      r_square <- summary(lm_h[[1]])$r.squared
      t1 <- h
    }
  }
  
  return(t1)
  
}

# calculate the linear model for a gene with expression per hours
calculate_lm <- function(gene) {
  tmp_lm <- lm(values ~ ind, gene)
  
  # if the slope is greater than 0 it returns NULL
  if (all(is.na(tmp_lm)) || (is.null(tmp_lm)) || 
      (!is.na(tmp_lm$coefficients[["ind"]]) && tmp_lm$coefficients[["ind"]] > (0.1 ^ 5))) {
    return(NULL)
  }
  
  if ((!is.na(tmp_lm$coefficients[["ind"]]) && tmp_lm$coefficients[["ind"]]) > 0 && 
      (!is.na(tmp_lm$coefficients[["ind"]]) && tmp_lm$coefficients[["ind"]]) <= (0.1 ^ 5)) {
    tmp_lm$coefficients[["ind"]] <- 0
  }
  
  return(list(tmp_lm,tmp_lm$fitted.values))
}

# get a gene with expression and hours
# calculate linear models for each start degradation time 
# compute r square for each model to determine the best t1
# return the gene with hours manipulated according to t1
find_t1_t2 <- function(gene, to_hour) {
  
  best_t <- choose_t1_t2(gene, last_hour = to_hour)
  gene$ind[gene$ind < best_t[[1]]] <- best_t[[1]]
  gene$ind[gene$ind > best_t[[2]]] <- best_t[[2]]
  
  return(gene)
  
}

# calculate all linear models for changing t1 (start checking from 0)
# take a gene with expression and hours
# returns the t1 for which the r square was the maximum
choose_t1_t2 <- function(gene, last_hour) {
  
  r_square <- 0
  t1 <- 0
  t2 <- gene$ind[(length(gene$ind))]
  
  # the sequence of samples from the first sample to the last_hour sample
  gene_hours_uniq <- unique(gene$ind)[1:which(unique(gene$ind) == last_hour)]
  
  for (h1 in gene_hours_uniq) { # gene$ind[nrow(gene)]
    gene_h <- gene
    gene_h$ind[gene_h$ind < h1] <- h1
    
    # the sequence of samples from h2 sample to the last sample (12 in Zhao et al. dataset)
    genes_h2_seq <- unique(gene$ind)[(which(unique(gene$ind) == h1)+2):length(unique(gene$ind))]
    # if this is the last sample, don't take na as an sample 
    if (h1 == gene_hours_uniq[length(gene_hours_uniq)]) {
      genes_h2_seq <- unique(gene$ind)[length(unique(gene$ind))]
    }
    
    # calculate t2 for the last few samples
    for (h2 in genes_h2_seq) {
      gene_h <- gene
      gene_h$ind[gene_h$ind < h1] <- h1
      gene_h$ind[gene_h$ind > h2] <- h2
      
      lm_h <- calculate_lm(gene_h)
      if (is.null(lm_h)) 
        next
      if (r_square < summary(lm_h[[1]])$r.squared) {
        r_square <- summary(lm_h[[1]])$r.squared
        t1 <- h1
        t2 <- h2
      }
      
    }
    
  }
  return(list(t1, t2))
}


# wilks test (LRT with chi square test) ---------------------------------------------------------------------
# 5.
# model examining with Wilks test (log likelihood ratio with chi square test) for the 3 nested models

# the statistics
calculate_LR <- function(gene, model_x, model_y) {
  # model_x[[2]] are the fitted values
  if (is.null(model_x) || is.null(model_y))
    return(NA)
  
  return(-2*(as.numeric(logLik(model_x)) - as.numeric(logLik(model_y))))
  
}

# choose the right model from M0 and M1
# and return the value of calculate_LR with the best model 
send_right_model <- function(gene, first_match, model0, model1, model2) {
  
  tmp_model <- if (first_match == "M0") model0 else model1
  return(calculate_LR(gene, tmp_model, model2))
}

# input p value and names of two models
# output the name of the accepted model
choose_better_model <- function(p_val, m_x,m_y) {
  first_match_gene <- if (p_val < 0.05) m_y else m_x
  return(first_match_gene)
}

dof_one_model <- function(model_num) {
  
  return(switch (model_num,
                 M0 = {1},
                 M1 = {2},
                 M2 = {3},
                 M3 = {4}
  ))
  
}

# returns the degrees of freedom for two models [0,2]
degrees_of_freedom <- function(m_x, m_y) {
  
  dof_x <- as.numeric(dof_one_model(m_x))
  dof_y <- as.numeric(dof_one_model(m_y))
  
  return(abs(dof_y - dof_x))
}

# get the gene, the two models while the x is the more primitive one
# get also the names of the models m_x, m_y
chi_test_1_comp <- function(gene, model_x, model_y, m_x, m_y) {
  #  dof <- degrees_of_freedom(m_x,m_y)
  #  print(model_x)
  p_val <- pchisq(calculate_LR(gene, model_x, model_y),degrees_of_freedom(m_x, m_y), lower.tail = FALSE)
  return(p_val)
}

# chi square test for the second comparison - 
# M2 against M0 or M1
chi_test_2_comp <- function(gene, p_val,model_x, model_y, model_z) {
  
  first_match <- choose_better_model(p_val, "M0", "M1")
  p_val <- pchisq(send_right_model(gene,first_match, model_x, model_y, model_z),
                  degrees_of_freedom(first_match, "M2"), lower.tail = FALSE)
  return(p_val)
}

model_by_name <- function(second_comp, model0,model1, model2, model3) {
  
  tmp_model <- if (second_comp == "M0") model0 else if (second_comp == "M1") model1 else model2
  return(tmp_model)
  
}

# chi square test for the second comparison -
# M3 against M0 or M1 or M2
chi_test_3_comp <- function(gene, second_comp,model_x, model_y, model_z, model_q) {
  
  previous_comp <- model_by_name(second_comp, model_x, model_y, model_z, model_q)
  p_val <- pchisq(calculate_LR(gene, previous_comp, model_q),
                  degrees_of_freedom(second_comp, "M3"), lower.tail = FALSE)
  return(p_val)
  
}


# get parameters ----------------------------------------------------------------------------
# 6. 
# I build the table of the parameters of the models

parameter_4_each_gene <- function(genes, fitted_models, models_list, parameter_type) {
  
  param_4_genes <- sapply(1:length(genes), function(x) switch_parameter(genes[x], fitted_models[x], models_list, parameter_type))
  
}

switch_parameter <- function(gene, model, models_list, parameter_type) {
  
  return( switch(parameter_type,
                 x0 = {find_x(gene, model, models_list, "t0")},
                 #               x1 = {find_x(gene, model, models_list, "t1")},
                 B = {find_B(gene, model, models_list)},
                 t0 = {find_t0(gene, model, models_list)},
                 t1 = {find_t1(gene, model, models_list)},
                 gene = {return(gene)},
                 mod = {return(model)},
                 {return(NULL)}
  ))
}

# onset time
find_t0 <- function(gene, gene_fitted_model, models_list) {
  
  model_number <- as.numeric(stri_sub(gene_fitted_model, -1))
  t0 <- models_list[[model_number+1]][[gene]][[1]]$model$ind[1]
  return(t0)
}

# offset time
find_t1 <- function(gene, gene_fitted_model, models_list) {
  
  model_number <- as.numeric(stri_sub(gene_fitted_model, -1))
  t1 <- models_list[[model_number+1]][[gene]][[1]]$model$ind[length(models_list[[model_number+1]][[gene]][[1]]$model$ind)]
  return(t1)
}

# expression state (constant)
find_x <- function(gene, gene_fitted_model, models_list, t0_or_1) {
  
  if ((gene_fitted_model == "M0" | gene_fitted_model == "M1") & t0_or_1 == "t0") {
    model_number <- as.numeric(stri_sub(gene_fitted_model, -1))
    return(models_list[[model_number+1]][[gene]][[1]]$coefficients[1])
  }
  
  # in order to find x0 i should find t0
  t <- switch(t0_or_1,
              t0 = find_t0(gene, gene_fitted_model, models_list),
              t1 = find_t1(gene, gene_fitted_model, models_list))
  
  
  if (gene_fitted_model == "M2" | gene_fitted_model == "M3") {
    # take the right model
    model_number <- as.numeric(stri_sub(gene_fitted_model, -1))
    model_details <- models_list[[model_number+1]][[gene]][[1]]$coefficients
    m <- as.numeric(model_details[2])
    b <- as.numeric(model_details[1])
    
    return(m*t + b)
    
  } else {
    return(NA)
  }
}

# degradation rate
find_B <- function(gene, gene_fitted_model, models_list) {
  
  model_number <- as.numeric(stri_sub(gene_fitted_model, -1))
  return(models_list[[model_number+1]][[gene]][[1]]$coefficients[2])
}


# plot model parameter distribution ----------------------------------------------------------
# 7.
# plot an histogram for each parameter of the models: X0, t0, t1, B

plot_parameter_hist <- function(param_vec, param_type) {
  
  ggplot(param_vec, aes(x)) +
    geom_histogram(aes(y=(..count..)/nrow(param_vec)), binwidth = 0.1) +
    labs(y = "% of genes", x = param_type, title = paste("Histogram of", param_type)) +
    theme(text = element_text(size = 18))
  
}

plot_params_list <- function(list_parameters_per_gene) {
  
  B <- data.frame(x = as.numeric(unlist(list_parameters_per_gene$B)))
  half_life <- data.frame(x = log(2)/(-B$x))
  plot_B <- plot_parameter_hist(half_life, "half life") 
  
  t0 <- data.frame(x = as.numeric(unlist(list_parameters_per_gene$t0)), mod = unlist(list_parameters_per_gene$mod)) %>%
    select(x)
  # t0 <- t0 %>%
  #   filter(mod == "M2" | mod == "M3") %>%
  #   select(x)
  plot_t0 <- plot_parameter_hist(t0, "t0") 
  
  t1 <- data.frame(x = as.numeric(unlist(list_parameters_per_gene$t1)), mod = unlist(list_parameters_per_gene$mod)) %>%
    select(x)
  # t1 <- t1 %>%
  #   filter(mod == "M3") %>%
  #   select(x)
  plot_t1 <- plot_parameter_hist(t1, "t1") 
  
  x0 <- data.frame(x = as.numeric(unlist(list_parameters_per_gene$x0)), mod = unlist(list_parameters_per_gene$mod)) %>%
    select(x)
  plot_x0 <- plot_parameter_hist(x0, "x0") 
  
  # x1 <- data.frame(x = as.numeric(unlist(list_parameters_per_gene$x1)), mod = unlist(list_parameters_per_gene$mod))
  # x1 <- x1 %>%
  #   filter(mod == "M3") %>%
  #   select(x)
  # plot_x1 <- plot_parameter_hist(x1, "x1") 
  
  pdf(file = "parameters_graphs_m1.pdf")
  par(mfrow=c(2,2))
  
  plot(plot_B)
  plot(plot_t0)
  plot(plot_t1)
  plot(plot_x0)
  #  plot(plot_x1)
  
  dev.off()
  
}


# goodness of fit and least square error -------------------------------------------------------
# 8.
# this part checking goodness of fit by chi square test of squared error

# calculate std for datasets with more than 1 repetition
# std for specific gene and given hour
std_specific_hour <- function(gene, h) {
  
  if (is_include_repetitions(gene$ind)) {
    indices <- which(gene$ind %in% h)
  }
  else {
    indices <- which(gene$ind %in% h)
    if (indices == length(gene$ind))
      indices <- c(indices - 1, indices)
    else
      indices <- c(indices, indices + 1)
  }
  
  return(sd(gene$values[indices]))
  
  #  return(sqrt((sum((abs(gene$values[indices] - mean(gene$values[indices])))^2) / (length(indices)-1))))
  
}

# calculate standard score for each xij
standard_score <- function(gene, index, vij_val_h) {
  
  std <- std_specific_hour(gene, gene$ind[index])
  if (std < 0.01) {
    std <- 0.01
  }  
  return((gene$values[index] - vij_val_h) / std)
  
}

# this function sums the standard score
sum_stand_score <- function(gene, model_fitted) {
  
  Xij <- unlist(lapply(1:nrow(gene), function(i) standard_score(gene, i, model_fitted[i])))
  return(sum((Xij)^2))
  
}

stat_chi_sqr <- function(gene, model_fitted, mod) {
  
  stati <- sum_stand_score(gene, model_fitted)
  return(chi_sqr_goodness_of_fit(stati, nrow(gene)))
  
}

chi_sqr_goodness_of_fit <- function(statistic, dof) {
  
  chi_sqr_p_val <- pchisq(statistic, dof, lower.tail = FALSE)
  return(chi_sqr_p_val)
  
}

apply_goodness_of_fit <- function(list_gene_mod) {
  p_vals <- lapply(list_gene_mod, function(x) stat_chi_sqr(x[[1]],
                                                           x[[2]][[2]], 
                                                           x[[3]]))
  return(p_vals)
}

calc_least_square_err <- function(gene,fitted_model) {
  r_sqr <- (fitted_model - gene$values)^2
  sum_squared_err <- sum(r_sqr)
  return(sum_squared_err)
}


# model choosing and ordering in data frame --------------------------------------------------------
# 9.
# order results in dataframe

calc_pval_mod <- function(m1_data, model_M0, model_M1, model_M2, model_M3) {
  
  # first comparison - model M0 vs. M1
  # for each gene, calculate the chi square distribution
  m0_vs_m1 <- unlist(lapply(names(m1_data), function(x) { chi_test_1_comp(m1_data[[x]], model_M0[[x]][[1]], model_M1[[x]][[1]], "M0", "M1")}))
  m0_vs_m1_corrected <- p.adjust(m0_vs_m1, "fdr")
  
  # compare M2 with the best model for each gene
  m_num_vs_m2 <- unlist(lapply(1:length(names(m1_data)), function(x) chi_test_2_comp(m1_data[[x]], m0_vs_m1_corrected[[x]], model_M0[[x]][[1]], model_M1[[x]][[1]], model_M2[[x]][[1]])))
  m_num_vs_m2_corrected <- p.adjust(m_num_vs_m2,"fdr")
  
  most_fitted_model <- data.frame(gene = names(m1_data), comp1 = m0_vs_m1_corrected, comp2 = m_num_vs_m2_corrected)
  
  most_fitted_model$first_comp_mod <- rep("M0", nrow(most_fitted_model))
  most_fitted_model$first_comp_mod[most_fitted_model$comp1 < 0.05] <- "M1"
  
  most_fitted_model$second_comp_mod <- rep("M0", nrow(most_fitted_model))
  most_fitted_model$second_comp_mod[most_fitted_model$comp2 < 0.05] <- "M1"
  
  # most_fitted_model <- most_fitted_model %>%
  #   mutate(first_comp_mod = derivedFactor("M1" = (comp1 < 0.05),"M0" = (comp1 >= 0.05),.method = "first",.default = NA)) %>%
  #   mutate(second_comp_mod = derivedFactor("M1" = (comp1 < 0.05),"M0" = (comp1 >= 0.05),.method = "first",.default = NA))
  
  most_fitted_model <- most_fitted_model %>%
    mutate(second_comp_mod = case_when((comp2 < 0.05) ~ "M2",
                                       (comp2 >= 0.05) ~ as.character(first_comp_mod)))
  
  # compare M3 with the best model for each gene 
  m_num_vs_m3 <- unlist(lapply(1:length(names(m1_data)), 
                               function(x) chi_test_3_comp(m1_data[[x]], 
                                                           most_fitted_model$second_comp_mod[x], 
                                                           model_M0[[x]][[1]], 
                                                           model_M1[[x]][[1]], 
                                                           model_M2[[x]][[1]], 
                                                           model_M3[[x]][[1]])))
  
  m_num_vs_m3_corrected <- p.adjust(m_num_vs_m3, "fdr")
  
  most_fitted_model$comp3 <- m_num_vs_m3_corrected 
  most_fitted_model <- most_fitted_model %>% 
    mutate(third_comp_mod = second_comp_mod) %>%
    mutate(third_comp_mod = case_when((comp3 < 0.05) ~ "M3",
                                      (comp3 >= 0.05) ~ as.character(second_comp_mod)))
  
  
  # cosmetical fix:
  # if the two models are the same and the more complicated model is chosen, 
  # than fix it.
  most_fitted_model <- cosmetical_model_fix(most_fitted_model, model_M0, model_M1, model_M2, model_M3)
  
  return(most_fitted_model)
  
}

# cosmetical fix:
# if the two models are the same and the more complicated model is chosen, 
# than fix it.
cosmetical_model_fix <- function(most_fitted_model, model_M0, model_M1, model_M2, model_M3) {
  
  for (g in most_fitted_model$gene) {
    
    if ((most_fitted_model[which(most_fitted_model$gene %in% g), c("first_comp_mod")] == "M1") & 
        all(model_M0[[g]][[1]]$fitted.values == model_M1[[g]][[1]]$fitted.values)) {
      
      most_fitted_model[which(most_fitted_model$gene %in% g), "first_comp_mod"] <- "M0"
    }
    else if ((most_fitted_model[which(most_fitted_model$gene %in% g), "second_comp_mod"] == "M2") & 
             all(model_M1[[g]][[1]]$fitted.values == model_M2[[g]][[1]]$fitted.values)) {
      
      most_fitted_model[which(most_fitted_model$gene %in% g), "second_comp_mod"] <- "M1"
    } 
    else if ((most_fitted_model[which(most_fitted_model$gene %in% g), "third_comp_mod"] == "M3") & 
             all(model_M2[[g]][[1]]$fitted.values == model_M3[[g]][[1]]$fitted.values)) {
      
      most_fitted_model[which(most_fitted_model$gene %in% g), "third_comp_mod"] <- "M2"
    }
    
  }
  return(most_fitted_model)
  
}

model_and_save <- function(m0_data, m1_data, m2_data, m3_data) {
  
  model_M0 <- lapply(m0_data, calculate_lm)
  model_M1 <- lapply(m1_data, calculate_lm)
  model_M2 <- lapply(m2_data, calculate_lm)
  model_M3 <- lapply(m3_data, calculate_lm)
  
  # remove from the list genes that may have a slope > 0 from models lists and from data lists
  # take the names of the genes that has a model null because of slope > 0
  null_model_1 <- names(m1_data[unlist(lapply(1:length(model_M1), function(i) if(is.null(model_M1[[i]])) return(i)))])
  null_model_2 <- names(m2_data[unlist(lapply(1:length(model_M2), function(i) if(is.null(model_M2[[i]])) return(i)))])
  null_model_3 <- names(m3_data[unlist(lapply(1:length(model_M3), function(i) if(is.null(model_M3[[i]])) return(i)))])
  
  if ((length(null_model_1) > 0) || (length(null_model_2) > 0) || (length(null_model_3) > 0)) {
    
    # intersect the names and remove from the models list and the data
    gene_names_to_drop <- union(union(null_model_1, null_model_2), null_model_3)
    if (is.null(gene_names_to_drop) == FALSE) {
      model_M0 <- model_M0[names(model_M0) %in% gene_names_to_drop == FALSE]
      model_M1 <- model_M1[names(model_M1) %in% gene_names_to_drop == FALSE]
      model_M2 <- model_M2[names(model_M2) %in% gene_names_to_drop == FALSE]
      model_M3 <- model_M3[names(model_M3) %in% gene_names_to_drop == FALSE]
      
      m0_data <- m0_data[names(m0_data) %in% gene_names_to_drop == FALSE]
      m1_data <- m1_data[names(m1_data) %in% gene_names_to_drop == FALSE]
      m2_data <- m2_data[names(m2_data) %in% gene_names_to_drop == FALSE]
      m3_data <- m3_data[names(m3_data) %in% gene_names_to_drop == FALSE]
      
    }
  }
  
  # save data for models
  save(m0_data, file = "data_M0_wo_slope_gt_0.rdata")
  save(m1_data, file = "data_M1_wo_slope_gt_0.rdata")
  save(m2_data, file = "data_M2_wo_slope_gt_0.rdata")
  save(m3_data, file = "data_M3_wo_slope_gt_0.rdata")
  
  save(model_M0, file = "model_M0.rdata")
  save(model_M1, file = "model_M1.rdata")
  save(model_M2, file = "model_M2.rdata")
  save(model_M3, file = "model_M3.rdata")
  
}

get_the_mod <- function(most_fitted_model, m1_data, model0, model1, model2, model3) {
  mod_fitted <- m1_data
  for (one_gene in names(m1_data)) {
    
    mod_tmp <- 
      switch (most_fitted_model$third_comp_mod[most_fitted_model$gene == one_gene],
              "M0" = {model0[[one_gene]]} , 
              "M1" = {model1[[one_gene]]} , 
              "M2" = {model2[[one_gene]]} , 
              "M3" = {model3[[one_gene]]})
    
    mod_fitted[[one_gene]] <- list(m1_data[[one_gene]], mod_tmp, 
                                   most_fitted_model$third_comp_mod[most_fitted_model$gene == one_gene])
    
  }
  return(mod_fitted)
}


# apply and call functions --------------------------------------------------------------------
# 10.

# this part take the classifications and the exp levels and 
# filter for mature transcript of maternal genes 
#' Title
#'
#' @return m0_matrix - 
#' @export
#'
#' @examples
load_data <- function() {
  
  class_criteria <- readRDS(file = "class_criteria.RDS")
  
  normalized_exp <- read.csv(file = "normalized_exp.csv", header = TRUE, sep = ",")[,-1]
  
  maternal_only <- normalized_exp %>%
    inner_join(class_criteria %>% 
                 select(gene_id, classification)) %>%
    filter(classification == "M") %>%
    select(gene_id, gene_name, starts_with("transcript_"))
  
  # take only exp levels for linear fit
  m0_matrix <- maternal_only %>% 
    select(starts_with("transcript_"))
  
  # move the name and id to row names
  gene_id_name <- unlist(sapply(1:nrow(maternal_only), function(x) paste0(maternal_only$gene_id[x], "_", maternal_only$gene_name[x])))
  row.names(m0_matrix) <- gene_id_name
  
  return(m0_matrix)
}

# this function is turning the hours, for stacked gene, into numeric
hours_2_numeric <- function(gene){
  if (is.factor(gene$ind)) {
    gene$ind <- as.numeric(levels(gene$ind))[gene$ind]
  }
  else {
    if (is.character(gene$ind)) {
      gene$ind <- as.numeric(gene$ind)
    }
  }
  return(gene)
}

#' Title
#'
#' @param m0_matrix 
#'
#' @return m_data - general data, basic for all the models.
#' @export
#'
#' @examples
prepare_m_data <- function(m0_matrix) {
  
  # change colnames to hours
  colnames(m0_matrix) <- gsub("transcript_|_\\d+", "", colnames(m0_matrix))
  
  m_not_numeric <- apply(m0_matrix, 1, stack)
  
  m_data <- lapply(m_not_numeric, hours_2_numeric)
  
  return(m_data)
}

# constant model
#' Title
#'
#' @param m0_matrix 
#'
#' @return m0_data - the hours all moved to zero so lm function 
#'  would make a constant model from it 
#' @export
#'
#' @examples
prepare_m0_data <- function(m0_matrix) {
  
  # change hours to 0 and create list of stacks for zeros
  m0_zeros <- m0_matrix
  colnames(m0_zeros) <- rep(0, times = ncol(m0_zeros))
  m0_data_zeros <- apply(m0_zeros,1, stack)
  m0_data <- lapply(m0_data_zeros, hours_2_numeric)
  
  return(m0_data)
}

# prepering data for M1 - simple linear model
#' Title
#'
#' @param m_data 
#'
#' @return m1_data - the data that lm can be fitted to simple linear model
#' @export
#'
#' @examples
prepare_m1_data <- function(m_data) {
  
  # expression per hour w/o changes for early degradation onset
  m1_data <- m_data
  
  return(m1_data)
}

# prepering data for M2 - complicated linear model
# first find the best onset for each hour
# second manipualate hours according to t1
#' Title
#'
#' @param m1_data 
#'
#' @return m2_data - the data for the complicated model
#'  with t_0 as onset time
#' @export
#'
#' @examples
prepare_m2_data <- function(m1_data) {
  
  if (length(unique(m1_data[[1]]$ind)) < 3) {
    latest_t1 <- unique(m1_data[[1]]$ind)[1]
  } else {
    latest_t1 <- unique(m1_data[[1]]$ind)[(length(unique(m1_data[[1]]$ind)))-2]
  }  
  # for each gene- calculate the t1
  m2_data <- lapply(m1_data, function(x) calc_t1(x, latest_t1))
  
  return(m2_data)
}

# prepering data for M3 - more complicated linear model
# first find the best onset and offset for each hour
# second manipualate hours according to t0 and t1 
#' Title
#'
#' @param m1_data 
#'
#' @return m3_data - the data for the more complicated model
#'  with t_0 as onset time and t_1 as offset
#' @export
#'
#' @examples
prepare_m3_data <- function(m1_data) {
  
  # limit onset time to 3 hours before the last timepoint
  latest_sample_hour <- unique(m1_data[[1]]$ind)[(length(unique(m1_data[[1]]$ind)))-2]
  
  # for each gene- calculate the t1
  m3_data <- lapply(m1_data, function(x) find_t1_t2(x, latest_sample_hour))
  
}

run_models_stat_tests <- function() {
  
  # call all m_data 
  m0_matrix <- load_data()
  m_data <- prepare_m_data(m0_matrix)
  m0_data <- prepare_m0_data(m0_matrix)
  m1_data <- prepare_m1_data(m_data)
  m2_data <- prepare_m2_data(m1_data)
  m3_data <- prepare_m3_data(m1_data)
  
  # model all m_data
  # and save
  model_and_save(m0_data, m1_data, m2_data, m3_data)
  
  load("model_M0.rdata")
  load("model_M1.rdata")
  load("model_M2.rdata")
  load("model_M3.rdata")
  
  load("data_M0_wo_slope_gt_0.rdata")
  load("data_M1_wo_slope_gt_0.rdata")
  load("data_M2_wo_slope_gt_0.rdata")
  load("data_M3_wo_slope_gt_0.rdata")
  
  most_fitted_model <- calc_pval_mod(m1_data, model_M0, model_M1, model_M2, model_M3)
  
  parameter_list <- c("gene","mod", "x0", "B", "t0", "t1")
  models_list <- list(model_M0, model_M1, model_M2, model_M3)
  
  list_genes_mod <- get_the_mod(most_fitted_model, m1_data, model_M0, model_M1, model_M2, model_M3)
  pvalues <- data.frame(p.v = unlist(apply_goodness_of_fit(list_genes_mod)))
  corrected_pv <- p.adjust(pvalues$p.v, "fdr")
  
  p_signif <- list(length(which(corrected_pv < 0.05)))
  saveRDS(file = "p_val_signif.RDS", corrected_pv)
  
  lse_all_genes <- unlist(lapply(list_genes_mod, function(x) calc_least_square_err(x[[1]],x[[2]][[2]])))
  mean_lse <- mean(lse_all_genes)
  
  saveRDS(file = "lse_list.RDS", lse_all_genes)
  
  # mse_all_genes <- unlist(lapply(list_genes_mod, function(x) MSE(x[[2]][[2]], x[[1]]$values)))
  # saveRDS(file = "mse_list.RDS", mse_all_genes)
  
  list_parameters_per_gene <- 
    as.data.frame(sapply(parameter_list, 
                         function(p) parameter_4_each_gene(most_fitted_model$gene,
                                                           most_fitted_model$third_comp_mod, 
                                                           models_list, p)))
  write.csv(list_parameters_per_gene, "model_parameters_per_gene.csv")
  
  plot_params_list(list_parameters_per_gene)
  saveRDS(list_parameters_per_gene, "list_parameters_lm.RDS")
  
}





