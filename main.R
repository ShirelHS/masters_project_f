rm(list = ls())
setwd("rabani lab/master_project/master_project")

source("first.R")
analyze_data()

source("second.R")
create_exp_class_tbl()

source("third.R")
run_models_stat_tests()