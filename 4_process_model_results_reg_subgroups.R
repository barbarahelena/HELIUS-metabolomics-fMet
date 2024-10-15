## Processing results subgroups

setwd('XGB')
source("functions.R")
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ggsci)
getwd()
# groups <- c("Afr", "NonAfr", "Fem", "Male", "Old", "Young", "DM", "All")
# 

# for(g in groups) {
#       li <- list.files(path = g)
# 
#       for (i in li) {
#         res <- ifelse(i=="DBPres"|i=="SBPres", TRUE, FALSE)
#         li2 <- list.files(path = str_c(g, i, sep = '/'))
#         a <- str_detect(li2, i)
#         b <- !(str_detect(li2, 'PERMUTED'))
#         c <- str_detect(li2, "input")
# 
#         path_true <- file.path(g, i, li2[which(a&b)])
#         path_permuted <- file.path(g, i, li2[which(a&!b)])
#         data_path <- file.path(g, i, li2[which(c)])
# 
#         compared_to_permuted_reg(path_true_outcome = path_true, path_permuted_outcome = path_permuted)
#         plot_features_tests_reg(data_path, path_true, top_n = 20, outcome_name = i, x_lab = 'Concentration', residuals = res)
#         plot_features_top_n_reg(data_path, path_true, outcome_name='SBP', x_lab='Concentration', residuals = res)
#         plot_feature_importance_microbiome(path_true, top_n = 20)
#         plot_feature_importance_class_microbiome(path_true, top_n = 20)
#       }
# }

# path_true_sbp <- 'All/SBP/output_XGB_reg_SBP_100_iterations_y_scaled_2020_11_06__16-41-19'
# pl1 <- plot_feature_importance_class_microbiome(path_true_sbp, top_n = 20) +
#   ggtitle("Best predictors SBP")
# pl1
# 
# path_true_dbp <- 'All/DBP/output_XGB_reg_DBP_100_iterations_y_scaled_2020_11_06__16-45-06'
# pl2 <- plot_feature_importance_class_microbiome(path_true_dbp, top_n = 20) +
#   ggtitle("Best predictors DBP")
# pl2
# 
# ggarrange(pl1, pl2, nrow = 2, ncol = 1, labels = c("A", "B"))
# ggsave('Relative_importance_All.pdf', width=9, height=10)


path_true <- file.path('BRS/output_XGB_reg_BRS_2020_12_16__12-27-15')
#path_permuted <- file.path(g, i, li2[which(a&!b)])
data_path <- file.path('BRS/input_data')

#compared_to_permuted_reg(path_true_outcome = path_true, path_permuted_outcome = path_permuted)
plot_features_tests_reg(data_path, path_true, top_n = 20, outcome_name = "BRS", x_lab = 'Concentration', residuals = FALSE)
plot_features_top_n_reg(data_path, path_true, outcome_name='BRS', x_lab='Concentration', residuals = FALSE)
plot_feature_importance_microbiome(path_true, top_n = 20)
plot_feature_importance_class_microbiome(path_true, top_n = 20)
