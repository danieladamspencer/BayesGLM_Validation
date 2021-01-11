# This is a script for running the group analysis on the HCP results
library(BayesfMRI)
results <- list.files("/Volumes/Macintosh HD/Users/Shared/HCP/5k_results/", full.names = T)
# results <- grep("left", results, value = T)
results <- grep("right", results, value = T)

library(INLA)
inla.setOption(pardiso.license = "~/licenses/pardiso.lic")

group_analysis <- BayesGLM2(results, excursion_type = '>')
# saveRDS(group_analysis,'~/Desktop/HCP_BayesGLM2_result_left_20210105.rds')
# This took 6:54:41 (the left side)
saveRDS(group_analysis,'~/Desktop/HCP_BayesGLM2_result_right_20210106.rds')

# Plot the results ----
group_analysis_left <- readRDS("~/Desktop/HCP_BayesGLM2_result_left_20210105.rds")
library(ciftiTools)
ciftiTools.setOption('wb_path',"/Users/danspen/workbench")

# >> Plot Estimates ----
indiv_results_left <- readRDS("/Volumes/Macintosh HD/Users/Shared/HCP/5k_results/103818_visit1_left_5k_20201120.rds")
cifti_left <- indiv_results_left$betas_Bayesian$avg
cifti_left$data$cortex_left <- group_analysis_left$estimates
plot(cifti_left, idx = 4, zlim = c(-1,1), color_mode = "diverging", title = "Tongue")

# >> Plot Activations ----
indiv_results_left <- readRDS("/Volumes/Macintosh HD/Users/Shared/HCP/5k_results/103818_visit1_left_5k_20201120.rds")
active_left <- indiv_results_left$betas_Bayesian$avg
active_left$data$cortex_left <- group_analysis_left$active
plot(active_left, idx = 4)
