# This is a script for running the group analysis on the HCP results
library(BayesfMRI)
result_dir <- "~/github/BayesGLM_Validation/HCP_results/5k_results/group"
# results <- list.files("/Volumes/Macintosh HD/Users/Shared/HCP/5k_results/", full.names = T)
results_all <- list.files("~/github/BayesGLM_Validation/HCP_results/5k_results/individual/", full.names = T)
library(INLA)
inla.setOption(pardiso.license = "~/licenses/pardiso.lic")
hem <- c("left","right")
gamma_0 <- c(0.5,1)

for(h in hem) {
  results <- grep(h, results_all, value = T)
  for(g in gamma_0) {
    start_time <- proc.time()[3]
    group_analysis <- BayesGLM2(results, excursion_type = '>', gamma = 0.5, no_cores = 4)
    group_analysis$total_time <- proc.time()[3] - start_time
    saveRDS(group_analysis,paste0(result_dir,'/HCP_BayesGLM2_result_',h,'_thresh',sub('\\.','',as.character(g)),'_',format(Sys.Date(),"%Y%m%d"),'.rds'))
  }
}




