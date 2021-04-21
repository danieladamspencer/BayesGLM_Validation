# This is a script for running the group analysis on the HCP results
library(BayesfMRI)
# result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group"
result_dir <- "~/Desktop"
# result_dir <- "/Users/Shared/HCP/5k_results/group"
data_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session"
# data_dir <- "~/github/BayesGLM_Validation/HCP_results/1k_results/PW/"
results_all <- list.files(data_dir, full.names = T)
# results_all <- list.files("~/github/BayesGLM_Validation/HCP_results/5k_results/individual/PW", full.names = T)
load("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/subjects.Rdata")
subjects <- subjects[1]
results_sub <- sapply(subjects, grep, x = results_all, value = T, simplify = F)
results_sub <- Reduce(c,results_sub)
# Check that all results are available for each subject
# per_subj <- vector('numeric',length = length(subjects))
# for(s in 1:length(per_subj)) {
#   per_subj[s] <- length(grep(subjects[s],results_sub))
# }
# if(any(per_subj != 4)) stop("Not all subjects have the same number of result files.")

library(INLA)
inla.setOption(pardiso.license = "~/licenses/pardiso.lic")
hem <- c("left")
gamma_0 <- c(0)

# 5k test, 4 cores, 1 subject, 2 sessions, visit 1 ----
for(h in hem) {
  results <- grep(h, results_sub, value = T)
  results <- grep("visit1", results, value = T)
  for(g in gamma_0) {
    start_time <- proc.time()[3]
    group_analysis <- BayesGLM2(results, excursion_type = '>', gamma = g, no_cores = NULL, alpha = 0.01)
    group_analysis$total_time <- proc.time()[3] - start_time
    saveRDS(group_analysis,paste0(result_dir,'/HCP_BayesGLM2_subj_',subjects,'_',h,'_thresh',sub('\\.','',as.character(g)),'_',format(Sys.Date(),"%Y%m%d"),'.rds'))
  }
}




