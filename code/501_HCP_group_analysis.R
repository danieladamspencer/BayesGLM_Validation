# This is a script for running the group analysis on the HCP results
library(BayesfMRI)
result_dir <- "~/Desktop"
data_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session"
results_all <- list.files(data_dir, full.names = T)
load("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/subjects.Rdata")
results_sub <- sapply(subjects, grep, x = results_all, value = T, simplify = F)
results_sub <- Reduce(c,results_sub)
library(INLA)
inla.setOption(pardiso.license = "~/licenses/pardiso.lic")
hem <- c("left","right")
gamma_0 <- c(0,0.5,1)

# 5k test, 4 cores, 1 subject, 2 sessions, visit 1 ----
for(subject in subjects) {
  for(visit in paste0("visit",1:2)) {
    for(h in hem) {
      results <- grep(subject, results_sub, value = T)
      results <- grep(visit,results, value=T)
      results <- grep(h, results, value = T)
      for(g in gamma_0) {
        start_time <- proc.time()[3]
        group_analysis <- BayesGLM2(results, excursion_type = '>', gamma = g, no_cores = 4, alpha = 0.01)
        group_analysis$total_time <- proc.time()[3] - start_time
        saveRDS(group_analysis,paste0(result_dir,'/HCP_BayesGLM2_subj_',subjects,'_',h,'_thresh',sub('\\.','',as.character(g)),'_',format(Sys.Date(),"%Y%m%d"),'.rds'))
      }
    }
  }
}




