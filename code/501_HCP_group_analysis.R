# # This is a script for running the group analysis on the HCP results
# # MULTIPLE SUBJECTS ----
# library(BayesfMRI)
# # result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group"
# result_dir <- "~/Desktop"
# # result_dir <- "/Users/Shared/HCP/5k_results/group"
# # data_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
# data_dir <- "~/github/BayesGLM_Validation/HCP_results/1k_results/PW/"
# results_all <- list.files(data_dir, full.names = T)
# # results_all <- list.files("~/github/BayesGLM_Validation/HCP_results/5k_results/individual/PW", full.names = T)
# load("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/subjects.Rdata")
# subjects <- subjects[seq(2)]
# results_sub <- sapply(subjects, grep, x = results_all, value = T, simplify = F)
# results_sub <- Reduce(c,results_sub)
# # Check that all results are available for each subject
# per_subj <- vector('numeric',length = length(subjects))
# for(s in 1:length(per_subj)) {
#   per_subj[s] <- length(grep(subjects[s],results_sub))
# }
# if(any(per_subj != 4)) stop("Not all subjects have the same number of result files.")
#
# library(INLA)
# inla.setOption(pardiso.license = "~/licenses/pardiso.lic")
# hem <- c("left")
# gamma_0 <- c(0)
#
# # 1k test, 1 cores, 2 subjects, 1k data ----
# for(h in hem) {
#   results <- grep(h, results_sub, value = T)
#   for(g in gamma_0) {
#     start_time <- proc.time()[3]
#     group_analysis <- BayesGLM2(results, excursion_type = '>', gamma = g, no_cores = NULL, alpha = 0.01)
#     group_analysis$total_time <- proc.time()[3] - start_time
#     saveRDS(group_analysis,paste0(result_dir,'/HCP_BayesGLM2_20subj_resultserial_',h,'_thresh',sub('\\.','',as.character(g)),'_',format(Sys.Date(),"%Y%m%d"),'.rds'))
#   }
# }

# SINGLE SUBJECT ----
library(INLA)
inla.setOption(pardiso.license = "~/licenses/pardiso.lic") # Dan's Mac Pro
library(BayesfMRI)
result_dir <- "/Users/Shared/Lab_Data/HCP_Motor_Task_Dan/5k_results/PW"
result_files <- list.files(result_dir, full.names = T)
load("/Users/Shared/Lab_Data/HCP_Motor_Task_Dan/subjects.Rdata")
for(subject in subjects) {
  for(hem in c("left","right")) {
    subject_files <- grep(subject,result_files, value = T)
    for(g in c(0,0.5,1)) {
      start_time <- proc.time()[3]
      subject_analysis <- BayesGLM2(grep(hem,subject_files,value = T),
                                    excursion_type = ">",
                                    gamma = g,
                                    no_cores = 4,
                                    alpha = 0.01)
      subject_analysis$total_time <- proc.time()[3] - start_time
      saveRDS(
        subject_analysis,
        file = paste0(
          "~/Desktop/501_BayesGLM2_subject_",
          subject,
          "_",
          hem,
          "_thr",
          g,
          ".rds"
        )
      )
    }
  }
}
