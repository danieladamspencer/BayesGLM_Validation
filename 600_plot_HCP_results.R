# SINGLE SUBJECT ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench/')
library(INLA)
inla.setOption(pardiso.lic = "~/pardiso.lic")
library(BayesfMRI)

# Estimates ----
result <- readRDS("HCP_results/5k_results/individual/PW/500_103818_visit1_left_5k_20210125.rds")
zlims <- c(-1,1)
plot(result$betas_Bayesian$LR, title = "Bayesian, LR, cue, 1k, PW", idx = 1, zlim = zlims)
plot(result2$betas_Bayesian$LR, title = "Bayesian, LR, cue, 1k", idx = 1, zlim = zlims)
plot(result$betas_Bayesian$RL, title = "Bayesian, RL, cue, 1k", idx = 1, zlim = zlims)
plot(result$betas_Bayesian$avg, title = "Bayesian, Average, cue, 1k", idx = 1, zlim = zlims)
plot(result$betas_Bayesian$LR, title = "Bayesian, LR, right foot, 1k", idx = 2, zlim = zlims)
plot(result$betas_Bayesian$RL, title = "Bayesian, RL, right foot, 1k", idx = 2, zlim = zlims)
plot(result$betas_Bayesian$avg, title = "Bayesian, Average, right foot, 1k", idx = 2, zlim = zlims)
plot(result$betas_Bayesian$LR, title = "Bayesian, LR, right hand, 1k", idx = 3, zlim = zlims)
plot(result$betas_Bayesian$RL, title = "Bayesian, RL, right hand, 1k", idx = 3, zlim = zlims)
plot(result$betas_Bayesian$avg, title = "Bayesian, Average, right hand, 1k", idx = 3, zlim = zlims)
plot(result$betas_Bayesian$LR, title = "Bayesian, LR, tongue, 1k", idx = 4, zlim = zlims)
plot(result$betas_Bayesian$RL, title = "Bayesian, RL, tongue, 1k", idx = 4, zlim = zlims)
plot(result$betas_Bayesian$avg, title = "Bayesian, Average, tongue, 1k", idx = 4, zlim = zlims)

# Activations ----
result <- readRDS("HCP_results/5k_results/individual/PW/500_103818_visit1_left_5k_20210125.rds")
act_single_subject <- id_activations.posterior(model_obj = result$GLMs_Bayesian$cortexL,field_name = NULL, threshold = 0,alpha = 0.01,area.limit = NULL)
cifti_temp <- result$betas_Bayesian$avg
cifti_temp$data$cortex_left <- act_single_subject$LR$active
plot(cifti_temp, title = "PW")
result_nonPW <- readRDS("HCP_results/5k_results/individual/103818_visit1_left_5k_20201120.rds")
act_single_subject_nonPW <- id_activations.posterior(model_obj = result_nonPW$GLMs_Bayesian$cortexL, field_name = NULL, threshold = 0, alpha = 0.01, area.limit = NULL)
cifti_temp2 <- result_nonPW$betas_Bayesian$avg
cifti_temp2$data$cortex_left <- act_single_subject_nonPW$LR$active
plot(cifti_temp2, title = "Non-PW")

# AR values ----
result <- readRDS("HCP_results/5k_results/individual/PW/500_103818_visit1_left_5k_20210125.rds")
cifti_AR <- result$betas_Bayesian$avg
cifti_AR$data$cortex_left <- result$prewhitening_info$left$AR_coeffs
for(p in 1:6) plot(cifti_AR, idx = p, title = paste("AR order",p), zlim = c(-.2,.2))
cifti_AR$data$cortex_left <- as.matrix(result$prewhitening_info$left$AR_var)
plot(cifti_AR,title = "AR Residual Variance")

# Estimates for multiple individual subjects ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench/')
library(INLA)
inla.setOption(pardiso.lic = "~/pardiso.lic")
library(BayesfMRI)
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
load("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
result_files <- list.files(result_dir, full.names = T)
result_files <- c(sapply(subjects,grep, x = result_files, value = T))
vis <- "visit1"
task_idx <- 4
task_name <- "Tongue"
result_files <- grep(vis, result_files, value = T)
# >> Bayesian ----
for(s in subjects) {
  cifti_template <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
  left_result <- readRDS(grep('left',grep(s,result_files, value = T), value = T))
  cifti_template$data$cortex_left <- left_result$betas_Bayesian$avg$data$cortex_left
  rm(left_result)
  right_result <- readRDS(grep('right',grep(s,result_files, value = T), value = T))
  cifti_template$data$cortex_right <- right_result$betas_Bayesian$avg$data$cortex_right
  rm(right_result)
  plot(cifti_template, idx = task_idx, #title = paste("Subject",s,task_name,"Task"),
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
       fname = paste0("plots/600_subject_",s,"_",tolower(task_name),"_estimates.png"),
       zlim = c(-1,1))
}
# >> Classical ----
for(s in subjects) {
  cifti_template <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
  left_result <- readRDS(grep('left',grep(s,result_files, value = T), value = T))
  cifti_template$data$cortex_left <- left_result$betas_classical$avg$data$cortex_left
  rm(left_result)
  right_result <- readRDS(grep('right',grep(s,result_files, value = T), value = T))
  cifti_template$data$cortex_right <- right_result$betas_classical$avg$data$cortex_right
  rm(right_result)
  plot(cifti_template, idx = task_idx, #title = paste("Subject",s,task_name,"Task"),
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
       fname = paste0("plots/600_subject_",s,"_",tolower(task_name),"_classical_estimates.png"),
       zlim = c(-1,1))
}

# GROUP RESULTS ----
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW"
result_files <- grep("_20210131.rds",list.files(result_dir), value = T)

# Estimate Plots ----
hem <- c("left","right")
cifti_obj <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
cifti_obj$data$cortex_left <- readRDS(file.path(result_dir,"HCP_BayesGLM2_45subj_result_left_thresh0_20210131.rds"))$estimates
cifti_obj$data$cortex_right <- readRDS(file.path(result_dir,"HCP_BayesGLM2_45subj_result_right_thresh0_20210131.rds"))$estimates
library(ciftiTools)
ciftiTools.setOption('wb_path', "/Applications/workbench")

# Grab analysis times from the result files ----
library(tidyverse)
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
result_files <- list.files(result_dir) %>%
  grep(pattern = "500_", x = ., value = T)
left_results <- grep("_left_", result_files, value = T)
right_results <- grep("_right_", result_files, value = T)
left_times <- sapply(left_results, function(file_n) {
  readRDS(file.path(result_dir, file_n))$total_time
})
right_times <- sapply(right_results, function(file_n) {
  readRDS(file.path(result_dir, file_n))$total_time
})
summary(left_times) / 3600
summary(right_times) / 3600
all_times <- c(rbind(left_times,right_times))
summary(all_times) / 60
sd(all_times / 3600)
hist(all_times/3600)
