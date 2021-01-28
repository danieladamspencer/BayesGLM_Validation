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
