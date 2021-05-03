library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')

# From the single-session results ----
est_ss <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
left_result_ss <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples/501_BayesGLM2_10subjset1_single_session_visit1_results__left_thresh0_20210426.rds")
right_result_ss <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples/501_BayesGLM2_10subjset1_single_session_visit1_results__right_thresh0_20210426.rds")
est_ss$data$cortex_left <- left_result_ss$estimates
est_ss$data$cortex_right <- right_result_ss$estimates
plot(est_ss, idx = 4, zlim = c(-1,1), title = "Tongue Estimate from Single Sessions",
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")

act_ss <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
left_result_ss$active[left_result_ss$active == 0] <- NA
right_result_ss$active[right_result_ss$active == 0] <- NA
act_ss$data$cortex_left <- left_result_ss$active
act_ss$data$cortex_right <- right_result_ss$active
plot(act_ss, idx = 4, title = "Tongue 0% Activations from Single Sessions",
     color_mode = "qualitative",
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")


# From the multi-session results ----
est_ms <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
left_result_ms <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples/501_HCP_BayesGLM2_10subj_sample1_visit1_result_left_thresh0_20210319.rds")
right_result_ms <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples/501_HCP_BayesGLM2_10subj_sample1_visit1_result_right_thresh0_20210319.rds")
est_ms$data$cortex_left <- left_result_ms$estimates
est_ms$data$cortex_right <- right_result_ms$estimates
plot(est_ms, idx = 4, zlim = c(-1,1), title = "Tongue Estimate from Multiple Sessions",
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")

act_ms <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
left_result_ms$active[left_result_ms$active == 0] <- NA
right_result_ms$active[right_result_ms$active == 0] <- NA
act_ms$data$cortex_left <- left_result_ms$active
act_ms$data$cortex_right <- right_result_ms$active
plot(act_ms, idx = 4, title = "Tongue 0% Activations from Multiple Sessions",
     color_mode = "qualitative",
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")
