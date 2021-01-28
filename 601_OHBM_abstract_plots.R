# This is a script to produce plots for the OHBM abstract submission
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench/')
library(INLA)
inla.setOption(pardiso.lic = "~/pardiso.lic")
library(BayesfMRI)

# HCP SINGLE SUBJECT ----
results_dir <- "~/github/BayesGLM_Validation/HCP_results/5k_results/individual"
result_files <- list.files(results_dir)
plot_subject <- "103818"
result_files <- grep(plot_subject,result_files, value = TRUE)
plot_hem <- "left"
result_files <- grep(plot_hem,result_files, value = TRUE)
plot_visit <- "visit1"
result_files <- grep(plot_visit,result_files,value = TRUE)

single_subject_result <- readRDS(file = file.path(results_dir,result_files))
num_betas <- ncol(single_subject_result$betas_Bayesian$LR$data$cortex_left)
task_names <- c("Visual Cue","Right Foot","Right Hand","Tongue")
file_task_names <- tolower(gsub(" ","_",task_names))
zlims <- list(c(-1,1),
              c(-1,1),
              c(-1,1),
              c(-1,1))

# >> Plot Estimates ----

# Bayesian plots
for(i in seq(num_betas)) {
  view_xifti_surface(
    xifti = single_subject_result$betas_Bayesian$avg,
    idx = i,
    hemisphere = plot_hem,
    view = "both",
    color_mode = "diverging",
    zlim = zlims[[i]],
    surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
    surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
    # title = task_names[i],
    fname = paste0("plots/601_single_subject_Bayes_", file_task_names[i]),
    save = TRUE
  )
}

# Classical_plots
for(i in seq(num_betas)) {
  view_xifti_surface(
    xifti = single_subject_result$betas_classical$avg,
    idx = i,
    hemisphere = plot_hem,
    view = "both",
    color_mode = "diverging",
    zlim = zlims[[i]],
    surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
    surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
    # title = task_names[i],
    fname = paste0("plots/601_single_subject_classical_", file_task_names[i]),
    save = TRUE
  )
}

# HCP GROUP PLOTS ----
group_analysis_left2 <- readRDS("~/github/BayesGLM_Validation/HCP_results/5k_results/group/HCP_BayesGLM2_result_left_thresh0_20210105.rds")
group_estimates_classical <- readRDS("HCP_results/5k_results/group/502_HCP_classical_group_estimates.rds")
library(ciftiTools)
# ciftiTools.setOption('wb_path',"/Users/danspen/workbench")
ciftiTools.setOption('wb_path',"/Applications/workbench")

# >> Plot Estimates ----
indiv_results_left <- readRDS("~/github/BayesGLM_Validation/HCP_results/5k_results/individual/103818_visit1_left_5k_20201120.rds")
# Bayesian
cifti_left <- indiv_results_left$betas_Bayesian$avg
cifti_left$data$cortex_left <- group_analysis_left$estimates
zlims <- list(c(-1.2,1.2),
              c(-1,1),
              c(-1,1),
              c(-1,1))

for(i in seq(num_betas)) {
  view_xifti_surface(
    xifti = cifti_left,
    idx = i,
    hemisphere = plot_hem,
    view = "both",
    color_mode = "diverging",
    zlim = zlims[[i]],
    surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
    surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
    # title = task_names[i],
    fname = paste0("plots/601_group_Bayes_", file_task_names[i]),
    save = TRUE
  )
}
# Classical
cifti_left_classical <- indiv_results_left$betas_Bayesian$avg
cifti_left_classical$data$cortex_left <- group_estimates_classical$left
for(i in seq(num_betas)) {
  view_xifti_surface(
    xifti = cifti_left_classical,
    idx = i,
    hemisphere = plot_hem,
    view = "both",
    color_mode = "diverging",
    zlim = zlims[[i]],
    surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
    surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
    # title = task_names[i],
    fname = paste0("plots/601_group_classical_", file_task_names[i]),
    save = TRUE
  )
}

# >> Plot Activations ----
# indiv_results_left <- readRDS("/Volumes/Macintosh HD/Users/Shared/HCP/5k_results/103818_visit1_left_5k_20201120.rds")
indiv_results_left <- readRDS("~/github/BayesGLM_Validation/HCP_results/5k_results/individual/103818_visit1_left_5k_20201120.rds")
active_left <- indiv_results_left$betas_Bayesian$avg
active_left$data$cortex_left <- group_analysis_left$active
plot(active_left, idx = 4)
