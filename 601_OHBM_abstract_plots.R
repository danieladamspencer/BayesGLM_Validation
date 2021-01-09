# This is a script to produce plots for the OHBM abstract submission
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench/')
library(INLA)
inla.setOption(pardiso.lic = "~/pardiso.lic")
library(BayesfMRI)

# HCP SINGLE SUBJECT ----
results_dir <- "~/github/BayesGLM_Validation/HCP_results"
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
zlims <- list(c(-1.5,1.5),c(-1,1),c(-2,2),c(-2.5,2.5))

# Bayesian plots
for(i in seq(num_betas)) {
  view_xifti_surface(
    xifti = single_subject_result$betas_Bayesian$avg,
    idx = i,
    hemisphere = plot_hem,
    view = "both",
    color_mode = "diverging",
    zlim = zlims[[i]],
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
    # title = task_names[i],
    fname = paste0("plots/601_single_subject_classical_", file_task_names[i]),
    save = TRUE
  )
}
