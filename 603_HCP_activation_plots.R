# This is a script meant to plot the activation areas at different thresholds
# for the HCP analysis results
# GROUP LEVEL ----
# Load libraries and grab result file names ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
result_dir <- "HCP_results/5k_results/group"
result_files <- list.files(result_dir, full.names = T)

# Make a cifti template at the 5k resolution ----
# This is being done as a general-purpose template to be able to add and plot
# data that aren't in a cifti (xifti) object.
# cifti_5k_template_left <-
#   readRDS("HCP_results/5k_results/individual/103818_visit1_left_5k_20201120.rds")$betas_Bayesian[[1]]
# cifti_5k_template_right <-
#   readRDS("HCP_results/5k_results/individual/103818_visit1_right_5k_20201120.rds")$betas_Bayesian[[1]]
# cifti_5k_template_whole <- cifti_5k_template_left
# cifti_5k_template_whole$data$cortex_right <-
#   cifti_5k_template_right$data$cortex_right
# cifti_5k_template_whole$meta$cortex$medial_wall_mask$right <-
#   cifti_5k_template_right$meta$cortex$medial_wall_mask$right
# saveRDS(cifti_5k_template_whole,"HCP_data/603_cifti_5k_template_whole.rds")

# Plot each task's activation ----
threshs <- c(0,0.5,1)
hem <- c('left','right')
thresh_chr <- sub("\\.","",as.character(threshs))
task_names <- c('visual cue','foot','hand','tongue')
activation_maps <- vector('list',length(task_names))
names(activation_maps) <- task_names
# Changing this to plot both hemispheres of the tongue task ONLY
for(task in c(4)) {
  # Load a template for the task
  activation_maps[[task]] <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
  for(h in hem) {
    # Set the data to match the dimensions needed for each task
    activation_maps[[task]]$data[[paste0('cortex_',h)]] <-
      matrix(NA,
             nrow = nrow(activation_maps[[task]]$data[[paste0('cortex_',h)]]),
             ncol = length(threshs))
    for(thr in seq(length(threshs))) {
      fname_hem_thr <-
        grep(paste0(h,"_thresh",thresh_chr[thr],"_"), result_files, value = T)
      group_results_hem_thr <- readRDS(fname_hem_thr)
      activation_maps[[task]]$data[[paste0('cortex_',h)]][,thr] <-
        ifelse(group_results_hem_thr$active[,task] == 1, threshs[thr], NA)
    }
    activation_maps[[task]]$data[[paste0('cortex_',h)]] <-
      as.matrix(
        apply(
          activation_maps[[task]]$data[[paste0('cortex_', h)]],
          1,
          function(x) {
            max_x <- max(x,na.rm = TRUE)
            max_x <- ifelse(max_x == -Inf, NA, max_x)
            return(max_x)
          }))
  }
}

for(task in task_names[4]) {
  plot(activation_maps[[task]], color_mode = 'qualitative', hemisphere = 'both', colors = c("yellow","orange"),
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
       save = TRUE, fname = paste0("plots/603_",sub(" ","_", task),"_activation_map"))
}

# SUBJECT LEVEL ----
# Load libraries and grab result file names ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
result_dir <- "HCP_results/5k_results/group"
result_files <- list.files(result_dir, full.names = T)
