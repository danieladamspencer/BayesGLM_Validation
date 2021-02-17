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
# >> Bayesian ----
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
result_files <- list.files(result_dir, full.names = T)
result_files <- grep("500_", result_files, value = T)
result_files <- grep("_alpha001_", result_files, value = T)
result_files <- grep("_thr01_", result_files, value = T, invert = T)
load("HCP_data/subjects.Rdata")
subjects <- subjects[seq(5)]
threshs <- paste0("_thr",c("0","05","1"),"_")
hems <- c("left","right")
library(viridisLite)
col_pal <- rev(plasma(5)[3:5])
for(subject in subjects) {
  subject_files <- grep(subject, result_files, value =T)
  subject_files_left <- grep("left",subject_files, value = T)
  subject_files_right <- grep("right",subject_files, value = T)
  left_act <- sapply(threshs, function(thr) {
    result <- readRDS(grep(thr,subject_files_left, value = T))
    return(result$active[,4])
  })
  left_act <- apply(left_act,1,sum)
  left_act <- as.matrix(ifelse(left_act == 0, NA, left_act))
  right_act <- sapply(threshs, function(thr) {
    result <- readRDS(grep(thr,subject_files_right, value = T))
    return(result$active[,4])
  })
  right_act <- apply(right_act,1,sum)
  right_act <- as.matrix(ifelse(right_act == 0, NA, right_act))
  cifti_obj <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
  cifti_obj$data$cortex_left <- left_act
  cifti_obj$data$cortex_right <- right_act
  plot(cifti_obj, color_mode = "qualitative", colors = col_pal,
       title = paste("Subject",subject, "Tongue Task"),
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
       save = paste0('plots/603_subject_',subject,"_tongue_task_activations.png"))
}

# >> Classical ----
# As a note, this would not be necessary in the case of the current (2021-02-10)
# state of the BayesfMRI package, but some retesting using the classical GLM
# is necessary in order to obtain the standard errors of the estimates.
results_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/classical"
load("HCP_data/subjects.Rdata")
subjects <- subjects[seq(5)]
result_files <- list.files(results_dir,full.names = T)
result_files <- c(sapply(subjects, grep, x = result_files, value = T))
result_files <- grep("visit1",result_files, value = T)
task_idx <- 4
library(BayesfMRI)
library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path',"/Applications/workbench")
subject <- subjects[1]
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
col_pal <- gg_color_hue(3)[2:3]
for(subject in subjects) {
  subject_files <- grep(subject,result_files, value = T)
  left_result <- readRDS(grep('left',subject_files, value=T))
  left_FDR <- id_activations.classical(left_result$GLMs_classical$cortexL, alpha = 0.1, correction = "FDR", field_inds = task_idx)
  left_FWER <- id_activations.classical(left_result$GLMs_classical$cortexL, alpha = 0.1, correction = "FWER", field_inds = task_idx)
  left_active <- left_FDR$active + left_FWER$active
  medial_left <- !is.na(left_active[,1])
  left_active <- as.matrix(left_active[,1][!is.na(left_active[,1])])
  left_active[,1] <- ifelse(left_active[,1] == 0,NA,left_active[,1])
  right_result <- readRDS(grep('right',subject_files, value=T))
  right_FDR <- id_activations.classical(right_result$GLMs_classical$cortexR, alpha = 0.1, correction = "FDR", field_inds = task_idx)
  right_FWER <- id_activations.classical(right_result$GLMs_classical$cortexR, alpha = 0.1, correction = "FWER", field_inds = task_idx)
  right_active <- right_FDR$active + right_FWER$active
  medial_right <- !is.na(right_active[,1])
  right_active <- as.matrix(right_active[,1][!is.na(right_active[,1])])
  right_active[,1] <- ifelse(right_active[,1] == 0,NA,right_active[,1])
  cifti_active <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
  cifti_active$data$cortex_left <- left_active
  cifti_active$data$cortex_right <- right_active
  plot(cifti_active, color_mode = "qualitative", colors = col_pal,
       title = paste("Subject",subject, "Tongue Task"),
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
       fname = paste0('plots/603_subject_',subject,"_tongue_task_activations_classical.png"))
}

