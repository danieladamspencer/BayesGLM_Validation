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
# result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
result_dir <- "HCP_results/5k_results"
result_files <- list.files(result_dir, full.names = T)
result_files <- grep("503_", result_files, value = T)
result_files <- grep("_alpha001_", result_files, value = T)
result_files <- grep("_thr01_", result_files, value = T, invert = T)
load("HCP_data/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
threshs <- paste0("_thr",c("0","05","1"),"_")
hems <- c("left","right")
# library(viridisLite)
# col_pal <- rev(plasma(5)[3:5])
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
col_pal <- c(light_orange,"red","purple")
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
       # title = paste("Subject",subject, "Tongue Task"),
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
       fname = paste0('plots/603_subject_',subject,"_tongue_task_activations.png"))
}

# >> Classical ----
# As a note, this would not be necessary in the case of the current (2021-02-10)
# state of the BayesfMRI package, but some retesting using the classical GLM
# is necessary in order to obtain the standard errors of the estimates.
results_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/classical"
load("HCP_data/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
result_files <- list.files(results_dir,full.names = T)
result_files <- c(unlist(sapply(subjects, grep, x = result_files, value = T)))
result_files <- grep("visit1",result_files, value = T)
# result_files <- grep("classical",result_files, value = T)
task_idx <- 4
library(BayesfMRI)
library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path',"/Applications/workbench")
# subject <- subjects[1]
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# col_pal <- gg_color_hue(3)[2:3]
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
col_pal <- c("yellow",light_orange,"red","purple")
for(subject in subjects) {
  subject_files <- grep(subject,result_files, value = T)
  left_result <- readRDS(grep('left',subject_files, value=T))
  left_FDR <-
    id_activations.classical(
      left_result$GLMs_classical$cortexL,
      alpha = 0.01,
      correction = "FDR",
      field_inds = task_idx
    )$active
  left_FWER <- sapply(c(0, 0.5, 1), function(thr) {
    out <-
      id_activations.classical(
        left_result$GLMs_classical$cortexL,
        alpha = 0.01,
        correction = "FWER",
        field_inds = task_idx,
        threshold = thr
      )$active
    return(out)
  }, simplify = T)
  left_FWER <- as.matrix(apply(left_FWER,1,sum))
  left_active <- left_FDR + left_FWER
  medial_left <- !is.na(left_active[,1])
  left_active <- as.matrix(left_active[,1][!is.na(left_active[,1])])
  left_active[,1] <- ifelse(left_active[,1] == 0,NA,left_active[,1])
  right_result <- readRDS(grep('right',subject_files, value=T))
  right_FDR <-
    id_activations.classical(
      right_result$GLMs_classical$cortexR,
      alpha = 0.01,
      correction = "FDR",
      field_inds = task_idx
    )$active
  right_FWER <- sapply(c(0, 0.5, 1), function(thr) {
    out <-
      id_activations.classical(
        right_result$GLMs_classical$cortexR,
        alpha = 0.01,
        correction = "FWER",
        field_inds = task_idx,
        threshold = thr
      )$active
    return(out)
  }, simplify = T)
  right_FWER <- as.matrix(apply(right_FWER,1,sum))
  right_active <- right_FDR + right_FWER
  medial_right <- !is.na(right_active[,1])
  right_active <- as.matrix(right_active[,1][!is.na(right_active[,1])])
  right_active[,1] <- ifelse(right_active[,1] == 0,NA,right_active[,1])
  cifti_active <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
  cifti_active$data$cortex_left <- left_active
  cifti_active$data$cortex_right <- right_active
  plot(cifti_active, color_mode = "qualitative", colors = col_pal,
       # title = paste("Subject",subject, "Tongue Task"),
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
       fname = paste0('plots/603_subject_',subject,"_tongue_task_activations_classical.png"))
}

# >> Add Bayesian and Classical activations to examine consistency ----
# >>>> Bayesian ----
result_dir <- "HCP_results/5k_results/individual/PW/activations"
load("HCP_data/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
task_idx <- 4 # This is for the tongue
thr <- 0.5
thr_chr <- paste0("_thr",sub("\\.","",thr))
alpha <- 0.01
alpha_chr <- paste0("_alpha",sub("\\.","",alpha))
# library(viridisLite)
# col_pal <- rev(plasma(5)[3:4])
# library(scales)
# col_pal <- alpha("red", alpha = c(0.1,0.2))
# col_pal <- sapply(col_pal, substring, first = 4)
col_pal <- c("pink","red")
subject_files <- grep("503_", list.files(result_dir, full.names = T), value = T)
subject_files <- unlist(sapply(subjects, grep, x = subject_files, value = T))
subject_files <- grep(thr_chr, subject_files,  value=T)
subject_files <- grep("_classical.rds",subject_files, value = T, invert = T)
for(subject in subjects) {
  for(h in c('left','right')){
    for(v in 1:2) {
      subject_file <- grep(paste0("_visit",v,"_subject_",subject,"_",h), subject_files, value = T)
      subject_act <- readRDS(subject_file)
      L_or_R <- toupper(substring(h,1,1))
      if(v == 1){
        if(h == 'left') {
          subject_cifti <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
        }
        subject_cifti$data[[paste0("cortex_",h)]] <- as.matrix(subject_act$active[,task_idx])
      }
      if(v == 2) {
        subject_cifti$data[[paste0("cortex_",h)]] <- as.matrix(subject_act$active[,task_idx]) +
          subject_cifti$data[[paste0("cortex_",h)]]
      }
    }
    subject_cifti$data[[paste0("cortex_",h)]] <-
      as.matrix(
        ifelse(subject_cifti$data[[paste0("cortex_", h)]] == 0,
               NA,
               subject_cifti$data[[paste0("cortex_", h)]]
               )
        )
  }
  plot(subject_cifti, fname = paste0("plots/603_subject_",subject,thr_chr,alpha_chr,"_num_active.png"),
       # title = paste("Subject",subject, "Tongue Task\nVisit Activations"),
       color_mode = "qualitative",colors = col_pal,
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")
}

# >>>> Classical ----
library(ciftiTools)
ciftiTools.setOption('wb_path',"/Applications/workbench")
result_dir <- "HCP_results/5k_results/individual/PW/activations"
load("HCP_data/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
task_idx <- 4 # This is for the tongue
thr <- 0.5
thr_chr <- paste0("_thr",sub("\\.","",thr))
alpha <- 0.01
alpha_chr <- paste0("_alpha",sub("\\.","",alpha))
# library(viridisLite)
# col_pal <- rev(plasma(5)[3:4])
# library(scales)
# col_pal <- alpha("red", alpha = c(0.1,0.2))
col_pal <- c("pink","red")
subject_files <- grep("503_", list.files(result_dir, full.names = T), value = T)
subject_files <- unlist(sapply(subjects, grep, x = subject_files, value = T))
subject_files <- grep(thr_chr, subject_files,  value=T)
subject_files <- grep("_classical.rds",subject_files, value = T, invert = F)
for(subject in subjects) {
  for(h in c('left','right')){
    for(v in 1:2) {
      subject_file <- grep(paste0("_visit",v,"_subject_",subject,"_",h), subject_files, value = T)
      subject_act <- readRDS(subject_file)
      L_or_R <- toupper(substring(h,1,1))
      if(v == 1){
        if(h == 'left') {
          subject_cifti <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
        }
        in_mask <- subject_cifti$meta$cortex$medial_wall_mask[[h]]
        subject_cifti$data[[paste0("cortex_",h)]] <-
          as.matrix(subject_act$active[in_mask,task_idx])
      }
      if(v == 2) {
        subject_cifti$data[[paste0("cortex_",h)]] <- as.matrix(subject_act$active[in_mask,task_idx]) +
          subject_cifti$data[[paste0("cortex_",h)]]
      }
    }
    subject_cifti$data[[paste0("cortex_",h)]] <-
      as.matrix(
        ifelse(subject_cifti$data[[paste0("cortex_", h)]] == 0,
               NA,
               subject_cifti$data[[paste0("cortex_", h)]]
        )
      )
  }
  plot(subject_cifti, fname = paste0("plots/603_subject_",subject,thr_chr,alpha_chr,"_num_active_classical.png"),
       # title = paste("Subject",subject, "Tongue Task\nVisit Activations"),
       color_mode = "qualitative",colors = col_pal,
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")
}

# Density plots for active locations ----

# >> Bayesian ----
library(INLA)
library(tidyverse)
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
result_files <- list.files(result_dir, full.names = T)
result_files <- grep("500_", result_files, value = T)
act_dir <- "HCP_results/5k_results"
act_files <- list.files(act_dir, full.names = T)
act_files <- grep("503_", act_files, value = T)
act_files <- grep("_alpha001_", act_files, value = T)
act_files <- grep("_thr01_", act_files, value = T, invert = T)
load("HCP_data/subjects.Rdata")
subjects <- subjects[seq(4)]
task_idx <- 4
my_cols <- c("orange","red","purple")
tnum <- c(0,0.5,1)
threshs <- paste0("_thr", sub("\\.","",tnum),"_")
result_files <- c(sapply(subjects, function(subj) {
  out <- grep(paste0("_",subj,"_"),result_files, value = T)
  return(out)
}, simplify = T))
result_files <- grep("_visit1_", result_files, value = T)
library(INLA)
# sapply(subjects, function(subj) {
for(subj in subjects){
  subject_files <- grep(paste0("_",subj,"_"), result_files, value = T)
  act_files_s <- grep(paste0("_",subj,"_"), act_files, value = T)
  act_files_left <- grep("left",act_files_s, value = T)
  act_files_right <- grep("right", act_files_s, value = T)
  subject_files_left <- grep("left",subject_files, value = T)
  subject_files_right <- grep("right",subject_files, value = T)
  left_act <- sapply(threshs, function(thr) {
    result <- readRDS(grep(thr,act_files_left, value = T))
    return(result$active[,task_idx])
  })
  left_act <- apply(left_act,1,sum)
  left_act <- as.matrix(ifelse(left_act == 0, NA, left_act))
  right_act <- sapply(threshs, function(thr) {
    result <- readRDS(grep(thr,act_files_right, value = T))
    return(result$active[,task_idx])
  })
  right_act <- apply(right_act,1,sum)
  right_act <- as.matrix(ifelse(right_act == 0, NA, right_act))
  active_idx <- list(
    left = which(!is.na(left_act)),
    right = which(!is.na(right_act))
  )
  left_result <- readRDS(grep("_left_", subject_files, value = T))
  marginal_obj_left <-
    left_result$GLMs_Bayesian$cortexL$INLA_result$marginals.random[[paste0("bbeta", task_idx)]][active_idx$left]
  right_result <- readRDS(grep("_right_", subject_files, value = T))
  marginal_obj_right <-
    right_result$GLMs_Bayesian$cortexR$INLA_result$marginals.random[[paste0("bbeta", task_idx)]][active_idx$right]
  active_densities <- list(
    left = sapply(marginal_obj_left, function(x) as.matrix(as.data.frame(inla.smarginal(x,factor = 2))), simplify = F),
    right = sapply(marginal_obj_right, function(x) as.matrix(as.data.frame(inla.smarginal(x,factor = 2))), simplify = F)
  )
  active_color <-
    list(
      left = left_act,
      right = right_act
    ) %>%
    reshape2::melt() %>%
    filter(!is.na(value)) %>%
    mutate(L2 = paste0("index.",Var1),
           color_act = paste0(tnum[value],"%"),
           color_act = factor(color_act,levels = c("0%","0.5%","1%"))) %>%
    select(L1,L2,color_act)
  dens_df <- reshape2::melt(active_densities) %>%
    tidyr::pivot_wider(names_from = Var2, values_from = value) %>%
    left_join(active_color)
  marg_dens_plot <- dens_df %>%
    filter(color_act != "0%") %>%
    ggplot() +
    geom_line(aes(x = x, y = y, color = color_act, group = L2)) +
    scale_color_manual(expression(paste(gamma," =")), values = my_cols[-1]) +
    facet_grid(~L1, scales = "free_x") +
    labs(
      x = "Activation Amplitude",
      y = "Marginal Posterior Density",
      title = paste("Subject", subj, "Tongue Task")
    ) +
    theme_bw()
  ggsave(
    paste0(
      "plots/603_subject_",
      subj,
      "_tongue_task_marginal_posterior_densities.png"
    ),
    plot = marg_dens_plot, width = 6, height = 2.5
  )
}
