# This is a script to produce plots for the OHBM abstract submission
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench/')
library(INLA)
inla.setOption(pardiso.lic = "~/pardiso.lic")
library(BayesfMRI)

# HCP SINGLE SUBJECT ----
# >> Estimates ----
# >>>> Bayes ----
## This section will not be used ##
# results_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/"
# result_files <- list.files(results_dir)
# plot_subject <- "103818"
# result_files <- grep(plot_subject,result_files, value = TRUE)
# plot_hem <- "left"
# result_files <- grep(plot_hem,result_files, value = TRUE)
# plot_visit <- "visit1"
# result_files <- grep(plot_visit,result_files,value = TRUE)
#
# single_subject_result <- readRDS(file = file.path(results_dir,result_files))
# num_betas <- ncol(single_subject_result$betas_Bayesian$LR$data$cortex_left)
# task_names <- c("Visual Cue","Right Foot","Right Hand","Tongue")
# file_task_names <- tolower(gsub(" ","_",task_names))
# zlims <- list(c(-1,1),
#               c(-1,1),
#               c(-1,1),
#               c(-1,1))
# Bayesian plots
# for(i in seq(num_betas)) {
#   view_xifti_surface(
#     xifti = single_subject_result$betas_Bayesian$avg,
#     idx = i,
#     hemisphere = plot_hem,
#     view = "both",
#     color_mode = "diverging",
#     zlim = zlims[[i]],
#     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
#     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
#     # title = task_names[i],
#     fname = paste0("plots/601_single_subject_Bayes_", file_task_names[i]),
#     save = TRUE
#   )
# }
#
#
#
# # Classical_plots
# for(i in seq(num_betas)) {
#   view_xifti_surface(
#     xifti = single_subject_result$betas_classical$avg,
#     idx = i,
#     hemisphere = plot_hem,
#     view = "both",
#     color_mode = "diverging",
#     zlim = zlims[[i]],
#     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
#     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
#     # title = task_names[i],
#     fname = paste0("plots/601_single_subject_classical_", file_task_names[i]),
#     save = TRUE
#   )
# }
## End ##

result_left <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_103818_visit1_left_5k_20210125.rds")
result_right <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_103818_visit1_right_5k_20210125.rds")

result_both <- result_left$betas_Bayesian$avg
result_both$data$cortex_right <- result_right$betas_Bayesian$avg$data$cortex_right
result_both$meta$cortex$medial_wall_mask$right <- result_right$betas_Bayesian$avg$meta$cortex$medial_wall_mask$right

rm(result_left,result_right)
plot(result_both, idx = 4, zlim = c(-1,1), save = TRUE, fname = "plots/601_subject_bayes_tongue_estimate.png",
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")

# >>>> Classical ----
result_left <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_103818_visit1_left_5k_20210125.rds")
result_right <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_103818_visit1_right_5k_20210125.rds")

result_both <- result_left$betas_classical$avg
result_both$data$cortex_right <- result_right$betas_classical$avg$data$cortex_right
result_both$meta$cortex$medial_wall_mask$right <- result_right$betas_classical$avg$meta$cortex$medial_wall_mask$right

rm(result_left,result_right)
plot(result_both, idx = 4, zlim = c(-1,1), save = TRUE, fname = "plots/601_subject_classical_tongue_estimate.png",
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")

# >> Activations ----
# >>>> Bayes ----
bayes_activation_files <- grep("503", list.files("HCP_results/5k_results", full.names = T), value = T)
thresholds <- c(0,0.5,1)
chr_thr <- paste0("_thr",sub("\\.","",as.character(thresholds)),"_")
bayes_activation_files <- c(sapply(chr_thr, grep, x = bayes_activation_files, value = T))
activations <- sapply(c('left','right'), function(hem) {
  sapply(thresholds, function(thr) {
    fn <-
      grep(paste0("thr", sub("\\.", "", as.character(thr)),"_"),
           grep(hem, bayes_activation_files, value = T) ,
           value = T)
    act_obj <- readRDS(fn)
    out_act <- ifelse(act_obj$active[,4] == 1, 1, NA)
    return(out_act)
  })
}, simplify = F)
activations$left <- apply(activations$left,1,sum,na.rm = T)
activations$left <- as.matrix(ifelse(activations$left == 0, NA, activations$left))
activations$right <- apply(activations$right,1,sum,na.rm = T)
activations$right <- as.matrix(ifelse(activations$right == 0, NA, activations$right))

subject_act_cifti <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
subject_act_cifti$data$cortex_left <- activations$left
subject_act_cifti$data$cortex_right <- activations$right

library(viridisLite)
col_pal <- rev(plasma(5)[3:5])

plot(subject_act_cifti, color_mode = "qualitative", colors = col_pal,
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
     save = TRUE, fname = "plots/601_subject_bayes_tongue_activation.png")

# >>>> Classical ----
class_act <- readRDS("HCP_results/5k_results/503_classical_tongue_activations.rds")
class_cifti <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
class_cifti$data$cortex_left <- class_act$left
class_cifti$data$cortex_right <- class_act$right

# library(viridisLite)
# col_pal <- viridis(2)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
col_pal <- gg_color_hue(3)[2:3]

plot(class_cifti, color_mode = "qualitative", colors = col_pal,
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
     save = TRUE, fname = "plots/601_subject_classical_tongue_activation.png")

# >> ICC plots ----
# ICC function definition ----
#' Intraclass Correlation
#'
#' @param X A matrix with N (# subjects) rows and K (# observations per subject)
#'   columns.
#'
#' @return The scalar value for the intraclass correlation
#' @export
#'
#' @examples
ICC <- function(X) {
  xbar <- mean(X)
  s2 <- sum((X - xbar)^2) / (prod(dim(X)))
  xbar_i <- apply(X,1,mean)
  r <- (ncol(X)*sum((xbar_i - xbar)^2) / (nrow(X)*s2) - 1) / (ncol(X) - 1)
  if(r < 0) r <- 0
  return(r)
}

# >>>> Bayes ----
avg_estimates <- readRDS("HCP_results/5k_results/602_avg_PW_estimates.rds")
library(abind)
abind_g <- function(...) abind(...,along = 4)
ICC_values_average <- sapply(avg_estimates, function(hem_est) {
  combined_visits <- Reduce(abind_g,hem_est)
  vertex_ICC <- apply(combined_visits, 1:2, ICC)
  return(vertex_ICC)
}, simplify = F)
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
plot_dir <- "~/github/BayesGLM_Validation/HCP_results/plots"
task_file_names <- c("visual_cue","foot","hand","tongue")
bayesian_icc_avg_cifti <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
bayesian_icc_avg_cifti$data$cortex_left <- ICC_values_average$left
bayesian_icc_avg_cifti$data$cortex_right <- ICC_values_average$right
library(viridisLite)
my_cols  <- viridis(4)
plot(bayesian_icc_avg_cifti, idx = 4, color_mode = "sequential",
     colors = my_cols, zlim = c(0,0.25,0.5,0.75), title = "ICC - Bayesian GLM",
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
     save = TRUE, fname = "plots/601_subject_bayes_tongue_icc.png")

# >>>> Classical ----
avg_estimates_classical <- readRDS("HCP_results/5k_results/602_avg_estimates_PW_classical.rds")
library(abind)
abind_g <- function(...) abind(...,along = 4)
ICC_values_average_classical <- sapply(avg_estimates_classical, function(hem_est) {
  combined_visits <- Reduce(abind_g,hem_est)
  vertex_ICC <- apply(combined_visits, 1:2, ICC)
  return(vertex_ICC)
}, simplify = F)
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
plot_dir <- "~/github/BayesGLM_Validation/HCP_results/plots"
task_file_names <- c("visual_cue","foot","hand","tongue")
classical_icc_avg_cifti <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
classical_icc_avg_cifti$data$cortex_left <- ICC_values_average_classical$left
classical_icc_avg_cifti$data$cortex_right <- ICC_values_average_classical$right
library(viridisLite)
my_cols  <- viridis(4)
plot(classical_icc_avg_cifti, idx = 4, color_mode = "sequential",
     colors = my_cols, zlim = c(0,0.25,0.5,0.75), title ="ICC - Classical GLM",
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
     save = TRUE, fname = "plots/601_subject_classical_tongue_icc.png")

# HCP GROUP PLOTS ----
group_analysis_left2 <- readRDS("~/github/BayesGLM_Validation/HCP_results/5k_results/group/HCP_BayesGLM2_result_left_thresh0_20210105.rds")
group_estimates_classical <- readRDS("HCP_results/5k_results/group/502_HCP_classical_group_estimates.rds")
library(ciftiTools)
# ciftiTools.setOption('wb_path',"/Users/danspen/workbench")
ciftiTools.setOption('wb_path',"/Applications/workbench")

# >> Plot Estimates ----
# >>>> Bayes ----
group_results_left <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/HCP_BayesGLM2_45subj_result_left_thresh0_20210129.rds")
group_results_right <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/HCP_BayesGLM2_45subj_result_right_thresh0_20210129.rds")
cifti_both <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
cifti_both$data$cortex_left <- as.matrix(group_results_left$estimates[,4])
cifti_both$data$cortex_right <- as.matrix(group_results_right$estimates[,4])
plot(cifti_both, zlim = c(-1,1),
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
     save = TRUE, fname = "plots/601_group_bayes_tongue_estimate.png")
# zlims <- list(c(-1.2,1.2),
#               c(-1,1),
#               c(-1,1),
#               c(-1,1))

# for(i in seq(num_betas)) {
#   view_xifti_surface(
#     xifti = cifti_left,
#     idx = i,
#     hemisphere = plot_hem,
#     view = "both",
#     color_mode = "diverging",
#     zlim = zlims[[i]],
#     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
#     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
#     # title = task_names[i],
#     fname = paste0("plots/601_group_Bayes_", file_task_names[i]),
#     save = TRUE
#   )
# }
# Classical
# cifti_left_classical <- indiv_results_left$betas_Bayesian$avg
# cifti_left_classical$data$cortex_left <- group_estimates_classical$left
# for(i in seq(num_betas)) {
#   view_xifti_surface(
#     xifti = cifti_left_classical,
#     idx = i,
#     hemisphere = plot_hem,
#     view = "both",
#     color_mode = "diverging",
#     zlim = zlims[[i]],
#     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
#     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
#     # title = task_names[i],
#     fname = paste0("plots/601_group_classical_", file_task_names[i]),
#     save = TRUE
#   )
# }

# >>>>  Classical ----
class_est <- readRDS("HCP_results/5k_results/502_HCP_classical_group_PW_estimates.rds")
cifti_both <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
cifti_both$data$cortex_left <- class_est$left[,4,drop=F]
cifti_both$data$cortex_right <- class_est$right[,4,drop=F]
plot(cifti_both, zlim = c(-1,1),
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
     save = TRUE, fname = "plots/601_group_classical_tongue_estimate.png")

# >> Plot Activations ----
# >>>>  Bayes ----
group_results_left0 <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/HCP_BayesGLM2_45subj_result_left_thresh0_20210131.rds")
group_results_left05 <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/HCP_BayesGLM2_45subj_result_left_thresh05_20210131.rds")
group_results_left1 <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/HCP_BayesGLM2_45subj_result_left_thresh1_20210131.rds")
group_active_left <- cbind(
  group_results_left0$active[,4],
  group_results_left05$active[,4],
  group_results_left1$active[,4]
)
group_active_left <- as.matrix(apply(group_active_left,1,sum))
group_active_left[,1] <- ifelse(group_active_left[,1] == 0, NA, group_active_left[,1])


group_results_right0 <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/HCP_BayesGLM2_45subj_result_right_thresh0_20210131.rds")
group_results_right05 <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/HCP_BayesGLM2_45subj_result_right_thresh05_20210131.rds")
group_results_right1 <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/HCP_BayesGLM2_45subj_result_right_thresh1_20210131.rds")
group_active_right <- cbind(
  group_results_right0$active[,4],
  group_results_right05$active[,4],
  group_results_right1$active[,4]
)
group_active_right <- as.matrix(apply(group_active_right,1,sum))
group_active_right[,1] <- ifelse(group_active_right[,1] == 0, NA, group_active_right[,1])

library(viridisLite)
my_cols <- viridis(3)

cifti_both <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
cifti_both$data$cortex_left <- group_active_left
cifti_both$data$cortex_right <- group_active_right
plot(cifti_both, color_mode = "qualitative", colors = my_cols,
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
     save = TRUE, fname = "plots/601_group_bayes_tongue_activation.png")

# >>>> Classical ----
combined_active_FDR <- readRDS("HCP_results/5k_results/502_HCP_classical_activations_PW_FDR.rds")
combined_active_FWER <- readRDS("HCP_results/5k_results/502_HCP_classical_activations_PW_FWER.rds")
combined_left <- cbind(
  combined_active_FDR$left$`0%`[,4],
  combined_active_FWER$left$`0%`[,4]
)
combined_left <- apply(combined_left,1,sum)
combined_left <- as.matrix(ifelse(combined_left == 0,NA, combined_left))
combined_right <- cbind(
  combined_active_FDR$right$`0%`[,4],
  combined_active_FWER$right$`0%`[,4]
)
combined_right <- apply(combined_right,1,sum)
combined_right <- as.matrix(ifelse(combined_right == 0,NA, combined_right))
cifti_both <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
cifti_both$data$cortex_left <- combined_left
cifti_both$data$cortex_right <- combined_right
library(viridisLite)
my_cols <- viridis(2)
plot(cifti_both, color_mode = "qualitative", colors = my_cols,
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
     save = TRUE, fname = "plots/601_group_classical_tongue_activation.png")
# >> Boxplots Dice Coefficients ----
task_names <- c("visual cue","tongue","right foot","right hand","left foot", "left hand")
model_names <- c("Bayesian GLM", "Classical GLM - FDR", "Classical GLM - FWER")
file_names <- c("HCP_results/5k_results/604_Dice_coefficient_PW_threshold0.rds",
                "HCP_results/5k_results/604_classical_FDR_Dice_coefficient_PW_threshold0.rds",
                "HCP_results/5k_results/604_classical_FWER_Dice_coefficient_PW_threshold0.rds")
dice_dfs <- mapply(function(fn,mn) {
  dice_mat <- readRDS(fn)
  colnames(dice_mat) <- task_names
  d_df <- as.data.frame(dice_mat)
  d_df <- tidyr::gather(d_df,key = "Task", value = "Dice")
  d_df$Model <- mn
  return(d_df)
}, fn = file_names, mn = model_names, SIMPLIFY = F)

dice_df <- Reduce(rbind,dice_dfs)

library(ggplot2)
dice_plot <- ggplot(dice_df) +
  geom_boxplot(aes(x = Task, y = Dice, fill = Model), position = "dodge") +
  scale_fill_discrete("") +
  labs(y = "Dice Coefficient") +
  theme_bw() +
  theme(legend.position = "top")
print(dice_plot)
ggsave(filename = "plots/601_dice_plot.png",width = 1600/227, height = 958/227, plot = dice_plot)

# >> Barplots wICC ----
# >>>> Weighted ICC Calculation ----
# The weights here all come from the averaged classical estimates, as they are unbiased
avg_estimates_classical <- readRDS("HCP_results/5k_results/602_avg_estimates_PW_classical.rds")
library(abind)
ICC_weights <- sapply(avg_estimates_classical, function(hem_est) {
  combined_est <- Reduce(abind,hem_est)
  avg_est <- apply(combined_est,1:2,mean)
  weight_est <- apply(avg_est,2,function(avg_task) abs(avg_task) / sum(abs(avg_task)))
  return(weight_est)
}, simplify = F)
ICC_weights_final <- list(
  left = matrix(0,4443,6),
  right = matrix(0,4444,6)
)
ICC_weights_final$left[,1:4] <- ICC_weights$left[,c(1,4,2,3)]
ICC_weights_final$right[,c(1,2,5,6)] <- ICC_weights$right[,c(1,4,2,3)]
# >>>> Bayesian ----
Bayes_ICC <- list(
  left = matrix(0,4443,6),
  right = matrix(0,4444,6)
)
Bayes_ICC$left[,1:4] <- ICC_values_average$left[,c(1,4,2,3)]
Bayes_ICC$right[,c(1,2,5,6)] <- ICC_values_average$right[,c(1,4,2,3)]
bayes_wICC_avg <- round(apply(mapply(function(wt,icc) {
  return(apply(wt*icc,2,sum))
}, wt = ICC_weights_final, icc = Bayes_ICC),1,sum),3)
# >>>> Classical ----
Class_ICC <- list(
  left = matrix(0,4443,6),
  right = matrix(0,4444,6)
)
Class_ICC$left[,1:4] <- ICC_values_average_classical$left[,c(1,4,2,3)]
Class_ICC$right[,c(1,2,5,6)] <- ICC_values_average_classical$right[,c(1,4,2,3)]
classical_wICC_avg <- round(apply(mapply(function(wt,icc) {
  return(apply(wt*icc,2,sum))
}, wt = ICC_weights_final, icc = Class_ICC),1,sum),3)
# >>>> Combined and plotted ----
wicc_df <- data.frame(
  model = c(rep("Bayesian GLM",6),rep("Classical GLM",6)),
  task = rep(c("visual cue","tongue","right foot","right hand","left foot","left hand"),2),
  wicc = c(bayes_wICC_avg,classical_wICC_avg)
)
library(ggplot2)
wicc_plot <- ggplot(wicc_df) +
  geom_point(aes(x = task, y = wicc, color = model), pch = 19,size  = 5) +
  labs(x = "Task",y = "weighted Intraclass Correlation\nCoefficient (wICC)") +
  scale_color_discrete("") +
  theme_bw() +
  theme(legend.position = "top")
print(wicc_plot)
ggsave("plots/601_wicc_plot.png",plot = wicc_plot,width = 1600/227, height = 796/227)
