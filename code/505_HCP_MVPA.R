# This is a script to try to implement the multivoxel pattern analysis that
# was described in "An fMRI-Based Neurologic Signature of Physical Pain" by
# Wager et al [2013]

# Things we have tried with MVPA:
# 1. Including all subjects and using the 0% group activation threshold and using all tasks combined for the PCA (doesn't really work, as subject 10 overpowers PC1 for the tongue and throws off the predictions. The LASSO regression chooses no regressors. R^2 = 0)
# 2. Removing subject 10 and using the 0% group activation threshold, all tasks combined (R^2 = 0.7)
# 3. All subjects, 0.5% group activation threshold combining all tasks (R^2 = 0.56)
# 4. Remove subject 10, 0.5% group activation combining all tasks (R^2 = 0)
# 5. Remove subject 10, use the 0% group activation, don't scale in prcomp (R^2 = 0.946) (GOLD STANDARD)
# 6. All subjects, 0% group activation masks, removing the temporal pole for the tongue task, don't scale prcomp (R^2 = 0.997)
# 7. All subjects, only the hand tasks, 0% group activation mask, don't scale prcomp. This one fails at the LASSO step, as zero predictors are found to be worth the penalty.

# 1) Make masks ----
# Using the group 0% active thresholds
group_active_left <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/501_HCP_BayesGLM2_45subj_result_left_thresh0_20210317.rds")$active
group_active_right <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/501_HCP_BayesGLM2_45subj_result_right_thresh0_20210318.rds")$active
# Using the Yeo 17 parcellation
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
yeo_17 <- read_cifti("/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/yeo_17_parcellation.dlabel.nii", resamp_res = 5000)
somatomotor_idx <- grep("Somatomotor", rownames(yeo_17$meta$cifti$labels$`#1`))
motor_mask <- list(
  left = as.matrix(yeo_17$data$cortex_left[,1] %in% somatomotor_idx),
  right = as.matrix(yeo_17$data$cortex_right[,1] %in% somatomotor_idx)
)

# Take out the regions around the temporal pole for the tongue task
cifti_bad <- readRDS("HCP_results/505_cifti_bad.rds")
group_active_left[cifti_bad$data$cortex_left == 1,4] <- 0
group_active_right[cifti_bad$data$cortex_right == 1,4] <- 0


# Plot the mask so we see which parts are included in the analysis
# cifti_group <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
# cifti_group$data$cortex_left <- group_active_left
# cifti_group$data$cortex_right <- group_active_right
# library(ciftiTools)
# ciftiTools.setOption('wb_path',"/Applications/workbench/")
# plot(cifti_group, idx = 3, color_mode = "qualitative",
#        surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
#        surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")

# This is for using the overlapping areas of activation for the mask
# all_session_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/all_sessions"
# all_session_files <- list.files(all_session_dir, full.names = T)
# all_session_files <- grep("_thr0.rds", all_session_files, value = T)
# load("HCP_data/subjects.Rdata")
# num_subjects <- length(subjects)
# task_names <- c("cue","foot","hand","tongue")
# task_idx <- 4
# percentile <- 75
# # for(percentile in 100*seq(0.5,1,by = 0.05)) {
#   activation_overlap <- sapply(c("left","right"), function(hem) {
#     all_hem_active <- sapply(subjects, function(subject) {
#       subject_file <- grep(subject, all_session_files, value = T)
#       subject_file <- grep(hem, subject_file, value = T)
#       obj <- readRDS(subject_file)
#       return(obj$active)
#     }, simplify = F)
#     overlap_out <- Reduce(`+`,all_hem_active)
#     overlap_out[overlap_out < percentile / 100 *num_subjects] <- 0
#     overlap_out[overlap_out > 0] <- 1
#     return(overlap_out)
#   },simplify = F)
#   # These need to be plotted for examination
#   library(ciftiTools)
#   ciftiTools.setOption('wb_path','/Applications/workbench')
#   sapply(activation_overlap, summary, simplify = F)
#   cifti_obj <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
#   cifti_obj$data$cortex_left <- activation_overlap$left
#   cifti_obj$data$cortex_right <- activation_overlap$right
#   plot(cifti_obj, idx = 4, color_mode = "qualitative",
#        title = paste0(task_names[task_idx]," - ",percentile,"th percentile"),
#   surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
#   surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
#   fname = paste0("~/Desktop/505_overlap_0.5activation_",task_names[task_idx],"_",percentile,"percentile.png"))
# }

# mask_left[mask_left == 0] <- NA # This is for making the ciftis easy to view
# mask_right[mask_right == 0] <- NA # This is for making the ciftis easy to view
# cifti_mask <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
# cifti_mask$data$cortex_left <- mask_left * 1
# cifti_mask$data$cortex_right <- mask_right * 1
# mask_left <- apply(activation_overlap$left,2,function(x) which(x == 1))
# mask_right <- apply(activation_overlap$right,2,function(x) which(x == 1))

mask_left <- apply(group_active_left,2,function(x) x == 1)
mask_right <- apply(group_active_right,2,function(x) x == 1)
masks <- list(left = mask_left, right = mask_right)

# 2) Extract individual amplitude estimates within the masks ----
subject_result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/all_sessions/"
subject_files <- list.files(subject_result_dir, full.names = T)
subject_files <- grep("thr0.rds", subject_files, value = T)
load("HCP_data/subjects.Rdata")
task_names <- c("cue","foot","hand","tongue")
all_mvpa_amplitudes <- sapply(subjects, function(subject) {
  all_hems <- sapply(c("left","right"), function(hem) {
    result_obj <- readRDS(grep(paste0(subject,"_",hem), subject_files, value = T))
    # all tasks
    all_tasks <- mapply(function(est, mas) {
      return(as.matrix(est[mas]))
    }, est = split(result_obj$estimates,col(result_obj$estimates)), mas = split(masks[[hem]], col(masks[[hem]])), SIMPLIFY = F)
    names(all_tasks) <- task_names
    return(all_tasks)
    # Just the hand tasks
    # hand_task <- result_obj$estimates[,3] * masks[[hem]]
    # return(hand_task)
  }, simplify = F)
}, simplify = F)

library(tidyverse)
demog_df <- read.csv("HCP_data/unrestricted_danspen_2_24_2021_12_55_28.csv") %>%
  mutate(subject = factor(Subject)) %>%
  filter(subject %in% subjects) %>%
  mutate(subject = droplevels(subject)) %>%
  select(subject,"Dexterity_Unadj") %>%
  rename(Dex = Dexterity_Unadj)

mvpa_df <-
  reshape2::melt(all_mvpa_amplitudes, value.name = "estimate") %>%
  select(-Var2) %>%
  # all tasks
  rename(task = L3,
         hem = L2,
         subject = L1) %>%
  pivot_wider(names_from = c(hem,task, Var1), values_from = estimate) %>%
  # hand tasks
  # rename(hem = L2,
  #        subject = L1) %>%
  # pivot_wider(names_from = c(hem, Var1), values_from = estimate) %>%
  left_join(demog_df) %>%
  select(-subject) %>%
  mutate(Dex = scale(Dex, scale = T))

# 3) Fit a final model using all 45 subjects to determine the shrinkage penalty ----

# Remove subject 10
# mvpa_df <- mvpa_df[-10,]

mvpa_pca <- prcomp(x = mvpa_df[,-ncol(mvpa_df)], scale. = F)
# heatmap(mvpa_pca$x)
library(glmnet)
n.times <- 100
r2 <- rep(0,n.times)
lam <- rep(0,n.times)
num_parms <- rep(0,n.times)
for(i in 1:n.times) {
  find_lambda <- cv.glmnet(mvpa_pca$x, mvpa_df$Dex,alpha = 1)
  # plot(find_lambda)
  lasso_lm <- glmnet(mvpa_pca$x,mvpa_df$Dex, alpha = 1, lambda = find_lambda$lambda.min)
  lasso_lm$beta
  y_hat <- predict(lasso_lm, newx = mvpa_pca$x)
  # plot(mvpa_df$Dex, y_hat)
  r2[i] <- cor(mvpa_df$Dex, y_hat)^2
  lam[i] <- find_lambda$lambda.min
  num_parms[i] <- sum(lasso_lm$beta != 0)
}
summary(r2)
summary(lam)
summary(num_parms)
table(num_parms)
plot(lam,r2)
plot(num_parms, r2)

final_model <- glmnet(mvpa_pca$x, mvpa_df$Dex, alpha = 1, lambda = median(lam))
final_model$beta
in_sample_PCS <- which(final_model$beta != 0)
y_hat <- predict(final_model, newx = mvpa_pca$x)
plot(mvpa_df$Dex, y_hat)
abline(0,1,col = 'red')
cor(mvpa_df$Dex, y_hat)^2
sqrt(mean((mvpa_df$Dex - y_hat)^2)) # Compared to SD = 1 of the response data

# 4) 45 LOO-CV models with OOS prediction using LASSO regression ----
# 1. First try: Using the same PCs and the same lambda as the in-subject models. OOS R^2 = 0.066
# 2. Calculate different PCs for each subject removal, use the same lambda as the in-subject model. OOS R^2 = 0.0262
# 3. Same PCs as in-subject model, new calculation of lambda using the LOO data. OOS R^2 = 0.003
# 4. Calculate new PCs for each subject removal, calculate a new lambda using the LOO data. OOS R^2 = 0.003
# 5. Use the same PCs as the in-sample model, and use the same PCs in an lm used in the final in-sample LASSO. OOS R^2 = 0.44, but the correlation is negative (r = -0.66)

subject_mvpa <- sapply(subjects, function(subject) {
  subject_idx <- which(subjects == subject)
  # 1.
  # loo_lasso <- glmnet(mvpa_pca$x[-subject_idx,], mvpa_df$Dex[-subject_idx], alpha = 1, lambda = median(lam))
  # yhat <- predict(loo_lasso, newx = mvpa_pca$x[subject_idx,,drop=F])
  # 2.
  # loo_pca <- prcomp(mvpa_df[-subject_idx,-ncol(mvpa_df)],scale. = F)
  # oos_pc <- as.matrix(mvpa_df[subject_idx,-ncol(mvpa_df)]) %*% loo_pca$rotation
  # loo_lasso <- glmnet(loo_pca$x, mvpa_df$Dex[-subject_idx], alpha = 1, lambda = median(lam))
  # yhat <- predict(loo_lasso, newx = oos_pc)
  # 3.
  # n.times <- 50
  # loo_lambda <- numeric(n.times)
  # for(i in 1:n.times) {
  #   loo_lambda[i] <- cv.glmnet(mvpa_pca$x[-subject_idx,], mvpa_df$Dex[-subject_idx], alpha = 1)$lambda.min
  # }
  # loo_lasso <- glmnet(mvpa_pca$x[-subject_idx,], mvpa_df$Dex[-subject_idx], alpha = 1, lambda = median(loo_lambda))
  # yhat <- predict(loo_lasso, newx = mvpa_pca$x[subject_idx,,drop=F])
  # 4.
  # loo_pca <- prcomp(mvpa_df[-subject_idx,-ncol(mvpa_df)],scale. = F)
  # oos_pc <- as.matrix(mvpa_df[subject_idx,-ncol(mvpa_df)]) %*% loo_pca$rotation
  # n.times <- 50
  # loo_lambda <- numeric(n.times)
  # for(i in 1:n.times) {
  #   loo_lambda[i] <- cv.glmnet(loo_pca$x, mvpa_df$Dex[-subject_idx], alpha = 1)$lambda.min
  # }
  # loo_lasso <- glmnet(loo_pca$x, mvpa_df$Dex[-subject_idx], alpha = 1, lambda = median(loo_lambda))
  # yhat <- predict(loo_lasso, newx = oos_pc)
  # 5.
  pc_df <- as.data.frame(mvpa_pca$x[-subject_idx,in_sample_PCS])
  pc_df$Dex <- mvpa_df$Dex[-subject_idx]
  oos_df <- as.data.frame(mvpa_pca$x[subject_idx,in_sample_PCS, drop = F])
  oos_df$Dex <- pc_df$Dex <- mvpa_df$Dex[subject_idx]
  loo_lm <- lm(Dex ~ -1 + ., data = pc_df)
  yhat <- predict(loo_lm, newdata = oos_df)
  return(yhat)
}, simplify = T)

plot(mvpa_df$Dex, subject_mvpa)
cor(mvpa_df$Dex, subject_mvpa)

# Plots of PCs used ----
pcs_used <- which(final_model$beta != 0)

pc_rotation <- mvpa_pca$rotation[,pcs_used]

for(pc_idx in pcs_used) {
  pcidx_idx <- which(pcs_used == pc_idx)
  for(task_idx in 1:4) {
    cifti_pc <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
    for(hem in c('left','right')) {
      pc_rows <- grep(paste0(hem,"_",task_names[task_idx]),rownames(pc_rotation))
      pc_loading <- masks[[hem]][,task_idx]
      pc_loading[pc_loading] <- pc_rotation[pc_rows,pcidx_idx] * sign(final_model$beta[pc_idx])
      pc_loading[pc_loading==0] <- NA
      cifti_pc$data[[paste0("cortex_",hem)]] <- as.matrix(pc_loading)
    }
    plot(cifti_pc, title = paste0("PC",pc_idx,": ", task_names[task_idx]," task"),
         zlim = c(-max(abs(pc_rotation[,pcidx_idx])),max(abs(pc_rotation[,pcidx_idx]))),
         fname = paste0("plots/505_mvpa_pc",pc_idx,"_",task_names[task_idx],".png"),
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")
  }
}

mdl_contrib <- pc_rotation %*% as.matrix(final_model$beta[pcs_used])
pc_contrib_cifti <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
for(task_idx in 1:4) {
  for(hem in c('left','right')) {
    pc_rows <- grep(paste0(hem,"_",task_names[task_idx]),rownames(pc_rotation))
    pc_loading <- masks[[hem]][,task_idx]
    pc_loading[pc_loading] <- mdl_contrib[pc_rows,1]
    pc_loading[pc_loading==0] <- NA
    pc_contrib_cifti$data[[paste0("cortex_",hem)]][,task_idx] <- pc_loading
  }
}

max_abs_contrib <- max(abs(mdl_contrib))
task_idx <- 4
for(task_idx in 1:4) {
  plot(pc_contrib_cifti, idx = task_idx,
       title = paste0("Model Contribution - ", task_names[task_idx]),
       zlim = c(-0.6*max_abs_contrib, 0.6*max_abs_contrib),
       fname = paste0("plots/505_mvpa_mdl_contrib_",task_names[task_idx],".png"),
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")
}

# Outlier work ----

#  >> Separate tasks PCA without subject 10 ----
# no_sub10 <- mvpa_df[-10,]
#
# pesel(no_sub10)
#
# cue_pesel <- pesel(select(no_sub10,starts_with("left_cue"), starts_with("right_cue")), npc.max = 45)
# plot(cue_pesel$vals)
# cue_mvpa <- prcomp(x = select(no_sub10,starts_with("left_cue"), starts_with("right_cue")), rank. = 40)
# cumsum(cue_mvpa$sdev^2 / sum(cue_mvpa$sdev^2))
# plot(cue_mvpa)
# heatmap(cue_mvpa$x)
# foot_mvpa <- prcomp(x = select(no_sub10,starts_with("left_foot"), starts_with("right_foot")), rank. = 40)
# heatmap(foot_mvpa$x)
# hand_mvpa <- prcomp(x = select(no_sub10,starts_with("left_hand"), starts_with("right_hand")), rank. = 40)
# heatmap(hand_mvpa$x)
# tongue_mvpa <- prcomp(x = select(no_sub10,starts_with("left_tongue"), starts_with("right_tongue")), rank. = 40)
# heatmap(tongue_mvpa$x)
#
# task_separate_pca <- cbind(cue_mvpa$x,foot_mvpa$x, hand_mvpa$x, tongue_mvpa$x)
# colnames(task_separate_pca) <- c(paste0("Cue_PC",1:40), paste0("Foot_PC",1:40),paste0("Hand_PC",1:40),paste0("Tongue_PC",1:40))
# # task_separate_pca <- task_separate_pca[,-136]
# heatmap(task_separate_pca)
# subject_10 <- subjects[10]
# left_result <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/all_sessions/501_BayesGLM2_subject_137128_left_thr0.rds")
# right_result <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/all_sessions/501_BayesGLM2_subject_137128_right_thr0.rds")
# library(ciftiTools)
# ciftiTools.setOption('wb_path',"/Applications/workbench/")
# cifti_obj <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
# cifti_obj$data$cortex_left <- left_result$estimates
# cifti_obj$data$cortex_right <- right_result$estimates
# cifti_obj <- cifti_obj * cifti_mask
# plot(cifti_obj, idx = 4, zlim = c(-1,1))

# cifti_pca_tongue <- cifti_obj
# sum(mask_left[,4], na.rm = T)
# # 1023
# sum(mask_right[,4], na.rm = T)
# # 1055
# cifti_pca_tongue$data$cortex_left[which(!is.na(cifti_pca_tongue$data$cortex_left[,4])),] <-
#   tongue_mvpa$rotation[seq(1023),1]
# cifti_pca_tongue$data$cortex_right[which(!is.na(cifti_pca_tongue$data$cortex_right[,4])),] <-
#   tongue_mvpa$rotation[-seq(1023),1]
# plot(cifti_pca_tongue,idx = 4)

# sub1_result_left <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/all_sessions/501_BayesGLM2_subject_103818_left_thr0.rds")
# sub1_result_right <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/all_sessions/501_BayesGLM2_subject_103818_right_thr0.rds")
# sub1_cifti <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
# sub1_cifti$data$cortex_left <- sub1_result_left$estimates
# sub1_cifti$data$cortex_right <- sub1_result_right$estimates
# sub1_cifti <- sub1_cifti * cifti_mask
# plot(sub1_cifti, idx = 4, zlim = c(-1,1))

# mvpa_df %>%
#   select(-Dex, starts_with("left_tongue"), starts_with("right_tongue")) %>%
#   pivot_longer(cols = -"subject", names_to = "location") %>%
#   mutate(sub10 = subject == subject_10) %>%
#   ggplot() +
#   geom_histogram(aes(x = value), alpha = 0.3, position = "dodge", bins = 100) +
#   facet_grid(sub10~., scales = "free") +
#   lims(x = c(-2,3))

# >> Session estimates for subject 10 ----
# visit1_left_result <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit1_left_5k_20210127.rds")
# visit1_right_result <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit1_right_5k_20210127.rds")
# visit2_left_result <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit2_left_5k_20210127.rds")
# visit2_right_result <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit2_right_5k_20210127.rds")

# subj10_results <- list(
#   visit1 = list(
#     left = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit1_left_5k_20210127.rds")$betas_Bayesian,
#     right = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit1_right_5k_20210127.rds")$betas_Bayesian
#   ),
#   visit2 = list(
#     left = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit2_left_5k_20210127.rds")$betas_Bayesian,
#     right = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit2_right_5k_20210127.rds")$betas_Bayesian
#   )
# )

# sub10_ciftis <- sapply(subj10_results, function(visit_results) {
#   session_ciftis <- sapply(c("LR","RL","avg"), function(session_name) {
#     cifti_out <- xifti_combine(visit_results$left[[session_name]], visit_results$right[[session_name]])
#     return(cifti_out)
#   }, simplify = F)
# }, simplify = F)

# library(ciftiTools)
# ciftiTools.setOption('wb_path','/Applications/workbench/')

# for(v in paste0("visit",1:2)) {
#   for(sess in c("LR","RL","avg")) {
#     plot(sub10_ciftis[[v]][[sess]], title = paste("Subject 10",v,sess),
#          surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
#          surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
#          idx = 4)
#   }
# }

# library(tidyverse)
# reshape2::melt(all_mvpa_amplitudes) %>%
#   filter(L3 == "tongue") %>%
#   mutate(is_sub10 = L1 == "137128") %>%
#   ggplot() +
#   geom_histogram(aes(x = value), bins = 100) +
#   facet_grid(is_sub10~., scales = "free")


# reshape2::melt(all_mvpa_amplitudes) %>%
#   filter(L3 == "tongue") %>%
#   mutate(is_sub10 = L1 == "137128") %>%
#   group_by(is_sub10) %>%
#   summarize(Min = min(value),
#             Mean = mean(value),
#             Median = median(value),
#             Max = max(value),
#             SD = sd(value))

# >> Look at data ----
# subj10_classical <- list(
#   visit1 = list(
#     left = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit1_left_5k_20210127.rds")$betas_classical,
#     right = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit1_right_5k_20210127.rds")$betas_classical
#   ),
#   visit2 = list(
#     left = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit2_left_5k_20210127.rds")$betas_classical,
#     right = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit2_right_5k_20210127.rds")$betas_classical
#   )
# )

# >> Examine the amplitude estimates ----
# sub10_est <- sapply(subj10_results, function(visit_results) {
#   session_ciftis <- sapply(c("LR","RL","avg"), function(session_name) {
#     est_out <- c(visit_results$left[[session_name]]$data$cortex_left[,4], visit_results$right[[session_name]]$data$cortex_right[,4])
#     return(est_out)
#   }, simplify = F)
# }, simplify = F)

# sapply(sub10_est, function(x) sapply(x, function(y) summary(y)), simplify = F)

# sub10_est_class <- sapply(subj10_classical, function(visit_results) {
#   session_ciftis <- sapply(c("LR","RL","avg"), function(session_name) {
#     est_out <- c(visit_results$left[[session_name]]$data$cortex_left[,4], visit_results$right[[session_name]]$data$cortex_right[,4])
#     return(est_out)
#   }, simplify = F)
# }, simplify = F)

# sapply(sub10_est_class, function(x) sapply(x, function(y) summary(y)), simplify = F)

# subj1_results <- list(
#   visit1 = list(
#     left = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_103818_visit1_left_5k_20210125.rds")$betas_Bayesian,
#     right = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_103818_visit1_right_5k_20210125.rds")$betas_Bayesian
#   ),
#   visit2 = list(
#     left = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_103818_visit2_left_5k_20210125.rds")$betas_Bayesian,
#     right = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_103818_visit2_right_5k_20210125.rds")$betas_Bayesian
#   )
# )
#
# sub1_est <- sapply(subj1_results, function(visit_results) {
#   session_ciftis <- sapply(c("LR","RL","avg"), function(session_name) {
#     est_out <- c(visit_results$left[[session_name]]$data$cortex_left[,4], visit_results$right[[session_name]]$data$cortex_right[,4])
#     return(est_out)
#   }, simplify = F)
# }, simplify = F)
#
# sapply(sub1_est, function(x) sapply(x, function(y) summary(y)), simplify = F)
#
#
#
# # >> Examine the single-session results ----
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session"
sub10_files <- grep("137128", list.files(result_dir, full.names = T), value = T)
ss_est <- sapply(sub10_files, function(x) {
  obj <- readRDS(x)
  out <- do.call(rbind, obj$betas_Bayesian[[1]]$data)
  return(out)
}, simplify = F)
summary(unlist(ss_est))
sapply(ss_est, summary)
paste(paste0("visit",1:2), c('left',"right"), c("LR,RL"))
# sapply(ss_est, function(x) summary(x[,4]))
#
# # Compare to subject 1
# sub1_files <- grep("103818", list.files(result_dir, full.names = T), value = T)
# ss1_est <- sapply(sub1_files, function(x) {
#   obj <- readRDS(x)
#   out <- do.call(rbind, obj$betas_Bayesian[[1]]$data)
#   return(out)
# }, simplify = F)
# summary(unlist(ss1_est))
# sapply(ss1_est, summary)
# sapply(ss1_est, function(x) summary(x[,3]))
#
# # >> Examine the design matrices ----
# sub10_designs <- sapply(sub10_files, function(x) {
#   readRDS(x)$design[[1]]
# }, simplify = F)
# names(sub10_designs) <- paste0("session",1:8)
#
# sub1_designs <- sapply(sub1_files, function(x) {
#   readRDS(x)$design[[1]][,4]
# }, simplify = F)
# names(sub1_designs) <- paste0("session",1:8)
#
# sapply(sub10_designs, summary)
# sapply(sub1_designs, summary)
#
# sub_files <- grep(".rds",list.files(result_dir, full.names = T), value = T)
# sub_files <- grep("_left_", sub_files, value = T)
# sapply(subjects, function(x) {
#   sfiles <- grep(x,sub_files, value = T)
#   length(sfiles)
# })
# grep("151526", sub_files, value = T)

# >> Calculate the VIF ----
# load("/Users/Shared/Lab_Data/HCP_Motor_Task_Dan/subjects.Rdata")
# main_dir <- "/Users/Shared/Lab_Data/HCP_Motor_Task_Dan"
# data_dir <- "/Users/Shared/Lab_Data/HCP_Motor_Task_Dan/visit1_data"
# ss_result_dir <- "HCP_results/5k_results/individual/PW/single_session"
# ss_result_files <- list.files(ss_result_dir, full.names = T)
# ss_result_files <- grep("500_", ss_result_files, value = T)
# library(clever)
# library(BayesfMRI)
# task_names <- c("cue","foot","hand","tongue")
# VIF_list <- list()
# for(subject in subjects) {
#   VIF_list[[subject]] <- list()
#   for(visit in 1:2) {
#     VIF_list[[subject]][[paste0("visit",visit)]] <- list()
#     #analyze hemispheres separately due to different in set of tasks
#     for(h in 1:2){
#       #h=1 -- left hemisphere
#       #h=2 -- right hemisphere
#       hem <- c('left','right')[h]
#       VIF_list[[subject]][[paste0("visit",visit)]][[hem]] <- list()
#
#       # find the VIF at the session level
#       for (s in 1:2) {
#         sess <- c("LR","RL")[s]
#         motion <-
#           as.matrix(
#             read.table(
#               file.path(
#                 main_dir,
#                 paste0("visit", visit, "_data"),
#                 subject,
#                 "MNINonLinear",
#                 "Results",
#                 paste0("tfMRI_MOTOR_",sess),
#                 "Movement_Regressors.txt"
#               )
#             , header = FALSE)
#           )
#         design <- readRDS(
#           grep(
#             paste0(subject,"_visit",visit,"_",hem,"_5k_session",sess),
#             ss_result_files,
#             value = T
#           )
#         )$design[[1]]
#         drift1 <- (1:284)/284
#         drift <- cbind(drift1, drift1^2)
#         # design2 <- BayesfMRI:::nuisance_regress(design, cbind(motion, drift))
#         # myFD <- FD(X = motion[,1:6])
#         VIF_session <- matrix(NA,4,1)
#         for(k in 1:4){
#           VIF_session[k,1] <- 1 / (1 - summary(lm(design[,k] ~ as.matrix(cbind(motion, drift))))$r.squared)
#         }
#         rownames(VIF_session) <- task_names
#         VIF_list[[subject]][[paste0("visit",visit)]][[hem]][[sess]] <- VIF_session
#       }
#     }
#   }
# }
# saveRDS(VIF_list, "~/Desktop/505_VIF_list.rds")

# VIF_list <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session/505_VIF_list.rds")
#
# library(tidyverse)
# task_names <- c("cue","foot","hand","tongue")
# # Individual VIFs
# reshape2::melt(VIF_list) %>%
#   rename(Task = Var1) %>%
#   ggplot() +
#   geom_jitter(aes(x = L1, y = value, color = Task), width = 0.25, height = 0) +
#   geom_vline(xintercept = 10) +
#   labs(x = "Subject",y = "VIF") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
#
# # Average VIFs
# reshape2::melt(VIF_list) %>%
#   rename(Task = Var1) %>%
#   group_by(L1,L2,L3,Task) %>%
#   summarize(value = min(value)) %>%
#   ggplot() +
#   geom_point(aes(x = L1, y = value, color = Task)) +
#   geom_vline(xintercept = 10) +
#   labs(x = "Subject", y = "Average VIF") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
#
# # >> Examine subject_10 results without nuisance regression or prewhitening ----
# result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/subject_10"
# result_files <- list.files(result_dir,full.names = T)
# result_files <- grep(".rds", result_files, value = T)
# sub1_ests <- sapply(paste0("visit",1:2), function(visit) {
#   sapply(c("left","right"), function(hem) {
#     sapply(c("LR","RL"), function(sess) {
#       filen <- grep(paste0("103818_",visit,"_",hem), result_files, value = T)
#       filen <- grep(sess, filen, value = T)
#       obj <- readRDS(filen)
#       ests <- obj$betas_classical[[sess]]$data[[paste0("cortex_",hem)]]
#       return(ests)
#     }, simplify = F)
#   }, simplify = F)
# }, simplify = F)
#
#
# library(tidyverse)
# sub10_ests <- sapply(paste0("visit",1:2), function(visit) {
#   sapply(c("left","right"), function(hem) {
#     sapply(c("LR","RL"), function(sess) {
#       filen <- grep(paste0("137128_",visit,"_",hem), result_files, value = T)
#       filen <- grep(sess, filen, value = T)
#       filen <- grep("nobold", filen, value = T)
#       obj <- readRDS(filen)
#       ests <- obj$betas_classical[[sess]]$data[[paste0("cortex_",hem)]]
#       return(ests)
#     }, simplify = F)
#   }, simplify = F)
# }, simplify = F)
#
# reshape2::melt(sub10_ests) %>%
#   filter(value > 10)
#
# task_names <- c("cue","foot","hand","tongue")
# reshape2::melt(sub10_ests) %>%
#   mutate(subject = "137128") %>%
#   full_join(
#     reshape2::melt(sub1_ests) %>%
#       mutate(subject = "103818")
#   ) %>%
#   mutate(task = task_names[Var2]) %>%
#   ggplot() +
#   geom_boxplot(aes(x = task, y = value, color = subject), position = "dodge") +
#   facet_grid(L1 + L2 ~ L3) +
#   labs(y = "Amplitude Estimates", title = "Subject 1 vs. Subject 10")
#
# sub1_visit2_left_LR <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/subject_10/500_103818_visit2_left_5k_classical_nonuisance_LR_20210415.rds")
# sub10_visit2_left_LR <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/subject_10/500_137128_visit2_left_5k_classical_nonuisance_LR_20210415.rds")
# sub10_visit2_left_RL <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/subject_10/500_137128_visit2_left_5k_classical_nonuisance_RL_20210415.rds")
# sub10_visit2_right_LR <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/subject_10/500_137128_visit2_right_5k_classical_nonuisance_LR_20210415.rds")
# sub10_visit2_right_RL <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/subject_10/500_137128_visit2_right_5k_classical_nonuisance_RL_20210415.rds")
#
# max(sub1_visit2_left_LR$design[[1]][,4])
# max(sub10_visit2_left_LR$design[[1]][,4])
#
# library(ciftiTools)
# ciftiTools.setOption('wb_path','/Applications/workbench/')
# plot(sub10_visit2_left_LR$betas_classical$LR, idx = 4, zlim = c(-20,20),
#      surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii")
#
# plot(sub10_visit2_left_RL$betas_classical$RL, idx = 4, zlim = c(-10,10),
#      surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii")
#
# plot(sub10_visit2_right_LR$betas_classical$LR, idx = 4, zlim = c(-10,10),
#      surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")
#
# plot(sub10_visit2_right_RL$betas_classical$RL, idx = 4, zlim = c(-20,20),
#      surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")
#
# # Try to find the tmeporal pole locations
# quantile(sub10_visit2_right_RL$betas_classical$RL$data$cortex_right[,4],probs = c(.95,.99))
# bad_right <- which(sub10_visit2_right_RL$betas_classical$RL$data$cortex_right[,4] > 7)
# quantile(sub10_visit2_left_LR$betas_classical$LR$data$cortex_left[,4],probs = c(.95,.99))
# bad_left <- which(sub10_visit2_left_LR$betas_classical$LR$data$cortex_left[,4] > 7)
#
# cifti_bad <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
# cifti_bad$data$cortex_left <- as.matrix(rep(0,4443))
# cifti_bad$data$cortex_left[bad_left,] <- 1
# cifti_bad$data$cortex_right <- as.matrix(rep(0,4444))
# cifti_bad$data$cortex_right[bad_right,] <- 1
# library(ciftiTools)
# ciftiTools.setOption('wb_path','/Applications/workbench')
# plot(cifti_bad, color_mode = "qualitative",
#      surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
#      surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")
# saveRDS(cifti_bad, "HCP_results/505_cifti_bad.rds")
#
# # Mean value map
# one_subj10_result <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session/500_137128_visit1_right_5k_sessionRL_20210406.rds")
# str(one_subj10_result$GLMs_Bayesian$cortexR$y)
# yy <- matrix(one_subj10_result$GLMs_Bayesian$cortexR$y, nrow = 284)
# plot(yy[,1], type='l')
# yy_avg <- apply(yy,2,mean)
# yy_avg_cifti <- sub10_visit2_right_RL$betas_classical$RL
# yy_avg_cifti$data$cortex_right <- as.matrix(yy_avg)
# plot(yy_avg_cifti,surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")
