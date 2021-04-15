# This is a script to try to implement the multivoxel pattern analysis that
# was described in "An fMRI-Based Neurologic Signature of Physical Pain" by
# Wager et al [2013]

# 1) Make masks based on 0.5% activation using all 45 subjects in the group analysis ----
group_active_left <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/501_HCP_BayesGLM2_45subj_result_left_thresh0_20210317.rds")$active
group_active_right <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/501_HCP_BayesGLM2_45subj_result_right_thresh0_20210318.rds")$active

# mask_left <- apply(group_active_left,2,function(x) x == 1)
# mask_left[mask_left == 0] <- NA
# mask_right <- apply(group_active_right,2,function(x) x == 1)
# mask_right[mask_right == 0] <- NA
# cifti_mask <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
# cifti_mask$data$cortex_left <- mask_left * 1
# cifti_mask$data$cortex_right <- mask_right * 1
mask_left <- apply(group_active_left,2,function(x) which(x == 1))
mask_right <- apply(group_active_right,2,function(x) which(x == 1))
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
    all_tasks <- mapply(function(est, mas) {
      return(as.matrix(est[mas]))
    }, est = split(result_obj$estimates,col(result_obj$estimates)), mas = masks[[hem]], SIMPLIFY = F)
    names(all_tasks) <- task_names
    return(all_tasks)
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
  rename(task = L3,
         hem = L2,
         subject = L1) %>%
  pivot_wider(names_from = c(hem,task, Var1), values_from = estimate) %>%
  left_join(demog_df) %>%
  select(-subject) %>%
  mutate(Dex = scale(Dex, scale = T))

heatmap(as.matrix(mvpa_df[,-ncol(mvpa_df)]))

# 3) Fit all 45 LOO-CV models using LASSO regression ----

no_sub10 <- mvpa_df[-10,]

pesel(no_sub10)

cue_pesel <- pesel(select(no_sub10,starts_with("left_cue"), starts_with("right_cue")), npc.max = 45)
plot(cue_pesel$vals)
cue_mvpa <- prcomp(x = select(no_sub10,starts_with("left_cue"), starts_with("right_cue")), rank. = 40)
cumsum(cue_mvpa$sdev^2 / sum(cue_mvpa$sdev^2))
plot(cue_mvpa)
heatmap(cue_mvpa$x)
foot_mvpa <- prcomp(x = select(no_sub10,starts_with("left_foot"), starts_with("right_foot")), rank. = 40)
heatmap(foot_mvpa$x)
hand_mvpa <- prcomp(x = select(no_sub10,starts_with("left_hand"), starts_with("right_hand")), rank. = 40)
heatmap(hand_mvpa$x)
tongue_mvpa <- prcomp(x = select(no_sub10,starts_with("left_tongue"), starts_with("right_tongue")), rank. = 40)
heatmap(tongue_mvpa$x)

task_separate_pca <- cbind(cue_mvpa$x,foot_mvpa$x, hand_mvpa$x, tongue_mvpa$x)
colnames(task_separate_pca) <- c(paste0("Cue_PC",1:40), paste0("Foot_PC",1:40),paste0("Hand_PC",1:40),paste0("Tongue_PC",1:40))
# task_separate_pca <- task_separate_pca[,-136]
heatmap(task_separate_pca)

mvpa_pca <- prcomp(x = no_sub10[,-ncol(no_sub10)], scale. = T, rank. = 43)
# heatmap(mvpa_pca$x)
library(glmnet)
# find_lambda <- cv.glmnet(task_separate_pca, no_sub10$Dex, alpha = 1)
find_lambda <- cv.glmnet(mvpa_pca$x, no_sub10$Dex,alpha = 1)
# find_lambda <- cv.glmnet(separate_pca2$x, no_sub10$Dex, alpha = 1)
plot(find_lambda)
lasso_lm <- glmnet(mvpa_pca$x,no_sub10$Dex, alpha = 1, lambda = find_lambda$lambda.min)
# separate_pca2 <- prcomp(x = task_separate_pca, rank. = 40)
# cumsum(separate_pca2$sdev^2 / sum(separate_pca2$sdev^2))
# heatmap(separate_pca2$x)
# find_lambda <- cv.glmnet(separate_pca2$x, no_sub10$Dex, alpha = 1)
# plot(find_lambda, ylim = c(0,1.5))
# lasso_lm <- glmnet(task_separate_pca,no_sub10$Dex, alpha = 1, lambda = find_lambda$lambda.min)
# lasso_lm <- glmnet(separate_pca2$x,no_sub10$Dex, alpha = 1, lambda = find_lambda$lambda.min)
lasso_lm$beta
y_hat <- predict(lasso_lm, newx = separate_pca2$x)
plot(mvpa_df$Dex[-10], y_hat)
cor(mvpa_df$Dex[-10], y_hat)
mvpa_pca_df <- as.data.frame(mvpa_pca$x)
mvpa_pca_df$Dex <- mvpa_df$Dex
cor(select(mvpa_pca_df,c("Dex",paste0("PC",1:6))))
pairs(select(mvpa_pca_df,c("Dex",paste0("PC",1:6))))
library(glmnet)
subject_mvpa <- sapply(subjects, function(subject) {
  subject_idx <- which(subjects == subject)
  loo_df <- mvpa_pca_df[-subject_idx,]
  lm_fit <- lm(Dex ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6,data = loo_df)
  # lasso_i <- cv.glmnet(x = as.matrix(loo_df[,-46]), y = loo_df[,46])
  # yhat <- predict(lasso_i, newx = as.matrix(mvpa_pca_df[subject_idx,-46,drop = F]), s = 'lambda.min')
  return(yhat)
}, simplify = T)

plot(mvpa_pca_df$Dex,subject_mvpa)


# Outlier work ----
subject_10 <- subjects[10]
left_result <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/all_sessions/501_BayesGLM2_subject_137128_left_thr0.rds")
right_result <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/all_sessions/501_BayesGLM2_subject_137128_right_thr0.rds")
library(ciftiTools)
ciftiTools.setOption('wb_path',"/Applications/workbench/")
cifti_obj <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
cifti_obj$data$cortex_left <- left_result$estimates
cifti_obj$data$cortex_right <- right_result$estimates
cifti_obj <- cifti_obj * cifti_mask
plot(cifti_obj, idx = 4, zlim = c(-1,1))

cifti_pca_tongue <- cifti_obj
sum(mask_left[,4], na.rm = T)
# 1023
sum(mask_right[,4], na.rm = T)
# 1055
cifti_pca_tongue$data$cortex_left[which(!is.na(cifti_pca_tongue$data$cortex_left[,4])),] <-
  tongue_mvpa$rotation[seq(1023),1]
cifti_pca_tongue$data$cortex_right[which(!is.na(cifti_pca_tongue$data$cortex_right[,4])),] <-
  tongue_mvpa$rotation[-seq(1023),1]
plot(cifti_pca_tongue,idx = 4)

sub1_result_left <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/all_sessions/501_BayesGLM2_subject_103818_left_thr0.rds")
sub1_result_right <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/all_sessions/501_BayesGLM2_subject_103818_right_thr0.rds")
sub1_cifti <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
sub1_cifti$data$cortex_left <- sub1_result_left$estimates
sub1_cifti$data$cortex_right <- sub1_result_right$estimates
sub1_cifti <- sub1_cifti * cifti_mask
plot(sub1_cifti, idx = 4, zlim = c(-1,1))

mvpa_df %>%
  select(-Dex, starts_with("left_tongue"), starts_with("right_tongue")) %>%
  pivot_longer(cols = -"subject", names_to = "location") %>%
  mutate(sub10 = subject == subject_10) %>%
  ggplot() +
  geom_histogram(aes(x = value), alpha = 0.3, position = "dodge", bins = 100) +
  facet_grid(sub10~., scales = "free") +
  lims(x = c(-2,3))

# >> Session estimates for subject 10 ----
visit1_left_result <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit1_left_5k_20210127.rds")
visit1_right_result <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit1_right_5k_20210127.rds")
visit2_left_result <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit2_left_5k_20210127.rds")
visit2_right_result <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit2_right_5k_20210127.rds")

subj10_results <- list(
  visit1 = list(
    left = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit1_left_5k_20210127.rds")$betas_Bayesian,
    right = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit1_right_5k_20210127.rds")$betas_Bayesian
  ),
  visit2 = list(
    left = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit2_left_5k_20210127.rds")$betas_Bayesian,
    right = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit2_right_5k_20210127.rds")$betas_Bayesian
  )
)

sub10_ciftis <- sapply(subj10_results, function(visit_results) {
  session_ciftis <- sapply(c("LR","RL","avg"), function(session_name) {
    cifti_out <- xifti_combine(visit_results$left[[session_name]], visit_results$right[[session_name]])
    return(cifti_out)
  }, simplify = F)
}, simplify = F)

library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench/')

for(v in paste0("visit",1:2)) {
  for(sess in c("LR","RL","avg")) {
    plot(sub10_ciftis[[v]][[sess]], title = paste("Subject 10",v,sess),
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         idx = 4)
  }
}

library(tidyverse)
reshape2::melt(all_mvpa_amplitudes) %>%
  filter(L3 == "tongue") %>%
  mutate(is_sub10 = L1 == "137128") %>%
  ggplot() +
  geom_histogram(aes(x = value), bins = 100) +
  facet_grid(is_sub10~., scales = "free")


reshape2::melt(all_mvpa_amplitudes) %>%
  filter(L3 == "tongue") %>%
  mutate(is_sub10 = L1 == "137128") %>%
  group_by(is_sub10) %>%
  summarize(Min = min(value),
            Mean = mean(value),
            Median = median(value),
            Max = max(value),
            SD = sd(value))

# >> Look at data ----
subj10_data <- list(
  visit1 = list(
    left = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit1_left_5k_20210127.rds")$GLMs_Bayesian$cortexL$,
    right = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit1_right_5k_20210127.rds")$betas_Bayesian
  ),
  visit2 = list(
    left = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit2_left_5k_20210127.rds")$betas_Bayesian,
    right = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_137128_visit2_right_5k_20210127.rds")$betas_Bayesian
  )
)

# >> Examine the amplitude estimates ----
sub10_est <- sapply(subj10_results, function(visit_results) {
  session_ciftis <- sapply(c("LR","RL","avg"), function(session_name) {
    est_out <- c(visit_results$left[[session_name]]$data$cortex_left[,4], visit_results$right[[session_name]]$data$cortex_right[,4])
    return(est_out)
  }, simplify = F)
}, simplify = F)

sapply(sub10_est, function(x) sapply(x, function(y) summary(y)), simplify = F)


subj1_results <- list(
  visit1 = list(
    left = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_103818_visit1_left_5k_20210125.rds")$betas_Bayesian,
    right = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_103818_visit1_right_5k_20210125.rds")$betas_Bayesian
  ),
  visit2 = list(
    left = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_103818_visit2_left_5k_20210125.rds")$betas_Bayesian,
    right = readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_103818_visit2_right_5k_20210125.rds")$betas_Bayesian
  )
)

sub1_est <- sapply(subj1_results, function(visit_results) {
  session_ciftis <- sapply(c("LR","RL","avg"), function(session_name) {
    est_out <- c(visit_results$left[[session_name]]$data$cortex_left[,4], visit_results$right[[session_name]]$data$cortex_right[,4])
    return(est_out)
  }, simplify = F)
}, simplify = F)

sapply(sub1_est, function(x) sapply(x, function(y) summary(y)), simplify = F)



# >> Examine the single-session results ----
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
sapply(ss_est, function(x) summary(x[,4]))

# Compare to subject 1
sub1_files <- grep("103818", list.files(result_dir, full.names = T), value = T)
ss1_est <- sapply(sub1_files, function(x) {
  obj <- readRDS(x)
  out <- do.call(rbind, obj$betas_Bayesian[[1]]$data)
  return(out)
}, simplify = F)
summary(unlist(ss1_est))
sapply(ss1_est, summary)
sapply(ss1_est, function(x) summary(x[,3]))

# >> Examine the design matrices ----
sub10_designs <- sapply(sub10_files, function(x) {
  readRDS(x)$design[[1]]
}, simplify = F)
names(sub10_designs) <- paste0("session",1:8)

sub1_designs <- sapply(sub1_files, function(x) {
  readRDS(x)$design[[1]][,4]
}, simplify = F)
names(sub1_designs) <- paste0("session",1:8)

sapply(sub10_designs, summary)
sapply(sub1_designs, summary)

sub_files <- grep(".rds",list.files(result_dir, full.names = T), value = T)
sub_files <- grep("_left_", sub_files, value = T)
sapply(subjects, function(x) {
  sfiles <- grep(x,sub_files, value = T)
  length(sfiles)
})
grep("151526", sub_files, value = T)

# >> Calculate the VIF ----
load("/Users/Shared/Lab_Data/HCP_Motor_Task_Dan/subjects.Rdata")
main_dir <- "/Users/Shared/Lab_Data/HCP_Motor_Task_Dan"
data_dir <- "/Users/Shared/Lab_Data/HCP_Motor_Task_Dan/visit1_data"
ss_result_dir <- "HCP_results/5k_results/individual/PW/single_session"
ss_result_files <- list.files(ss_result_dir, full.names = T)
ss_result_files <- grep("500_", ss_result_files, value = T)
library(clever)
library(BayesfMRI)
task_names <- c("cue","foot","hand","tongue")
VIF_list <- list()
for(subject in subjects) {
  VIF_list[[subject]] <- list()
  for(visit in 1:2) {
    VIF_list[[subject]][[paste0("visit",visit)]] <- list()
    #analyze hemispheres separately due to different in set of tasks
    for(h in 1:2){
      #h=1 -- left hemisphere
      #h=2 -- right hemisphere
      hem <- c('left','right')[h]
      VIF_list[[subject]][[paste0("visit",visit)]][[hem]] <- list()

      # find the VIF at the session level
      for (s in 1:2) {
        sess <- c("LR","RL")[s]
        motion <-
          as.matrix(
            read.table(
              file.path(
                main_dir,
                paste0("visit", visit, "_data"),
                subject,
                "MNINonLinear",
                "Results",
                paste0("tfMRI_MOTOR_",sess),
                "Movement_Regressors.txt"
              )
            , header = FALSE)
          )
        design <- readRDS(
          grep(
            paste0(subject,"_visit",visit,"_",hem,"_5k_session",sess),
            ss_result_files,
            value = T
          )
        )$design[[1]]
        drift1 <- (1:284)/284
        drift <- cbind(drift1, drift1^2)
        # design2 <- BayesfMRI:::nuisance_regress(design, cbind(motion, drift))
        # myFD <- FD(X = motion[,1:6])
        VIF_session <- matrix(NA,4,1)
        for(k in 1:4){
          VIF_session[k,1] <- 1 / (1 - summary(lm(design[,k] ~ as.matrix(cbind(motion, drift))))$r.squared)
        }
        rownames(VIF_session) <- task_names
        VIF_list[[subject]][[paste0("visit",visit)]][[hem]][[sess]] <- VIF_session
      }
    }
  }
}
saveRDS(VIF_list, "~/Desktop/505_VIF_list.rds")
