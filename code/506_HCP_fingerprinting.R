# This is a script for a fingerprinting analysis of the HCP retest data

# BAYESIAN ----
# 1) Make masks of the active regions for each task ----
mask_left <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/HCP_BayesGLM2_45subj_result_left_thresh0_20210131.rds")$active
mask_right <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/HCP_BayesGLM2_45subj_result_right_thresh0_20210131.rds")$active
masks <- list(
  left = mask_left,
  right = mask_right
)
# What if there is no masking?
masks <- list(
  left = matrix(1,4443,4),
  right = matrix(1,4444,4)
)
# What about masking out the limbic system and the DMN?


# 2) Grab the subject-level activation estimates ----
# load("~/github/BayesGLM_Validation/HCP_data/subjects.Rdata")
# result_files <- list.files("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW", full.names = T)
#
# visit_estimates <-
#   sapply(paste0("visit",1:2), function(visit) {
#     visit_files <- grep(visit, result_files, value = T)
#     subject_esitmates <- sapply(subjects, function(subject) {
#       subject_files <- grep(subject, visit_files, value = T)
#       subject_mask_est <- mapply(function(hem,mas) {
#         result_obj <- readRDS(grep(hem,subject_files, value = T))
#         result_avg <- result_obj$betas_Bayesian$avg$data[[paste0("cortex_",hem)]]
#         return(result_avg * mas)
#       }, hem = c("left","right"), mas = masks, SIMPLIFY = FALSE)
#     }, simplify = FALSE)
#   }, simplify = FALSE)
# saveRDS(visit_estimates, "~/github/BayesGLM_Validation/HCP_results/5k_results/506_HCP_masked_estimates.rds")

visit_estimates <- readRDS("~/github/BayesGLM_Validation/HCP_results/5k_results/506_HCP_masked_estimates.rds")
task_names <- c("cue","foot","hand","tongue")
library(tidyverse)
est_df <-
  reshape2::melt(visit_estimates) %>%
  filter(value != 0) %>%
  mutate(Var2 = task_names[Var2],
         Var1 = ifelse(L3 == 'left', Var1 + 4444, Var1),
         side = ifelse(L3 == 'left','right','left'),
         task = ifelse(
           Var2 %in% c('foot','hand'),
           paste0(side,"_",Var2),
           Var2)) %>%
  rename(subject = L2,
         visit = L1,
         location = Var1) %>%
  select(subject, visit, task, location, value)

# 3) Perform the fingerprinting ----
#  >> All tasks together ----
visit1_subj <-
  filter(est_df, visit == "visit1") %>%
  pivot_wider(id_cols = c(task,location), names_from = subject, values_from = value) %>%
  select(-task,-location) %>%
  as.matrix()

visit2_subj <-
  filter(est_df, visit == "visit2") %>%
  pivot_wider(id_cols = c(task,location), names_from = subject, values_from = value) %>%
  select(-task,-location) %>%
  as.matrix()

fingerprint_corr <- matrix(NA,45,45)
for(i in 1:45) {
  for(j in 1:45) {
    fingerprint_corr[i,j] <- cor(visit1_subj[,i], visit2_subj[,j])
  }
}

visit1ref <- apply(fingerprint_corr,1,which.max)
sum(apply(cbind(1:45,visit1ref),1,function(x) x[1] == x[2]))/45

visit2ref <- apply(fingerprint_corr,2,which.max)
sum(apply(cbind(1:45,visit2ref),1,function(x) x[1] == x[2]))/45

reshape2::melt(fingerprint_corr) %>%
  ggplot() +
  geom_raster(aes(x = Var1, y = Var2, fill = value)) +
  scale_fill_distiller("Correlation",palette = "Greys") +
  labs(x = "Visit 1 Subject", y = "Visit 2 Subject") +
  theme_classic()

# >> Tasks separately ----
for(each_task in unique(est_df$task)) {
  visit1_subj <-
    filter(est_df, visit == "visit1", task == each_task) %>%
    pivot_wider(id_cols = c(task,location), names_from = subject, values_from = value) %>%
    select(-task,-location) %>%
    as.matrix()

  visit2_subj <-
    filter(est_df, visit == "visit2", task == each_task) %>%
    pivot_wider(id_cols = c(task,location), names_from = subject, values_from = value) %>%
    select(-task,-location) %>%
    as.matrix()

  fingerprint_corr <- matrix(NA,45,45)
  for(i in 1:45) {
    for(j in 1:45) {
      fingerprint_corr[i,j] <- cor(visit1_subj[,i], visit2_subj[,j])
    }
  }

  visit1ref <- apply(fingerprint_corr,1,which.max)
  cat("Fingerprinting success with a visit 1 reference for",each_task,"is",sum(apply(cbind(1:45,visit1ref),1,function(x) x[1] == x[2]))/45,"\n")

  visit2ref <- apply(fingerprint_corr,2,which.max)
  cat("Fingerprinting success with a visit 2 reference for",each_task,"is",sum(apply(cbind(1:45,visit2ref),1,function(x) x[1] == x[2]))/45,"\n")
}

# Fingerprinting success with a visit 1 reference for cue is 0.9333333
# Fingerprinting success with a visit 2 reference for cue is 0.9777778
# Fingerprinting success with a visit 1 reference for right_foot is 0.4444444
# Fingerprinting success with a visit 2 reference for right_foot is 0.4666667
# Fingerprinting success with a visit 1 reference for right_hand is 0.4888889
# Fingerprinting success with a visit 2 reference for right_hand is 0.4666667
# Fingerprinting success with a visit 1 reference for tongue is 0.5111111
# Fingerprinting success with a visit 2 reference for tongue is 0.4888889
# Fingerprinting success with a visit 1 reference for left_foot is 0.5111111
# Fingerprinting success with a visit 2 reference for left_foot is 0.4
# Fingerprinting success with a visit 1 reference for left_hand is 0.6222222
# Fingerprinting success with a visit 2 reference for left_hand is 0.5555556

# This is obviously very successful for the cue, so just do that one
# >> Fingerprinting with cue ----
visit1_subj <-
  filter(est_df, visit == "visit1", task =="cue") %>%
  pivot_wider(id_cols = c(task,location), names_from = subject, values_from = value) %>%
  select(-task,-location) %>%
  as.matrix()

visit2_subj <-
  filter(est_df, visit == "visit2", task == "cue") %>%
  pivot_wider(id_cols = c(task,location), names_from = subject, values_from = value) %>%
  select(-task,-location) %>%
  as.matrix()

fingerprint_corr <- matrix(NA,45,45)
for(i in 1:45) {
  for(j in 1:45) {
    fingerprint_corr[i,j] <- cor(visit1_subj[,i], visit2_subj[,j])
  }
}

visit1ref <- apply(fingerprint_corr,1,which.max)
sum(apply(cbind(1:45,visit1ref),1,function(x) x[1] == x[2]))/45

visit2ref <- apply(fingerprint_corr,2,which.max)
sum(apply(cbind(1:45,visit2ref),1,function(x) x[1] == x[2]))/45

reshape2::melt(fingerprint_corr) %>%
  ggplot() +
  geom_raster(aes(x = Var1, y = Var2, fill = value)) +
  scale_fill_distiller("Correlation",palette = "Greys") +
  labs(x = "Visit 1 Subject", y = "Visit 2 Subject") +
  theme_classic()

# CLASSICAL ----
# Fingerprinting using the Bayesian masks ----
# 1) Grab the Bayesian masks ----
mask_left <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/HCP_BayesGLM2_45subj_result_left_thresh0_20210131.rds")$active
mask_right <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/HCP_BayesGLM2_45subj_result_right_thresh0_20210131.rds")$active
masks <- list(
  left = mask_left,
  right = mask_right
)

# 2) Grab the classical estimates ----
classical_files <- list.files("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/classical", full.names = T)
load("~/github/BayesGLM_Validation/HCP_data/subjects.Rdata")

visit_estimates_classical <-
  sapply(paste0("visit",1:2), function(visit) {
    visit_files <- grep(visit, classical_files, value = T)
    subject_esitmates <- sapply(subjects, function(subject) {
      subject_files <- grep(subject, visit_files, value = T)
      subject_mask_est <- mapply(function(hem,mas) {
        result_obj <- readRDS(grep(hem,subject_files, value = T))
        result_avg <- result_obj$betas_classical$avg$data[[paste0("cortex_",hem)]]
        return(result_avg * mas)
      }, hem = c("left","right"), mas = masks, SIMPLIFY = FALSE)
    }, simplify = FALSE)
  }, simplify = FALSE)

task_names <- c("cue","foot","hand","tongue")
library(tidyverse)
est_df <-
  reshape2::melt(visit_estimates_classical) %>%
  filter(value != 0) %>%
  mutate(Var2 = task_names[Var2],
         Var1 = ifelse(L3 == 'left', Var1 + 4443, Var1),
         side = ifelse(L3 == 'left','right','left'),
         task = ifelse(
           Var2 %in% c('foot','hand'),
           paste0(side,"_",Var2),
           Var2)) %>%
  rename(subject = L2,
         visit = L1,
         location = Var1) %>%
  select(subject, visit, task, location, value)

# 3) Perform fingerprinting with cue ----
visit1_subj <-
  filter(est_df, visit == "visit1", task =="cue") %>%
  pivot_wider(id_cols = c(task,location), names_from = subject, values_from = value) %>%
  select(-task,-location) %>%
  as.matrix()

visit2_subj <-
  filter(est_df, visit == "visit2", task == "cue") %>%
  pivot_wider(id_cols = c(task,location), names_from = subject, values_from = value) %>%
  select(-task,-location) %>%
  as.matrix()

fingerprint_corr <- matrix(NA,45,45)
for(i in 1:45) {
  for(j in 1:45) {
    fingerprint_corr[i,j] <- cor(visit1_subj[,i], visit2_subj[,j])
  }
}

visit1ref <- apply(fingerprint_corr,1,which.max)
sum(apply(cbind(1:45,visit1ref),1,function(x) x[1] == x[2]))/45

visit2ref <- apply(fingerprint_corr,2,which.max)
sum(apply(cbind(1:45,visit2ref),1,function(x) x[1] == x[2]))/45

reshape2::melt(fingerprint_corr) %>%
  ggplot() +
  geom_raster(aes(x = Var1, y = Var2, fill = value)) +
  scale_fill_distiller("Correlation",palette = "Greys") +
  labs(x = "Visit 1 Subject", y = "Visit 2 Subject") +
  theme_classic()

# Fingerprinting using the classical masks ----
# 1) Grab the classical masks ----
classical_activations <- readRDS("HCP_results/5k_results/group/502_HCP_classical_activations_FWER.rds")
masks <- list(
  left = classical_activations$left$`0%`,
  right = classical_activations$right$`0%`
)
# What if there is no masking?
masks <- list(
  left = matrix(1,4443,4),
  right = matrix(1,4444,4)
)


# 2) Grab the classical estimates ----
classical_files <- list.files("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/classical", full.names = T)
load("~/github/BayesGLM_Validation/HCP_data/subjects.Rdata")

visit_estimates_classical <-
  sapply(paste0("visit",1:2), function(visit) {
    visit_files <- grep(visit, classical_files, value = T)
    subject_esitmates <- sapply(subjects, function(subject) {
      subject_files <- grep(subject, visit_files, value = T)
      subject_mask_est <- mapply(function(hem,mas) {
        result_obj <- readRDS(grep(hem,subject_files, value = T))
        result_avg <- result_obj$betas_classical$avg$data[[paste0("cortex_",hem)]]
        return(result_avg * mas)
      }, hem = c("left","right"), mas = masks, SIMPLIFY = FALSE)
    }, simplify = FALSE)
  }, simplify = FALSE)

task_names <- c("cue","foot","hand","tongue")
library(tidyverse)
est_df <-
  reshape2::melt(visit_estimates_classical) %>%
  filter(value != 0) %>%
  mutate(Var2 = task_names[Var2],
         Var1 = ifelse(L3 == 'left', Var1 + 4444, Var1),
         side = ifelse(L3 == 'left','right','left'),
         task = ifelse(
           Var2 %in% c('foot','hand'),
           paste0(side,"_",Var2),
           Var2)) %>%
  rename(subject = L2,
         visit = L1,
         location = Var1) %>%
  select(subject, visit, task, location, value)

# 3) Perform fingerprinting ----
# >> Cue only ----
visit1_subj <-
  filter(est_df, visit == "visit1", task =="cue") %>%
  pivot_wider(id_cols = c(task,location), names_from = subject, values_from = value) %>%
  select(-task,-location) %>%
  as.matrix()

visit2_subj <-
  filter(est_df, visit == "visit2", task == "cue") %>%
  pivot_wider(id_cols = c(task,location), names_from = subject, values_from = value) %>%
  select(-task,-location) %>%
  as.matrix()

fingerprint_corr <- matrix(NA,45,45)
for(i in 1:45) {
  for(j in 1:45) {
    fingerprint_corr[i,j] <- cor(visit1_subj[,i], visit2_subj[,j])
  }
}

visit1ref <- apply(fingerprint_corr,1,which.max)
sum(apply(cbind(1:45,visit1ref),1,function(x) x[1] == x[2]))/45

visit2ref <- apply(fingerprint_corr,2,which.max)
sum(apply(cbind(1:45,visit2ref),1,function(x) x[1] == x[2]))/45

reshape2::melt(fingerprint_corr) %>%
  ggplot() +
  geom_raster(aes(x = Var1, y = Var2, fill = value)) +
  scale_fill_distiller("Correlation",palette = "Greys") +
  labs(x = "Visit 1 Subject", y = "Visit 2 Subject") +
  theme_classic()

# >> Tasks separately ----
for(each_task in unique(est_df$task)) {
  visit1_subj <-
    filter(est_df, visit == "visit1", task == each_task) %>%
    pivot_wider(id_cols = c(task,location), names_from = subject, values_from = value) %>%
    select(-task,-location) %>%
    as.matrix()

  visit2_subj <-
    filter(est_df, visit == "visit2", task == each_task) %>%
    pivot_wider(id_cols = c(task,location), names_from = subject, values_from = value) %>%
    select(-task,-location) %>%
    as.matrix()

  fingerprint_corr <- matrix(NA,45,45)
  for(i in 1:45) {
    for(j in 1:45) {
      fingerprint_corr[i,j] <- cor(visit1_subj[,i], visit2_subj[,j])
    }
  }

  visit1ref <- apply(fingerprint_corr,1,which.max)
  cat("Fingerprinting success with a visit 1 reference for",each_task,"is",sum(apply(cbind(1:45,visit1ref),1,function(x) x[1] == x[2]))/45,"\n")

  visit2ref <- apply(fingerprint_corr,2,which.max)
  cat("Fingerprinting success with a visit 2 reference for",each_task,"is",sum(apply(cbind(1:45,visit2ref),1,function(x) x[1] == x[2]))/45,"\n")
}
