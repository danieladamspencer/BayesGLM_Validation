# This is a script for performing some basic analyses using the demographics
# data for the Human Connectome Project

# Compile area of activation sizes for each subject ----
# Access the sizes of the areas of activation at different thresholds for each
# subject
# library(BayesfMRI)
# result_dir <-
#   "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
# result_files <-
#   grep("500_", list.files(result_dir, full.names = F), value = T)
# load("HCP_data/subjects.Rdata")
# hem <- c("left", 'right')
# visit <- paste0("visit", 1:2)
# thresholds <- c(0, 0.5, 1)
# task_names <- c("visual cue", "foot", "hand", "tongue")
# activation_sizes <-
#   expand.grid(
#     subject = subjects,
#     visit = visit,
#     hem = hem,
#     thr = thresholds,
#     task = task_names
#   )
# activation_sizes$area <- NA
# activation_sizes <-
#   dplyr::arrange(activation_sizes, subject, visit, hem, thr)
# for (subj in subjects) {
#   for (v in visit) {
#     for (h in hem) {
#       subject_file <-
#         grep(paste0("500_", subj, "_", v, "_", h, "_5k_"),
#              result_files,
#              value = T)
#       result_obj <- readRDS(file.path(result_dir, subject_file))
#       L_or_R <- toupper(substring(h, 1, 1))
#       for (thr in thresholds) {
#         exc_obj <-
#           id_activations.posterior(result_obj$GLMs_Bayesian[[paste0("cortex", L_or_R)]],
#                                    alpha = 0.01,
#                                    threshold = thr)
#         areas <- matrixStats::colSums2(exc_obj$active)
#         df_idx <- which(
#           activation_sizes$subject == subj &
#             activation_sizes$visit == v &
#             activation_sizes$hem == h &
#             activation_sizes$thr == thr
#         )
#         activation_sizes$area[df_idx] <- areas
#       }
#     }
#   }
# }
# saveRDS(activation_sizes,
#         "HCP_results/5k_results/504_activation_sizes_df.rds")

# EDA on activation_sizes ----
active_df <- readRDS("HCP_results/5k_results/504_activation_sizes_df.rds")

library(tidyverse)
active_df %>%
  mutate(thr = factor(thr, levels = c(0,0.5,1))) %>%
  filter(visit == "visit1") %>%
  ggplot() +
  geom_boxplot(aes(y = area, fill = subject)) +
  facet_grid(thr~task, scales = "free") +
  labs(y = "# Active Locations", x = "") +
  theme_bw() +
  theme(legend.position = 'none',
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())


active_df %>%
  mutate(thr = factor(thr, levels = c(0,0.5,1))) %>%
  filter(hem == "left") %>%
  ggplot() +
  geom_point(aes(x = subject, y = area, color = subject, shape = visit)) +
  facet_grid(thr~task, scales = "free") +
  labs(x = "", y = "# Active Locations", title = "Left Hemisphere") +
  theme_bw() +
  theme(legend.position = 'none',
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

active_df %>%
  pivot_wider(names_from = visit, values_from = area) %>%
  group_by(hem, thr, task) %>%
  summarize(correlation = cor(visit1, visit2)) %>%
  # pivot_wider(names_from = task, values_from = correlation)
  ggplot(aes(x = thr, y = correlation)) +
  geom_point() +
  geom_smooth( method = "loess") +
  facet_grid(task~hem, scales = "free") +
  labs(x = "Activation Threshold", y = "Inter-visit Correlation\n(# of active locations)")

# MODELS USING DEMOGRAPHICS ----

# Combining active areas by visit ----
# Note: This is an important step that should effectively combine the # of
# active data locations by only using the # of overlapping active regions. In
# order to complete this calculation, the Dice coefficient is used as a shortcut
# to circumvent the need to calculate the excursions all over again.
# dice_dir <- "HCP_results/5k_results"
# overlap_files <- grep("604_overlap_size_PW_", list.files(dice_dir), value = T)[c(2,1,3)]
# task_names <- c("visual cue","tongue","right foot","right hand","left foot","left hand")
# load("HCP_data/subjects.Rdata")
# overlap_df <- sapply(c(0,0.5,1), function(thr){
#   out_df <- readRDS(file.path(dice_dir,paste0("604_overlap_size_PW_threshold",thr,".rds")))
#   colnames(out_df) <- task_names
#   out_df <- as.data.frame(out_df)
#   out_df$thr <- thr
#   out_df$subject <- subjects
#   return(out_df)
# }, simplify = F)
# overlap_df <- Reduce(rbind,overlap_df)
# saveRDS(overlap_df, "HCP_results/5k_results/504_overlap_sizes_df.rds")

# Load and read the data ----
library(tidyverse)
# active_df <- readRDS("HCP_results/5k_results/504_activation_sizes_df.rds")
# demog_df <- read.csv("HCP_data/unrestricted_danspen_2_24_2021_12_55_28.csv") %>%
#   mutate(Subject = factor(Subject))
# hcp_df <-
#   active_df %>%
#   left_join(demog_df, by = c("subject" = "Subject")) %>%
#   select(subject, visit, hem, thr, task, area, Endurance_Unadj,
#          Strength_Unadj, ProcSpeed_Unadj, Dexterity_Unadj) %>%
#   pivot_wider(names_from = c(visit,hem,task,thr), values_from = area)

overlap_df <- readRDS("HCP_results/5k_results/504_overlap_sizes_df.rds")
demog_df <- read.csv("HCP_data/unrestricted_danspen_2_24_2021_12_55_28.csv") %>%
  mutate(Subject = factor(Subject))
hcp_df <-
  overlap_df %>%
  pivot_longer(c(`visual cue`, tongue, `right foot`, `right hand`, `left foot`, `left hand`), names_to="task",values_to = "area") %>%
  mutate(area = 100* area / (4443+4444)) %>% # This line changes area from absolute to proportional
  left_join(demog_df, by = c("subject" = "Subject")) %>%
  select(subject, thr, task, area, Endurance_Unadj,
         Strength_Unadj, ProcSpeed_Unadj, Dexterity_Unadj) %>%
  pivot_wider(names_from = c(task,thr), values_from = area)

hcp_df %>%
  select(-ends_with("Unadj"), - subject) %>%
  `*`((4443+4444)/100) %>% summary

# Endurance ----
# endurance_df <- select(hcp_df, Endurance_Unadj, starts_with("visit"))
# plot(density(endurance_df$Endurance_Unadj))
# pairs(select(endurance_df, Endurance_Unadj, starts_with("visit1")))
# pairs(select(endurance_df, Endurance_Unadj, starts_with("visit2")))
plot(density(hcp_df$Endurance_Unadj))
pairs(Endurance_Unadj ~ ., data = hcp_df[,-1])
endurance_df <- select(hcp_df, -Strength_Unadj, -ProcSpeed_Unadj, -Dexterity_Unadj, -subject) %>%
  mutate(Endurance_Unadj = scale(Endurance_Unadj))
end_cor <- cor(endurance_df)
sort(end_cor[-1,1])
endur_big_lm <- lm(Endurance_Unadj~ ., data = endurance_df)
summary(endur_big_lm)
endur_step <- step(endur_big_lm)
summary(endur_step)
plot(endur_step)
plot(endur_step$model$Endurance_Unadj,endur_step$fitted.values,
     xlab = "Actual", ylab = "Fitted", main = "Endurance")
abline(a=0,b=1, col='red')

# Speed ----
plot(density(hcp_df$ProcSpeed_Unadj))
pairs(ProcSpeed_Unadj ~ ., data = hcp_df[,-1])
speed_df <- select(hcp_df, -Strength_Unadj, -Endurance_Unadj, -Dexterity_Unadj, -subject) %>%
  # mutate(ProcSpeed_Unadj = scale(ProcSpeed_Unadj))
  mutate(hand_0 = `left hand_0` + `right hand_0`,
         hand_0.5 = `left hand_0.5` + `right hand_0.5`,
         hand_1 = `left hand_1` + `right hand_1`,
         foot_0 = `left foot_0` + `right foot_0`,
         foot_0.5 = `left foot_0.5` + `right foot_0.5`,
         foot_1 = `left foot_1` + `right foot_1`) %>%
  select(-starts_with("right"), -starts_with("left"))
spd_cor <- cor(speed_df)
sort(spd_cor[-1,1])
speed_big_lm <- lm(ProcSpeed_Unadj~ ., data = speed_df)
summary(speed_big_lm)
speed_step <- step(speed_big_lm)
summary(speed_step)
plot(speed_step)
plot(speed_step$model$ProcSpeed_Unadj,speed_step$fitted.values,
     xlab = "Actual", ylab = "Fitted", main = "ProcSpeed")
abline(a=0,b=1, col='red')

spd_lm0 <- lm(ProcSpeed_Unadj~., data = select(speed_df, ProcSpeed_Unadj, ends_with("_0")))
summary(spd_lm0)
spd_lm0.5 <- lm(ProcSpeed_Unadj~., data = select(speed_df, ProcSpeed_Unadj, ends_with("_0.5")))
summary(spd_lm0.5)
spd_lm0_0.5 <- lm(ProcSpeed_Unadj~., data = select(speed_df, ProcSpeed_Unadj, ends_with("_0"),ends_with("_0.5")))
summary(spd_lm0_0.5)
spd_lm1 <- lm(ProcSpeed_Unadj~., data = select(speed_df, ProcSpeed_Unadj, ends_with("_1")))
summary(spd_lm1)
spd_lm0_1 <- lm(ProcSpeed_Unadj~., data = select(speed_df, ProcSpeed_Unadj, ends_with("_0"),ends_with("_1")))
summary(spd_lm0_1)

# >>  Elastic Net Regression ----
spd_01_df <- select(speed_df, ProcSpeed_Unadj, ends_with("_0"), ends_with("_1"))
library(glmnet)
spd_cv_lasso <-
  cv.glmnet(
    x = as.matrix(spd_01_df[, -1]),
    y = as.matrix(spd_01_df[, 1]),
    family = "gaussian",
    alpha = 1
  )
spd_lambda_lasso <- spd_cv_lasso$lambda[which.min(spd_cv_lasso$cvm)]
spd_lasso <-
  glmnet(
    x = as.matrix(spd_01_df[, -1]),
    y = as.matrix(spd_01_df[, 1]),
    family = "gaussian",
    lambda = spd_lambda,
    alpha = 1
  )
spd_cv_ridge <-
  cv.glmnet(
    x = as.matrix(spd_01_df[, -1]),
    y = as.matrix(spd_01_df[, 1]),
    family = "gaussian",
    alpha = 0
  )
spd_lambda_ridge <- spd_cv_ridge$lambda[which.min(spd_cv_ridge$cvm)]
spd_ridge <-
  glmnet(
    x = as.matrix(spd_01_df[, -1]),
    y = as.matrix(spd_01_df[, 1]),
    family = "gaussian",
    lambda = spd_lambda_ridge,
    alpha = 0
  )
cbind(lasso = spd_lasso$beta, ridge =  spd_ridge$beta)
min(spd_cv_lasso$cvm)
min(spd_cv_ridge$cvm)

# Strength ----
plot(density(hcp_df$Strength_Unadj))
pairs(Strength_Unadj ~ ., data = hcp_df[,-1])
strength_df <- select(hcp_df, -ProcSpeed_Unadj, -Endurance_Unadj, -Dexterity_Unadj, -subject) %>%
  mutate(Strength_Unadj = scale(Strength_Unadj))
str_cor <- cor(strength_df)
sort(str_cor[-1,1])
strength_big_lm <- lm(Strength_Unadj~ ., data = strength_df)
summary(strength_big_lm)
strength_step <- step(strength_big_lm)
summary(strength_step)
plot(strength_step)
plot(strength_step$model$Strength_Unadj,strength_step$fitted.values,
     xlab = "Actual", ylab = "Fitted", main = "Strength")
abline(a=0,b=1, col='red')

# Dexterity ----
plot(density(hcp_df$Dexterity_Unadj))
pairs(Dexterity_Unadj ~ ., data = hcp_df[,-1])
dexterity_df <- select(hcp_df, -ProcSpeed_Unadj, -Endurance_Unadj, -Strength_Unadj, -subject) %>%
  # mutate(Dexterity_Unadj = scale(Dexterity_Unadj))
  mutate(hand_0 = `left hand_0` + `right hand_0`,
         hand_0.5 = `left hand_0.5` + `right hand_0.5`,
         hand_1 = `left hand_1` + `right hand_1`,
         foot_0 = `left foot_0` + `right foot_0`,
         foot_0.5 = `left foot_0.5` + `right foot_0.5`,
         foot_1 = `left foot_1` + `right foot_1`) %>%
  select(-starts_with("right"), -starts_with("left"))
dex_cor <- cor(dexterity_df)
sort(dex_cor[-1,1])
dexterity_big_lm <- lm(Dexterity_Unadj~., data = dexterity_df)
summary(dexterity_big_lm)
dexterity_step <- step(dexterity_big_lm)
summary(dexterity_step)
plot(dexterity_step)
plot(dexterity_step$model$Dexterity_Unadj,dexterity_step$fitted.values,
     xlab = "Actual", ylab = "Fitted", main = "Dexterity")
cor(dexterity_step$model$Dexterity_Unadj,dexterity_step$fitted.values)
abline(a=0,b=1, col='red')

dex_lm0 <- lm(Dexterity_Unadj~., data = select(dexterity_df, Dexterity_Unadj, ends_with("_0")))
summary(dex_lm0)
dex_lm0.5 <- lm(Dexterity_Unadj~., data = select(dexterity_df, Dexterity_Unadj, ends_with("_0.5")))
summary(dex_lm0.5)
dex_lm0_0.5 <- lm(Dexterity_Unadj~., data = select(dexterity_df, Dexterity_Unadj, ends_with("_0"),ends_with("_0.5")))
summary(dex_lm0_0.5)
dex_lm1 <- lm(Dexterity_Unadj~., data = select(dexterity_df, Dexterity_Unadj, ends_with("_1")))
summary(dex_lm1)
dex_lm0_1 <- lm(Dexterity_Unadj~., data = select(dexterity_df, Dexterity_Unadj, ends_with("_0"),ends_with("_1")))
summary(dex_lm0_1)
dex_final <- lm(Dexterity_Unadj~foot_1 + hand_1 + tongue_0 + `visual cue_1`, data = dexterity_df)
summary(dex_final)
pairs(cbind(dex_final$residuals, dexterity_df))
sort(cor(cbind(dex_final$residuals, dexterity_df))[1,-1])

# >> Elastic Net regression ----
library(glmnet)
dex_cv_lasso <-
  cv.glmnet(
    x = as.matrix(dexterity_df[, -1]),
    y = as.matrix(dexterity_df[, 1]),
    family = "gaussian",
    alpha = 1
  )
dex_lambda_lasso <- dex_cv_lasso$lambda[which.min(dex_cv_lasso$cvm)]
dex_lasso <-
  glmnet(
    x = as.matrix(dexterity_df[, -1]),
    y = as.matrix(dexterity_df[, 1]),
    family = "gaussian",
    lambda = dex_lambda_lasso,
    alpha = 1
  )
dex_cv_ridge <-
  cv.glmnet(
    x = as.matrix(dexterity_df[, -1]),
    y = as.matrix(dexterity_df[, 1]),
    family = "gaussian",
    alpha = 0
  )
dex_lambda_ridge <- dex_cv_ridge$lambda[which.min(dex_cv_ridge$cvm)]
dex_ridge <-
  glmnet(
    x = as.matrix(dexterity_df[, -1]),
    y = as.matrix(dexterity_df[, 1]),
    family = "gaussian",
    lambda = dex_lambda_ridge,
    alpha = 0
  )
cbind(lasso = dex_lasso$beta, ridge =  dex_ridge$beta)
dex_cv_lasso$cvm[which.min(dex_cv_lasso$cvm)]
dex_cv_ridge$cvm[which.min(dex_cv_ridge$cvm)]

dex_01_df <- select(dexterity_df, Dexterity_Unadj, ends_with("_0"), ends_with("_1"))
dex_cv_lasso <-
  cv.glmnet(
    x = as.matrix(dex_01_df[, -1]),
    y = as.matrix(dex_01_df[, 1]),
    family = "gaussian",
    alpha = 1
  )
dex_lambda_lasso <- dex_cv_lasso$lambda[which.min(dex_cv_lasso$cvm)]
dex_lasso <-
  glmnet(
    x = as.matrix(dex_01_df[, -1]),
    y = as.matrix(dex_01_df[, 1]),
    family = "gaussian",
    lambda = dex_lambda_lasso,
    alpha = 1
  )
dex_cv_ridge <-
  cv.glmnet(
    x = as.matrix(dex_01_df[, -1]),
    y = as.matrix(dex_01_df[, 1]),
    family = "gaussian",
    alpha = 0
  )
dex_lambda_ridge <- dex_cv_ridge$lambda[which.min(dex_cv_ridge$cvm)]
dex_ridge <-
  glmnet(
    x = as.matrix(dex_01_df[, -1]),
    y = as.matrix(dex_01_df[, 1]),
    family = "gaussian",
    lambda = dex_lambda_ridge,
    alpha = 0
  )
cbind(lasso = dex_lasso$beta, ridge =  dex_ridge$beta)
min(dex_cv_lasso$cvm)
min(dex_cv_ridge$cvm)
cor(predict(dex_ridge, as.matrix(dex_01_df[,-1])), dex_01_df$Dexterity_Unadj)^2
cor(predict(dex_lasso, as.matrix(dex_01_df[,-1])), dex_01_df$Dexterity_Unadj)^2
