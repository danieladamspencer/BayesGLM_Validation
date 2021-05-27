# This is a script meant to provide reproducible code for each figure in the
# manuscript

# FIGURE 1: Graphical Abstract (not used in figure) ----
# >> Prior Precision ----
tau2 <- 2
kappa2 <- 2
C <- diag(11)
G <- matrix(
  c(
    1,1,0,1,0,0,0,0,0,0,0,
    1,1,1,1,1,0,0,0,0,0,0,
    0,1,1,0,1,1,0,0,0,0,0,
    1,1,0,1,1,0,1,0,0,0,0,
    0,1,1,1,1,1,1,1,0,0,0,
    0,0,1,0,1,1,0,1,1,0,0,
    0,0,0,1,1,0,1,1,0,1,0,
    0,0,0,0,1,1,1,1,1,1,0,
    0,0,0,0,0,1,0,1,1,1,1,
    0,0,0,0,0,0,0,1,1,1,1,
    0,0,0,0,0,0,0,0,1,1,1
  ),
  nrow = 11,
  ncol = 11,
  byrow = T
)
reshape2::melt(G) %>%
  ggplot() +
  geom_raster(aes(x = Var1, y = Var2, fill = value)) +
  scale_fill_gradient2(low = "blue",high = "black") +
  guides(fill = F) +
  theme_void()

Q <- tau2*(kappa2^2*C + 2*kappa2*G + G%*%solve(C)%*%G)
Q_inv <- solve(Q)
reshape2::melt(Q) %>%
  ggplot() +
  geom_raster(aes(x = Var1, y = Var2, fill = value)) +
  scale_fill_gradient2(low = "blue",high = "red") +
  guides(fill = F) +
  theme_void()

# FIGURE 2: Single-subject estimates ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench/')
library(INLA)
inla.setOption(pardiso.lic = "~/pardiso.lic")
library(BayesfMRI)
load("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
# >> Single run ----
# >>>> Bayesian ----
library(ciftiTools)
ciftiTools.setOption('wb_path',"/Applications/workbench")
library(BayesfMRI)
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session"
result_files <- list.files(result_dir, full.names = T)
result_files <- grep(".rds", result_files, value= T)
load("HCP_data/subjects.Rdata")
subjects <- subjects[1]
sessions <- c("LR")
task_idx <- 3
task_names <- c("visual_cue","foot","hand","tongue")
for(subject in subjects) {
  for(v in 1) {
    for(sess in sessions) {
      cifti_obj <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
      for(h in c("left","right")) {
        result_obj <-
          readRDS(grep(
            paste0("500_", subject, "_visit", v, "_", h, "_5k_session", sess),
            result_files,
            value = T
          ))
        cifti_obj$data[[paste0("cortex_",h)]] <- result_obj$betas_Bayesian[[sess]]$data[[paste0("cortex_",h)]]
      }
      if(task_idx %in% c(1,4)) {
        plot(
          cifti_obj,
          idx = task_idx,zlim = c(-1,1), legend_embed = F,
          fname = paste0(
            "plots/607_bayes_",
            subject,
            "_visit",
            v,
            "_session",
            sess,
            "_",
            task_names[task_idx],
            "_estimate.png"
          ),
          surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
          surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii"
        )
      }
      if(task_idx %in% c(2,3)) {
        for(hem in c('left','right')) {
          plot(
            cifti_obj, hemisphere = hem,
            idx = task_idx,zlim = c(-1,1), legend_embed = F,
            fname = paste0(
              "plots/607_bayes_",
              subject,
              "_visit",
              v,
              "_session",
              sess,
              "_",
              hem,
              "_",
              task_names[task_idx],
              "_estimate.png"
            ),
            surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
            surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii"
          )
        }
      }
    }
  }
}

# >>>> Classical ----
library(ciftiTools)
ciftiTools.setOption('wb_path',"/Applications/workbench")
library(BayesfMRI)
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session"
result_files <- list.files(result_dir, full.names = T)
result_files <- grep(".rds", result_files, value= T)
load("HCP_data/subjects.Rdata")
subjects <- subjects[1]
sessions <- c("LR")
task_idx <- 3
task_names <- c("visual_cue","foot","hand","tongue")
for(subject in subjects) {
  for(v in 1) {
    for(sess in sessions) {
      cifti_obj <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
      for(h in c("left","right")) {
        result_obj <-
          readRDS(grep(
            paste0("500_", subject, "_visit", v, "_", h, "_5k_session", sess),
            result_files,
            value = T
          ))
        cifti_obj$data[[paste0("cortex_",h)]] <- result_obj$betas_classical[[sess]]$data[[paste0("cortex_",h)]]
      }
      if(task_idx %in% c(1,4)) {
        plot(
          cifti_obj,
          idx = task_idx,zlim = c(-1,1), legend_embed = F,
          fname = paste0(
            "plots/607_classical_",
            subject,
            "_visit",
            v,
            "_session",
            sess,
            "_",
            task_names[task_idx],
            "_estimate.png"
          ),
          surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
          surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii"
        )
      }
      if(task_idx %in% c(2,3)) {
        for(hem in c("left","right")) {
          plot(
            cifti_obj, hemisphere = hem,
            idx = task_idx,zlim = c(-1,1), legend_embed = F,
            fname = paste0(
              "plots/607_classical_",
              subject,
              "_visit",
              v,
              "_session",
              sess,
              "_",
              hem,
              "_",
              task_names[task_idx],
              "_estimate.png"
            ),
            surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
            surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii"
          )
        }
      }
    }
  }
}

# >> Multirun ----
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
result_files <- list.files(result_dir, full.names = T)
result_files <- c(sapply(subjects,grep, x = result_files, value = T))
vis <- "visit1"
task_idx <- 3
# task_name <- "Tongue"
task_names <- c("visual_cue","foot","hand","tongue")
result_files <- grep(vis, result_files, value = T)
# >>>> Bayesian ----
for(s in subjects) {
  cifti_template <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
  left_result <- readRDS(grep('left',grep(s,result_files, value = T), value = T))
  cifti_template$data$cortex_left <- left_result$betas_Bayesian$avg$data$cortex_left
  rm(left_result)
  right_result <- readRDS(grep('right',grep(s,result_files, value = T), value = T))
  cifti_template$data$cortex_right <- right_result$betas_Bayesian$avg$data$cortex_right
  rm(right_result)
  if(task_idx %in% c(1,4)) {
    plot(cifti_template, idx = task_idx, #title = paste("Subject",s,task_name,"Task"),
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         fname = paste0("plots/600_subject_",s,"_",task_names[task_idx],"_estimates.png"),
         zlim = c(-1,1), legend_embed = F)
  }
  if(task_idx %in% c(2,3)) {
    for(hem in c("left","right")) {
      plot(cifti_template, idx = task_idx, hemisphere = hem,
           surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
           surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
           fname = paste0("plots/600_subject_",s,"_",hem,"_",task_names[task_idx],"_estimates.png"),
           zlim = c(-1,1), legend_embed = F)
    }
  }
}
# >>>> Classical ----
# task_idx <- 1
# task_names <- c("visual_cue","foot","hand","tongue")
for(s in subjects) {
  cifti_template <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
  left_result <- readRDS(grep('left',grep(s,result_files, value = T), value = T))
  cifti_template$data$cortex_left <- left_result$betas_classical$avg$data$cortex_left
  rm(left_result)
  right_result <- readRDS(grep('right',grep(s,result_files, value = T), value = T))
  cifti_template$data$cortex_right <- right_result$betas_classical$avg$data$cortex_right
  rm(right_result)
  if(task_idx %in% c(1,4)) {
    plot(cifti_template, idx = task_idx, #title = paste("Subject",s,task_name,"Task"),
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         fname = paste0("plots/600_subject_",s,"_",task_names[task_idx],"_classical_estimates.png"),
         zlim = c(-1,1), legend_embed = F)
  }
  if(task_idx %in% c(2,3)) {
    for(hem in c("left","right")) {
      plot(cifti_template, idx = task_idx, hemisphere = hem,
           surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
           surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
           fname = paste0("plots/600_subject_",s,"_",hem,"_",task_names[task_idx],"_classical_estimates.png"),
           zlim = c(-1,1), legend_embed = F)
    }
  }
}

# FIGURE 3: ICC quality barplots ----
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
  # xbar <- mean(X)
  # s2 <- sum((X - xbar)^2) / (prod(dim(X)))
  # xbar_i <- apply(X,1,mean)
  sig2_t <- mean(apply(X,2,var))
  sig2_w <- var(Reduce(`-`,split(X,col(X)))) / 2
  sig2_b <- sig2_t - sig2_w
  sig2_b[sig2_b < 0] <- 0
  r <- sig2_b / sig2_t
  # r <- (ncol(X)*sum((xbar_i - xbar)^2) / (nrow(X)*s2) - 1) / (ncol(X) - 1)
  # if(r < 0) r <- 0
  return(r)
}

# Grab the indices of the data locations that are active for Bayes group at 0%
group_left <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples/501_HCP_BayesGLM2_45subj_visit1_result_left_thresh0_20210318.rds")$active
group_right <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples/501_HCP_BayesGLM2_45subj_visit1_result_right_thresh0_20210319.rds")$active
group_active <- list(
  left = group_left,
  right = group_right
)

# >> Bayesian ----
result_dir <- "HCP_results/5k_results"
avg_estimates <- readRDS(file.path(result_dir,"602_avg_estimates.rds"))
library(abind)
abind_g <- function(...) abind(...,along = 4)
ICC_values_average <- sapply(avg_estimates, function(hem_est) {
  combined_visits <- Reduce(abind_g,hem_est)
  vertex_ICC <- apply(combined_visits, 1:2, ICC)
  return(vertex_ICC)
}, simplify = F)

ICC_quality <- sapply(ICC_values_average, function(hem_icc) {
  apply(hem_icc,2,function(task_icc) {
    icc_rating <-
      cut(
        task_icc,
        breaks = c(-Inf, 0.4, 0.6, 0.75, Inf),
        labels = c(1, 2, 3, 4)
      )
    return(as.numeric(icc_rating))
  })
})

# >> Classical ----
avg_estimates_classical <- readRDS(file.path(result_dir,"602_avg_estimates_classical.rds"))
library(abind)
abind_g <- function(...) abind(...,along = 4)
ICC_values_average_classical <- sapply(avg_estimates_classical, function(hem_est) {
  combined_visits <- Reduce(abind_g,hem_est)
  vertex_ICC <- apply(combined_visits, 1:2, ICC)
  return(vertex_ICC)
}, simplify = F)

ICC_quality_classical <- sapply(ICC_values_average_classical, function(hem_icc) {
  icc_out <- apply(hem_icc,2,function(task_icc) {
    icc_rating <-
      cut(
        task_icc,
        breaks = c(-Inf, 0.4, 0.6, 0.75, Inf),
        labels = c(1, 2, 3, 4)
      )
    return(as.numeric(icc_rating))
  })
  colnames(icc_out) <- NULL
  return(as.matrix(icc_out))
})

task_file_names <- c("visual_cue","foot","hand","tongue")
task_names <- c("Visual Cue","Foot","Hand","Tongue")

library(RColorBrewer)
my_pal <- brewer.pal(6, "Accent")[c(3,5,6)]

library(tidyverse)
ICC_df <-
  reshape2::melt(ICC_quality) %>%
  mutate(Model = "Bayesian") %>%
  full_join(
    reshape2::melt(ICC_quality_classical) %>%
      mutate(Model = "Classical")
  ) %>%
  left_join(
    reshape2::melt(group_active, value.name = "active")
  ) %>%
  filter(active == 1) %>%
  select(-active) %>%
  mutate(task = task_names[Var2],
         quality = factor(value, labels = c("Poor","Fair","Good","Excellent")))

ICC_df %>% filter(task == "Visual Cue") %>% group_by(Model) %>%  summarize(mean(value == 4))

ICC_quality_df <-
  ICC_df %>%
  group_by(L1, Model, task, quality) %>%
  tally %>%
  pivot_wider(names_from = L1, values_from = n) %>%
  mutate(both = left+right) %>%
  pivot_longer(cols = c(left,right,both), names_to = "hem") %>%
  filter((task %in% c("Foot", "Hand") & hem %in% c("left","right")) |
           (task %in% c("Visual Cue","Tongue") & hem == "both")) %>%
  mutate(hem = str_to_title(hem),
         task = paste(hem,task),
         task = sub("Both ","", task)) %>%
  select(-hem)

# Facet by task
icc_quality_plot <- group_by(ICC_quality_df,Model, task) %>%
  mutate(total = sum(value),
         prop = value / total) %>%
  filter(quality != "Poor") %>%
  ggplot(aes(x = Model, y = prop, fill = quality)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~task, scales = "free_x") +
  # scale_y_continuous(breaks = c(0,1), labels = c(0,1))+
  scale_fill_manual("ICC Quality",values = my_pal) +
  labs(x = "", y = "Proportion of Data Locations") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_flip()

icc_quality_plot

ggsave("plots/607_icc_quality_plot.png",plot = icc_quality_plot, height = 2.5, width = 10)

# FIGURE 4: Subject Test-Retest MSE and Correlation ----
# >> MSE ----
bayes_result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
classical_result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/classical"

load("HCP_data/subjects.Rdata")

classical_visit2 <- sapply(subjects, function(subject) {
  classical_file_left <- grep(paste0(subject,"_visit2_left"), list.files(classical_result_dir, full.names = T), value = T)
  classical_file_right <- grep(paste0(subject,"_visit2_right"), list.files(classical_result_dir, full.names = T), value = T)
  class_visit2_est <- list(
    left = readRDS(classical_file_left)$betas_classical$avg$data$cortex_left,
    right = readRDS(classical_file_right)$betas_classical$avg$data$cortex_right
  )
  return(class_visit2_est)
}, simplify = FALSE)

# bayes_visit1_est <- list()
# for(subject in subjects) {
#   result_file_left <- grep(paste0(subject,"_visit1_left"), list.files(bayes_result_dir, full.names = T), value = T)
#   result_file_right <- grep(paste0(subject,"_visit1_right"), list.files(bayes_result_dir, full.names = T), value = T)
#   ests <- list(
#     left = readRDS(result_file_left)$betas_Bayesian$avg$data$cortex_left,
#     right = readRDS(result_file_right)$betas_Bayesian$avg$data$cortex_right
#   )
#   bayes_visit1_est[[subject]] <- ests
# }
# saveRDS(bayes_visit1_est, "HCP_results/5k_results/607_bayes_visit1_est.rds")
bayes_visit1 <- readRDS("HCP_results/5k_results/607_bayes_visit1_est.rds")

classical_visit1 <- sapply(subjects, function(subject) {
  classical_file_left <- grep(paste0(subject,"_visit1_left"), list.files(classical_result_dir, full.names = T), value = T)
  classical_file_right <- grep(paste0(subject,"_visit1_right"), list.files(classical_result_dir, full.names = T), value = T)
  class_visit1_est <- list(
    left = readRDS(classical_file_left)$betas_classical$avg$data$cortex_left,
    right = readRDS(classical_file_right)$betas_classical$avg$data$cortex_right
  )
  return(class_visit1_est)
}, simplify = FALSE)

classical_mse <- mapply(function(est,tru) {
  subj_sqerr <- mapply(function(e,t){
    err <- (e - t)^2
    return(err)
  }, e = est, t = tru, SIMPLIFY = F)
  tongue <- mean(c(subj_sqerr[[1]][,1],subj_sqerr[[1]][,1]))
  cue <- mean(c(subj_sqerr[[1]][,4],subj_sqerr[[1]][,4]))
  right_foot <- mean(subj_sqerr[[1]][,2])
  left_foot <- mean(subj_sqerr[[2]][,2])
  right_hand <- mean(subj_sqerr[[1]][,3])
  left_hand <- mean(subj_sqerr[[2]][,3])
  subj_mse_df <- data.frame(
    cue = cue,
    right_foot = right_foot,
    right_hand = right_hand,
    left_foot = left_foot,
    left_hand = left_hand,
    tongue = tongue
  )
  return(subj_mse_df)
}, est = classical_visit1, tru = classical_visit2, SIMPLIFY = F)

bayes_mse <- mapply(function(est,tru) {
  subj_sqerr <- mapply(function(e,t){
    err <- (e - t)^2
    return(err)
  }, e = est, t = tru, SIMPLIFY = F)
  tongue <- mean(c(subj_sqerr[[1]][,1],subj_sqerr[[1]][,1]))
  cue <- mean(c(subj_sqerr[[1]][,4],subj_sqerr[[1]][,4]))
  right_foot <- mean(subj_sqerr[[1]][,2])
  left_foot <- mean(subj_sqerr[[2]][,2])
  right_hand <- mean(subj_sqerr[[1]][,3])
  left_hand <- mean(subj_sqerr[[2]][,3])
  subj_mse_df <- data.frame(
    cue = cue,
    right_foot = right_foot,
    right_hand = right_hand,
    left_foot = left_foot,
    left_hand = left_hand,
    tongue = tongue
  )
  return(subj_mse_df)
}, est = bayes_visit1, tru = classical_visit2, SIMPLIFY = F)


library(tidyverse)
max_vals <- reshape2::melt(bayes_mse) %>%
  mutate(Model = "Bayesian") %>%
  full_join(
    reshape2::melt(classical_mse) %>%
      mutate(Model = "Classical")
  ) %>%
  mutate(variable = sub("_"," ", variable)) %>%
  group_by(variable) %>%
  summarize(max_value = max(value)) %>%
  filter(variable == "tongue")

subj_mse_plot <-
  reshape2::melt(bayes_mse) %>%
  mutate(Model = "Bayesian") %>%
  full_join(
    reshape2::melt(classical_mse) %>%
      mutate(Model = "Classical")
  ) %>%
  mutate(variable = sub("_"," ", variable)) %>%
  pivot_wider(names_from = Model, values_from = value) %>%
  mutate(diff = abs(Bayesian - Classical)) %>%
  group_by(variable) %>%
  mutate(diff = diff / max(diff),
         diff = ifelse(diff > 0.8,0.8,diff)) %>%
  # filter(variable == "tongue") %>%
  ggplot() +
  geom_point(aes(x = max_value, y = max_value), color = "white", data = max_vals) +
  geom_segment(aes(x = Classical, y = Bayesian, xend = Classical, yend = Classical, color = diff)) +
  # scale_color_gradientn("",colors = c('yellow','turquoise2','blue2'), limits = c(0,0.8)) +
  scale_color_gradientn("",
                        colors = c(rgb(242,238,73,maxColorValue = 255),
                                   rgb(79,195,239,maxColorValue = 255),
                                   rgb(71, 109, 174, maxColorValue = 255)),
                        limits = c(0, 0.8))+
  geom_point(aes(y = Bayesian, x = Classical)) +
  geom_abline(intercept = 0, slope = 1, color = 'grey70') +
  labs(
    y = "Bayesian GLM",
    x = "Classical GLM"
  ) +
  facet_wrap(~variable, scales = "free") +
  guides(color = F) +
  theme_classic() +
  theme(text=element_text(size = 14))
  # theme(plot.title = element_text(hjust = 0.5,size = 7),
  #       axis.title = element_text(size = 9))

subj_mse_plot

ggsave("plots/607_subject_mse_plot.png", plot = subj_mse_plot, width = 7.5, height = 5)
# ggsave("~/Desktop/SMI 2021 poster images/607_subject_mse_plot.png", plot = subj_mse_plot, width = 2.5, height = 3)

# >> Correlation ----
bayes_result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
classical_result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/classical"

load("HCP_data/subjects.Rdata")

classical_visit2 <- sapply(subjects, function(subject) {
  classical_file_left <- grep(paste0(subject,"_visit2_left"), list.files(classical_result_dir, full.names = T), value = T)
  classical_file_right <- grep(paste0(subject,"_visit2_right"), list.files(classical_result_dir, full.names = T), value = T)
  class_visit2_est <- list(
    left = readRDS(classical_file_left)$betas_classical$avg$data$cortex_left,
    right = readRDS(classical_file_right)$betas_classical$avg$data$cortex_right
  )
  return(class_visit2_est)
}, simplify = FALSE)

# bayes_visit1_est <- list()
# for(subject in subjects) {
#   result_file_left <- grep(paste0(subject,"_visit1_left"), list.files(bayes_result_dir, full.names = T), value = T)
#   result_file_right <- grep(paste0(subject,"_visit1_right"), list.files(bayes_result_dir, full.names = T), value = T)
#   ests <- list(
#     left = readRDS(result_file_left)$betas_Bayesian$avg$data$cortex_left,
#     right = readRDS(result_file_right)$betas_Bayesian$avg$data$cortex_right
#   )
#   bayes_visit1_est[[subject]] <- ests
# }
# saveRDS(bayes_visit1_est, "HCP_results/5k_results/607_bayes_visit1_est.rds")
bayes_visit1 <- readRDS("HCP_results/5k_results/607_bayes_visit1_est.rds")

classical_visit1 <- sapply(subjects, function(subject) {
  classical_file_left <- grep(paste0(subject,"_visit1_left"), list.files(classical_result_dir, full.names = T), value = T)
  classical_file_right <- grep(paste0(subject,"_visit1_right"), list.files(classical_result_dir, full.names = T), value = T)
  class_visit1_est <- list(
    left = readRDS(classical_file_left)$betas_classical$avg$data$cortex_left,
    right = readRDS(classical_file_right)$betas_classical$avg$data$cortex_right
  )
  return(class_visit1_est)
}, simplify = FALSE)

classical_cor <- mapply(function(est,tru) {
  tongue <- cor(c(est[[1]][,1],est[[2]][,1]),c(tru[[1]][,1],tru[[2]][,1]))
  cue <- cor(c(est[[1]][,4],est[[2]][,4]),c(tru[[1]][,4],tru[[2]][,4]))
  right_foot <- cor(est[[1]][,2], tru[[1]][,2])
  left_foot <- cor(est[[2]][,2], tru[[2]][,2])
  right_hand <- cor(est[[1]][,3], tru[[1]][,3])
  left_hand <- cor(est[[2]][,3], tru[[2]][,3])
  subj_cor_df <- data.frame(
    cue = cue,
    right_foot = right_foot,
    right_hand = right_hand,
    left_foot = left_foot,
    left_hand = left_hand,
    tongue = tongue
  )
  return(subj_cor_df)
}, est = classical_visit1, tru = classical_visit2, SIMPLIFY = F)

bayes_cor <- mapply(function(est,tru) {
  tongue <- cor(c(est[[1]][,1],est[[2]][,1]),c(tru[[1]][,1],tru[[2]][,1]))
  cue <- cor(c(est[[1]][,4],est[[2]][,4]),c(tru[[1]][,4],tru[[2]][,4]))
  right_foot <- cor(est[[1]][,2], tru[[1]][,2])
  left_foot <- cor(est[[2]][,2], tru[[2]][,2])
  right_hand <- cor(est[[1]][,3], tru[[1]][,3])
  left_hand <- cor(est[[2]][,3], tru[[2]][,3])
  subj_cor_df <- data.frame(
    cue = cue,
    right_foot = right_foot,
    right_hand = right_hand,
    left_foot = left_foot,
    left_hand = left_hand,
    tongue = tongue
  )
  return(subj_cor_df)
}, est = bayes_visit1, tru = classical_visit2, SIMPLIFY = F)


library(tidyverse)
max_cors <- reshape2::melt(bayes_cor) %>%
  mutate(Model = "Bayesian") %>%
  full_join(
    reshape2::melt(classical_cor) %>%
      mutate(Model = "Classical")
  ) %>%
  mutate(variable = sub("_"," ", variable)) %>%
  group_by(variable) %>%
  summarize(max_value = max(value))

library(viridis)
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
col_pal <- colorRampPalette(c('red','orange','yellow','green','blue','violet'))

subj_cor_plot <-
  reshape2::melt(bayes_cor) %>%
  mutate(Model = "Bayesian") %>%
  full_join(
    reshape2::melt(classical_cor) %>%
      mutate(Model = "Classical")
  ) %>%
  mutate(variable = sub("_"," ", variable)) %>%
  pivot_wider(names_from = Model, values_from = value) %>%
  mutate(diff = Bayesian - Classical,
         diff = ifelse(diff < 0, 0, diff),
         diff = ifelse(diff > 0.1, 0.1, diff)) %>%
  ggplot() +
  geom_point(aes(x = max_value, y = max_value), color = "white", data = max_cors) +
  geom_segment(aes(x = Classical, y = Bayesian, xend = Classical, yend = Classical, color = diff), alpha = 0.8) +
  # scale_color_gradientn("",colors = c('red','orange','yellow','green','cornflowerblue'), limits = c(0,0.1)) +
  scale_color_gradientn("",colors = c(
    'red',
    'orange',
    rgb(242,238,73,maxColorValue = 255),
    rgb(79,195,239,maxColorValue = 255),
    rgb(71, 109, 174, maxColorValue = 255)
    ), limits = c(0,0.1)) +
  # scale_color_gradientn("",colors = rev(c('red','orange','yellow','green','blue')), limits = c(-.18,.18)) +
  # scale_color_manual("", values = col_pal) +
  # scale_color_viridis(option = "C", limits = c(0,.1)) +
  # scale_color_gradient2("",low = "blue",high = "red", limits = c(-.2,.2)) +
  geom_point(aes(y = Bayesian, x = Classical)) +
  geom_abline(intercept = 0, slope = 1, color = 'grey70') +
  labs(y = "Bayesian GLM", x = "Classical GLM") +
  facet_wrap(~variable, scales = "free") +
  guides(color = FALSE) +
  theme_classic() +
  theme(text = element_text(size = 14))

subj_cor_plot

ggsave("plots/607_subject_cor_plot.png", plot = subj_cor_plot, width = 7.5, height = 5)

# FIGURE 5: Multi-run subject activations ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
# >> Bayesian ----
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/single_subject"
result_files <- list.files(result_dir, full.names = T)
result_files <- grep("501_", result_files, value = T)
result_files <- grep("visit1", result_files, value = T)
# result_files <- grep("_alpha001_", result_files, value = T)
# result_files <- grep("_thr01_", result_files, value = T, invert = T)
load("HCP_data/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
threshs <- paste0("_thr",c(0,0.5,1),".rds")
hems <- c("left","right")
# library(viridisLite)
# col_pal <- rev(plasma(5)[3:5])
task_names <- c("cue","foot","hand","tongue")
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
col_pal <- c(light_orange,"red","purple")
for(subject in subjects) {
  for(task_idx in 1:4) {
    subject_files <- grep(subject, result_files, value =T)
    subject_files_left <- grep("left",subject_files, value = T)
    subject_files_right <- grep("right",subject_files, value = T)
    left_act <- sapply(threshs, function(thr) {
      result <- readRDS(grep(thr,subject_files_left, value = T))
      return(result$active[,task_idx])
    })
    left_act <- apply(left_act,1,sum)
    left_act <- as.matrix(ifelse(left_act == 0, NA, left_act))
    right_act <- sapply(threshs, function(thr) {
      result <- readRDS(grep(thr,subject_files_right, value = T))
      return(result$active[,task_idx])
    })
    right_act <- apply(right_act,1,sum)
    right_act <- as.matrix(ifelse(right_act == 0, NA, right_act))
    cifti_obj <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
    cifti_obj$data$cortex_left <- left_act
    cifti_obj$data$cortex_right <- right_act
    if(task_idx %in% c(1,4)) {
      plot(cifti_obj, color_mode = "qualitative", colors = col_pal,
           # title = paste("Subject",subject, "Tongue Task"),
           surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
           surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
           fname = paste0('plots/607_bayes_',
                          subject,
                          "_visit1_",
                          task_names[task_idx],
                          "_activations.png"))
    }
    if(task_idx %in% c(2,3)) {
      for(hem in c('left','right')) {
        plot(cifti_obj, hemisphere = hem,
             color_mode = "qualitative", colors = col_pal,
             surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
             surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
             fname = paste0('plots/607_bayes_',
                            subject,
                            "_visit1_",
                            hem,
                            "_",
                            task_names[task_idx],
                            "_activations.png"))
      }
    }
  }
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
# task_idx <- 4
library(BayesfMRI)
library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path',"/Applications/workbench")
task_names <- c("cue",'foot','hand','tongue')
# subject <- subjects[1]
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# col_pal <- gg_color_hue(3)[2:3]
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
# col_pal <- c("yellow",light_orange,"red","purple")
col_pal <- c(light_orange,"red","purple")
for(subject in subjects) {
  for(task_idx in 1:4) {
    subject_files <- grep(subject,result_files, value = T)
    left_result <- readRDS(grep('left',subject_files, value=T))
    # left_FDR <-
    #   id_activations.classical(
    #     left_result$GLMs_classical$cortexL,
    #     alpha = 0.01,
    #     correction = "FDR",
    #     field_inds = task_idx
    #   )$active
    left_FDR <- sapply(c(0, 0.5, 1), function(thr) {
      out <-
        id_activations.classical(
          left_result$GLMs_classical$cortexL,
          alpha = 0.01,
          correction = "FDR",
          field_inds = task_idx,
          threshold = thr
        )$active
      return(out)
    }, simplify = T)
    left_FWER <- left_FDR
    # left_FWER <- sapply(c(0, 0.5, 1), function(thr) {
    #   out <-
    #     id_activations.classical(
    #       left_result$GLMs_classical$cortexL,
    #       alpha = 0.01,
    #       correction = "FWER",
    #       field_inds = task_idx,
    #       threshold = thr
    #     )$active
    #   return(out)
    # }, simplify = T)
    left_FWER <- as.matrix(apply(left_FWER,1,sum))
    # left_active <- left_FDR + left_FWER
    left_active <- left_FWER
    medial_left <- !is.na(left_active[,1])
    left_active <- as.matrix(left_active[,1][!is.na(left_active[,1])])
    left_active[,1] <- ifelse(left_active[,1] == 0,NA,left_active[,1])
    right_result <- readRDS(grep('right',subject_files, value=T))
    # right_FDR <-
    #   id_activations.classical(
    #     right_result$GLMs_classical$cortexR,
    #     alpha = 0.01,
    #     correction = "FDR",
    #     field_inds = task_idx
    #   )$active
    right_FDR <- sapply(c(0, 0.5, 1), function(thr) {
      out <-
        id_activations.classical(
          right_result$GLMs_classical$cortexR,
          alpha = 0.01,
          correction = "FDR",
          field_inds = task_idx,
          threshold = thr
        )$active
      return(out)
    }, simplify = T)
    right_FWER <- right_FDR
    # right_FWER <- sapply(c(0, 0.5, 1), function(thr) {
    #   out <-
    #     id_activations.classical(
    #       right_result$GLMs_classical$cortexR,
    #       alpha = 0.01,
    #       correction = "FWER",
    #       field_inds = task_idx,
    #       threshold = thr
    #     )$active
    #   return(out)
    # }, simplify = T)
    right_FWER <- as.matrix(apply(right_FWER,1,sum))
    # right_active <- right_FDR + right_FWER
    right_active <- right_FWER
    medial_right <- !is.na(right_active[,1])
    right_active <- as.matrix(right_active[,1][!is.na(right_active[,1])])
    right_active[,1] <- ifelse(right_active[,1] == 0,NA,right_active[,1])
    cifti_active <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
    cifti_active$data$cortex_left <- left_active
    cifti_active$data$cortex_right <- right_active
    if(task_idx %in% c(1,4)) {
      plot(cifti_active, color_mode = "qualitative", colors = col_pal,
           surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
           surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
           fname = paste0('plots/607_classical_',
                          subject,
                          "_visit1_",
                          task_names[task_idx],
                          "_FDR_activations.png"))
    }
    if(task_idx %in% c(2,3)) {
      for(hem in c('left','right')) {
        plot(cifti_active, hemisphere = hem,
             color_mode = "qualitative", colors = col_pal,
             surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
             surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
             fname = paste0('plots/607_classical_',
                            subject,
                            "_visit1_",
                            hem,
                            "_",
                            task_names[task_idx],
                            "_FDR_activations.png"))
      }
    }
  }
}

# FIGURE 6: Subject test-retest activation overlap ----
# First, make the palette
pale_red <- colorRampPalette(c("white","red"), space = "rgb")(4)[2]
pale_blue <- colorRampPalette(c("white","blue"), space = "rgb")(4)[2]
purple_overlap <- colorRampPalette(c("red","blue"), space = "rgb")(3)[2]
col_pal <- c(pale_red,pale_blue,purple_overlap)

# >> Bayesian ----
library(ciftiTools)
ciftiTools.setOption('wb_path',"/Applications/workbench/")
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/single_subject"
load("HCP_data/subjects.Rdata")
subjects <- subjects[1]
task_names <- c("cue","foot","hand","tongue")
# task_idx <- 4 # This is for the tongue
thr <- c(0.5,1)
# thr_chr <- paste0("_thr",sub("\\.","",thr),"_")
alpha <- 0.01
alpha_chr <- paste0("alpha",sub("\\.","",alpha))
# col_pal is now defined above
# col_pal <- c("pink","red")
# col_pal <- grDevices::colorRampPalette(c("white","red"))(3)[-1]
subject_files <- list.files(result_dir, full.names = T)
subject_files <- grep("501_", subject_files, value = T)
subject_files <- unlist(sapply(subjects, grep, x = subject_files, value = T))
# subject_files <- unlist(sapply(thr_chr, grep, x = subject_files,  value=T))
# subject_files <- grep("_classical.rds",subject_files, value = T, invert = T)
for(task_idx in 2:3) {
  for(subject in subjects) {
    for(thresh in thr) {
      for(h in c('left','right')){
        for(v in 1:2) {
          # subject_file <- grep(paste0("_visit",v,"_subject_",subject,"_", h, thresh), subject_files, value = T)
          subject_file <- grep(paste0("_subject_",subject,"_visit",v, "_", h, "_thr", thresh,".rds"), subject_files, value = T)
          subject_act <- readRDS(subject_file)
          L_or_R <- toupper(substring(h,1,1))
          if(v == 1){
            if(h == 'left') {
              subject_cifti <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
            }
            subject_cifti$data[[paste0("cortex_",h)]] <- as.matrix(subject_act$active[,task_idx]*v)
          }
          if(v == 2) {
            subject_cifti$data[[paste0("cortex_",h)]] <- as.matrix(subject_act$active[,task_idx]*v) +
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
      if(task_idx %in% c(1,4)) {
        needs_colors <- unique(c(subject_cifti$data$cortex_left[,1],subject_cifti$data$cortex_right[,1]))
        needs_colors <- sort(na.omit(needs_colors))
        plot(subject_cifti,
             fname = paste0("plots/607_bayes_",subject,"_",task_names[task_idx],"_thr",sub("\\.","",thresh),"_",alpha_chr,"_num_active.png"),
             # title = paste("Subject",subject, "Tongue Task\nVisit Activations"),
             color_mode = "qualitative",colors = col_pal[needs_colors],
             surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
             surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")
      }
      if(task_idx %in% 2:3) {
        for(hem in c('left','right')) {
          needs_colors <- unique(subject_cifti$data[[paste0("cortex_",hem)]][,1])
          needs_colors <- sort(na.omit(needs_colors))
          plot(subject_cifti, hemisphere = hem,
               fname = paste0("plots/607_bayes_",subject,"_",hem,"_",task_names[task_idx],"_thr",sub("\\.","",thresh),"_",alpha_chr,"_num_active.png"),
               # title = paste("Subject",subject, "Tongue Task\nVisit Activations"),
               color_mode = "qualitative",colors = col_pal[needs_colors],
               surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
               surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")
        }
      }
    }
  }
}


# >> Classical ----
library(ciftiTools)
ciftiTools.setOption('wb_path',"/Applications/workbench")
result_dir <- "HCP_results/5k_results/individual/PW/activations"
load("HCP_data/subjects.Rdata")
subjects <- subjects[c(1,2)]
task_idx <- 4 # This is for the tongue
thr <- c(0,0.5)
thr_chr <- paste0("_thr",sub("\\.","",thr),"_")
alpha <- 0.01
alpha_chr <- paste0("alpha",sub("\\.","",alpha))
# col_pal <- c("pink","red") # col_pal is now defined above
subject_files <- grep("503_", list.files(result_dir, full.names = T), value = T)
subject_files <- unlist(sapply(subjects, grep, x = subject_files, value = T))
subject_files <- unlist(sapply(thr_chr,grep,x = subject_files,  value=T))
subject_files <- grep("_classical.rds",subject_files, value = T, invert = F)
for(subject in subjects) {
  for(thresh in thr_chr) {
    for(h in c('left','right')){
      for(v in 1:2) {
        subject_file <-
          grep(paste0("_visit", v, "_subject_", subject, "_", h,thresh),
               subject_files,
               value = T)
        subject_act <- readRDS(subject_file)
        L_or_R <- toupper(substring(h,1,1))
        if(v == 1){
          if(h == 'left') {
            subject_cifti <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
          }
          in_mask <- subject_cifti$meta$cortex$medial_wall_mask[[h]]
          subject_cifti$data[[paste0("cortex_",h)]] <-
            apply(subject_act$active,2,as.numeric)[in_mask,task_idx]
          # as.matrix(subject_act$active[in_mask,task_idx])
        }
        if(v == 2) {
          subject_cifti$data[[paste0("cortex_", h)]] <-
            apply(subject_act$active,2,as.numeric)[in_mask,task_idx]*v +
            # as.matrix(subject_act$active[in_mask, task_idx]) +
            subject_cifti$data[[paste0("cortex_", h)]]
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
    plot(
      subject_cifti,
      fname = paste0(
        "plots/603_subject_",
        subject,
        thresh,
        alpha_chr,
        "_num_active_classical.png"
      ),
      color_mode = "qualitative",
      colors = col_pal,
      surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
      surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii"
    )
  }
}


# FIGURE 7: Spherical vs midthickness ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
xifti_combine <- function(left,right) {
  # Checks
  if(!is.matrix(left$data$cortex_left))
    stop("The left xifti object has no data for the left cortex.")
  if (!is.matrix(right$data$cortex_right))
    stop("The right xifti object has no data for the right cortex.")
  if(is.null(left$meta$cortex$medial_wall_mask$left) |
     is.null(right$meta$cortex$medial_wall_mask$right))
    stop("One or both of the hemisphere xiftis is missing a medial wall mask.")

  cifti_out <- left
  cifti_out$data$cortex_right <- right$data$cortex_right
  cifti_out$meta$cortex$medial_wall_mask$right <-
    right$meta$cortex$medial_wall_mask$right
  return(cifti_out)
}
sphere_dir <- "HCP_results/5k_results/individual/PW/spherical_analyses"
mid_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
diff_dirs <- c(mid_dir, sphere_dir)
load("HCP_data/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
task_names <- c("cue","foot","hand","tongue")
# >> Estimates ----
zlims <- list(
  c(-1,1),
  c(-1,1),
  NULL
)
for(subject in subjects[1]) {
  for(visit in paste0("visit",1:2)[1]) {
    surf_ciftis <- sapply(1:2, function(surf) {
      dir_files <- list.files(diff_dirs[surf], full.names = T)
      subj_files <- grep(subject, dir_files, value = T)
      for(hem in c('left','right')) {
        filen <- grep(paste0(visit,"_",hem), subj_files, value = T)
        result_obj <- readRDS(filen)
        if(hem == 'left') cifti_obj <- result_obj$betas_Bayesian$avg
        if(hem == 'right') cifti_obj <- xifti_combine(cifti_obj, result_obj$betas_Bayesian$avg)
      }
      return(cifti_obj)
    }, simplify = F)
    names(surf_ciftis) <- c("midthickness","spherical")
    surf_ciftis$diff <- surf_ciftis$midthickness - surf_ciftis$spherical
    for(surf in 1:3) {
      for(task_idx in c(1,4)) {
        plot(surf_ciftis[[surf]], idx = task_idx, zlim = zlims[[surf]], legend_embed = F,
             surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
             surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
             fname = paste0("plots/607_bayes_",c('midthickness','spherical','diff')[surf],"_",subject,"_",visit,"_",task_names[task_idx],"_estimate.png"))
      }
      for(task_idx in 2:3) {
        for(hem in c('left','right')) {
          plot(surf_ciftis[[surf]], idx = task_idx, zlim = zlims[[surf]], hemisphere = hem,
               legend_embed = F,
               surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
               surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
               fname = paste0("plots/607_bayes_",c('midthickness','spherical','diff')[surf],"_",subject,"_",visit,"_",hem,"_",task_names[task_idx],"_estimate.png"))
        }
      }
    }
  }
}

# >> Activations ----
pale_red <- colorRampPalette(c("white","red"), space = "rgb")(4)[2]
pale_blue <- colorRampPalette(c("white","blue"), space = "rgb")(4)[2]
purple_overlap <- colorRampPalette(c("red","blue"), space = "rgb")(3)[2]
col_pal <- c(pale_red,pale_blue,purple_overlap)
library(INLA)
inla.setOption(pardiso.license = "~/licenses/pardiso.lic") # Dan's Macbook Pro
for(subject in subjects[1]) {
  for(visit in paste0("visit",1:2)[1]) {
    job::job({surf_ciftis <- sapply(1:2, function(surf) {
      dir_files <- list.files(diff_dirs[surf], full.names = T)
      subj_files <- grep(subject, dir_files, value = T)
      for(hem in c('left','right')) {
        filen <- grep(paste0(visit,"_",hem), subj_files, value = T)
        result_obj <- readRDS(filen)
        L_or_R <- toupper(substring(hem,1,1))
        exc_obj <-
          BayesGLM2(
            list(result_obj$GLMs_Bayesian[[paste0("cortex", L_or_R)]]),
            excursion_type = ">",
            gamma = 0.5,
            alpha = 0.01,
            no_cores = 4
          )
        if(hem == 'left') {
          cifti_left <- result_obj$betas_Bayesian$avg
          cifti_left$data$cortex_left <- exc_obj$active
        }
        if(hem == 'right') {
          cifti_right <- result_obj$betas_Bayesian$avg
          cifti_right$data$cortex_right <- exc_obj$active
          cifti_both <- xifti_combine(cifti_left, cifti_right)
        }
      }
      return(cifti_both)
    }, simplify = F)},import = "auto")
    names(surf_ciftis) <- c("midthickness","spherical")
    surf_ciftis$diff <- surf_ciftis$midthickness + (2*surf_ciftis$spherical)
    surf_ciftis$diff$data$cortex_left[surf_ciftis$diff$data$cortex_left == 0] <- NA
    surf_ciftis$diff$data$cortex_right[surf_ciftis$diff$data$cortex_right == 0] <- NA
    for(surf in 3) {
      for(task_idx in c(1,4)) {
        use_colors <- sort(unique(c(surf_ciftis[[surf]]$data$cortex_left[,task_idx],
                                    surf_ciftis[[surf]]$data$cortex_right[,task_idx])))
        plot(surf_ciftis[[surf]], idx = task_idx, color_mode = "qualitative",
             legend_embed = F, colors = col_pal[use_colors],
             surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
             surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
             fname = paste0("plots/607_bayes_",c('midthickness','spherical','diff')[surf],"_",subject,"_",visit,"_",task_names[task_idx],"_activation.png"))
      }
      for(task_idx in 2:3) {
        for(hem in c('left','right')) {
          use_colors <-
            sort(unique(surf_ciftis[[surf]]$data[[paste0("cortex_", hem)]][, task_idx]))
          plot(surf_ciftis[[surf]], idx = task_idx, color_mode = "qualitative",
               hemisphere = hem, legend_embed = F, colors = col_pal[use_colors],
               surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
               surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
               fname = paste0("plots/607_bayes_",c('midthickness','spherical','diff')[surf],"_",subject,"_",visit,"_",hem,"_",task_names[task_idx],"_activation.png"))
        }
      }
    }
  }
}

# FIGURE 7: Dice Overlap CIs and Dice x Area scatterplot ----
# >> Dice Overlap CIs ----
bstrap_CIs <- function(x, num_resamp = 5000, alpha = 0.05) {
  bsamp_means <- sapply(seq(num_resamp), function(s){
    bstrap_sample <- sample(x, length(x), replace = T)
    mean_out <- mean(bstrap_sample, na.rm = T)
    return(mean_out)
  })
  lo_mean_hi <-
    c(
      Low = as.numeric(quantile(bsamp_means, probs = alpha / 2)),
      Mean = mean(bsamp_means),
      High = as.numeric(quantile(bsamp_means, probs = 1 - alpha / 2))
    )
  return(lo_mean_hi)
}

dice_result_files <- list.files("HCP_results/5k_results", full.names = T)
dice_result_files <- grep("605_", dice_result_files, value = T)
dice_result_files <- grep("_PW_", dice_result_files, value = T)
dice_result_files <- grep("_Dice_", dice_result_files, value = T)
dice_files <- list(Bayesian = grep("classical",dice_result_files, invert = T, value = T)[c(2,1,3)],
                   Classical = grep("FWER",dice_result_files, value = T)[c(2,1,3)])
dice_results <- sapply(dice_files, function(df) sapply(df,readRDS, simplify = F), simplify = F)
task_names <- c("visual cue","tongue","right foot","right hand","left foot","left hand")
names(dice_results$Bayesian) <- c("0%","0.5%","1%")
names(dice_results$Classical) <- c("0%","0.5%","1%")

bstrap_ci_results <- sapply(dice_results, function(method_results) {
  sapply(method_results, function(thr_results) {
    apply(thr_results,2,bstrap_CIs)
  },simplify = F)
}, simplify = F)

classical_standard <- reshape2::melt(bstrap_ci_results) %>%
  mutate(task = task_names[Var2],
         threshold = factor(L2, levels = c("0%","0.5%","1%"))) %>%
  pivot_wider(names_from = Var1, values_from = value) %>%
  filter(L1 == "Classical", threshold == "0%")

library(tidyverse)
dice_ci_plot <- reshape2::melt(bstrap_ci_results) %>%
  mutate(task = task_names[Var2],
         threshold = factor(L2, levels = c("0%","0.5%","1%"))) %>%
  pivot_wider(names_from = Var1, values_from = value) %>%
  ggplot() +
  geom_vline(aes(xintercept = Mean, color = L1), data = classical_standard, show.legend = FALSE,lty = 2) +
  geom_errorbar(aes(xmin = Low, xmax = High, y = threshold, color = L1), position = position_dodge(width = 0.2), width = 0.2) +
  geom_point(aes(x = Mean, y = threshold, color = L1), position = position_dodge(width = 0.2)) +
  # facet_wrap(vars(task), scales = "fixed", nrow = 3, ncol = 2) +
  facet_grid(task~., scales = "fixed") +
  labs(y = "Activation Threshold", x = "Dice Coefficient") +
  scale_color_discrete("") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "top")

dice_ci_plot

ggsave("plots/607_dice_ci_plot.png", width = 3, height = 7, plot = dice_ci_plot)

# >> Dice x Area Scatterplot ----

area_result_files <- list.files("HCP_results/5k_results", full.names = T)
area_result_files <- grep("605_", area_result_files, value = T)
area_result_files <- grep("_PW_", area_result_files, value = T)
area_result_files <- grep("_overlap_", area_result_files, value = T)
area_files <- list(Bayesian = grep("classical",area_result_files, invert = T, value = T)[c(2,1,3)],
                   Classical = grep("FWER",area_result_files, value = T)[c(2,1,3)])
area_results <- sapply(area_files, function(df) sapply(df,readRDS, simplify = F), simplify = F)
names(area_results$Bayesian) <- c("0%","0.5%","1%")
names(area_results$Classical) <- c("0%","0.5%","1%")

count <- 0
breaks_fun <- function(x) {
  count <<- count %% 3 + 1L
  switch(
    count,
    c(0, 750, 1500),
    c(0, 50, 100),
    c(0, 100, 200)
  )
}

area_by_dice <- reshape2::melt(dice_results, value.name = "dice") %>%
  cbind(overlap = reshape2::melt(area_results, value.name = "overlap")$overlap) %>%
  mutate(Var2 = task_names[Var2],
         Var2 = as.factor(Var2),
         L2 = factor(L2, levels = c("0%","0.5%","1%"))) %>%
  ggplot() +
  geom_point(aes(x = overlap, y = dice, color = L1)) +
  facet_grid(Var2 ~ L2, scales = "free_x") +
  scale_color_discrete("") +
  scale_x_continuous(breaks = breaks_fun, limits = c(0,NA))+
  labs(y = "Dice Coefficient", x = "Size of Overlap (# of vertices)") +
  theme_bw() +
  theme(legend.position = "top", panel.grid = element_blank())

area_by_dice

ggsave("plots/607_area_by_dice.png", width = 4.5, height = 7, plot = area_by_dice)

# >> Poster Image ----
count <- 0
breaks_fun <- function(x) {
  count <<- count %% 3 + 1L
  switch(
    count,
    c(0, 350, 700),
    c(0, 50, 100),
    c(0, 100, 200)
  )
}

area_by_dice <-
  reshape2::melt(dice_results, value.name = "dice") %>%
  cbind(overlap = reshape2::melt(area_results, value.name = "overlap")$overlap) %>%
  mutate(Var2 = task_names[Var2],
         Var2 = as.factor(Var2),
         L2 = factor(L2, levels = c("0%","0.5%","1%")),
         L2_num = as.numeric(sub("%","",as.character(L2)))) %>%
  filter(Var2 == "tongue") %>%
  ggplot() +
  geom_point(aes(x = overlap, y = dice, color = L1)) +
  facet_grid(. ~ L2_num, scales = "free_x", labeller = label_bquote(cols = gamma == .(L2_num)~'%')) +
  scale_color_discrete("") +
  scale_x_continuous(breaks = breaks_fun, limits = c(0,NA))+
  labs(y = "Dice Coefficient", x = "Size of Overlap (# of vertices)") +
  theme_bw() +
  theme(legend.position = "top", panel.grid = element_blank())

abd <- area_by_dice + theme(text = element_text(size = 18))
ggsave("~/Desktop/SMI 2021 poster images/607_area_by_dice.png", plot = abd, width = 5, height = 6)

# FIGURE 8: Group Estimates and Activations ----
# >> Estimates ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench/')
task_file_names <- c("visual_cue","foot","hand","tongue")
# >>>> Bayesian ----
group_results_left <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples/501_HCP_BayesGLM2_45subj_visit1_result_left_thresh0_20210318.rds")
group_results_right <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples/501_HCP_BayesGLM2_45subj_visit1_result_right_thresh0_20210319.rds")
cifti_both <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
cifti_both$data$cortex_left <- as.matrix(group_results_left$estimates)
cifti_both$data$cortex_right <- as.matrix(group_results_right$estimates)
# >>>>>>  Visual Cue and Tongue (idx = c(1,4)) ----
for(i in c(1,4)) {
  plot(cifti_both, zlim = c(-1,1), idx = i, hemisphere = "both",
       legend_embed = F,
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
       fname = paste0("plots/607_group_bayes_",task_file_names[i],"_estimate.png"))
}

# >>>>>> Foot and Hand Tasks (idx = c(2,3)) ----
for(i in c(2,3)) {
  for(hem in c("left","right")) {
    plot(cifti_both, zlim = c(-1,1), idx = i, hemisphere = hem,
         legend_embed = F,
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         fname = paste0("plots/607_group_bayes_",hem,"_",task_file_names[i],"_estimate.png"))
  }
}

# >>>> Classical ----
group_results_classical <- readRDS("HCP_results/5k_results/502_HCP_classical_group_PW_estimates_visit1.rds")
cifti_classical <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
cifti_classical$data$cortex_left <- as.matrix(group_results_classical$left)
cifti_classical$data$cortex_right <- as.matrix(group_results_classical$right)
# >>>>>>  Visual Cue and Tongue (idx = c(1,4)) ----
for(i in c(1,4)) {
  plot(cifti_classical, zlim = c(-1,1), idx = i, hemisphere = "both",
       legend_embed = F,
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
       fname = paste0("plots/607_group_classical_",task_file_names[i],"_estimate.png"))
}

# >>>>>> Foot and Hand Tasks (idx = c(2,3)) ----
for(i in c(2,3)) {
  for(hem in c("left","right")) {
    plot(cifti_classical, zlim = c(-1,1), idx = i, hemisphere = hem,
         legend_embed = F,
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         fname = paste0("plots/607_group_classical_",hem,"_",task_file_names[i],"_estimate.png"))
  }
}

# >> Activations ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
# >>>> Bayesian ----
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples"
# result_dir <- "HCP_results/5k_results/individual/PW"
result_files <- list.files(result_dir, full.names = T)
result_files <- grep("501_HCP_BayesGLM2_45subj_visit1", result_files, value = T)
left_files <- grep("left",result_files, value = T)
right_files <- grep("right",result_files, value = T)
threshs <- paste0("_thresh",c("0","05","1"),"_")
hems <- c("left","right")
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
col_pal <- c(light_orange,"red","purple")
task_names <- c("visual_cue","foot","hand","tongue")
for(task_idx in 1:4) {
  left_act <- sapply(threshs, function(thr) {
    result <- readRDS(grep(thr,left_files, value = T))
    return(result$active[,task_idx])
  })
  left_act <- apply(left_act,1,sum)
  left_act <- as.matrix(ifelse(left_act == 0, NA, left_act))
  right_act <- sapply(threshs, function(thr) {
    result <- readRDS(grep(thr,right_files, value = T))
    return(result$active[,task_idx])
  })
  right_act <- apply(right_act,1,sum)
  right_act <- as.matrix(ifelse(right_act == 0, NA, right_act))
  cifti_obj <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
  cifti_obj$data$cortex_left <- left_act
  cifti_obj$data$cortex_right <- right_act
  if(task_idx %in% c(1,4)) {
    plot(cifti_obj, color_mode = "qualitative", colors = col_pal, hemisphere = "both",
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         fname = paste0('plots/607_group_bayes_',task_names[task_idx],"_activations.png"))
  }
  if(task_idx %in% c(2,3)) {
    for(hem in c("left","right")) {
      plot(cifti_obj, color_mode = "qualitative", colors = col_pal, hemisphere = hem,
           surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
           surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
           fname = paste0('plots/607_group_bayes_',hem,"_",task_names[task_idx],"_activations.png"))
    }
  }
}


## >>>>>> Poster ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
# >>>>>>>> Bayesian ----
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples"
result_files <- list.files(result_dir, full.names = T)
result_files <- grep("501_HCP_BayesGLM2_45subj_visit1", result_files, value = T)
left_files <- grep("left",result_files, value = T)
grep("10subj_visit1", list.files(result_dir), value = T)
threshs <- paste0("_thresh",c("0","05","1"),"_")
hems <- c("left","right")
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
col_pal <- c(light_orange,"red","purple")
task_names <- c("visual_cue","foot","hand","tongue")
task_idx <- 4
left_act <- sapply(threshs, function(thr) {
  result <- readRDS(grep(thr,left_files, value = T))
  return(result$active[,task_idx])
})
left_act <- apply(left_act,1,sum)
left_act <- as.matrix(ifelse(left_act == 0, NA, left_act))
cifti_obj <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
cifti_obj$data$cortex_left <- left_act
plot(cifti_obj, color_mode = "qualitative", colors = col_pal, hemisphere = "left",
     view = "lateral",
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
     fname = paste0('plots/607_group_bayes_',task_names[task_idx],"_activations.png"))


# >>>> Classical ----
# As a note, this would not be necessary in the case of the current (2021-02-10)
# state of the BayesfMRI package, but some retesting using the classical GLM
# is necessary in order to obtain the standard errors of the estimates.
# results_dir <- "HCP_results/5k_results/group"
# result_files <- list.files(results_dir,full.names = T)
# FDR_result <- readRDS(grep("classical_activations_PW_FDR", result_files, value= T))
FWER_result <- readRDS("HCP_results/5k_results/group/502_HCP_classical_activations_PW_visit1_FWER.rds")
library(BayesfMRI)
library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path',"/Applications/workbench")
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
# col_pal <- c("yellow",light_orange,"red","purple")
col_pal <- c(light_orange,"red","purple")

# left_active <- FDR_result$left$`0%`
# left_active <- left_active + Reduce(`+`, FWER_result$left)
left_active <- Reduce(`+`, FWER_result$left)
left_active[left_active == 0] <- NA
# right_active <- FDR_result$right$`0%`
# right_active <- right_active + Reduce(`+`, FWER_result$right)
right_active <- Reduce(`+`, FWER_result$right)
right_active[right_active == 0] <- NA
cifti_classical <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
cifti_classical$data$cortex_left <- left_active
cifti_classical$data$cortex_right <- right_active
task_names <- c("visual_cue","foot","hand","tongue")
for(task_idx in c(1,4)) {
  plot(cifti_classical, hemisphere = "both", idx = task_idx,
       color_mode = "qualitative", colors = col_pal,
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
       fname = paste0('plots/607_group_classical_',task_names[task_idx],"_activations.png"))
}

for(task_idx in c(2,3)) {
  for(hem in c("left","right")) {
    plot(cifti_classical, hemisphere = hem, idx = task_idx,
         color_mode = "qualitative", colors = col_pal,
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         fname = paste0('plots/607_group_classical_',hem,"_",task_names[task_idx],"_activations.png"))
  }
}

# FIGURE 9: Group Test-Retest MSE and Correlation ----
# >> MSE ----
library(tidyverse)
truth <- readRDS("HCP_results/5k_results/502_HCP_classical_group_PW_estimates_visit2.rds")
classical_est <- readRDS("HCP_results/5k_results/502_HCP_classical_group_PW_estimates_visit1.rds")
classical_subsample_est <- readRDS("HCP_results/5k_results/group/502_HCP_classical_subsample_estimates.rds")
bayes_est <- list(
  left = readRDS("HCP_results/5k_results/group/501_HCP_BayesGLM2_45subj_result_left_thresh0_20210317.rds")$estimates,
  right = readRDS("HCP_results/5k_results/group/501_HCP_BayesGLM2_45subj_result_right_thresh0_20210318.rds")$estimates
)
bayes_sub_files <- list.files("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples", full.names = T) %>%
  grep(pattern = "subj_sample",x = ., value = T) %>%
  grep(pattern = "_thresh0_", x = ., value = T)
bayes_subsample_est <- sapply(paste0(c(10,20,30),"subj"), function(sample_size) {
  sapply(paste0("sample",1:10), function(sample_num) {
    sapply(c("left","right"), function(hem) {
      readRDS(grep(paste0(
        "_",sample_size,"_",sample_num,"_visit1_result_",hem
      ),bayes_sub_files,value = T))$estimates
    }, simplify = F)
  }, simplify = F)
}, simplify = F)

mse_bayes <- mapply(function(est_hem,tru_hem){
  mapply(function(est,tru) {
    mean((est - tru)^2)
  }, est = split(est_hem, col(est_hem)), tru = split(tru_hem, col(tru_hem)))
}, est_hem = bayes_est, tru_hem = truth)

mse_classical <- mapply(function(est_hem,tru_hem){
  mapply(function(est,tru) {
    mean((est - tru)^2)
  }, est = split(est_hem, col(est_hem)), tru = split(tru_hem, col(tru_hem)))
}, est_hem = classical_est, tru_hem = truth)

rownames(mse_bayes) <- rownames(mse_classical) <- c("cue","foot","hand","tongue")

mse_bayes
mse_classical

mse_classical_subsample <- sapply(classical_subsample_est, function(size_ests) {
  sapply(size_ests, function(sampl_est) {
    out <- mapply(function(est_hem,tru_hem){
      mapply(function(est,tru) {
        mean((est - tru)^2)
      }, est = split(est_hem, col(est_hem)), tru = split(tru_hem, col(tru_hem)))
    }, est_hem = sampl_est, tru_hem = truth)
    rownames(out) <- c("cue","foot","hand","tongue")
    return(out)
  }, simplify = F)
}, simplify = F)

mse_bayes_subsample <- sapply(bayes_subsample_est, function(size_ests) {
  sapply(size_ests, function(sampl_est) {
    out <- mapply(function(est_hem,tru_hem){
      mapply(function(est,tru) {
        mean((est - tru)^2)
      }, est = split(est_hem, col(est_hem)), tru = split(tru_hem, col(tru_hem)))
    }, est_hem = sampl_est, tru_hem = truth)
    rownames(out) <- c("cue","foot","hand","tongue")
    return(out)
  }, simplify = F)
}, simplify = F)

library(tidyverse)
mse_classical_df <-
  reshape2::melt(mse_classical) %>%
  pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
  mutate(tongue = (left_tongue + right_tongue)/2,
         cue = (left_cue + right_cue)/2) %>%
  select(-ends_with("_tongue"), -ends_with("_cue")) %>%
  pivot_longer(cols = everything() ,names_to = "task", values_to = "MSE") %>%
  mutate(task = sub("_"," ", task),
         model = "Classical",
         num_subjects = "45")

mse_bayes_df <-
  reshape2::melt(mse_bayes) %>%
  pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
  mutate(tongue = (left_tongue + right_tongue)/2,
         cue = (left_cue + right_cue)/2) %>%
  select(-ends_with("_tongue"), -ends_with("_cue")) %>%
  pivot_longer(cols = everything() ,names_to = "task", values_to = "MSE") %>%
  mutate(task = sub("_"," ", task),
         model = "Bayes",
         num_subjects = "45")
classical_mse_subsample_df <- reshape2::melt(mse_classical_subsample) %>%
  mutate(num_subjects = sub("subjects","", L1),
         sample_num = sub("sample","",L2)) %>%
  select(-L1, -L2) %>%
  pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
  mutate(tongue = (left_tongue + right_tongue)/2,
         cue = (left_cue + right_cue)/2) %>%
  select(-ends_with("_tongue"), -ends_with("_cue")) %>%
  pivot_longer(-c(num_subjects,sample_num),names_to = "task", values_to = "MSE") %>%
  mutate(task = sub("_"," ", task),
         model = "Classical")

bayes_mse_subsample_df <- reshape2::melt(mse_bayes_subsample) %>%
  mutate(num_subjects = sub("subj","", L1),
         sample_num = sub("sample","",L2)) %>%
  select(-L1, -L2) %>%
  pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
  mutate(tongue = (left_tongue + right_tongue)/2,
         cue = (left_cue + right_cue)/2) %>%
  select(-ends_with("_tongue"), -ends_with("_cue")) %>%
  pivot_longer(-c(num_subjects,sample_num),names_to = "task", values_to = "MSE") %>%
  mutate(task = sub("_"," ", task),
         model = "Bayes")

median_mse_df <- full_join(classical_mse_subsample_df,bayes_mse_subsample_df) %>%
  group_by(model,task,num_subjects) %>%
  summarize(MSE = median(MSE)) %>%
  full_join(mse_bayes_df) %>%
  full_join(mse_classical_df)

mse_plot <- full_join(classical_mse_subsample_df,bayes_mse_subsample_df) %>%
  ggplot() +
  # geom_boxplot(aes(x = task, y = MSE, fill = num_subjects, color = num_subjects)) +
  geom_jitter(aes(x = num_subjects, y = MSE, color = model), width = 0.1, height = 0, alpha = 0.3) +
  geom_point(aes(x = num_subjects, y = MSE, color = model), data = mse_classical_df, size = 4) +
  geom_point(aes(x = num_subjects, y = MSE, color = model), data = mse_bayes_df, size = 4) +
  geom_line(aes(x = as.numeric(factor(num_subjects)), y = MSE, color = model), data = median_mse_df) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  # geom_segment(aes(x = as.numeric(factor(num_subjects))-0.3,xend = as.numeric(factor(num_subjects))+0.3, y = MSE,yend = MSE, color = model), data = median_mse_df) +
  labs(x = "Subsample Size") +
  scale_color_discrete("") +
  # scale_fill_manual(name = "Subsample Size", values = my_pal) +
  # scale_color_manual(name = "Subsample Size", values = my_pal) +
  facet_wrap(~task, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "top",
        text = element_text(size = 14))

# mse_plot <-
#   reshape2::melt(mse_bayes) %>%
#   mutate(model = "Bayes") %>%
#   full_join(
#     reshape2::melt(mse_classical) %>%
#       mutate(model = "Classical")
#   ) %>%
#   pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
#   mutate(tongue = (left_tongue + right_tongue)/2,
#          cue = (left_cue + right_cue)/2) %>%
#   select(-ends_with("_cue"), -ends_with("_tongue")) %>%
#   pivot_longer(-model,names_to = "task", values_to = "MSE") %>%
#   mutate(task = sub("_"," ", task)) %>%
#   ggplot() +
#   geom_point(aes(x = task, y = MSE, color = model)) +
#   labs(x = "Task", y = "Mean Square Error") +
#   scale_color_discrete("") +
#   theme_classic()

mse_plot
ggsave("plots/607_group_mse.png",plot = mse_plot, width = 6, height = 4.5)

# >>>> Poster image ----
mse_classical_df <-
  reshape2::melt(mse_classical) %>%
  pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
  mutate(tongue = (left_tongue + right_tongue)/2,
         cue = (left_cue + right_cue)/2) %>%
  select(-ends_with("_tongue"), -ends_with("_cue")) %>%
  pivot_longer(cols = everything() ,names_to = "task", values_to = "MSE") %>%
  mutate(task = sub("_"," ", task),
         model = "Classical",
         num_subjects = "45") %>%
  filter(task == "tongue")

mse_bayes_df <-
  reshape2::melt(mse_bayes) %>%
  pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
  mutate(tongue = (left_tongue + right_tongue)/2,
         cue = (left_cue + right_cue)/2) %>%
  select(-ends_with("_tongue"), -ends_with("_cue")) %>%
  pivot_longer(cols = everything() ,names_to = "task", values_to = "MSE") %>%
  mutate(task = sub("_"," ", task),
         model = "Bayes",
         num_subjects = "45") %>%
  filter(task == "tongue")

median_mse_df <- full_join(classical_mse_subsample_df,bayes_mse_subsample_df) %>%
  group_by(model,task,num_subjects) %>%
  summarize(MSE = median(MSE)) %>%
  full_join(mse_bayes_df) %>%
  full_join(mse_classical_df) %>%
  filter(task == "tongue")

mse_plot <-
  full_join(classical_mse_subsample_df,bayes_mse_subsample_df) %>%
  filter(task == "tongue") %>%
  ggplot() +
  # geom_boxplot(aes(x = task, y = MSE, fill = num_subjects, color = num_subjects)) +
  geom_jitter(aes(x = num_subjects, y = MSE, color = model), width = 0.1, height = 0, alpha = 0.3) +
  geom_point(aes(x = num_subjects, y = MSE, color = model), data = mse_classical_df, size = 4) +
  geom_point(aes(x = num_subjects, y = MSE, color = model), data = mse_bayes_df, size = 4) +
  geom_line(aes(x = as.numeric(factor(num_subjects)), y = MSE, color = model), data = median_mse_df) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  # geom_segment(aes(x = as.numeric(factor(num_subjects))-0.3,xend = as.numeric(factor(num_subjects))+0.3, y = MSE,yend = MSE, color = model), data = median_mse_df) +
  labs(x = "Subsample Size") +
  scale_color_discrete("") +
  # scale_fill_manual(name = "Subsample Size", values = my_pal) +
  # scale_color_manual(name = "Subsample Size", values = my_pal) +
  # facet_wrap(~task, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "top",
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9))

mse_plot
ggsave("~/Desktop/SMI 2021 poster images/607_group_mse_plot.png", plot = mse_plot, width = 3, height = 4)

# >> Correlation ----
library(tidyverse)
# Calculate the correlation for the Group estimates
truth <- readRDS("HCP_results/5k_results/502_HCP_classical_group_PW_estimates_visit2.rds")
classical_est <- readRDS("HCP_results/5k_results/502_HCP_classical_group_PW_estimates_visit1.rds")
classical_subsample_est <- readRDS("HCP_results/5k_results/group/502_HCP_classical_subsample_estimates.rds")
bayes_est <- list(
  left = readRDS("HCP_results/5k_results/group/501_HCP_BayesGLM2_45subj_result_left_thresh0_20210317.rds")$estimates,
  right = readRDS("HCP_results/5k_results/group/501_HCP_BayesGLM2_45subj_result_right_thresh0_20210318.rds")$estimates
)

bayes_sub_files <- list.files("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples", full.names = T) %>%
  grep(pattern = "subj_sample",x = ., value = T) %>%
  grep(pattern = "_thresh0_", x = ., value = T)
bayes_subsample_est <- sapply(paste0(c(10,20,30),"subj"), function(sample_size) {
  sapply(paste0("sample",1:10), function(sample_num) {
    sapply(c("left","right"), function(hem) {
      readRDS(grep(paste0(
        "_",sample_size,"_",sample_num,"_visit1_result_",hem
      ),bayes_sub_files,value = T))$estimates
    }, simplify = F)
  }, simplify = F)
}, simplify = F)

cor_bayes <- mapply(function(est_hem,tru_hem){
  mapply(function(est,tru) {
    cor(est,tru)
  }, est = split(est_hem, col(est_hem)), tru = split(tru_hem, col(tru_hem)))
}, est_hem = bayes_est, tru_hem = truth)

cor_classical <- mapply(function(est_hem,tru_hem){
  mapply(function(est,tru) {
    cor(est,tru)
  }, est = split(est_hem, col(est_hem)), tru = split(tru_hem, col(tru_hem)))
}, est_hem = classical_est, tru_hem = truth)

cor_bayes
cor_classical

rownames(cor_bayes) <- rownames(cor_classical) <- c("cue","foot","hand","tongue")

cor_classical_subsample <- sapply(classical_subsample_est, function(size_ests) {
  sapply(size_ests, function(sampl_est) {
    out <- mapply(function(est_hem,tru_hem){
      mapply(function(est,tru) {
        cor(est,tru)
      }, est = split(est_hem, col(est_hem)), tru = split(tru_hem, col(tru_hem)))
    }, est_hem = sampl_est, tru_hem = truth)
    rownames(out) <- c("cue","foot","hand","tongue")
    return(out)
  }, simplify = F)
}, simplify = F)

cor_bayes_subsample <- sapply(bayes_subsample_est, function(size_ests) {
  sapply(size_ests, function(sampl_est) {
    out <- mapply(function(est_hem,tru_hem){
      mapply(function(est,tru) {
        cor(est, tru)
      }, est = split(est_hem, col(est_hem)), tru = split(tru_hem, col(tru_hem)))
    }, est_hem = sampl_est, tru_hem = truth)
    rownames(out) <- c("cue","foot","hand","tongue")
    return(out)
  }, simplify = F)
}, simplify = F)

cor_classical_df <-
  reshape2::melt(cor_classical) %>%
  pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
  mutate(tongue = (left_tongue + right_tongue)/2,
         cue = (left_cue + right_cue)/2) %>%
  select(-ends_with("_tongue"), -ends_with("_cue")) %>%
  pivot_longer(cols = everything() ,names_to = "task", values_to = "corr") %>%
  mutate(task = sub("_"," ", task),
         model = "Classical",
         num_subjects = "45")

cor_bayes_df <-
  reshape2::melt(cor_bayes) %>%
  pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
  mutate(tongue = (left_tongue + right_tongue)/2,
         cue = (left_cue + right_cue)/2) %>%
  select(-ends_with("_tongue"), -ends_with("_cue")) %>%
  pivot_longer(cols = everything() ,names_to = "task", values_to = "corr") %>%
  mutate(task = sub("_"," ", task),
         model = "Bayes",
         num_subjects = "45")
classical_cor_subsample_df <- reshape2::melt(cor_classical_subsample) %>%
  mutate(num_subjects = sub("subjects","", L1),
         sample_num = sub("sample","",L2)) %>%
  select(-L1, -L2) %>%
  pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
  mutate(tongue = (left_tongue + right_tongue)/2,
         cue = (left_cue + right_cue)/2) %>%
  select(-ends_with("_tongue"), -ends_with("_cue")) %>%
  pivot_longer(-c(num_subjects,sample_num),names_to = "task", values_to = "corr") %>%
  mutate(task = sub("_"," ", task),
         model = "Classical")

bayes_cor_subsample_df <- reshape2::melt(cor_bayes_subsample) %>%
  mutate(num_subjects = sub("subj","", L1),
         sample_num = sub("sample","",L2)) %>%
  select(-L1, -L2) %>%
  pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
  mutate(tongue = (left_tongue + right_tongue)/2,
         cue = (left_cue + right_cue)/2) %>%
  select(-ends_with("_tongue"), -ends_with("_cue")) %>%
  pivot_longer(-c(num_subjects,sample_num),names_to = "task", values_to = "corr") %>%
  mutate(task = sub("_"," ", task),
         model = "Bayes")

median_cor_df <- full_join(classical_cor_subsample_df,bayes_cor_subsample_df) %>%
  group_by(model,task,num_subjects) %>%
  summarize(corr = median(corr)) %>%
  full_join(cor_bayes_df) %>%
  full_join(cor_classical_df)

cor_plot <- full_join(classical_cor_subsample_df,bayes_cor_subsample_df) %>%
  ggplot() +
  # geom_boxplot(aes(x = task, y = MSE, fill = num_subjects, color = num_subjects)) +
  geom_jitter(aes(x = num_subjects, y = corr, color = model), width = 0.1, height = 0, alpha = 0.3) +
  geom_point(aes(x = num_subjects, y = corr, color = model), data = cor_classical_df, size = 4) +
  geom_point(aes(x = num_subjects, y = corr, color = model), data = cor_bayes_df, size = 4) +
  geom_line(aes(x = as.numeric(factor(num_subjects)), y = corr, color = model), data = median_cor_df) +
  geom_hline(aes(yintercept = 1), lty = 2) +
  # geom_segment(aes(x = as.numeric(factor(num_subjects))-0.3,xend = as.numeric(factor(num_subjects))+0.3, y = MSE,yend = MSE, color = model), data = median_mse_df) +
  labs(x = "Subsample Size", y = "Correlation") +
  scale_color_discrete("") +
  # scale_fill_manual(name = "Subsample Size", values = my_pal) +
  # scale_color_manual(name = "Subsample Size", values = my_pal) +
  facet_wrap(~task, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "top",
        text = element_text(size = 14))
cor_plot

# library(tidyverse)
# cor_plot <-
#   reshape2::melt(cor_bayes) %>%
#   mutate(model = "Bayes") %>%
#   full_join(
#     reshape2::melt(cor_classical) %>%
#       mutate(model = "Classical")
#   ) %>%
#   pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
#   mutate(tongue = (left_tongue + right_tongue)/2,
#          cue = (left_cue + right_cue)/2) %>%
#   select(-ends_with("_cue"), -ends_with("_tongue")) %>%
#   pivot_longer(-model,names_to = "task", values_to = "cor") %>%
#   mutate(task = sub("_"," ", task)) %>%
#   ggplot() +
#   geom_point(aes(x = task, y = cor, color = model)) +
#   labs(x = "Task", y = "Correlation") +
#   scale_color_discrete("") +
#   theme_classic()
#
# cor_plot
ggsave("plots/607_group_cor.png",plot = cor_plot, width = 6, height = 4.5)

#  FIGURES 10 & 11: Number of subjects x Number of activations (Group) ----
# >> Bayesian ----
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples"
result_files <- list.files(result_dir, full.names = T)
result_files <- grep("_visit1_", result_files, value = T)
result_files <- grep("_single_session", result_files, value = T, invert = T)
# result_files <- grep("_sample", result_files, value = T)

bayes_num_act <-
  sapply(as.character(c(10,20,30)), function(ssize) {
    sapply(as.character(1:10), function(samp_num) {
      sapply(c('left','right'), function(hem) {
        sapply(as.character(c(0,0.5,1)), function(thr) {
          thr_chr <- sub("\\.","",thr)
          filen <- grep(paste0(ssize,"subj_sample",samp_num,"_"), result_files, value = T)
          filen <- grep(hem, filen, value = T)
          filen <- grep(paste0("thresh",thr_chr,"_"), filen, value = T)
          result_obj <- readRDS(filen)
          return(as.matrix(apply(result_obj$active,2,sum)))
        }, simplify = F)
      }, simplify = F)
    }, simplify = F)
  }, simplify = F)

bayes_num_act_45 <-
  sapply(as.character(45), function(ssize) {
    sapply(as.character(1), function(samp_num) {
      sapply(c('left','right'), function(hem) {
        sapply(as.character(c(0,0.5,1)), function(thr) {
          thr_chr <- sub("\\.","",thr)
          filen <- grep(paste0(ssize,"subj_"), result_files, value = T)
          filen <- grep(hem, filen, value = T)
          filen <- grep(paste0("thresh",thr_chr,"_"), filen, value = T)
          result_obj <- readRDS(filen)
          return(as.matrix(apply(result_obj$active,2,sum)))
        }, simplify = F)
      }, simplify = F)
    }, simplify = F)
  }, simplify = F)

# >> Classical ----
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/classical"
library(abind)
result_files <- list.files(result_dir, full.names = T)
result_files <- grep("_visit1_", result_files, value = T)
load("HCP_data/subjects.Rdata")
classical_estimates <- sapply(c('left','right'), function(hem) {
  sapply(subjects, function(subject) {
    filen <- grep(paste0(subject, "_visit1_",hem), result_files, value = T)
    result_obj <- readRDS(filen)
    out <-
      abind(result_obj$betas_classical$LR$data[[paste0("cortex_", hem)]],
            result_obj$betas_classical$RL$data[[paste0("cortex_", hem)]],
            along = 3)
    return(out)
  }, simplify = "array")
}, simplify = F)

subsample_subjects <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples/501_sample_subjects_20210319.rds")
load("HCP_data/subjects.Rdata")
subsample_idx <- sapply(subsample_subjects, function(sam){
  out <- apply(sam,1:2, function(x) which(subjects == x))
  return(out)
}, simplify = F)
threshs <- c(0,0.5,1)
class_num_act <- sapply(subsample_idx, function(ssize) {
  sample_num_active <- sapply(1:10,function(samp) {
    sapply(c('left','right'), function(hem) {
      samp_est <- classical_estimates[[hem]][,,,ssize[,samp]]
      num_locs <- dim(samp_est)[1]
      bonferroni_cutoff <- 0.01 / num_locs # alpha = 0.01
      data_df <- reshape2::melt(samp_est)
      vertex_lists <- split(data_df,data_df$Var1)
      thresh_active <- sapply(c(0,0.5,1), function(thr) {
        active_all <- sapply(vertex_lists, function(v_df){
          vt_df <- split(v_df, v_df$Var2)
          active_vertex <- sapply(vt_df, function(vt) {
            ttest_res <- t.test(x = vt$value, mu = thr, alternative = "greater")
            active_v <- ifelse(ttest_res$p.value < bonferroni_cutoff, 1, 0)
            return(active_v)
          }, simplify = "array")
          return(active_vertex)
        }, simplify = 'array')
        out <- apply(active_all,1,sum)
        return(as.matrix(out))
      }, simplify = F)
      names(thresh_active) <- paste0(threshs,"%")
      return(thresh_active)
    }, simplify = F)
  }, simplify = F)
  names(sample_num_active) <- as.character(1:10)
  return(sample_num_active)
}, simplify = F)
names(class_num_act) <- c("10","20","30")

class_num_act_45 <- sapply("45", function(ssize) {
  sample_num_active <- sapply(1,function(samp) {
    sapply(c('left','right'), function(hem) {
      samp_est <- classical_estimates[[hem]]
      num_locs <- dim(samp_est)[1]
      bonferroni_cutoff <- 0.01 / num_locs # alpha = 0.01
      data_df <- reshape2::melt(samp_est)
      vertex_lists <- split(data_df,data_df$Var1)
      thresh_active <- sapply(c(0,0.5,1), function(thr) {
        active_all <- sapply(vertex_lists, function(v_df){
          vt_df <- split(v_df, v_df$Var2)
          active_vertex <- sapply(vt_df, function(vt) {
            ttest_res <- t.test(x = vt$value, mu = thr, alternative = "greater")
            active_v <- ifelse(ttest_res$p.value < bonferroni_cutoff, 1, 0)
            return(active_v)
          }, simplify = "array")
          return(active_vertex)
        }, simplify = 'array')
        out <- apply(active_all,1,sum)
        return(as.matrix(out))
      }, simplify = F)
      names(thresh_active) <- paste0(threshs,"%")
      return(thresh_active)
    }, simplify = F)
  }, simplify = F)
  names(sample_num_active) <- "1"
  return(sample_num_active)
}, simplify = F)

library(tidyverse)
classical_num_act <- reshape2::melt(class_num_act) %>%
  full_join(reshape2::melt(class_num_act_45)) %>%
  mutate(L4 = sub("%","", L4),
         model = "Classical GLM")

# >> Plot ----
col_pal <- scales::hue_pal()(5)
library(viridis)
col_pal <- viridis(7,option = "C")[1:5]
col_pal <- rgb(c(230, 86, 0, 0,204),c(159,180,158,114,121),c(0,233,115,178,167), maxColorValue = 255)
library(tidyverse)
task_names <- c("cue","foot","hand","tongue")
num_act_plots <- sapply(c(0,0.5), function(thr) {
  # sapply(c("One Hemisphere","Two Hemispheres"), function(hemi) {
  #   if(hemi == "One Hemisphere") {
  #     my_colors <- col_pal[1:4]
  #   } else {
  #     my_colors <- col_pal[5]
  #   }
  # sapply(c("Bayesian","classical"), function(mdl) {
  num_act_df <-
    reshape2::melt(bayes_num_act) %>%
    full_join(reshape2::melt(bayes_num_act_45)) %>%
    mutate(model = "Bayesian GLM") %>%
    full_join(classical_num_act) %>%
    filter(L4 == as.character(thr)) %>%
    mutate(task = task_names[Var1],
           not_hem = ifelse(L3 == 'left', 'right','left'),
           task = paste(not_hem,task)) %>%
    select(-Var1,-Var2,-L3, -not_hem) %>%
    pivot_wider(names_from = task, values_from = value) %>%
    mutate(cue = `left cue` + `right cue`,
           tongue = `left tongue` + `right tongue`) %>%
    select(-ends_with(" cue"), -ends_with(" tongue")) %>%
    mutate(threshold = as.numeric(L4),
           num_subjects = as.numeric(L1)) %>%
    select(-L1,-L2,-L4) %>%
    pivot_longer(cols = -c(threshold,num_subjects,model), names_to = "task", values_to = "num_active") %>%
    mutate(both_hems = ifelse(task %in% c("cue","tongue"),"Two Hemispheres","One Hemisphere")) %>%
    filter(#both_hems == hemi,
      task != "cue",
      # model == mdl
    )
  sub_df <- filter(num_act_df, num_subjects < 45)
  means_df <- group_by(num_act_df,num_subjects, task, threshold, model) %>%
    summarize(num_active = mean(num_active))
  num_act_plot <- ggplot(sub_df, aes(x = num_subjects, y = num_active, color = task)) +
    geom_jitter(width = 1, height = 0, alpha = 0.3) +
    geom_line(data = means_df,lwd = 1.6) +
    geom_point(data = filter(means_df,num_subjects == 45), size = 4) +
    # facet_grid(threshold~model, scales = "free", labeller = label_bquote(rows = gamma == .(threshold) ~'%')) +
    facet_wrap(~model, scales = "free") +
    # labs(x = "Number of Subjects", y = "Number of Active Locations", title = expression(paste0("Activation Threshold (",gamma,"=",thr,"%)"))) +
    labs(x = "Number of Subjects", y = "Number of Active Locations") +
    scale_color_manual("", values = col_pal) +
    # geom_hline(yintercept = 0) +
    theme_classic() +
    theme(legend.position = "right",
          text = element_text(size = 14))
  # num_act_plot
  return(num_act_plot)
  # },simplify = F)
}, simplify = F)

num_act_plots[[2]]
num_act_plots[[1]]

library(ggpubr)
ggarrange(plotlist = num_act_plots, ncol = 2,nrow = 1)
ggarrange(plotlist = num_act_plots[[2]], ncol = 2,nrow = 1)
ggarrange(plotlist = num_act_plots[[1]], ncol = 2,nrow = 1)

ggsave("plots/607_num_act_plots_thr0.png", plot = num_act_plots[[1]], width = 7, height = 6)
ggsave("plots/607_num_act_plots_thr0.5.png", plot = num_act_plots[[2]], width = 7, height = 6)

# FIGURE 12: Subject vs. group number of activations ----
library(tidyverse)
# >> Bayesian ----
# >>>> Group ----
bayes_group_files <- list.files(
  "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples",
  full.names = T) %>%
  grep(pattern = "501_",x = .,value = T) %>%
  grep(pattern = "_visit1",x = .,value = T, invert = F) %>%
  grep(pattern = "_45subj_", x = ., value = T)

bayes_group_num_active <- sapply(paste0("_thresh",c("0","05","1"),"_"), function(thr_chr) {
  sapply(c("left","right"), function(hem) {
    obj <- readRDS(grep(paste0(hem,thr_chr), bayes_group_files, value = T))
    return(as.matrix(apply(obj$active,2,sum))) # Have to do as.matrix for reshape2::melt
  }, simplify = F)
}, simplify = F)

names(bayes_group_num_active) <- c(0,0.5,1)

reshape2::melt(bayes_group_num_active)

# >>>> Subject ----
load("HCP_data/subjects.Rdata")
bayes_subject_files <- list.files(
  "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/single_subject",
  full.names = T) %>%
  grep(pattern = "501_", x = ., value = T)
# bayes_subject_files <- list.files(
#   "HCP_results/5k_results/individual/PW/activations",
#   full.names = T) %>%
#   grep(pattern = "503_",x = ., value = T) %>%
#   grep(pattern = "classical",x = ., value = T, invert = T)

# subjects <- subjects[-seq(3)]
bayes_subject_num_active <- sapply(subjects, function(subject) {
  sapply(paste0("visit",1), function(vis) {
    sapply(c('left','right'), function(hem) {
      # L_or_R <- toupper(substring(hem,1,1))
      all_thr <- sapply(c(0,0.5,1), function(thr) {
        # thr_chr <- paste0("thr",sub("\\.","", as.character(thr)),"_")
        # result_file <- grep(
        #   paste0(vis,"_subject_",subject,"_",hem,"_",thr_chr),
        #   bayes_subject_files, value = T
        # )
        result_file <- grep(
          paste0(subject,"_",vis,"_",hem,"_thr",thr,".rds"),
          bayes_subject_files, value = T
        )
        result_obj <- readRDS(result_file)
        out <- apply(result_obj$active,2,sum)
        return(as.matrix(out))
      }, simplify = F)
      names(all_thr) <- c(0,0.5,1)
      return(all_thr)
    }, simplify = F)
  }, simplify = F)
}, simplify = F)

reshape2::melt(bayes_subject_num_active) %>% str

# >> Classical ----
# >>>> Group ----
classical_group_num_active <- readRDS("HCP_results/5k_results/group/502_HCP_classical_activations_PW_FWER.rds")
classical_group_num_active <- sapply(classical_group_num_active, function(hem_act) {
  out <- sapply(hem_act, function(thr_act) {
    as.matrix(apply(thr_act,2,sum)) # Have to do as.matrix for reshape2::melt
  }, simplify = F)
  names(out) <- c(0,0.5,1)
  return(out)
}, simplify = F)

reshape2::melt(classical_group_num_active)

# >>>> Subject ----
load("HCP_data/subjects.Rdata")
classical_subject_files <- list.files(
  "HCP_results/5k_results/individual/PW/activations",
  full.names = T) %>%
  grep(pattern = "503_",x = ., value = T) %>%
  grep(pattern = "classical",x = ., value = T, invert = F)

classical_subject_num_active <- sapply(subjects, function(subject) {
  sapply(paste0("visit",1), function(vis) {
    sapply(c('left','right'), function(hem) {
      L_or_R <- toupper(substring(hem,1,1))
      all_thr <- sapply(c(0,0.5,1), function(thr) {
        thr_chr <- paste0("thr",sub("\\.","", as.character(thr)),"_")
        result_file <- grep(
          paste0(vis,"_subject_",subject,"_",hem,"_",thr_chr),
          classical_subject_files, value = T
        )
        result_obj <- readRDS(result_file)
        out <- apply(result_obj$active,2,function(x) sum(as.numeric(x)))
        return(as.matrix(out))
      }, simplify = F)
      names(all_thr) <- c(0,0.5,1)
      return(all_thr)
    }, simplify = F)
  }, simplify = F)
}, simplify = F)

reshape2::melt(classical_subject_num_active) %>% str

# >> Plot ----
task_names <- c("cue","foot","hand","tongue")
subject_activations_df <- reshape2::melt(bayes_subject_num_active) %>%
  mutate(model = "Bayes") %>%
  full_join(
    reshape2::melt(classical_subject_num_active) %>%
      mutate(model = "Classical")
  ) %>%
  mutate(task = task_names[Var1]) %>%
  select(-starts_with("Var")) %>%
  rename(thresh = L4,
         hem = L3,
         visit = L2,
         subject = L1) %>%
  pivot_wider(names_from = c(hem,task), values_from = value) %>%
  mutate(cue = left_cue + right_cue,
         tongue = left_tongue + right_tongue) %>%
  select(-ends_with("_cue"), -ends_with("_tongue")) %>%
  pivot_longer(
    cols = -c(thresh, visit, subject, model),
    names_to = "task",
    values_to = "activations"
  ) %>%
  mutate(task = sub("_"," ",task),
         thresh = as.numeric((thresh)))

group_activations_df <-
  reshape2::melt(bayes_group_num_active) %>%
  mutate(model = "Bayes") %>%
  rename(thresh = L1,
         hem = L2) %>%
  full_join(
    reshape2::melt(classical_group_num_active) %>%
      mutate(model = "Classical") %>%
      rename(thresh = L2,
             hem = L1)
  ) %>%
  mutate(task = task_names[Var1]) %>%
  select(-starts_with("Var")) %>%
  pivot_wider(names_from = c(hem,task), values_from = value) %>%
  mutate(cue = left_cue + right_cue,
         tongue = left_tongue + right_tongue) %>%
  select(-ends_with("_cue"), -ends_with("_tongue")) %>%
  pivot_longer(cols = -c(thresh,model), names_to = "task", values_to = "activations") %>%
  mutate(task = sub("_"," ", task),
         thresh = as.numeric(thresh)) %>%
  filter(task != "cue")

my_shapes <- c("B","C")

num_activations_plots <- map(c(0,0.5,1), function(thr) {
  group_df <- filter(group_activations_df, thresh ==thr)
  out <- subject_activations_df %>%
    filter(thresh == thr,
           task != "cue") %>%
    ggplot() +
    # geom_boxplot(aes(x = model, y = activations + 1, color = model)) +
    geom_jitter(aes(x = task, y = activations, color = model), width = 0.1, height = 0, alpha = 0.3) +
    geom_point(aes(x = task, y = activations, shape = model), data = group_df, size = 4, show.legend = FALSE) +
    scale_shape_manual("",values = c("B","C")) +
    # facet_wrap(~task) +
    facet_grid(~thresh, scales = "free", labeller = label_bquote(cols = gamma == .(thresh)~'%')) +
    geom_hline(yintercept = 0) +
    labs(x = "",y = "Number of Active Locations") +
    # scale_y_log10() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5), legend.position = "none")
  out
  return(out)
})

group_bayes_df <- filter(group_activations_df,
                         model == "Bayes" & thresh == 0.5)
group_classical_df <- filter(group_activations_df,
                             model == "Classical" & thresh == 0)
group_df <- full_join(group_bayes_df, group_classical_df) %>%
  mutate(model = ifelse(model == "Bayes", "Bayesian GLM","Classical GLM"),
         model_thresh = paste0(model," (",thresh,"%)"))
subj_bayes_df <- filter(subject_activations_df,
                        model == "Bayes", thresh == 0.5)
subj_classical_df <- filter(subject_activations_df,
                            model == "Classical", thresh == 0)
subj_df <- full_join(subj_bayes_df, subj_classical_df)

col_pal <- rgb(c(230, 86, 0, 0,204),c(159,180,158,114,121),c(0,233,115,178,167), maxColorValue = 255)

compare_num_plot <-
  subj_df %>%
  filter(task != "cue") %>%
  mutate(model = ifelse(model == "Bayes", "Bayesian GLM","Classical GLM"),
         model_thresh = paste0(model," (",thresh,"%)")) %>%
  ggplot() +
  # geom_boxplot(aes(x = model, y = activations + 1, color = model)) +
  geom_boxplot(aes(x = task, y = activations), outlier.alpha = 0, width = 0.4) +
  geom_jitter(aes(x = task, y = activations, color = task), width = 0.2, height = 0, alpha = 0.3) +
  geom_point(aes(x = task, y = activations, color = task), shape = 19, data = group_df, size = 4, show.legend = FALSE) +
  geom_point(aes(x = task, y = activations), shape = 1, data = group_df, size = 4, show.legend = FALSE) +
  scale_color_manual("", values = col_pal) +
  # scale_shape_manual("",values = c("B","C")) +
  # facet_wrap(~task) +
  # facet_grid(~thresh, scales = "free", labeller = label_bquote(cols = gamma == .(thresh)~'%')) +
  facet_grid(model_thresh~., scales = "free") +
  geom_hline(yintercept = 0) +
  labs(x = "",y = "Number of Active Locations") +
  # scale_y_log10() +
  theme_classic() +
  # theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5), legend.position = "none")
  theme(legend.position = "none", text = element_text(size = 10))

compare_num_plot

ggsave("plots/607_compare_num_plot.png", width = 5, height = 4)

# library(gridExtra)
# marrangeGrob(grobs = num_activations_plots, nrow = 1, ncol = 3, top ="")

library(ggpubr)
num_activations_plot <- ggarrange(plotlist = num_activations_plots, nrow = 1, ncol = 3)
num_activations_plot
ggsave("plots/607_num_activations_plot.png", width = 10, height = 7)

# FIGURE 13: Full-res classical ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
xifti_combine <- function(left,right) {
  # Checks
  if(!is.matrix(left$data$cortex_left))
    stop("The left xifti object has no data for the left cortex.")
  if (!is.matrix(right$data$cortex_right))
    stop("The right xifti object has no data for the right cortex.")
  if(is.null(left$meta$cortex$medial_wall_mask$left) |
     is.null(right$meta$cortex$medial_wall_mask$right))
    stop("One or both of the hemisphere xiftis is missing a medial wall mask.")

  cifti_out <- left
  cifti_out$data$cortex_right <- right$data$cortex_right
  cifti_out$meta$cortex$medial_wall_mask$right <-
    right$meta$cortex$medial_wall_mask$right
  return(cifti_out)
}
# >> Unsmoothed ----
result_dir <- "HCP_results/32k_results"
result_files <- list.files(result_dir,full.names = T)
load("HCP_data/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
task_names <- c("cue","foot","hand","tongue")
for(subject in subjects) {
  for(hem in c('left','right')) {
    filen <- grep(paste0(subject,"_visit1_",hem), result_files, value = T)
    result_object <- readRDS(filen)
    if(hem == 'left') cifti_obj <- result_object$betas_classical$avg
    if(hem == 'right') cifti_obj <- xifti_combine(cifti_obj,
                                                  result_object$betas_classical$avg)
  }
  for(task_idx in c(1,4)) {
    plot(cifti_obj, idx = task_idx, zlim = c(-1,1), legend_embed = F,
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         fname = paste0("plots/607_classical_",
                        subject,"_",task_names[task_idx],
                        "_32k_unsmoothed_estimate.png"))
  }
  for(task_idx in 2:3) {
    for(hem in c('left','right')) {
      plot(cifti_obj, idx = task_idx, zlim = c(-1,1), legend_embed = F,
           hemisphere = hem,
           surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
           surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
           fname = paste0("plots/607_classical_",
                          subject,"_",hem,"_",task_names[task_idx],
                          "_32k_unsmoothed_estimate.png"))
    }
  }
}

# >> Smoothing the estimates ----
result_dir <- "HCP_results/32k_results"
result_files <- list.files(result_dir,full.names = T)
load("HCP_data/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
task_names <- c("cue","foot","hand","tongue")
for(subject in subjects) {
  for(hem in c('left','right')) {
    filen <- grep(paste0(subject,"_visit1_",hem), result_files, value = T)
    result_object <- readRDS(filen)
    if(hem == 'left') cifti_obj <- result_object$betas_classical$avg
    if(hem == 'right') cifti_obj <- xifti_combine(cifti_obj,
                                                  result_object$betas_classical$avg)
  }
  cifti_obj_smoothed <- smooth_cifti(
    x = cifti_obj,
    cifti_target_fname = paste0("HCP_results/32k_results/607_",
                                subject,"_visit1_32k_smoothed_estimates.dscalar.nii"),
    surf_FWHM = 6,
    surfL_fname = paste0("HCP_data/subject_surfaces/",subject,
                         ".L.midthickness.32k_fs_LR.surf.gii"),
    surfR_fname = paste0("HCP_data/subject_surfaces/",subject,
                         ".R.midthickness.32k_fs_LR.surf.gii")
  )
  cifti_obj <- read_cifti(paste0("HCP_results/32k_results/607_", subject,
                              "_visit1_32k_smoothed_estimates.dscalar.nii"))
  for(task_idx in c(1,4)) {
    plot(cifti_obj, idx = task_idx, zlim = c(-1,1), legend_embed = F,
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         fname = paste0("plots/607_classical_",
                        subject,"_",task_names[task_idx],
                        "_32k_smoothed_estimate.png"))
  }
  for(task_idx in 2:3) {
    for(hem in c('left','right')) {
      plot(cifti_obj, idx = task_idx, zlim = c(-1,1), legend_embed = F,
           hemisphere = hem,
           surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
           surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
           fname = paste0("plots/607_classical_",
                          subject,"_",hem,"_",task_names[task_idx],
                          "_32k_smoothed_estimate.png"))
    }
  }
}

# >> Estimates from smoothed data ----
result_dir <- "HCP_results/32k_results/smoothed"
result_files <- list.files(result_dir,full.names = T)
load("HCP_data/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
task_names <- c("cue","foot","hand","tongue")
for(subject in subjects) {
  for(hem in c('left','right')) {
    filen <- grep(paste0(subject,"_visit1_",hem), result_files, value = T)
    result_object <- readRDS(filen)
    if(hem == 'left') cifti_obj <- result_object$betas_classical$avg
    if(hem == 'right') cifti_obj <- xifti_combine(cifti_obj,
                                                  result_object$betas_classical$avg)
  }
  for(task_idx in c(1,4)) {
    plot(cifti_obj, idx = task_idx, zlim = c(-1,1), legend_embed = F,
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         fname = paste0("plots/607_smoothed_classical_",
                        subject,"_",task_names[task_idx],
                        "_32k_estimate.png"))
  }
  for(task_idx in 2:3) {
    for(hem in c('left','right')) {
      plot(cifti_obj, idx = task_idx, zlim = c(-1,1), legend_embed = F,
           hemisphere = hem,
           surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
           surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
           fname = paste0("plots/607_smoothed_classical_",
                          subject,"_",hem,"_",task_names[task_idx],
                          "_32k_estimate.png"))
    }
  }
}

# FIGURE 14: Subject test-retest metrics comparison to full-res classical analyses -----

bayes_result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
classical_result_dir_unsmoothed <- "HCP_results/32k_results"
classical_result_dir_smoothed <- "HCP_results/32k_results/smoothed"
classical_result_dir_5k <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/classical"
sphere_5k_result_dir <- "HCP_results/5k_results/individual/PW/spherical_analyses"

load("HCP_data/subjects.Rdata")

classical_32k_unsmoothed_files <- list.files(classical_result_dir_unsmoothed, full.names = T)

keep_subjects <- sapply(subjects, function(subject) {
  return(length(grep(paste0("500_",subject), classical_32k_unsmoothed_files)) == 4)
})
subjects <- subjects[keep_subjects]

classical_visit2_32k <- sapply(subjects, function(subject) {
  classical_file_left <- grep(paste0(subject,"_visit2_left"), list.files(classical_result_dir_unsmoothed, full.names = T), value = T)
  classical_file_right <- grep(paste0(subject,"_visit2_right"), list.files(classical_result_dir_unsmoothed, full.names = T), value = T)
  class_visit2_est <- list(
    left = readRDS(classical_file_left)$betas_classical$avg$data$cortex_left,
    right = readRDS(classical_file_right)$betas_classical$avg$data$cortex_right
  )
  return(class_visit2_est)
}, simplify = FALSE)

classical_visit2_5k <- sapply(subjects, function(subject) {
  classical_file_left <- grep(paste0(subject,"_visit2_left"), list.files(classical_result_dir_5k, full.names = T), value = T)
  classical_file_right <- grep(paste0(subject,"_visit2_right"), list.files(classical_result_dir_5k, full.names = T), value = T)
  class_visit2_est <- list(
    left = readRDS(classical_file_left)$betas_classical$avg$data$cortex_left,
    right = readRDS(classical_file_right)$betas_classical$avg$data$cortex_right
  )
  return(class_visit2_est)
}, simplify = FALSE)

avail_files <- sapply(subjects, function(subject) {
  subject_visit2_files <- grep(subject, list.files(sphere_5k_result_dir))
  return(length(subject_visit2_files))
})

# This is unnecessary, as the classical doesn't depend on the surface
# classical_visit2_5k_sphere <- list()
#   # sapply(subjects, function(subject) {
#   for(subject in subjects[-seq(25)]) {
#   classical_file_left <- grep(paste0(subject,"_visit2_left"), list.files(sphere_5k_result_dir, full.names = T), value = T)
#   classical_file_right <- grep(paste0(subject,"_visit2_right"), list.files(sphere_5k_result_dir, full.names = T), value = T)
#   ests <- list(
#     left = readRDS(classical_file_left)$betas_classical$avg$data$cortex_left,
#     right = readRDS(classical_file_right)$betas_classical$avg$data$cortex_right
#   )
#   classical_visit2_5k_sphere[[subject]] <- ests
#   }
#   # return(class_visit2_est)
# # }, simplify = F)
# saveRDS(classical_visit2_5k_sphere,"HCP_results/5k_results/607_classical_visit2_sphere_est.rds")
# classical_visit2_5k_sphere <- readRDS("HCP_results/5k_results/607_classical_visit2_sphere_est.rds")

# bayes_visit1_est <- list()
# for(subject in subjects) {
#   result_file_left <- grep(paste0(subject,"_visit1_left"), list.files(bayes_result_dir, full.names = T), value = T)
#   result_file_right <- grep(paste0(subject,"_visit1_right"), list.files(bayes_result_dir, full.names = T), value = T)
#   ests <- list(
#     left = readRDS(result_file_left)$betas_Bayesian$avg$data$cortex_left,
#     right = readRDS(result_file_right)$betas_Bayesian$avg$data$cortex_right
#   )
#   bayes_visit1_est[[subject]] <- ests
# }
# saveRDS(bayes_visit1_est, "HCP_results/5k_results/607_bayes_visit1_est.rds")
bayes_visit1 <- readRDS("HCP_results/5k_results/607_bayes_visit1_est.rds")
bayes_visit1 <- bayes_visit1[keep_subjects]

bayes_visit1_sphere <- list()
  # sapply(subjects, function(subject) {
  for(subject in subjects) {
  bayes_file_left <- grep(paste0(subject,"_visit1_left"), list.files(sphere_5k_result_dir, full.names = T), value = T)
  bayes_file_right <- grep(paste0(subject,"_visit1_right"), list.files(sphere_5k_result_dir, full.names = T), value = T)
  ests <- list(
    left = readRDS(bayes_file_left)$betas_Bayesian$avg$data$cortex_left,
    right = readRDS(bayes_file_right)$betas_Bayesian$avg$data$cortex_right
  )
  bayes_visit1_sphere[[subject]] <- ests
}
#   return(bayes_visit1_est)
# }, simplify = F)
saveRDS(bayes_visit1_sphere,"HCP_results/5k_results/607_bayes_visit1_sphere_est.rds")
bayes_visit1_sphere <- readRDS("HCP_results/5k_results/607_classical_visit2_sphere_est.rds")

classical_visit1_5k <- sapply(subjects, function(subject) {
  classical_file_left <- grep(paste0(subject,"_visit1_left"), list.files(classical_result_dir_5k, full.names = T), value = T)
  classical_file_right <- grep(paste0(subject,"_visit1_right"), list.files(classical_result_dir_5k, full.names = T), value = T)
  class_visit2_est <- list(
    left = readRDS(classical_file_left)$betas_classical$avg$data$cortex_left,
    right = readRDS(classical_file_right)$betas_classical$avg$data$cortex_right
  )
  return(class_visit2_est)
}, simplify = FALSE)

classical_visit1_unsmoothed <- sapply(subjects, function(subject) {
  classical_file_left <- grep(paste0(subject,"_visit1_left"), list.files(classical_result_dir_unsmoothed, full.names = T), value = T)
  classical_file_right <- grep(paste0(subject,"_visit1_right"), list.files(classical_result_dir_unsmoothed, full.names = T), value = T)
  class_visit1_est <- list(
    left = readRDS(classical_file_left)$betas_classical$avg$data$cortex_left,
    right = readRDS(classical_file_right)$betas_classical$avg$data$cortex_right
  )
  return(class_visit1_est)
}, simplify = FALSE)

classical_visit1_smoothed <- sapply(subjects, function(subject) {
  classical_file_left <- grep(paste0(subject,"_visit1_left"), list.files(classical_result_dir_smoothed, full.names = T), value = T)
  classical_file_right <- grep(paste0(subject,"_visit1_right"), list.files(classical_result_dir_smoothed, full.names = T), value = T)
  class_visit1_est <- list(
    left = readRDS(classical_file_left)$betas_classical$avg$data$cortex_left,
    right = readRDS(classical_file_right)$betas_classical$avg$data$cortex_right
  )
  return(class_visit1_est)
}, simplify = FALSE)

# >> MSE ----
classical_mse_32k_unsmoothed <- mapply(function(est,tru) {
  subj_sqerr <- mapply(function(e,t){
    err <- (e - t)^2
    return(err)
  }, e = est, t = tru, SIMPLIFY = F)
  tongue <- mean(c(subj_sqerr[[1]][,1],subj_sqerr[[1]][,1]))
  cue <- mean(c(subj_sqerr[[1]][,4],subj_sqerr[[1]][,4]))
  right_foot <- mean(subj_sqerr[[1]][,2])
  left_foot <- mean(subj_sqerr[[2]][,2])
  right_hand <- mean(subj_sqerr[[1]][,3])
  left_hand <- mean(subj_sqerr[[2]][,3])
  subj_mse_df <- data.frame(
    cue = cue,
    right_foot = right_foot,
    right_hand = right_hand,
    left_foot = left_foot,
    left_hand = left_hand,
    tongue = tongue
  )
  return(subj_mse_df)
}, est = classical_visit1_unsmoothed, tru = classical_visit2_32k, SIMPLIFY = F)

classical_mse_32k_smoothed <- mapply(function(est,tru) {
  subj_sqerr <- mapply(function(e,t){
    err <- (e - t)^2
    return(err)
  }, e = est, t = tru, SIMPLIFY = F)
  tongue <- mean(c(subj_sqerr[[1]][,1],subj_sqerr[[1]][,1]))
  cue <- mean(c(subj_sqerr[[1]][,4],subj_sqerr[[1]][,4]))
  right_foot <- mean(subj_sqerr[[1]][,2])
  left_foot <- mean(subj_sqerr[[2]][,2])
  right_hand <- mean(subj_sqerr[[1]][,3])
  left_hand <- mean(subj_sqerr[[2]][,3])
  subj_mse_df <- data.frame(
    cue = cue,
    right_foot = right_foot,
    right_hand = right_hand,
    left_foot = left_foot,
    left_hand = left_hand,
    tongue = tongue
  )
  return(subj_mse_df)
}, est = classical_visit1_smoothed, tru = classical_visit2_32k, SIMPLIFY = F)

classical_mse_5k <- mapply(function(est,tru) {
  subj_sqerr <- mapply(function(e,t){
    err <- (e - t)^2
    return(err)
  }, e = est, t = tru, SIMPLIFY = F)
  tongue <- mean(c(subj_sqerr[[1]][,1],subj_sqerr[[1]][,1]))
  cue <- mean(c(subj_sqerr[[1]][,4],subj_sqerr[[1]][,4]))
  right_foot <- mean(subj_sqerr[[1]][,2])
  left_foot <- mean(subj_sqerr[[2]][,2])
  right_hand <- mean(subj_sqerr[[1]][,3])
  left_hand <- mean(subj_sqerr[[2]][,3])
  subj_mse_df <- data.frame(
    cue = cue,
    right_foot = right_foot,
    right_hand = right_hand,
    left_foot = left_foot,
    left_hand = left_hand,
    tongue = tongue
  )
  return(subj_mse_df)
}, est = classical_visit1_5k, tru = classical_visit2_5k, SIMPLIFY = F)

bayes_mse <- mapply(function(est,tru) {
  subj_sqerr <- mapply(function(e,t){
    err <- (e - t)^2
    return(err)
  }, e = est, t = tru, SIMPLIFY = F)
  tongue <- mean(c(subj_sqerr[[1]][,1],subj_sqerr[[1]][,1]))
  cue <- mean(c(subj_sqerr[[1]][,4],subj_sqerr[[1]][,4]))
  right_foot <- mean(subj_sqerr[[1]][,2])
  left_foot <- mean(subj_sqerr[[2]][,2])
  right_hand <- mean(subj_sqerr[[1]][,3])
  left_hand <- mean(subj_sqerr[[2]][,3])
  subj_mse_df <- data.frame(
    cue = cue,
    right_foot = right_foot,
    right_hand = right_hand,
    left_foot = left_foot,
    left_hand = left_hand,
    tongue = tongue
  )
  return(subj_mse_df)
}, est = bayes_visit1, tru = classical_visit2_5k, SIMPLIFY = F)

library(tidyverse)
max_vals <- reshape2::melt(bayes_mse) %>%
  mutate(Model = "Bayesian, 5k") %>%
  full_join(
    reshape2::melt(classical_mse_5k) %>%
      mutate(Model = "Classical, 5k, unsmoothed")
  ) %>%
  full_join(
    reshape2::melt(classical_mse_32k_unsmoothed) %>%
      mutate(Model = "Classical, 32k, unsmoothed")
  ) %>%
  full_join(
    reshape2::melt(classical_mse_32k_smoothed) %>%
      mutate(Model = "Classical, 32k, smoothed")
  ) %>%
  mutate(variable = sub("_"," ", variable)) %>%
  group_by(variable) %>%
  summarize(max_value = max(value)) #%>%
  # filter(variable == "tongue")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

color_pal <- c(gg_color_hue(2)[1], gg_color_hue(2)[2], "turquoise2","blue3")

subj_mse_plot <-
  reshape2::melt(bayes_mse) %>%
  mutate(Model = "Bayesian, 5k") %>%
  # full_join(
  #   reshape2::melt(classical_mse) %>%
  #     mutate(Model = "Classical")
  # ) %>%
  full_join(
    reshape2::melt(classical_mse_5k) %>%
      mutate(Model = "Classical, 5k")
  ) %>%
  full_join(
    reshape2::melt(classical_mse_32k_unsmoothed) %>%
      mutate(Model = "Classical, 32k")
  ) %>%
  full_join(
    reshape2::melt(classical_mse_32k_smoothed) %>%
      mutate(Model = "Classical, 32k (smoothed)")
  ) %>%
  mutate(variable = sub("_"," ", variable),
         Model = fct_relevel(Model,
                             "Bayesian, 5k",
                             "Classical, 5k",
                             "Classical, 32k",
                             "Classical, 32k (smoothed)")) %>%
  # pivot_wider(names_from = Model, values_from = value) %>%
  # mutate(diff = abs(Bayesian - Classical)) %>%
  # group_by(variable) %>%
  # mutate(diff = diff / max(diff),
         # diff = ifelse(diff > 0.8,0.8,diff)) %>%
  # filter(variable == "tongue") %>%
  ggplot() +
  geom_boxplot(aes(x = Model, y = value, color = Model)) +
  # geom_point(aes(x = max_value, y = max_value), color = "white", data = max_vals) +
  # geom_segment(aes(x = Classical, y = Bayesian, xend = Classical, yend = Classical, color = diff)) +
  # scale_color_gradientn("",colors = c('yellow','turquoise2','blue2'), limits = c(0,0.8)) +
  # scale_color_gradientn("",
  #                       colors = c(rgb(242,238,73,maxColorValue = 255),
  #                                  rgb(79,195,239,maxColorValue = 255),
  #                                  rgb(71, 109, 174, maxColorValue = 255)),
  #                       limits = c(0, 0.8))+
  # geom_point(aes(y = Bayesian, x = Classical)) +
  # geom_abline(intercept = 0, slope = 1, color = 'grey70') +
  labs(
    y = "Mean Squared Error",
    x = ""
  ) +
  scale_color_manual("",values = color_pal) +
  facet_wrap(~variable, scales = "free",nrow = 1) +
  # guides(color = T) +
  theme_classic() +
  theme(text=element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.margin = margin(0,1,1,1,"pt"))
# theme(plot.title = element_text(hjust = 0.5,size = 7),
#       axis.title = element_text(size = 9))

subj_mse_plot #+ ggtitle("MSE (Bayes found at 5k, Classical found at 32k)")

ggsave("plots/607_mse_vs_full_res_classical.png",plot = subj_mse_plot, width = 7, height = 5)

# >> Correlation ----
classical_cor_5k <- mapply(function(est,tru) {
  tongue <- cor(c(est[[1]][,1],est[[2]][,1]),c(tru[[1]][,1],tru[[2]][,1]))
  cue <- cor(c(est[[1]][,4],est[[2]][,4]),c(tru[[1]][,4],tru[[2]][,4]))
  right_foot <- cor(est[[1]][,2], tru[[1]][,2])
  left_foot <- cor(est[[2]][,2], tru[[2]][,2])
  right_hand <- cor(est[[1]][,3], tru[[1]][,3])
  left_hand <- cor(est[[2]][,3], tru[[2]][,3])
  subj_cor_df <- data.frame(
    cue = cue,
    right_foot = right_foot,
    right_hand = right_hand,
    left_foot = left_foot,
    left_hand = left_hand,
    tongue = tongue
  )
  return(subj_cor_df)
}, est = classical_visit1_5k, tru = classical_visit2_5k, SIMPLIFY = F)

classical_cor_32k_unsmoothed <- mapply(function(est,tru) {
  tongue <- cor(c(est[[1]][,1],est[[2]][,1]),c(tru[[1]][,1],tru[[2]][,1]))
  cue <- cor(c(est[[1]][,4],est[[2]][,4]),c(tru[[1]][,4],tru[[2]][,4]))
  right_foot <- cor(est[[1]][,2], tru[[1]][,2])
  left_foot <- cor(est[[2]][,2], tru[[2]][,2])
  right_hand <- cor(est[[1]][,3], tru[[1]][,3])
  left_hand <- cor(est[[2]][,3], tru[[2]][,3])
  subj_cor_df <- data.frame(
    cue = cue,
    right_foot = right_foot,
    right_hand = right_hand,
    left_foot = left_foot,
    left_hand = left_hand,
    tongue = tongue
  )
  return(subj_cor_df)
}, est = classical_visit1_unsmoothed, tru = classical_visit2_32k, SIMPLIFY = F)

classical_cor_32k_smoothed <- mapply(function(est,tru) {
  tongue <- cor(c(est[[1]][,1],est[[2]][,1]),c(tru[[1]][,1],tru[[2]][,1]))
  cue <- cor(c(est[[1]][,4],est[[2]][,4]),c(tru[[1]][,4],tru[[2]][,4]))
  right_foot <- cor(est[[1]][,2], tru[[1]][,2])
  left_foot <- cor(est[[2]][,2], tru[[2]][,2])
  right_hand <- cor(est[[1]][,3], tru[[1]][,3])
  left_hand <- cor(est[[2]][,3], tru[[2]][,3])
  subj_cor_df <- data.frame(
    cue = cue,
    right_foot = right_foot,
    right_hand = right_hand,
    left_foot = left_foot,
    left_hand = left_hand,
    tongue = tongue
  )
  return(subj_cor_df)
}, est = classical_visit1_smoothed, tru = classical_visit2_32k, SIMPLIFY = F)

bayes_cor <- mapply(function(est,tru) {
  tongue <- cor(c(est[[1]][,1],est[[2]][,1]),c(tru[[1]][,1],tru[[2]][,1]))
  cue <- cor(c(est[[1]][,4],est[[2]][,4]),c(tru[[1]][,4],tru[[2]][,4]))
  right_foot <- cor(est[[1]][,2], tru[[1]][,2])
  left_foot <- cor(est[[2]][,2], tru[[2]][,2])
  right_hand <- cor(est[[1]][,3], tru[[1]][,3])
  left_hand <- cor(est[[2]][,3], tru[[2]][,3])
  subj_cor_df <- data.frame(
    cue = cue,
    right_foot = right_foot,
    right_hand = right_hand,
    left_foot = left_foot,
    left_hand = left_hand,
    tongue = tongue
  )
  return(subj_cor_df)
}, est = bayes_visit1, tru = classical_visit2_5k, SIMPLIFY = F)

bayes_cor_sphere <- mapply(function(est,tru) {
  tongue <- cor(c(est[[1]][,1],est[[2]][,1]),c(tru[[1]][,1],tru[[2]][,1]))
  cue <- cor(c(est[[1]][,4],est[[2]][,4]),c(tru[[1]][,4],tru[[2]][,4]))
  right_foot <- cor(est[[1]][,2], tru[[1]][,2])
  left_foot <- cor(est[[2]][,2], tru[[2]][,2])
  right_hand <- cor(est[[1]][,3], tru[[1]][,3])
  left_hand <- cor(est[[2]][,3], tru[[2]][,3])
  subj_cor_df <- data.frame(
    cue = cue,
    right_foot = right_foot,
    right_hand = right_hand,
    left_foot = left_foot,
    left_hand = left_hand,
    tongue = tongue
  )
  return(subj_cor_df)
}, est = bayes_visit1_sphere, tru = classical_visit2_5k, SIMPLIFY = F)

library(tidyverse)
# max_cors <- reshape2::melt(bayes_cor) %>%
#   mutate(Model = "Bayesian") %>%
#   full_join(
#     reshape2::melt(classical_cor) %>%
#       mutate(Model = "Classical")
#   ) %>%
#   mutate(variable = sub("_"," ", variable)) %>%
#   group_by(variable) %>%
#   summarize(max_value = max(value))
#
# library(viridis)
# library(ciftiTools)
# ciftiTools.setOption('wb_path','/Applications/workbench')
# col_pal <- colorRampPalette(c('red','orange','yellow','green','blue','violet'))
#
# subj_cor_plot <-
#   reshape2::melt(bayes_cor) %>%
#   mutate(Model = "Bayesian") %>%
#   full_join(
#     reshape2::melt(classical_cor) %>%
#       mutate(Model = "Classical")
#   ) %>%
#   mutate(variable = sub("_"," ", variable)) %>%
#   pivot_wider(names_from = Model, values_from = value) %>%
#   mutate(diff = Bayesian - Classical,
#          diff = ifelse(diff < 0, 0, diff),
#          diff = ifelse(diff > 0.1, 0.1, diff)) %>%
#   ggplot() +
#   geom_point(aes(x = max_value, y = max_value), color = "white", data = max_cors) +
#   geom_segment(aes(x = Classical, y = Bayesian, xend = Classical, yend = Classical, color = diff), alpha = 0.8) +
#   # scale_color_gradientn("",colors = c('red','orange','yellow','green','cornflowerblue'), limits = c(0,0.1)) +
#   scale_color_gradientn("",colors = c(
#     'red',
#     'orange',
#     rgb(242,238,73,maxColorValue = 255),
#     rgb(79,195,239,maxColorValue = 255),
#     rgb(71, 109, 174, maxColorValue = 255)
#   ), limits = c(0,0.1)) +
#   # scale_color_gradientn("",colors = rev(c('red','orange','yellow','green','blue')), limits = c(-.18,.18)) +
#   # scale_color_manual("", values = col_pal) +
#   # scale_color_viridis(option = "C", limits = c(0,.1)) +
#   # scale_color_gradient2("",low = "blue",high = "red", limits = c(-.2,.2)) +
#   geom_point(aes(y = Bayesian, x = Classical)) +
#   geom_abline(intercept = 0, slope = 1, color = 'grey70') +
#   labs(y = "Bayesian GLM", x = "Classical GLM") +
#   facet_wrap(~variable, scales = "free") +
#   guides(color = FALSE) +
#   theme_classic() +
#   theme(text = element_text(size = 14))
#
# subj_cor_plot + ggtitle("Correlation (Bayes found at 5k, Classical found at 32k)")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

color_pal <- c(gg_color_hue(2)[1], "red", gg_color_hue(2)[2], "turquoise2","blue3")

subj_cor_plot <-
  reshape2::melt(bayes_cor) %>%
  mutate(Model = "Bayesian, 5k") %>%
  full_join(
    reshape2::melt(bayes_cor_sphere) %>%
      mutate(Model = "Bayesian, 5k (spherical)")
  ) %>%
  full_join(
    reshape2::melt(classical_cor_5k) %>%
      mutate(Model = "Classical, 5k")
  ) %>%
  full_join(
    reshape2::melt(classical_cor_32k_unsmoothed) %>%
      mutate(Model = "Classical, 32k")
  ) %>%
  full_join(
    reshape2::melt(classical_cor_32k_smoothed) %>%
      mutate(Model = "Classical, 32k (smoothed)")
  ) %>%
  mutate(variable = sub("_"," ", variable),
         Model = fct_relevel(Model,
                             "Bayesian, 5k",
                             "Bayesian, 5k (spherical)",
                             "Classical, 5k",
                             "Classical, 32k",
                             "Classical, 32k (smoothed)")) %>%
  # pivot_wider(names_from = Model, values_from = value) %>%
  # mutate(diff = abs(Bayesian - Classical)) %>%
  # group_by(variable) %>%
  # mutate(diff = diff / max(diff),
  # diff = ifelse(diff > 0.8,0.8,diff)) %>%
  # filter(variable == "tongue") %>%
  ggplot() +
  geom_boxplot(aes(x = Model, y = value, color = Model)) +
  # geom_point(aes(x = max_value, y = max_value), color = "white", data = max_vals) +
  # geom_segment(aes(x = Classical, y = Bayesian, xend = Classical, yend = Classical, color = diff)) +
  # scale_color_gradientn("",colors = c('yellow','turquoise2','blue2'), limits = c(0,0.8)) +
  # scale_color_gradientn("",
  #                       colors = c(rgb(242,238,73,maxColorValue = 255),
  #                                  rgb(79,195,239,maxColorValue = 255),
  #                                  rgb(71, 109, 174, maxColorValue = 255)),
  #                       limits = c(0, 0.8))+
  # geom_point(aes(y = Bayesian, x = Classical)) +
  # geom_abline(intercept = 0, slope = 1, color = 'grey70') +
  labs(
    y = "Correlation",
    x = ""
  ) +
  scale_color_manual("",values = color_pal) +
  facet_wrap(~variable, scales = "fixed",nrow = 1) +
  guides(color = guide_legend(nrow = 2)) +
  theme_classic() +
  theme(text=element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.box.spacing = unit(-15,"points"),
        legend.title = element_blank(),
        plot.margin = margin(0,1,1,1,"pt"))
# theme(plot.title = element_text(hjust = 0.5,size = 7),
#       axis.title = element_text(size = 9))

subj_cor_plot #+ ggtitle("MSE (Bayes found at 5k, Classical found at 32k)")

ggsave("plots/607_cor_vs_full_res_classical.png",plot = subj_cor_plot, width = 7, height = 4)

# <<<< NOT USED IN MANUSCRIPT >>>>----
# FIGURE 4: ICC per task ----
# >> ICC function definition ----
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
  # xbar <- mean(X)
  # s2 <- sum((X - xbar)^2) / (prod(dim(X)))
  # xbar_i <- apply(X,1,mean)
  sig2_t <- mean(apply(X,2,var))
  sig2_w <- var(Reduce(`-`,split(X,col(X)))) / 2
  sig2_b <- sig2_t - sig2_w
  sig2_b[sig2_b < 0] <- 0
  r <- sig2_b / sig2_t
  # r <- (ncol(X)*sum((xbar_i - xbar)^2) / (nrow(X)*s2) - 1) / (ncol(X) - 1)
  # if(r < 0) r <- 0
  return(r)
}

# >> weighted ICC (I2C2) function ----
I2C2 <- function(img_X,img_wt) {
  vars_bt <- sapply(asplit(img_X,1),function(X){
    sig2_t <- mean(apply(X,2,var))
    sig2_w <- var(Reduce(`-`,asplit(X,2))) / 2
    sig2_b <- sig2_t - sig2_w
    sig2_b[sig2_b < 0] <- 0
    return(c(sig2_b,sig2_t))
  })
  weighted_vars <- t(vars_bt) * cbind(img_wt, img_wt)
  sum_weighted_vars <- apply(weighted_vars,2,sum)
  out <- sum_weighted_vars[1] / sum_weighted_vars[2]
  return(out)
}

# >> Calculate the weights ----
avg_estimates_classical <- readRDS(file.path(result_dir,"602_avg_estimates_classical.rds"))
ICC_weights <- sapply(avg_estimates_classical, function(hem_est) {
  combined_est <- Reduce(abind,hem_est)
  avg_est <- apply(combined_est,1:2,mean)
  weight_est <- apply(avg_est,2,function(avg_task) abs(avg_task) / sum(abs(avg_task)))
  return(weight_est)
}, simplify = F)

# >> Bayesian ----
result_dir <- "HCP_results/5k_results"
avg_estimates <- readRDS(file.path(result_dir,"602_avg_estimates.rds"))
library(abind)
abind_g <- function(...) abind(...,along = 4)
ICC_values_average <- sapply(avg_estimates, function(hem_est) {
  combined_visits <- Reduce(abind_g,hem_est)
  vertex_ICC <- apply(combined_visits, 1:2, ICC)
  return(vertex_ICC)
}, simplify = F)

library(abind)
abind_g <- function(...) abind(...,along = 4)
bayes_wICC <- round(mapply(function(hem_est, wts) {
  combined_visits <- Reduce(abind_g,hem_est)
  vertex_wICC <- mapply(I2C2, img_X = asplit(combined_visits,2),img_wt = asplit(wts,2))
  return(vertex_wICC)
}, hem_est = avg_estimates,wts = ICC_weights,SIMPLIFY = T),3)

task_file_names <- c("visual_cue","foot","hand","tongue")
bayesian_icc_avg_cifti <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
bayesian_icc_avg_cifti$data$cortex_left <- ICC_values_average$left
bayesian_icc_avg_cifti$data$cortex_right <- ICC_values_average$right
library(viridisLite)
my_pal <- c("black",viridis(4),"white")
# >>>> Visual Cue and Tongue (idx = c(1,4)) ----
# task_idx <- 1
for(task_idx in c(1,4)) {
  plot(
    bayesian_icc_avg_cifti,
    idx = task_idx,
    hemisphere = 'both',
    title = paste("wICC (left, right) =", paste(bayes_wICC[task_idx,], collapse = ",")),
    cex.title = 2,
    color_mode = 'sequential',
    colors = my_pal,
    legend_embed = F,
    zlim = seq(0,0.6,length.out = 6),
    surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
    surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
    fname = paste0("plots/602_",task_file_names[task_idx],"_icc_bayesian.png")
  )
}

# >>>> Foot and Hand Tasks (idx = c(2,3)) ----
task_idx <- 2:3
hem <- c("left","right")
for(ti in task_idx) {
  for(h in 1:2) {
    plot(
      bayesian_icc_avg_cifti,
      idx = ti,
      hemisphere = hem[h],
      title = paste("wICC =", bayes_wICC[ti,h]),
      color_mode = 'sequential',
      colors = my_pal,
      legend_embed = F,
      zlim = seq(0,0.6,length.out = 6),
      surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
      surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
      fname = paste0("plots/602_",hem[h],"_",task_file_names[ti],"_icc_bayesian.png")
    )
  }
}

# >> Classical ----
avg_estimates_classical <- readRDS(file.path(result_dir,"602_avg_estimates_classical.rds"))
library(abind)
abind_g <- function(...) abind(...,along = 4)
ICC_values_average_classical <- sapply(avg_estimates_classical, function(hem_est) {
  combined_visits <- Reduce(abind_g,hem_est)
  vertex_ICC <- apply(combined_visits, 1:2, ICC)
  return(vertex_ICC)
}, simplify = F)

classical_wICC <- round(mapply(function(hem_est, wts) {
  combined_visits <- Reduce(abind_g,hem_est)
  vertex_wICC <- mapply(I2C2, img_X = asplit(combined_visits,2),img_wt = asplit(wts,2))
  # vertex_ICC <- apply(combined_visits, 2, I2C2, img_wt = wts)
  return(vertex_wICC)
}, hem_est = avg_estimates_classical,wts = ICC_weights,SIMPLIFY = T),3)

classical_icc_avg_cifti <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
classical_icc_avg_cifti$data$cortex_left <- ICC_values_average_classical$left
classical_icc_avg_cifti$data$cortex_right <- ICC_values_average_classical$right
library(viridisLite)
my_pal <- c("black",viridis(4),"white")
# >>>> Visual Cue and Tongue (idx = c(1,4)) ----
# task_idx <- 1
for(task_idx in c(1,4)) {
  plot(
    classical_icc_avg_cifti,
    idx = task_idx,
    hemisphere = 'both',
    title = paste("wICC (left, right) =", paste(classical_wICC[task_idx,], collapse = ",")),
    cex.title = 2,
    color_mode = 'sequential',
    colors = my_pal,
    legend_embed = F,
    zlim = seq(0,0.6,length.out = 6),
    surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
    surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
    fname = paste0("plots/602_",task_file_names[task_idx],"_icc_classical")
  )
}

# >>>> Foot and Hand Tasks (idx = c(2,3)) ----
task_idx <- 2:3
hem <- c("left","right")
for(ti in task_idx) {
  for(h in 1:2) {
    plot(
      classical_icc_avg_cifti,
      idx = ti,
      hemisphere = hem[h],
      title = paste("wICC =", classical_wICC[ti,h]),
      color_mode = 'sequential',
      colors = my_pal,
      legend_embed = F,
      zlim = seq(0,0.6,length.out = 6),
      surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
      surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
      fname = paste0("plots/602_",hem[h],"_",task_file_names[ti],"_icc_classical.png")
    )
  }
}

# FIGURE 5: weighted ICC by task ----

# >> weighted ICC (I2C2) function ----
I2C2 <- function(img_X,img_wt) {
  vars_bt <- sapply(asplit(img_X,1),function(X){
    sig2_t <- mean(apply(X,2,var))
    sig2_w <- var(Reduce(`-`,asplit(X,2))) / 2
    sig2_b <- sig2_t - sig2_w
    sig2_b[sig2_b < 0] <- 0
    return(c(sig2_b,sig2_t))
  })
  weighted_vars <- t(vars_bt) * cbind(img_wt, img_wt)
  sum_weighted_vars <- apply(weighted_vars,2,sum)
  out <- sum_weighted_vars[1] / sum_weighted_vars[2]
  return(out)
}

# >> Calculate the weights ----
result_dir <- "HCP_results/5k_results"
avg_estimates_classical <- readRDS(file.path(result_dir,"602_avg_estimates_classical.rds"))
ICC_weights <- sapply(avg_estimates_classical, function(hem_est) {
  combined_est <- Reduce(abind,hem_est)
  avg_est <- apply(combined_est,1:2,mean)
  weight_est <- apply(avg_est,2,function(avg_task) abs(avg_task) / sum(abs(avg_task)))
  return(weight_est)
}, simplify = F)


# >> Bayesian ----
result_dir <- "HCP_results/5k_results"
session_estimates <- readRDS(file.path(result_dir,"602_session_estimates.rds"))
avg_estimates <- readRDS(file.path(result_dir,"602_avg_estimates.rds"))
Bayesian_estimates <- sapply(1:2, function(sess_num) {
  sapply(c("left","right"), function(hem) {
    sapply(paste0("visit",1:2), function(visit_num) {
      session_estimates[[hem]][[visit_num]][,,sess_num,]
    }, simplify = F)
  }, simplify = F)
}, simplify = F)
names(Bayesian_estimates) <- c("LR","RL")
Bayesian_estimates$Average <- avg_estimates
library(abind)
abind_g <- function(...) abind(...,along = 4)
bayes_wICC <- sapply(Bayesian_estimates, function(sess) {
  sess_wICC <- round(mapply(function(hem_est, wts) {
    combined_visits <- Reduce(abind_g,hem_est)
    vertex_wICC <- mapply(I2C2, img_X = asplit(combined_visits,2),img_wt = asplit(wts,2))
    return(vertex_wICC)
  }, hem_est = sess,wts = ICC_weights,SIMPLIFY = T),3)
  rownames(sess_wICC) <- c("cue","foot","hand","tongue")
  return(sess_wICC)
}, simplify = F)

# >> Classical ----
result_dir <- "HCP_results/5k_results"
avg_estimates_classical <- readRDS(file.path(result_dir,"602_avg_estimates_classical.rds"))
session_estimates_classical <- readRDS(file.path(result_dir,"602_session_estimates_classical.rds"))
Classical_estimates <- sapply(1:2, function(sess_num) {
  sapply(c("left","right"), function(hem) {
    sapply(paste0("visit",1:2), function(visit_num) {
      session_estimates_classical[[hem]][[visit_num]][,,sess_num,]
    }, simplify = F)
  }, simplify = F)
}, simplify = F)
names(Classical_estimates) <- c("LR","RL")
Classical_estimates$Average <- avg_estimates_classical
library(abind)
abind_g <- function(...) abind(...,along = 4)
classical_wICC <- sapply(Classical_estimates, function(sess) {
  sess_wICC <- round(mapply(function(hem_est, wts) {
    combined_visits <- Reduce(abind_g,hem_est)
    vertex_wICC <- mapply(I2C2, img_X = asplit(combined_visits,2),img_wt = asplit(wts,2))
    return(vertex_wICC)
  }, hem_est = sess,wts = ICC_weights,SIMPLIFY = T),3)
  rownames(sess_wICC) <- c("cue","foot","hand","tongue")
  return(sess_wICC)
}, simplify = F)

# >> Plotting ----
library(tidyverse)
wicc_scatter <- reshape2::melt(bayes_wICC) %>%
  mutate(model = "Bayesian") %>%
  full_join(
    reshape2::melt(classical_wICC) %>%
      mutate(model = "Classical")
  ) %>%
  mutate(sess_or_avg = ifelse(L1 == "Average","Average","LR/RL")) %>%
  ggplot() +
  geom_point(aes(x = Var2, y = value, color = model, shape = L1), size = 5) +
  facet_grid(sess_or_avg~Var1) +
  scale_color_discrete("Model") +
  scale_shape_discrete("Session") +
  labs(y = "weighted ICC", x = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5))

wicc_scatter

ggsave("plots/607_wicc_task_scatterplot.png", plot = wicc_scatter, width = 5, height = 4)


# FIGURE 7: Group Estimates and Activations----
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW"
group_result_files <- list.files(result_dir)

# >> Estimates ----
# >>>> Bayesian ----
library(tidyverse)
bayes_files <- grep("BayesGLM2",group_result_files, value = T) %>%
  grep(pattern = "_45subj_", x = ., value = T) %>%
  grep(pattern = "_thresh0_20210131.rds", x = ., value = T)
left_result <- readRDS(file.path(result_dir,grep("_left_",bayes_files,value = T)))
right_result <- readRDS(file.path(result_dir,grep("_right_",bayes_files,value = T)))

library(ciftiTools)
ciftiTools.setOption('wb_path', "/Applications/workbench/")

bayes_estimates <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
bayes_estimates$data$cortex_left <- left_result$estimates
bayes_estimates$data$cortex_right <- right_result$estimates

plot(bayes_estimates, idx = 4, zlim = c(-1,1),
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
     fname = "plots/607_group_tongue_estimates.png")

# >>>> Classical ----
class_group <- readRDS("HCP_results/5k_results/502_HCP_classical_group_PW_estimates.rds")
class_estimates <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
class_estimates$data$cortex_left <- class_group$left
class_estimates$data$cortex_right <- class_group$right
plot(class_estimates, zlim = c(-1,1), idx = 4,
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
     fname = "plots/607_group_tongue_estimates_classical.png")

# >> Poster images ----
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW"
group_result_files <- list.files(result_dir)
library(tidyverse)
bayes_10 <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples/501_HCP_BayesGLM2_10subj_sample1_visit1_result_left_thresh0_20210319.rds")$estimates
bayes_45 <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/501_HCP_BayesGLM2_45subj_visit1_result_left_thresh0_20210318.rds")$estimates

classical_subs <- readRDS("HCP_results/5k_results/group/502_HCP_classical_subsample_estimates.rds")
classical_10 <- classical_subs$`10subjects`$sample1$left
classical_45 <- readRDS("HCP_results/5k_results/502_HCP_classical_group_PW_estimates.rds")$left

library(ciftiTools)
ciftiTools.setOption('wb_path',"/Applications/workbench/")
cifti_obj <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
cifti_obj <- remove_xifti(cifti_obj, "cortex_right")
task_idx <- 4
task_names <- c("cue","foot","hand","tongue")
cifti_obj$data$cortex_left <-
  cbind(
    classical_10[,task_idx],
    classical_45[,task_idx],
    bayes_10[,task_idx],
    bayes_45[,task_idx]
  )
model_names <- c("classical_10","classical_45","bayes_10","bayes_45")
for(i in 1:4) {
  plot(cifti_obj, hemisphere = "left", idx = i,
       fname = paste0(
         "~/Desktop/SMI 2021 poster images/607_",
         model_names[i],
         "_",task_names[task_idx],
         ".png"
         ),
       view = "lateral", legend_embed = F, zlim = c(-1,1),
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")
}

# FIGURE 8: ICC Thresholded Images ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench/')


task_file_names <- c("visual_cue","foot","hand","tongue")
bayesian_icc_avg_cifti <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
bayesian_icc_avg_cifti$data$cortex_left <- ICC_quality$left
bayesian_icc_avg_cifti$data$cortex_right <- ICC_quality$right
library(viridisLite)
# my_pal <- c("black",viridis(4),"white")
my_pal <- c(viridis(4))
my_pal[1] <- "white"
# >>>> Visual Cue and Tongue (idx = c(1,4)) ----
# task_idx <- 1
for(task_idx in c(1,4)) {
  plot(
    bayesian_icc_avg_cifti,
    idx = task_idx,
    hemisphere = 'both',
    # title = paste("wICC (left, right) =", paste(bayes_wICC[task_idx,], collapse = ",")),
    # cex.title = 2,
    color_mode = 'qualitative',
    colors = my_pal,
    legend_embed = F,
    # zlim = seq(0,0.6,length.out = 6),
    surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
    surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
    fname = paste0("plots/607_",task_file_names[task_idx],"_icc_quality_bayesian.png")
  )
}

# >>>> Foot and Hand Tasks (idx = c(2,3)) ----
task_idx <- 2:3
hem <- c("left","right")
for(ti in task_idx) {
  for(h in 1:2) {
    plot(
      bayesian_icc_avg_cifti,
      idx = ti,
      hemisphere = hem[h],
      # title = paste("wICC =", bayes_wICC[ti,h]),
      color_mode = 'qualitative',
      colors = my_pal,
      legend_embed = F,
      # zlim = seq(0,0.6,length.out = 6),
      surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
      surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
      fname = paste0("plots/607_",hem[h],"_",task_file_names[ti],"_icc_quality_bayesian.png")
    )
  }
}

# >> Classical ----
result_dir <- "HCP_results/5k_results"
avg_estimates <- readRDS(file.path(result_dir,"602_average_estimates_classical.rds"))
library(abind)
abind_g <- function(...) abind(...,along = 4)
ICC_values_average <- sapply(avg_estimates, function(hem_est) {
  combined_visits <- Reduce(abind_g,hem_est)
  vertex_ICC <- apply(combined_visits, 1:2, ICC)
  return(vertex_ICC)
}, simplify = F)

ICC_quality_classical <- sapply(ICC_values_average, function(hem_icc) {
  apply(hem_icc,2,function(task_icc) {
    icc_rating <-
      cut(
        task_icc,
        breaks = c(-Inf, 0.4, 0.6, 0.75, Inf),
        labels = c(1, 2, 3, 4)
      )
    return(as.numeric(icc_rating))
  })
})

task_file_names <- c("visual_cue","foot","hand","tongue")
classical_icc_avg_cifti <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
classical_icc_avg_cifti$data$cortex_left <- ICC_quality_classical$left
classical_icc_avg_cifti$data$cortex_right <- ICC_quality_classical$right
library(viridisLite)
# my_pal <- c("black",viridis(4),"white")
my_pal <- c(viridis(4))
my_pal[1] <- "white"
# >>>> Visual Cue and Tongue (idx = c(1,4)) ----
# task_idx <- 1
for(task_idx in c(1,4)) {
  plot(
    classical_icc_avg_cifti,
    idx = task_idx,
    hemisphere = 'both',
    # title = paste("wICC (left, right) =", paste(bayes_wICC[task_idx,], collapse = ",")),
    # cex.title = 2,
    color_mode = 'qualitative',
    colors = my_pal,
    legend_embed = F,
    # zlim = seq(0,0.6,length.out = 6),
    surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
    surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
    fname = paste0("plots/607_",task_file_names[task_idx],"_icc_quality_classical.png")
  )
}

# >>>> Foot and Hand Tasks (idx = c(2,3)) ----
task_idx <- 2:3
hem <- c("left","right")
for(ti in task_idx) {
  for(h in 1:2) {
    plot(
      classical_icc_avg_cifti,
      idx = ti,
      hemisphere = hem[h],
      # title = paste("wICC =", bayes_wICC[ti,h]),
      color_mode = 'qualitative',
      colors = my_pal,
      legend_embed = F,
      # zlim = seq(0,0.6,length.out = 6),
      surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
      surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
      fname = paste0("plots/607_",hem[h],"_",task_file_names[ti],"_icc_quality_classical.png")
    )
  }
}





# FIGURE 11: Group areas of activation ----


# FIGURE 12: Single-session estimates ----

# FIGURE 13: Single-session activations ----
# >> Bayesian ----
library(ciftiTools)
ciftiTools.setOption('wb_path',"/Applications/workbench")
library(BayesfMRI)
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session"
result_files <- list.files(result_dir, full.names = T)
result_files <- grep(".rds", result_files, value= T)
load("HCP_data/subjects.Rdata")
subjects <- subjects[1]
sessions <- c("LR")
# task_idx <- 1
task_names <- c("cue","foot","hand","tongue")
thresholds <- c(0,0.5,1)
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
col_pal <- c(light_orange,"red","purple")
for(subject in subjects) {
  for(v in 1) {
    for(sess in sessions) {
      cifti_obj <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
      for(task_idx in 1:3) {
        for(h in c("left","right")) {
          not_h <- grep(h, c("left","right"), value = T, invert = T)
          if(task_idx %in% c(1,4)) fld_nm <- task_names[task_idx]
          if(task_idx %in% c(2,3)) fld_nm <- paste0(not_h,"_",task_names[task_idx])
          L_or_R <- toupper(substring(h,1,1))
          result_obj <-
            readRDS(grep(
              paste0("500_", subject, "_visit", v, "_", h, "_5k_session", sess),
              result_files,
              value = T
            ))
          active <- vector('numeric', length = nrow(result_obj$betas_Bayesian[[sess]]$data[[paste0("cortex_",h)]]))
          for(thr in thresholds){
            act_obj <-
              id_activations_cifti(
                model_obj = result_obj,
                field_names = fld_nm,
                alpha = 0.01,
                method = "Bayesian",
                threshold = thr
              )
            active <- active + act_obj$activations[[paste0("cortex",L_or_R)]]$active[,1]
          }
          active[active == 0] <- NA
          cifti_obj$data[[paste0("cortex_",h)]] <- as.matrix(active)
        }
        if(task_idx %in% c(1,4)) {
          plot(
            cifti_obj,
            zlim = c(-1,1), legend_embed = F,
            color_mode = "qualitative", colors = col_pal,
            fname = paste0(
              "plots/607_bayes_",
              subject,
              "_visit",
              v,
              "_session",
              sess,
              "_",
              task_names[task_idx],
              "_activations.png"
            ),
            surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
            surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii"
          )
        }
        if(task_idx %in% c(2,3)) {
          for(hem in c("left","right")) {
            plot(
              cifti_obj, hemisphere = hem,
              zlim = c(-1,1), legend_embed = F,
              color_mode = "qualitative", colors = col_pal,
              fname = paste0(
                "plots/607_bayes_",
                subject,
                "_visit",
                v,
                "_session",
                sess,
                "_",
                hem,
                "_",
                task_names[task_idx],
                "_activations.png"
              ),
              surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
              surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii"
            )
          }
        }
      }
    }
  }
}
# >> Classical ----
library(ciftiTools)
ciftiTools.setOption('wb_path',"/Applications/workbench")
library(BayesfMRI)
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session"
result_files <- list.files(result_dir, full.names = T)
result_files <- grep(".rds", result_files, value= T)
load("HCP_data/subjects.Rdata")
subjects <- subjects[1]
sessions <- c("LR")
# task_idx <- 4
task_names <- c("cue","foot","hand","tongue")
thresholds <- c(0,0.5,1)
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
col_pal <- c(light_orange,"red","purple")
for(subject in subjects) {
  for(v in 1) {
    for(sess in sessions) {
      for(task_idx in 1:4) {
        cifti_obj <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
        for(h in c("left","right")) {
          not_h <- grep(h, c("left","right"), value = T, invert = T)
          if(task_idx %in% c(1,4)) fld_nm <- task_names[task_idx]
          if(task_idx %in% c(2,3)) fld_nm <- paste0(not_h,"_",task_names[task_idx])
          L_or_R <- toupper(substring(h,1,1))
          result_obj <-
            readRDS(grep(
              paste0("500_", subject, "_visit", v, "_", h, "_5k_session", sess),
              result_files,
              value = T
            ))
          active <- vector('numeric', length = nrow(result_obj$betas_Bayesian[[sess]]$data[[paste0("cortex_",h)]]))
          for(thr in thresholds){
            act_obj <-
              id_activations_cifti(
                model_obj = result_obj,
                field_names = fld_nm,
                alpha = 0.01,
                method = "classical",
                threshold = thr
              )
            active <- active + act_obj$activations[[paste0("cortex",L_or_R)]]$active[,1]
          }
          active[active == 0] <- NA
          cifti_obj$data[[paste0("cortex_",h)]] <- as.matrix(active)
        }
        if(task_idx %in% c(1,4)) {
          plot(
            cifti_obj,
            zlim = c(-1,1), legend_embed = F,
            color_mode = "qualitative", colors = col_pal,
            fname = paste0(
              "plots/607_classical_",
              subject,
              "_visit",
              v,
              "_session",
              sess,
              "_",
              task_names[task_idx],
              "_activations.png"
            ),
            surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
            surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii"
          )
        }
        if(task_idx %in% c(2,3)) {
          for(hem in c('left','right')){
            plot(
              cifti_obj, hemisphere = hem,
              zlim = c(-1,1), legend_embed = F,
              color_mode = "qualitative", colors = col_pal,
              fname = paste0(
                "plots/607_classical_",
                subject,
                "_visit",
                v,
                "_session",
                sess,
                "_",
                hem,
                "_",
                task_names[task_idx],
                "_activations.png"
              ),
              surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
              surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii"
            )
          }
        }
      }
    }
  }
}

# FIGURE 16: Group Dice Coefficients and activation areas ----
library(tidyverse)
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples"
group_activation_files <- list.files(result_dir, full.names = T)
for(num_subj in paste0(c(10,20,30),"subj")) {

}


# FIGURE 18: Single- vs. Multi-session estimates ----
# This uses plots from figures on single session estimates and single-subject estimates

# FIGURE 19: Single- vs. Multi-session activation ----
# This uses plots from figures on single session activations and single-subject areas of activation

#' Combine left and right xifti objects into a single bilateral xifti object
#'
#' @param left an object of class \code{xifti} that has data and metadata for the left cortex
#' @param right an object of class \code{xifti} that has data and metadata for the right cortex
#'
#' @export
#' @return an object of class \code{xifti} with data and metadata for both the left and right hemispheres
#'
#' @examples
xifti_combine <- function(left,right) {
  # Checks
  if(!is.matrix(left$data$cortex_left))
    stop("The left xifti object has no data for the left cortex.")
  if (!is.matrix(right$data$cortex_right))
    stop("The right xifti object has no data for the right cortex.")
  if(is.null(left$meta$cortex$medial_wall_mask$left) |
     is.null(right$meta$cortex$medial_wall_mask$right))
    stop("One or both of the hemisphere xiftis is missing a medial wall mask.")

  cifti_out <- left
  cifti_out$data$cortex_right <- right$data$cortex_right
  cifti_out$meta$cortex$medial_wall_mask$right <-
    right$meta$cortex$medial_wall_mask$right
  return(cifti_out)
}

# FIGURE 20: Dice x area for single vs. multisession data ----
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session"
result_files <- list.files(result_dir, full.names = T)
load("HCP_data/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
all_activations <-
  sapply(subjects, function(subject) {
    sapply(c("LR","RL"), function(session) {
      session_activations <- sapply(paste0("visit",1:2), function(visit) {
        sapply(c("left","right"), function(hem) {
          L_or_R <- toupper(substring(hem,1,1))
          result_obj <-
            readRDS(grep(
              paste0(subject, "_", visit, "_", hem, "_5k_session", session),
              result_files,
              value = T
            ))
          threshold_activations <- sapply(c(0,0.5,1), function(threshold) {
            bayes_activations <-
              id_activations_cifti(
                model_obj = result_obj,
                alpha = 0.01,
                method = "Bayesian",
                threshold = threshold
              )
            classical_activations <-
              id_activations_cifti(
                model_obj = result_obj,
                alpha = 0.01,
                method = "classical",
                threshold = threshold
              )
            return(
              list(
                Bayes = bayes_activations$activations[[paste0("cortex",L_or_R)]]$active,
                Classical = classical_activations$activations[[paste0("cortex",L_or_R)]]$active
              )
            )
          }, simplify = F)
        })
      })
    })
  })

# FIGURE 21: Motor Cortex Parcellation ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench/')
library(BayesfMRI)
load("HCP_data/subjects.Rdata")
subjects <- subjects[c(1,2,3)]
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/all_sessions"
result_files <- list.files(result_dir, full.names = T)
result_files <- unlist(sapply(subjects, grep, x = result_files, value = T))
result_files <- grep("_thr0.5.rds", result_files, value = T)
out <- list()
for(subject in subjects) {
  out[[subject]] <- list()
  for(hem in c("left","right")) {
    out[[subject]][[hem]] <- list()
    L_or_R <- toupper(substring(hem,1,1))
    # result_obj <- readRDS(grep(paste0(subject,"_visit1_",hem), result_files, value = T))
    result_obj <- readRDS(grep(paste0(subject,"_",hem), result_files, value = T))
    # act_result <- id_activations_cifti(result_obj,alpha = 0.01,threshold = 1)
    # act_mat <- act_result$activations[[paste0("cortex",L_or_R)]]$active
    # overlap_idx <- which(apply(act_mat,1,sum) > 1)
    # View(act_mat[overlap_idx,])
    # out[[subject]][[hem]] <- table(apply(act_mat,1,sum))
    act_mat <- result_obj$active
    out[[subject]][[hem]] <- act_mat
  }
}

sapply(out, function(x) sapply(x, function(y) table(apply(y,1,sum))))

task_names <- c("cue","foot","hand","tongue")

# subj_1_labels <- mapply(function(act, nm) {
#   out <- act
#   out[out == 1] <- nm
#   return(out)
# }, act = split(out$`103818`$left, col(out$`103818`$left)), nm = task_names)


# sub1_label_vec <- unlist(apply(subj_1_labels,1,function(x) {
#   out <- NULL
#   u_x <- unique(x)
#   if(length(u_x) == 1) out <- NA
#   if(length(u_x) > 1) out <- paste(u_x[which(u_x != "0")], collapse = " & ")
#   return(out)
# }))

cifti_parcellation <- sapply(out, function(act_i) {
  subject_labels <- sapply(act_i, function(act_ih) {
    labels_mat <- mapply(function(act, nm) {
      out <- act
      out[out == 1] <- nm
      return(out)
    }, act = split(act_ih[,-1], col(act_ih[,-1])), nm = task_names[-1])
    labels_vec <- unlist(apply(labels_mat,1,function(x) {
      out <- NULL
      u_x <- unique(x)
      if(length(u_x) == 1) out <- NA
      if(length(u_x) > 1) out <- paste(u_x[which(u_x != "0")], collapse = " & ")
      return(out)
    }))
    return(as.matrix(as.numeric(factor(labels_vec))))
  }, simplify = F)
  cifti_out <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
  cifti_out$data$cortex_left <- subject_labels$left
  cifti_out$data$cortex_right <- subject_labels$right
  return(cifti_out)
}, simplify = F)

view_xifti_surface(cifti_parcellation$`103818`,
                   color_mode = "qualitative",
                   title = "Subject A")
view_xifti_surface(cifti_parcellation$`105923`,
                   color_mode = "qualitative",
                   title = "Subject B")
view_xifti_surface(cifti_parcellation$`114823`,
                   color_mode = "qualitative",
                   title = "Subject C")

# FIGURE 22: AR coefficient map ----
left_obj <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session/500_103818_visit1_left_5k_sessionLR_20210322.rds")
right_obj <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session/500_103818_visit1_right_5k_sessionLR_20210322.rds")
cifti_obj <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
cifti_obj$data$cortex_left <- left_obj$prewhitening_info$left$AR_coeffs
cifti_obj$data$cortex_right <- right_obj$prewhitening_info$right$AR_coeffs
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
plot(cifti_obj, hemisphere = 'left',idx = 1,
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")

# FIGURE 27: Mesh structure ----
library(ciftiTools)
ciftiTools.setOption('wb_path',"/Applications/workbench")

# >> midthickness ----
cifti_obj <- read_cifti(cifti_fname = cifti_fname <- c("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/103818/MNINonLinear/Results/tfMRI_MOTOR_LR/tfMRI_MOTOR_LR_Atlas.dtseries.nii"),
                        surfL_fname = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/103818/MNINonLinear/fsaverage_LR32k/103818.L.midthickness.32k_fs_LR.surf.gii",
                        surfR_fname = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/103818/MNINonLinear/fsaverage_LR32k/103818.R.midthickness.32k_fs_LR.surf.gii",resamp_res = 5000)

cifti_obj$data$cortex_left <- array(0, dim = dim(cifti_obj$data$cortex_left))
cifti_obj$data$cortex_right <- array(0, dim = dim(cifti_obj$data$cortex_right))

view_xifti_surface(
  xifti = cifti_obj, color_mode = "qualitative", colors = "white",
  surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/103818/MNINonLinear/fsaverage_LR32k/103818.L.midthickness.32k_fs_LR.surf.gii",
  surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/103818/MNINonLinear/fsaverage_LR32k/103818.R.midthickness.32k_fs_LR.surf.gii",borders = F, view = "lateral", hemisphere = "left", edge_color = "black",legend_embed = F,
  fname = "plots/607_mesh_cifti.png")

# >> spherical ----
cifti_obj <- read_cifti(cifti_fname = cifti_fname <- c("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/103818/MNINonLinear/Results/tfMRI_MOTOR_LR/tfMRI_MOTOR_LR_Atlas.dtseries.nii"),
                        surfL_fname = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Sphere.32k.L.surf.gii",
                        surfR_fname = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Sphere.32k.R.surf.gii",resamp_res = 5000)

cifti_obj$data$cortex_left <- array(0, dim = dim(cifti_obj$data$cortex_left))
cifti_obj$data$cortex_right <- array(0, dim = dim(cifti_obj$data$cortex_right))

view_xifti_surface(
  xifti = cifti_obj, color_mode = "qualitative", colors = "white",
  surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Sphere.32k.L.surf.gii",
  surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Sphere.32k.R.surf.gii",borders = F, view = "lateral", hemisphere = "left", edge_color = "black",legend_embed = F,
  fname = "plots/607_spherical_mesh_cifti.png")

# FIGURE 28: Full-resolution classical results ----
library(ciftiTools)
ciftiTools.setOption('wb_path',"/Applications/workbench")
result_dir <- "HCP_results/32k_results"
result_files <- list.files(result_dir, full.names = T)
load("HCP_data/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
task_names <- c("cue",'foot','hand','tongue')
xifti_combine <- function(left,right) {
  # Checks
  if(!is.matrix(left$data$cortex_left))
    stop("The left xifti object has no data for the left cortex.")
  if (!is.matrix(right$data$cortex_right))
    stop("The right xifti object has no data for the right cortex.")
  if(is.null(left$meta$cortex$medial_wall_mask$left) |
     is.null(right$meta$cortex$medial_wall_mask$right))
    stop("One or both of the hemisphere xiftis is missing a medial wall mask.")

  cifti_out <- left
  cifti_out$data$cortex_right <- right$data$cortex_right
  cifti_out$meta$cortex$medial_wall_mask$right <-
    right$meta$cortex$medial_wall_mask$right
  return(cifti_out)
}
# >> Estimates ----
for(subject in subjects) {
  for(visit in paste0('visit',1:2)) {
      filen_L <- grep(
        paste0(subject,"_",visit,"_left"),
        result_files, value = T)
      filen_R <- grep(
        paste0(subject,"_",visit,"_right"),
        result_files, value = T)
      result_obj_L <- readRDS(filen_L)
      result_obj_R <- readRDS(filen_R)
      avg_cifti <- xifti_combine(result_obj_L$betas_classical$avg,result_obj_R$betas_classical$avg)
      for(task_idx in 1:4) {
        if(task_idx %in% c(1,4)) {
          plot(avg_cifti, idx = task_idx, zlim = c(-1,1), legend_embed =F,
               surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
               surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
               fname = paste0("plots/607_15k_classical_",
                              subject,
                              "_",
                              visit,
                              "_",
                              task_names[task_idx],
                              "_estimate.png"))
        }
        if(task_idx %in% c(2,3)) {
          for(hem in c('left','right')) {
            plot(avg_cifti, idx = task_idx, zlim = c(-1,1), legend_embed =F, hemisphere = hem,
                 surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
                 surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
                 fname = paste0("plots/607_15k_classical_",
                                subject,
                                "_",
                                visit,
                                "_",
                                hem,
                                "_",
                                task_names[task_idx],
                                "_estimate.png"))
          }
        }
      }
  }
}

# >> Activations ----
library(BayesfMRI)
for(subject in subjects) {
  for(visit in paste0('visit',1:2)) {
    # for(hem in c('left','right')) {
    # stop()
    filen_L <- grep(
      paste0(subject,"_",visit,"_left"),
      result_files, value = T)
    filen_R <- grep(
      paste0(subject,"_",visit,"_right"),
      result_files, value = T)
    result_obj_L <- readRDS(filen_L)
    act_L <- id_activations_cifti()
    result_obj_R <- readRDS(filen_R)
    avg_cifti <- xifti_combine(result_obj_L$betas_classical$avg,result_obj_R$betas_classical$avg)
    plot(avg_cifti, idx = 1, zlim = c(-1,1), legend_embed =F,
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         fname = paste0("plots/607_classical_",subject,"_",visit,"_15k_estimate.png"))
    # }
  }
}

#  FIGURES 7 & 8: Number of runs x Number of activations (Subject) ----
load("HCP_data/subjects.Rdata")

multi_result_dir <- "HCP_results/5k_results/individual/PW/activations"
multi_result_files <- list.files(multi_result_dir, full.names = T)
multi_result_files <- grep("_visit1_", multi_result_files, value = T)
multi_result_files <- grep("_single_session", multi_result_files, value = T, invert = T)
# result_files <- grep("_sample", result_files, value = T)

bayes_multi_num_act <-
  sapply(subjects, function(subject) {
    sapply(c('left','right'), function(hem) {
      sapply(as.character(c(0,0.5,1)), function(thr) {
        thr_chr <- sub("\\.","",thr)
        L_or_R <- toupper(substring(hem,1,1))
        filen <- grep(subject, multi_result_files, value = T)
        filen <- grep(hem, filen, value = T)
        filen <- grep(paste0("thr",thr_chr,"_"), filen, value = T)
        filen <- grep("classical", filen, invert = T, value = T)
        result_obj <- readRDS(filen)
        return(as.matrix(apply(result_obj$activations[[paste0("cortex",L_or_R)]]$active,2,sum)))
      }, simplify = F)
    }, simplify = F)
  }, simplify = F)

classical_multi_num_act <-
  sapply(subjects, function(subject) {
    sapply(c('left','right'), function(hem) {
      sapply(as.character(c(0,0.5,1)), function(thr) {
        thr_chr <- sub("\\.","",thr)
        L_or_R <- toupper(substring(hem,1,1))
        filen <- grep(subject, multi_result_files, value = T)
        filen <- grep(hem, filen, value = T)
        filen <- grep(paste0("thr",thr_chr,"_"), filen, value = T)
        filen <- grep("classical", filen, invert = F, value = T)
        result_obj <- readRDS(filen)
        return(as.matrix(apply(result_obj$active,2,sum)))
      }, simplify = F)
    }, simplify = F)
  }, simplify = F)

# single_result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session"
# single_result_files <- list.files(single_result_dir, full.names = T)
# single_result_files <- grep("visit1", single_result_files, value = T)
# single_result_files <- grep("sessionLR", single_result_files, value = T)
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
library(INLA)
inla.setOption(pardiso.license = "~/licenses/pardiso.lic") # Dan's Macbook Pro
library(BayesfMRI)

# single_num_act <-
#   sapply(subjects, function(subject) {
#     sapply(c('left','right'), function(hem) {
#       L_or_R <- toupper(substring(hem,1,1))
#       filen <- grep(paste0(subject,"_visit1_",hem), single_result_files, value = T)
#       result_obj <- readRDS(filen)
#       hem_output <- sapply(c('Bayesian','classical'), function(mthd) {
#         mthd_output <- sapply(c(0,0.5,1), function(thr) {
#           active_obj <- id_activations_cifti(
#             model_obj = result_obj,
#             field_names = NULL,
#             session_name = "LR",
#             alpha = 0.01,
#             method = mthd,
#             threshold = thr
#           )
#           if(mthd == "Bayesian") {
#             colnames(active_obj$activations[[paste0("cortex",L_or_R)]]$active) <- NULL
#             return(as.matrix(apply(active_obj$activations[[paste0("cortex",L_or_R)]]$active,2,sum)))
#           }
#           if(mthd == 'classical') {
#             return(as.matrix(apply(active_obj$activations[[paste0("cortex",L_or_R)]]$active,2,sum)))
#           }
#         }, simplify = F)
#         return(mthd_output)
#       }, simplify = F)
#       return(hem_output)
#     }, simplify = F)
#   }, simplify = F)
single_num_act <- readRDS("HCP_results/5k_results/individual/PW/activations/all_single_run_activations.rds")

# >> Plot ----
col_pal <- scales::hue_pal()(5)
library(viridis)
col_pal <- viridis(7,option = "C")[1:5]
col_pal <- rgb(c(230, 86, 0, 0,204),c(159,180,158,114,121),c(0,233,115,178,167), maxColorValue = 255)
library(tidyverse)
task_names <- c("cue","foot","hand","tongue")
num_act_plots <- sapply(c(0,0.5,1), function(thr) {
  single_num_act_df <-
    reshape2::melt(single_num_act, value.name = "number_active") %>%
    mutate(threshold = c(0,0.5,1)[L4],
           task = task_names[Var1],
           runs = "Single") %>%
    rename(subject = L1,
           hem = L2,
           Model = L3) %>%
    select(-starts_with("Var"), -L4) %>%
    filter(threshold == thr)
  multi_bayes_num_act_df <-
    reshape2::melt(bayes_multi_num_act, value.name = "number_active") %>%
    mutate(Model = "Bayesian",
           task = task_names[Var1],
           runs = "Multiple",
           L3 = as.numeric(L3)) %>%
    rename(subject = L1,
           hem = L2,
           threshold = L3) %>%
    select(-Var1,-Var2) %>%
    filter(threshold == thr)
  multi_classical_num_act_df <-
    reshape2::melt(classical_multi_num_act, value.name = "number_active") %>%
    mutate(Model = "classical",
           task = task_names[Var1],
           runs = "Multiple",
           L3 = as.numeric(L3)) %>%
    rename(subject = L1,
           hem = L2,
           threshold = L3) %>%
    select(-Var1,-Var2) %>%
    filter(threshold == thr)
  combined_df <-
    full_join(single_num_act_df, multi_bayes_num_act_df) %>%
    full_join(multi_classical_num_act_df) %>%
    mutate(not_hem = ifelse(hem == 'left', 'right','left'),
           task = paste(not_hem,task)) %>%
    select(-hem, -not_hem) %>%
    pivot_wider(names_from = task, values_from = number_active) %>%
    mutate(cue = `left cue` + `right cue`,
           tongue = `left tongue` + `right tongue`,
           runs = ifelse(runs == "Single", 1, 2)) %>%
    select(-ends_with(" cue"), -ends_with(" tongue")) %>%
    pivot_longer(cols = -c(threshold,subject,Model,runs), names_to = "task", values_to = "num_active") %>%
    filter(task != "cue")
  means_df <-
    group_by(combined_df, Model, threshold, runs, task) %>%
    summarize(num_active = median(num_active)) %>%
    filter(task != "cue")
  num_act_plot <- ggplot(combined_df, aes(x = runs, y = num_active, color = Model)) +
    geom_line(aes(group = interaction(subject,Model)),alpha = 0.1) +
    geom_line(aes(group = Model),data = means_df,lwd = 2) +
    facet_wrap(~task, scales = "free",ncol = 5) +
    labs(x = "Number of Runs", y = "Number of Active Locations") +
    # scale_color_manual("", values = col_pal) +
    scale_color_discrete("") +
    scale_x_continuous(breaks=c(1,2), limits = c(0.9,2.1)) +
    theme_classic() +
    theme(legend.position = "right",
          text = element_text(size = 14))
  # num_act_plot
  return(num_act_plot)
  # },simplify = F)
}, simplify = F)

num_act_plots[[3]] # 1% threshold
num_act_plots[[2]] # 0.5% threshold
num_act_plots[[1]] # 0% threshold

ggsave("plots/607_num_act_plots_thr0.png", plot = num_act_plots[[1]], width = 7, height = 6)
ggsave("plots/607_num_act_plots_thr0.5.png", plot = num_act_plots[[2]], width = 7, height = 6)
