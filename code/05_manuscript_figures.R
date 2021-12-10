# This is a script meant to provide reproducible code for each figure in the
# manuscript. If opening this in RStudio, use Ctrl + O to view the script
# outline, which will make this script much more navigable.

# TABLE 1: Computation times ----
bayes_single_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/experimental_inla"

bayes_files <- list.files(bayes_single_dir, full.names = T) |> grep(pattern = ".rds", value = TRUE) |> grep(pattern = "500_", value = T)

single_subject_times <- sapply(bayes_files, function(bf) {
  read_obj <- readRDS(bf)
  subject <- substring(bf, 91,96)
  visit <- substring(bf, 98,103)
  hem <- substring(bf, 105,109)
  hem <- sub("_","",hem)
  L_or_R <- toupper(substring(hem,1,1))
  out <- data.frame(
    subject = subject,
    visit = visit,
    hem = hem,
    Bayesian_time = read_obj$GLMs_Bayesian[[paste0("cortex",L_or_R)]]$total_time,
    classical_time = read_obj$GLMs_classical[[paste0("cortex",L_or_R)]]$total_time,
    # em_time = read_obj$GLMs_EM[[paste0("cortex",L_or_R)]]$total_time,
    overall_time = read_obj$total_time
  )
  return(out)
}, simplify = F)

saveRDS(single_subject_times,file = "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/05_single_subject_times_PW.rds")

# Classical smoothed times
result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
result_files <- list.files(result_dir,full.names = TRUE) |>
  grep(pattern = "FWHM6", value = T)
result_files <- grep(".rds", result_files, value = T)

classical_times <- sapply(result_files, function(x) {
  y <- readRDS(x)
  return(y$GLMs_classical[[which(!sapply(y$GLMs_classical,is.null))]]$total_time)
})
mean(classical_times/60) #0.1190365
sd(classical_times/60) #0.05329523

# Load the Bayesian subject times
single_subject_times <- readRDS(file.path(bayes_single_dir,"05_single_subject_times_PW.rds"))

single_subject_times <- Reduce(rbind,single_subject_times)
names(single_subject_times)[4] <- "Bayesian_time"

library(tidyverse)

# Find out which subjects don't have both hemispheres of data
single_subject_times %>%
  group_by(subject) %>%
  tally %>%
  filter(n != 2)
# A tibble: 1 × 2
# subject     n
# <chr>   <int>
#   1 250427      1

# single_subj_df <-
single_subject_times %>%
  filter(subject != "250427") %>%
  pivot_longer(cols = ends_with("_time"), names_to = "model", values_to = "Time") %>%
  group_by(subject, visit, model) %>%
  summarize(Time = sum(Time)) %>%
  pivot_wider(names_from = model, values_from = Time) %>%
  mutate(preprocessing_time = overall_time - Bayesian_time - classical_time) %>%
  pivot_longer(cols = -c("subject","visit"),names_to = "Category", values_to = "Time") %>%
  mutate(Category = sub("_time","", Category)) %>%
  filter(Category != "em", Category != "overall") %>%
  group_by(Category) %>%
  summarize(Mean = round(mean(Time)/60,3), SD = round(sd(Time)/60,3)) %>%
  ungroup() %>%
  mutate(Mean_SD = paste0(Mean," (",SD,")")) %>%
  select(-Mean,-SD) %>%
  pivot_wider(names_from = Category, values_from = Mean_SD)

group_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples"
group_files <- list.files(group_dir, full.names = T) |>
  grep(pattern = "501_HCP_", value = TRUE) #|>
# grep(pattern = "_visit1_", value = TRUE) |>
# grep(pattern = "_sample", value = TRUE) #|>
# grep(pattern = "_thresh0_", value = TRUE)

subsample_times <- sapply(paste0(c(10,20,30),"subj"), function(N) {
  sapply(paste0("sample",1:10), function(sample_num) {
    sapply(paste0("visit",1:2), function(visit_num) {
      sapply(c("left","right"), function(hem) {
        sapply(paste0("thresh",c("0","05","1")), function(thr) {
          gf <- grep(paste0(N,"_",sample_num,"_",visit_num,"_result_",hem,"_", thr),
                     group_files, value = T)
          if(length(gf) > 1) gf <- gf[1]
          return(readRDS(gf)$total_time)
        }, simplify = F)
      }, simplify = F)
    }, simplify = F)
  }, simplify = F)
}, simplify = F)

subsample_df <- reshape2::melt(subsample_times, value.name = "Time") %>%
  rename(threshold = L5,
         hem = L4,
         visit_num = L3,
         sample_num = L2,
         num_subjects = L1)

group_times <- sapply(paste0("visit",1:2), function(visit_num) {
  sapply(c("left","right"), function(hem) {
    sapply(paste0("thresh",c("0","05","1")), function(thr) {
      gf <- grep(paste0("45subj_result_",hem,"_", thr),
                 group_files, value = T)
      if(length(gf) > 1) gf <- gf[1]
      return(readRDS(gf)$total_time)
    }, simplify = F)
  }, simplify = F)
}, simplify = F)

group_df <- reshape2::melt(group_times, value.name = "Time") %>%
  rename(
    threshold = L3,
    hem = L2,
    visit_num = L1
  ) %>%
  mutate(num_subjects = "45subj", sample_num = "sample1")

group_times_df <- full_join(subsample_df, group_df)

group_time_summary <- group_times_df %>%
  group_by(sample_num, num_subjects, visit_num, threshold) %>%
  summarize(Time = sum(Time)) %>%
  group_by(num_subjects) %>%
  summarize(Mean = round(mean(Time/60),2),
            SD = round(sd(Time/60), 2))

paste(group_time_summary$Mean, collapse = " & ")
paste(group_time_summary$SD, collapse = " & ")

group_times <- sapply(group_files, function(gf) {
  num_subjects <- substring(gf,100,101)
  sample_num <- substring(gf, 113, 114)
  sample_num <- sub("_","", sample_num)
  hem <- substring(gf, 129, 133)
  hem <- sub("_","", hem)
  if(hem == "righ") hem <- "right"
  total_time <- readRDS(gf)$total_time
  out <- data.frame(N = num_subjects, Time = total_time, hem = hem, sample_num = sample_num)
  return(out)
}, simplify = F)

group_times <- Reduce(rbind,group_times)

group_times %>%
  filter(N == "30", hem == "right", sample_num == "10")
group_by(N, hem, sample_num) %>%
  tally %>%
  filter(n != 6)

group_df <-
  group_times %>%
  mutate(N = paste0("N = ",N)) %>%
  group_by(N) %>%
  summarize(Mean = round(mean(Time)/60,2), SD = round(sd(Time)/60,2),
            Mean_SD = paste0(Mean," (",SD,")")) %>%
  select(-Mean,-SD) %>%
  pivot_wider(names_from = N, values_from = Mean_SD)

group45_files <- list.files(group_dir, full.names = TRUE) |>
  grep(pattern = "45subj", value = T)
# grep(pattern = "_202104", value = T)

group45_times <- sapply(group45_files, function(g45) {
  out <- data.frame(
    N = 45,
    Time = readRDS(g45)$total_time
  )
  return(out)
}, simplify = F)

group45_times <- Reduce(rbind, group45_times)

group45_df <-
  group45_times %>%
  mutate(N = paste0("N = ",N)) %>%
  group_by(N) %>%
  summarize(Mean = round(mean(Time)/60,2), SD = round(sd(Time)/60,2),
            Mean_SD = paste0(Mean," (",SD,")")) %>%
  select(-Mean,-SD) %>%
  pivot_wider(names_from = N, values_from = Mean_SD)


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
load("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/subjects.Rdata")
# load("HCP_data/subjects.Rdata")
subjects <- subjects[1]
sessions <- c("LR")
task_idx <- 3
task_names <- c("visual_cue","foot","hand","tongue")
# subject <- subjects[1]; v <- 1; sess <- sessions[1]; h <- 'left'
for(subject in subjects) {
  for(v in 1) {
    for(sess in sessions) {
      cifti_obj <- readRDS("HCP_data/cifti_5k_template_whole.rds")
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
            "plots/5_bayes_",
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
              "plots/5_bayes_",
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
result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
result_files <- list.files(result_dir, full.names = T) |>
  grep(pattern = "FWHM6", value= T)
load("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
task_names <- c("visual_cue","foot","hand","tongue")
# subject <- subjects[1]; visit <- 1; run <- "LR"; hem <- "left"
for(subject in subjects) {
  for(visit in 1) {
    for(run in c("avg")) {
      cifti_obj <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/603_cifti_5k_template_whole.rds")
      for(hem in c("left","right")) {
        result_obj <- readRDS(
          grep(pattern = subject,result_files,value = T) |>
            grep(pattern = paste0("visit",visit),value = TRUE) |>
            grep(pattern = hem,value = T)
        )
        cifti_obj$data[[paste0("cortex_",hem)]] <- result_obj$betas_classical[[run]]$data[[paste0("cortex_",hem)]]
      }
        for(task_idx in 1:4) {
          # cifti_obj <- result_obj$betas_classical[[run]]
          if(task_idx %in% c(1,4)) {
            plot(
              cifti_obj,
              idx = task_idx,zlim = c(-1,1), legend_embed = F,
              fname = paste0("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/600_subject_",subject,"_",task_names[task_idx],"_classical_estimates.png"
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
                fname = paste0("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/600_subject_",subject,"_",hem,"_",task_names[task_idx],"_classical_estimates.png"
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
  cifti_template <- readRDS("HCP_data/cifti_5k_template_whole.rds")
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
         fname = paste0("plots/5_subject_",s,"_",task_names[task_idx],"_estimates.png"),
         zlim = c(-1,1), legend_embed = F)
  }
  if(task_idx %in% c(2,3)) {
    for(hem in c("left","right")) {
      plot(cifti_template, idx = task_idx, hemisphere = hem,
           surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
           surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
           fname = paste0("plots/5_subject_",s,"_",hem,"_",task_names[task_idx],"_estimates.png"),
           zlim = c(-1,1), legend_embed = F)
    }
  }
}
# >>>> Classical ----
# task_idx <- 1
# task_names <- c("visual_cue","foot","hand","tongue")
for(s in subjects) {
  cifti_template <- readRDS("HCP_data/cifti_5k_template_whole.rds")
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
         fname = paste0("plots/5_subject_",s,"_",task_names[task_idx],"_classical_estimates.png"),
         zlim = c(-1,1), legend_embed = F)
  }
  if(task_idx %in% c(2,3)) {
    for(hem in c("left","right")) {
      plot(cifti_template, idx = task_idx, hemisphere = hem,
           surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
           surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
           fname = paste0("plots/5_subject_",s,"_",hem,"_",task_names[task_idx],"_classical_estimates.png"),
           zlim = c(-1,1), legend_embed = F)
    }
  }
}

# >>>> Classical (unsmoothed) ----
# task_idx <- 1
# task_names <- c("visual_cue","foot","hand","tongue")
unsmoothed_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/PW"
result_files <- list.files(unsmoothed_dir, full.names = T)
load("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
task_names <- c("visual_cue","foot","hand","tongue")
task_idx <- 4
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
for(s in subjects) {
  cifti_template <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/603_cifti_5k_template_whole.rds")
  left_file <- grep('left',grep(s,result_files, value = T), value = T)
  if(length(left_file) > 1) left_file <- left_file[1]
  left_result <- readRDS(left_file)
  cifti_template$data$cortex_left <- left_result$betas_classical$avg$data$cortex_left
  rm(left_result)
  right_file <- grep('right',grep(s,result_files, value = T), value = T)
  if(length(right_file) > 1) right_file <- right_file[1]
  right_result <- readRDS(right_file)
  cifti_template$data$cortex_right <- right_result$betas_classical$avg$data$cortex_right
  rm(right_result)
  if(task_idx %in% c(1,4)) {
    plot(cifti_template, idx = task_idx, #title = paste("Subject",s,task_name,"Task"),
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         fname = file.path(plot_dir,paste0("5_subject_",s,"_",task_names[task_idx],"_classical_unsmoothed_estimates.png")),
         zlim = c(-1,1), legend_embed = F)
  }
  if(task_idx %in% c(2,3)) {
    for(hem in c("left","right")) {
      plot(cifti_template, idx = task_idx, hemisphere = hem,
           surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
           surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
           fname = file.path(plot_dir,paste0("plots/5_subject_",s,"_",hem,"_",task_names[task_idx],"_classical_unsmoothed_estimates.png")),
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

# Grab estimates ----
# >> Bayesian ----
result_dir <- "HCP_results/5k_results"
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
result_files <- list.files(result_dir,full.names = TRUE)
result_files <- grep(".rds", result_files, value = T)
# Individual sessions
session_estimates <- sapply(c('left','right'), function(hem) {
  out_list_l1 <- sapply(c('visit1','visit2'), function(visit) {
    file_locs <- grep(hem, grep(visit, result_files,value = TRUE), value = TRUE)
    session_estimates <-
      sapply(file_locs,
             function(file_n){
               result_obj <- readRDS(file_n)
               which_cortex <- paste0('cortex_',hem)
               dim_est <- dim(result_obj$betas_Bayesian$LR$data[[which_cortex]])
               out_obj <- array(NA,dim = c(dim_est,2)) # Third index corresponds to session (1-LR, 2-RL)
               out_obj[,,1] <- result_obj$betas_Bayesian$LR$data[[which_cortex]]
               out_obj[,,2] <- result_obj$betas_Bayesian$RL$data[[which_cortex]]
               return(out_obj)
             }, simplify = 'array')
  }, simplify = F)
  return(out_list_l1)
}, simplify = F)
saveRDS(session_estimates, file.path(result_dir,"5_session_estimates.rds"))
# Averages
avg_estimates <- sapply(c('left','right'), function(hem) {
  out_list_l1 <- sapply(c('visit1','visit2'), function(visit) {
    file_locs <- grep(hem, grep(visit, result_files,value = TRUE), value = TRUE)
    session_estimates <-
      sapply(file_locs,
             function(file_n){
               result_obj <- readRDS(file_n)
               which_cortex <- paste0('cortex_',hem)
               out_obj <- result_obj$betas_Bayesian$avg$data[[which_cortex]]
               return(out_obj)
             }, simplify = 'array')
  }, simplify = F)
  return(out_list_l1)
}, simplify = F)
saveRDS(avg_estimates, file.path(result_dir,"5_avg_estimates.rds"))
avg_estimates <- readRDS(file.path(result_dir,"5_avg_estimates.rds"))
# >>  Classical ----
# Individual sessions
# session_estimates_classical <- sapply(c('left','right'), function(hem) {
#   out_list_l1 <- sapply(c('visit1','visit2'), function(visit) {
#     file_locs <- grep(hem, grep(visit, result_files,value = TRUE), value = TRUE)
#     session_estimates <-
#       sapply(file_locs,
#              function(file_n){
#                result_obj <- readRDS(file_n)
#                which_cortex <- paste0('cortex_',hem)
#                dim_est <- dim(result_obj$betas_classical$LR$data[[which_cortex]])
#                out_obj <- array(NA,dim = c(dim_est,2)) # Third index corresponds to session (1-LR, 2-RL)
#                out_obj[,,1] <- result_obj$betas_classical$LR$data[[which_cortex]]
#                out_obj[,,2] <- result_obj$betas_classical$RL$data[[which_cortex]]
#                return(out_obj)
#              }, simplify = 'array')
#   }, simplify = F)
#   return(out_list_l1)
# }, simplify = F)
# saveRDS(session_estimates_classical, file.path(result_dir,"5_session_estimates_classical.rds"))
# Averages
result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
result_files <- list.files(result_dir,full.names = TRUE) |>
  grep(pattern = "FWHM6", value = T)
result_files <- grep(".rds", result_files, value = T)

load("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/subjects.Rdata")
avg_estimates_classical <- list(
  left = list(
    visit1 = array(NA, dim = c(4443,4,length(subjects))),
    visit2 = array(NA, dim = c(4443,4,length(subjects)))
  ),
  right = list(
    visit1 = array(NA, dim = c(4444,4,length(subjects))),
    visit2 = array(NA, dim = c(4444,4,length(subjects)))
  )
)
for(subject in subjects) {
  subject_idx <- which(subjects == subject)
  for(visit in paste0("visit",1:2)) {
    for(hem in c('left','right')) {
      subj_hem_file <- grep(pattern = subject, result_files, value = T) |>
        grep(pattern = visit, value = T) |>
        grep(pattern = hem, value = T)
      if(length(subj_hem_file) > 1) subj_hem_file <- subj_hem_file[1]
      result_obj <- readRDS(subj_hem_file)
      avg_estimates_classical[[hem]][[visit]][,,subject_idx]<-
        result_obj$betas_classical$avg$data[[paste0("cortex_",hem)]]
    }
  }
}
saveRDS(avg_estimates_classical, file.path(result_dir,"05_avg_classical_estimates.rds"))
avg_estimates_classical <- readRDS(file.path(result_dir,"05_avg_classical_estimates.rds"))
# >> Classical (unsmoothed) ----
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/classical"
result_files <- list.files(result_dir,full.names = TRUE) |>
  grep(pattern = "500_", value = T)

load("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/subjects.Rdata")
avg_estimates_classical <- list(
  left = list(
    visit1 = array(NA, dim = c(4443,4,length(subjects))),
    visit2 = array(NA, dim = c(4443,4,length(subjects)))
  ),
  right = list(
    visit1 = array(NA, dim = c(4444,4,length(subjects))),
    visit2 = array(NA, dim = c(4444,4,length(subjects)))
  )
)
for(subject in subjects) {
  subject_idx <- which(subjects == subject)
  for(visit in paste0("visit",1:2)) {
    for(hem in c('left','right')) {
      subj_hem_file <- grep(pattern = subject, result_files, value = T) |>
        grep(pattern = visit, value = T) |>
        grep(pattern = hem, value = T)
      if(length(subj_hem_file) > 1) subj_hem_file <- subj_hem_file[1]
      result_obj <- readRDS(subj_hem_file)
      avg_estimates_classical[[hem]][[visit]][,,subject_idx]<-
        result_obj$betas_classical$avg$data[[paste0("cortex_",hem)]]
    }
  }
}
saveRDS(avg_estimates_classical, "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/unsmoothed/05_avg_classical_unsmoothed_estimates.rds")

# Calculate the ICC ----
# >> Bayesian ----
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
avg_estimates <- readRDS(file.path(result_dir,"5_avg_estimates.rds"))
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
result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
avg_estimates_classical <- readRDS(file.path(result_dir,"05_avg_classical_estimates.rds"))
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

# >> Classical (unsmoothed) ----
result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/unsmoothed"
avg_estimates_classical_unsmoothed <-
  readRDS(file.path(result_dir,"05_avg_classical_unsmoothed_estimates.rds"))

library(abind)
abind_g <- function(...) abind(...,along = 4)
ICC_values_average_classical_unsmoothed <- sapply(avg_estimates_classical_unsmoothed, function(hem_est) {
  combined_visits <- Reduce(abind_g,hem_est)
  vertex_ICC <- apply(combined_visits, 1:2, ICC)
  return(vertex_ICC)
}, simplify = F)

ICC_quality_classical_unsmoothed <-
  sapply(ICC_values_average_classical_unsmoothed, function(hem_icc) {
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

# >> Make the plot ----
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
      mutate(Model = "Classical (6mm)")
  ) %>%
  full_join(
    reshape2::melt(ICC_quality_classical_unsmoothed) %>%
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

ICC_df %>%
  group_by(L1, Model, task, quality) %>%
  tally %>%
  filter(Model == "Classical", L1 == "right")

ICC_df %>%
  group_by(L1, Model, task, quality) %>%
  tally %>%
  pivot_wider(names_from = L1, values_from = n) %>%
  mutate(both = left+right) %>%
  pivot_longer(cols = c(left,right,both), names_to = "hem") %>%
  filter((task %in% c("Foot", "Hand") & hem %in% c("left","right")) |
           (task %in% c("Visual Cue","Tongue") & hem == "both")) %>%
  mutate(hem = ifelse(hem == 'left','right',ifelse(hem == 'right','left',hem)),
         hem = str_to_title(hem),
         task = paste(hem,task),
         task = sub("Both ","", task),
         value = ifelse(is.na(value),0,value)) %>%
  select(-hem) %>%
  ungroup() %>%
  group_by(task, Model) %>%
  summarize(total = sum(value),
            prop = value / total,
            quality = unique(quality)) %>%
  pivot_wider(names_from = Model, values_from = prop) %>%
  mutate(B_C6_pval = prop.test(x = c(Bayesian, `Classical (6mm)`)*total,
                              n = rep(total,2))$p.value,
         C6_C_pval = prop.test(x = c(`Classical (6mm)`, Classical)*total,
                               n = rep(total,2))$p.value) %>%
  filter(quality != "Poor")

ICC_quality_df <-
  ICC_df %>%
  group_by(L1, Model, task, quality) %>%
  tally %>%
  pivot_wider(names_from = L1, values_from = n) %>%
  mutate(both = left+right) %>%
  pivot_longer(cols = c(left,right,both), names_to = "hem") %>%
  filter((task %in% c("Foot", "Hand") & hem %in% c("left","right")) |
           (task %in% c("Visual Cue","Tongue") & hem == "both")) %>%
  mutate(hem = ifelse(hem == 'left','right',ifelse(hem == 'right','left',hem)),
         hem = str_to_title(hem),
         task = paste(hem,task),
         task = sub("Both ","", task),
         value = ifelse(is.na(value),0,value),
         Model = factor(Model, levels = c("Classical","Classical (6mm)","Bayesian"))) %>%
  select(-hem)

# Facet by task
icc_quality_plot <- group_by(ICC_quality_df,Model, task) %>%
  mutate(total = sum(value),
         prop = value / total,
         Model = factor(Model, levels = c("Classical","Classical (6mm)","Bayesian"),
                        labels = c("C","C6mm","B"))) %>%
  filter(quality != "Poor") %>%
  ggplot(aes(x = Model, y = prop, fill = quality)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~task, scales = "fixed", nrow = 1) +
  # scale_y_continuous(breaks = c(0,1), labels = c(0,1))+
  scale_fill_manual("ICC Quality",values = my_pal) +
  labs(x = "", y = "Proportion of Data Locations") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5))
  # coord_flip()

icc_quality_plot

ggsave("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/607_icc_quality_plot.png",plot = icc_quality_plot, height = 3.5, width = 10)

# FIGURE 4: Subject Test-Retest MSE and Correlation ----
# >> MSE ----
bayes_result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
classical_result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
classical_result_files <- grep("FWHM6", list.files(classical_result_dir, full.names = T), value = T) |>
  grep(pattern = "classical", value = T)
classical_visit2_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/classical"
classical_visit2_files <- list.files(classical_visit2_dir, full.names = T) |>
  grep(pattern = "visit2", value = T)

load("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/subjects.Rdata")
# subject <- subjects[1]
classical_visit2 <- sapply(subjects, function(subject) {
  classical_file_left <- grep(paste0(subject,"_visit2_left"), classical_visit2_files, value = T)
  classical_file_right <- grep(paste0(subject,"_visit2_right"), classical_visit2_files, value = T)
  if(length(classical_file_left) > 1) classical_file_left <- classical_file_left[1]
  if(length(classical_file_right) > 1) classical_file_right <- classical_file_right[1]
  class_visit2_est <- list(
    left = readRDS(classical_file_left)$betas_classical$avg$data$cortex_left,
    right = readRDS(classical_file_right)$betas_classical$avg$data$cortex_right
  )
  return(class_visit2_est)
}, simplify = FALSE)

# bayes_visit1_est <- sapply(subjects, function(x) out = list(left = matrix(NA,4443,4), right = matrix(NA,4444,4)), simplify = F)
# for(subject in subjects) {
#   result_file_left <- grep(paste0(subject,"_visit1_left"), list.files(bayes_result_dir, full.names = T), value = T)
#   result_file_right <- grep(paste0(subject,"_visit1_right"), list.files(bayes_result_dir, full.names = T), value = T)
#   bayes_visit1_est[[subject]]$left <- readRDS(result_file_left)$betas_Bayesian$LR$data$cortex_left
#   bayes_visit1_est[[subject]]$right <- readRDS(result_file_right)$betas_Bayesian$LR$data$cortex_right
# }
# saveRDS(bayes_visit1_est, "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/5_bayes_visit1_LR_est.rds")
# bayes_visit1 <- readRDS("HCP_results/5k_results/5_bayes_visit1_est.rds")
bayes_visit1 <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/5_bayes_visit1_est.rds")
# bayes_visit1 <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/5_bayes_visit1_LR_est.rds")

classical_visit1 <- sapply(subjects, function(subject) {
  classical_file_left <- grep(paste0(subject,"_visit1_left"), classical_result_files, value = T)
  classical_file_right <- grep(paste0(subject,"_visit1_right"), classical_result_files, value = T)
  if(length(classical_file_left) > 1) classical_file_left <- classical_file_left[1]
  if(length(classical_file_right) > 1) classical_file_right <- classical_file_right[1]
  class_visit1_est <- list(
    left = readRDS(classical_file_left)$betas_classical$avg$data$cortex_left,
    right = readRDS(classical_file_right)$betas_classical$avg$data$cortex_right
  )
  return(class_visit1_est)
}, simplify = FALSE)

# If we want to mask to just the active voxels, we can do something like this:
classical_group_activations <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/group/502_HCP_classical_activations_PW_visit2_FWER.rds")
classical_group_0 <- sapply(classical_group_activations,`[[`,"0%", simplify = F)
keep_left <- apply(classical_group_0$left, 2, function(x) which(x != 0))
keep_right <- apply(classical_group_0$right, 2, function(x) which(x != 0))

# If no masking is desired, use these
# keep_left <- sapply(1:4, function(x) seq(4443), simplify = F)
# keep_right <- sapply(1:4, function(x) seq(4444), simplify = F)

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

subj_mse_df <-
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
         diff = ifelse(diff > 0.8,0.8,diff)) #%>%
  # filter(L1 != "562345") %>%
  # filter(variable == "tongue") %>%

# subj_mse_test <-
  subj_mse_df %>%
  select(variable, Bayesian, Classical) %>%
  group_by(variable) %>%
  summarize(estimate = t.test(x = Bayesian, y = Classical, paired = TRUE)$estimate,
          p.value = t.test(x = Bayesian, y = Classical, paired = TRUE)$p.value)

  subj_mse_df %>%
    select(variable, Bayesian, Classical) %>%
    mutate(diff = Bayesian - Classical) %>%
    # pivot_longer(cols = -variable, names_to = "Model", values_to = "MSE") %>%
    ggplot() +
    # geom_boxplot(aes(x = Model,
    #                  color = Model,
    #                  y = MSE),outlier.alpha = 0) +
    # geom_jitter(aes(x = Model, color = Model, y = MSE)) +
    geom_boxplot(aes(y = diff), outlier.alpha = 0) +
    geom_jitter(aes(x = 0,y = diff)) +
    geom_hline(aes(yintercept = 0), lty = 2) +
    labs(y = "Difference in MSE\n(Bayesian - Classical)", x = "") +
    facet_wrap(~variable, scales = "free_y") +
    theme_classic() +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank())

subj_mse_plot <- ggplot(subj_mse_df) +
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
    # title = "Test-Retest MSE"
  ) +
  facet_wrap(~variable, scales = "free") +
  guides(color = 'none') +
  theme_classic() +
  theme(text=element_text(size = 14))
  # theme(plot.title = element_text(hjust = 0.5,size = 14),
  #       axis.title = element_text(size = 14))

subj_mse_plot



ggsave("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/05_subject_mse_plot.png", plot = subj_mse_plot, width = 7.5, height = 5)
# ggsave("~/Desktop/SMI 2021 poster images/5_subject_mse_plot.png", plot = subj_mse_plot, width = 3.5, height = 4)

# >>>> Now try to see if a paired t-test for the differences can be performed ----
reshape2::melt(bayes_mse) %>%
  mutate(Model = "Bayesian") %>%
  full_join(
    reshape2::melt(classical_mse) %>%
      mutate(Model = "Classical")
  ) %>%
  mutate(variable = sub("_"," ", variable)) %>%
  pivot_wider(names_from = Model, values_from = value) %>%
  mutate(diff = Bayesian - Classical) %>%
  ggplot(aes(x = diff)) +
  geom_histogram(bins = 100) +
  facet_wrap(~variable)

reshape2::melt(bayes_mse) %>%
  mutate(Model = "Bayesian") %>%
  full_join(
    reshape2::melt(classical_mse) %>%
      mutate(Model = "Classical")
  ) %>%
  mutate(variable = sub("_"," ", variable)) %>%
  pivot_wider(names_from = Model, values_from = value) %>%
  mutate(diff = Bayesian - Classical) %>%
  group_by(variable) %>%
  summarize(p_values = t.test(Bayesian,Classical, paired = T)$p.value,
            Bayes_minus_Classical = t.test(Bayesian,Classical, paired = T)$estimate)

# >> Correlation ----
bayes_result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
classical_result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
classical_result_files <- grep("FWHM6", list.files(classical_result_dir, full.names = T), value = T)
classical_visit2_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/classical"
classical_visit2_files <- list.files(classical_visit2_dir, full.names = T) |>
  grep(pattern = "visit2", value = T)

load("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/subjects.Rdata")

classical_visit2 <- sapply(subjects, function(subject) {
  classical_file_left <- grep(paste0(subject,"_visit2_left"), classical_result_files, value = T)
  classical_file_right <- grep(paste0(subject,"_visit2_right"), classical_result_files, value = T)
  if(length(classical_file_left) > 1) classical_file_left <- classical_file_left[1]
  if(length(classical_file_right) > 1) classical_file_right <- classical_file_right[1]
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
# saveRDS(bayes_visit1_est, "HCP_results/5k_results/5_bayes_visit1_est.rds")
bayes_visit1 <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/5_bayes_visit1_est.rds")

classical_visit1 <- sapply(subjects, function(subject) {
  classical_file_left <- grep(paste0(subject,"_visit1_left"), classical_result_files, value = T)
  classical_file_right <- grep(paste0(subject,"_visit1_right"), classical_result_files, value = T)
  if(length(classical_file_left) > 1) classical_file_left <- classical_file_left[1]
  if(length(classical_file_right) > 1) classical_file_right <- classical_file_right[1]
  class_visit1_est <- list(
    left = readRDS(classical_file_left)$betas_classical$avg$data$cortex_left,
    right = readRDS(classical_file_right)$betas_classical$avg$data$cortex_right
  )
  return(class_visit1_est)
}, simplify = FALSE)

# If we want to mask to just the active voxels, we can do something like this:
classical_group_activations <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/group/502_HCP_classical_activations_PW_visit2_FWER.rds")
classical_group_0 <- sapply(classical_group_activations,`[[`,"0%", simplify = F)
keep_left <- apply(classical_group_0$left, 2, function(x) which(x != 0))
keep_right <- apply(classical_group_0$right, 2, function(x) which(x != 0))

# If no masking is desired, use these
# keep_left <- sapply(1:4, function(x) seq(4443), simplify = F)
# keep_right <- sapply(1:4, function(x) seq(4444), simplify = F)

classical_cor <- mapply(function(est,tru) {
  tongue <-
    cor(c(est[[1]][keep_left[[4]], 4], est[[2]][keep_right[[4]], 4]),
        c(tru[[1]][keep_left[[4]], 4], tru[[2]][keep_right[[4]], 4]))
  cue <- cor(c(est[[1]][keep_left[[1]],1],est[[2]][keep_right[[1]],1]),
             c(tru[[1]][keep_left[[1]],1],tru[[2]][keep_right[[1]],1]))
  right_foot <- cor(est[[1]][keep_left[[2]],2], tru[[1]][keep_left[[2]],2])
  left_foot <- cor(est[[2]][keep_right[[2]],2], tru[[2]][keep_right[[2]],2])
  right_hand <- cor(est[[1]][keep_left[[3]],3], tru[[1]][keep_left[[3]],3])
  left_hand <- cor(est[[2]][keep_right[[3]],3], tru[[2]][keep_right[[3]],3])
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
  tongue <- cor(c(est[[1]][keep_left[[4]],4],est[[2]][keep_right[[4]],4]),
                c(tru[[1]][keep_left[[4]],4],tru[[2]][keep_right[[4]],4]))
  cue <- cor(c(est[[1]][keep_left[[1]],1],est[[2]][keep_right[[1]],1]),
             c(tru[[1]][keep_left[[1]],1],tru[[2]][keep_right[[1]],1]))
  right_foot <- cor(est[[1]][keep_left[[2]],2], tru[[1]][keep_left[[2]],2])
  left_foot <- cor(est[[2]][keep_right[[2]],2], tru[[2]][keep_right[[2]],2])
  right_hand <- cor(est[[1]][keep_left[[3]],3], tru[[1]][keep_left[[3]],3])
  left_hand <- cor(est[[2]][keep_right[[3]],3], tru[[2]][keep_right[[3]],3])
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

subj_cor_df <-
  reshape2::melt(bayes_cor) %>%
  mutate(Model = "Bayesian") %>%
  full_join(
    reshape2::melt(classical_cor) %>%
      mutate(Model = "Classical")
  ) %>%
  mutate(variable = sub("_"," ", variable)) %>%
  pivot_wider(names_from = Model, values_from = value) %>%
  mutate(diff = Bayesian - Classical) #%>% #,
         # diff = ifelse(diff < 0, 0, diff),
         # diff = ifelse(diff > 0.1, 0.1, diff))
  # filter(variable == "tongue") %>%

subj_cor_plot <- ggplot(subj_cor_df) +
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
  labs(y = "Bayesian GLM",
       x = "Classical GLM"
       # title = "Test-Retest Correlation"
       ) +
  facet_wrap(~variable, scales = "free") +
  guides(color = 'none') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.title = element_text(size = 10))
  # theme(text = element_text(size = 14),
  #       plot.title = element_text(hjust = 0.5))

subj_cor_plot

ggsave("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/05_subject_cor_plot.png", plot = subj_cor_plot, width = 7.5, height = 5)
# ggsave("~/Desktop/SMI 2021 poster images/5_subject_cor_plot.png", plot = subj_cor_plot, width = 3.5, height = 4)

# >>>> Let's see if a hypothesis test would work here ----
# Let's see if the differences look normal
reshape2::melt(bayes_cor) %>%
  mutate(Model = "Bayesian") %>%
  full_join(
    reshape2::melt(classical_cor) %>%
      mutate(Model = "Classical")
  ) %>%
  mutate(variable = sub("_"," ", variable)) %>%
  pivot_wider(names_from = Model, values_from = value) %>%
  mutate(diff = Bayesian - Classical) %>%
  ggplot(aes(x = diff)) +
  geom_histogram(bins = 100) +
  facet_wrap(~variable)
# The differences **do** look relatively normal

reshape2::melt(bayes_cor) %>%
  mutate(Model = "Bayesian") %>%
  full_join(
    reshape2::melt(classical_cor) %>%
      mutate(Model = "Classical")
  ) %>%
  mutate(variable = sub("_"," ", variable)) %>%
  pivot_wider(names_from = Model, values_from = value) %>%
  mutate(diff = Bayesian - Classical) %>%
  group_by(variable) %>%
  summarize(p_values = t.test(Bayesian,Classical, paired = T)$p.value,
            Bayes_minus_Classical = t.test(Bayesian,Classical, paired = T)$estimate)

# >> Combine MSE and Correlation in one plot ----
mse_cor_df <-
  mutate(subj_mse_df, metric = "MSE") %>%
  full_join(
    mutate(subj_cor_df, metric = "Correlation")
  ) %>%
  mutate(diff = Bayesian - Classical) %>%
  group_by(metric, variable) %>%
  mutate(p_value = t.test(x = Bayesian, y = Classical, paired = TRUE)$p.value,
         significant = ifelse(6*p_value < 0.05,"Yes","No"))

mse_cor_plot <- ggplot(mse_cor_df) +
  geom_boxplot(aes(x = variable, y = diff, fill = significant), outlier.alpha = 0) +
  geom_jitter(aes(x = variable, y = diff)) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  labs(x = "", y = "Difference\n(Bayesian - Classical)") +
  scale_fill_manual("Sig. Diff.", values = c("white","forestgreen")) +
  facet_wrap(~metric, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
mse_cor_plot
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
ggsave(filename = file.path(plot_dir,"05_subject_mse_cor_plot.png"),
       plot = mse_cor_plot, width = 8, height = 5)

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
    cifti_obj <- readRDS("HCP_data/cifti_5k_template_whole.rds")
    cifti_obj$data$cortex_left <- left_act
    cifti_obj$data$cortex_right <- right_act
    if(task_idx %in% c(1,4)) {
      plot(cifti_obj, color_mode = "qualitative", colors = col_pal,
           # title = paste("Subject",subject, "Tongue Task"),
           surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
           surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
           fname = paste0('plots/5_bayes_',
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
             fname = paste0('plots/5_bayes_',
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
results_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
load("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
result_files <- list.files(results_dir,full.names = T) |> grep(pattern = "FWHM6", value = T)
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
    # left_FDR <- sapply(c(0, 0.5, 1), function(thr) {
    #   out <-
    #     id_activations.classical(
    #       left_result$GLMs_classical$cortexL,
    #       alpha = 0.01,
    #       correction = "FDR",
    #       field_inds = task_idx,
    #       threshold = thr
    #     )$active
    #   return(out)
    # }, simplify = T)
    # left_FWER <- left_FDR
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
    # right_FDR <- sapply(c(0, 0.5, 1), function(thr) {
    #   out <-
    #     id_activations.classical(
    #       right_result$GLMs_classical$cortexR,
    #       alpha = 0.01,
    #       correction = "FDR",
    #       field_inds = task_idx,
    #       threshold = thr
    #     )$active
    #   return(out)
    # }, simplify = T)
    # right_FWER <- right_FDR
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
    # right_active <- right_FDR + right_FWER
    right_active <- right_FWER
    medial_right <- !is.na(right_active[,1])
    right_active <- as.matrix(right_active[,1][!is.na(right_active[,1])])
    right_active[,1] <- ifelse(right_active[,1] == 0,NA,right_active[,1])
    cifti_active <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/603_cifti_5k_template_whole.rds")
    cifti_active$data$cortex_left <- left_active
    cifti_active$data$cortex_right <- right_active
    if(task_idx %in% c(1,4)) {
      plot(cifti_active, color_mode = "qualitative", colors = col_pal,
           surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
           surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
           fname = paste0('/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/607_classical_',
                          subject,
                          "_visit1_",
                          task_names[task_idx],
                          "_activations.png"))
    }
    if(task_idx %in% c(2,3)) {
      for(hem in c('left','right')) {
        plot(cifti_active, hemisphere = hem,
             color_mode = "qualitative", colors = col_pal,
             surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
             surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
             fname = paste0('/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/607_classical_',
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

# >> Classical (unsmoothed) ----
library(ciftiTools)
ciftiTools.setOption("wb_path","/Applications/workbench")
library(BayesfMRI)
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results_EM"
result_files <- list.files(result_dir, full.names = T) |>
  grep(pattern = "visit1", value = T) |>
  grep(pattern = "left", value = T)
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
col_pal <- c(light_orange,"red","purple")
load("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
for(subject in subjects) {
  subject_file <- grep(subject, result_files, value = T)
  result_obj <- readRDS(subject_file)
  result_active <- list()
  for(thr in c(0,0.5,1)) {
    thr_chr <- as.character(thr)
    result_active[[thr_chr]] <- id_activations_cifti(result_obj,alpha = 0.01,method = "classical",correction = "FWER",threshold = thr)$activations_xifti
  }
  result_active <- Reduce(`+`,result_active)
  result_active$data$cortex_left[result_active$data$cortex_left == 0] <- NA
  plot(result_active, idx = 4, hemisphere = "left", view = "lateral",
       colors = col_pal, color_mode = "qualitative",
       fname = file.path(plot_dir, paste0("607_classical_unsmoothed_", subject, "_visit1_tongue_activations.png")))
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
              subject_cifti <- readRDS("HCP_data/cifti_5k_template_whole.rds")
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
             fname = paste0("plots/5_bayes_",subject,"_",task_names[task_idx],"_thr",sub("\\.","",thresh),"_",alpha_chr,"_num_active.png"),
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
               fname = paste0("plots/5_bayes_",subject,"_",hem,"_",task_names[task_idx],"_thr",sub("\\.","",thresh),"_",alpha_chr,"_num_active.png"),
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
            subject_cifti <- readRDS("HCP_data/cifti_5k_template_whole.rds")
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
        "plots/subject_",
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


# FIGURE 7: Dice Overlap CIs and Dice x Area scatterplot ----
library(tidyverse)
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

dice_result_files <- list.files("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results", full.names = T)
dice_result_files <- grep("6_", dice_result_files, value = T)
dice_result_files <- grep("_PW_", dice_result_files, value = T)
dice_result_files <- grep("_Dice_", dice_result_files, value = T)
dice_files <- list(Bayesian = grep("6_classical",dice_result_files, invert = T, value = T)[c(2,1,3)],
                   Classical = grep("06_classical_FWER",dice_result_files, value = T)[c(2,1,3)])
dice_results <- sapply(dice_files, function(df) sapply(df,readRDS, simplify = F), simplify = F)
task_names <- c("visual cue","tongue","right foot","right hand","left foot","left hand")
names(dice_results$Bayesian) <- c("0%","0.5%","1%")
names(dice_results$Classical) <- c("0%","0.5%","1%")

# This is a comparison of the Dice values for the Bayesian 0.5% and Classical 0%
dice_ls <-
  list(
    Bayesian = dice_results$Bayesian$`0.5%`,
    Classical = dice_results$Classical$`0%`
  )

dice_tests <- mapply(function(B,C) {
  test_out <- t.test(x = B, y = C, paired = TRUE)
  return(list(
    Estimate = test_out$estimate,
    P = test_out$p.value
  ))
}, B = split(dice_ls$Bayesian, col(dice_ls$Bayesian)),
C = split(dice_ls$Classical, col(dice_ls$Classical)), SIMPLIFY = F)
# dice_tests <- as.matrix(dice_tests)
# colnames(dice_tests) <- task_names
# dice_tests
dt_df <-
  reshape2::melt(dice_tests) %>%
  rename(Var2 = L1) %>%
  mutate(Var2 = task_names[as.numeric(Var2)]) %>%
  pivot_wider(names_from = L2, values_from = value) %>%
  mutate(corrected_P = 6*P,
         greater_than = corrected_P < 0.05,
         P = round(P,3)) %>%
  select(-Estimate)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
light_teal <- colorRampPalette(c("white",gg_color_hue(2)[2]))(5)[2]

col_pal <- c("grey",light_teal,gg_color_hue(2)[2])

dimnames(dice_ls[[1]])[[1]] <- seq(45)
dimnames(dice_ls[[2]])[[1]] <- seq(45)

subject_dice_plot <- reshape2::melt(dice_ls) %>%
  mutate(Var2 = task_names[Var2]) %>%
  pivot_wider(names_from = L1, values_from = value) %>%
  mutate(diff = Bayesian - Classical) %>%
  left_join(dt_df) %>%
  mutate(p_cat = cut(P, breaks = c(-Inf,0.01,0.05,1)),
         p_lab = ifelse(P == 0, "p < 0.001",paste0("p = ",P))) %>%
  ggplot() +
  geom_boxplot(aes(x = Var2, y = diff, fill = p_cat), outlier.alpha = 0) +
  geom_jitter(aes(x = Var2, y = diff), alpha = 0.3) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  geom_text(aes(x = Var2, y = -0.35, label = p_lab)) +
  scale_fill_manual("Paired t-test", values = rev(col_pal), labels = rev(c("p > 0.05","0.01 < p < 0.05","p < 0.01"))) +
  labs(x = "", y = "Difference in Dice Coefficient\n(Bayesian - Classical)") +
  theme_classic() +
  theme(text = element_text(size = 16))
subject_dice_plot
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
ggsave(filename = file.path(plot_dir,"05_dice_paired_difference.png"), width = 9, height = 6)

bstrap_ci_results <- sapply(dice_results, function(method_results) {
  sapply(method_results, function(thr_results) {
    apply(thr_results,2,bstrap_CIs)
  },simplify = F)
}, simplify = F)

library(tidyverse)
classical_standard <- reshape2::melt(bstrap_ci_results) %>%
  mutate(task = task_names[Var2],
         threshold = factor(L2, levels = c("0%","0.5%","1%"))) %>%
  pivot_wider(names_from = Var1, values_from = value) %>%
  filter(L1 == "Classical", threshold == "0%")

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
  guides(color = guide_legend(nrow = 2)) +
  scale_color_manual("", values = rev(gg_color_hue(2)), labels = c("Bayesian (unsmoothed)","Classical (FWHM = 6mm)")) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,-10,-10,-10))

dice_ci_plot

ggsave("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/05_dice_ci_plot.png", width = 3, height = 7, plot = dice_ci_plot)

# >> Dice x Area Scatterplot ----

area_result_files <- list.files("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results", full.names = T)
area_result_files <- grep("6_", area_result_files, value = T)
area_result_files <- grep("_PW_", area_result_files, value = T)
area_result_files <- grep("_overlap_", area_result_files, value = T)
area_files <- list(Bayesian = grep("classical",area_result_files, invert = T, value = T)[c(2,1,3)],
                   Classical = grep("06_classical_FWER",area_result_files, value = T)[c(2,1,3)])
area_results <- sapply(area_files, function(df) sapply(df,readRDS, simplify = F), simplify = F)
names(area_results$Bayesian) <- c("0%","0.5%","1%")
names(area_results$Classical) <- c("0%","0.5%","1%")

# Used for the OHBM 2021 poster
color_pal <- c(
  rgb(39,86,197, maxColorValue = 255),
  rgb(169,99,36, maxColorValue = 255)
)

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

# This is just for the tongue task
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
# End

area_by_dice <-
  reshape2::melt(dice_results, value.name = "dice") %>%
  cbind(overlap = reshape2::melt(area_results, value.name = "overlap")$overlap) %>%
  mutate(Var2 = task_names[Var2],
         Var2 = as.factor(Var2),
         label_gamma = as.numeric(sub("%","",L2)),
         L2 = factor(L2, levels = c("0%","0.5%","1%"))) %>%
  # filter(Var2 == "tongue") %>%
  ggplot() +
  geom_point(aes(x = overlap, y = dice, color = L1)) +
  facet_grid(Var2 ~ L2, scales = "free_x") +
  # facet_grid(. ~ label_gamma, scales = "free_x", labeller = label_bquote(cols = gamma == .(label_gamma)~"%")) +
  # scale_color_discrete("") +
  scale_color_manual("", values = rev(gg_color_hue(2)),
                     labels = c("Bayesian (unsmoothed)","Classical (FWHM = 6mm)")) +
  scale_x_continuous(breaks = breaks_fun, limits = c(0,NA))+
  labs(y = "Dice Coefficient", x = "Size of Overlap (# of vertices)") +
  guides(color=guide_legend(nrow = 2)) +
  theme_bw() +
  theme(legend.position = "top", panel.grid = element_blank(),
        # legend.text = element_text(size = 12),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,-10,-10,-10)
        )

area_by_dice

ggsave("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/05_area_by_dice.png", width = 4.5, height = 7, plot = area_by_dice)
# ggsave("~/Desktop/SMI 2021 poster images/5_area_by_dice.png", width = 6, height = 2.5, plot = area_by_dice)

# FIGURE 8: Group Estimates and Activations ----
# >> Estimates ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench/')
task_file_names <- c("visual_cue","foot","hand","tongue")
# >>>> Bayesian ----
group_results_left <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples/501_HCP_BayesGLM2_45subj_visit1_result_left_thresh0_20210318.rds")
group_results_right <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples/501_HCP_BayesGLM2_45subj_visit1_result_right_thresh0_20210319.rds")
cifti_both <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/603_cifti_5k_template_whole.rds")
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
# group_results_classical <- readRDS("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed/group/03_HCP_classical_group_PW_estimates_visit1.rds")
classical_result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed/"
classical_result_files <- list.files(classical_result_dir, full.names = T) |>
  grep(pattern = "classical", value = T) |>
  grep(pattern = "visit1", value = T) |>
  grep(pattern = "sessionLR", value = T, invert = T) |>
  grep(pattern = "FWHM6", value = T)
left_files <- grep(pattern = "left", classical_result_files, value = T)
right_files <- grep(pattern = "right", classical_result_files, value = T)
library(BayesfMRI)
group_left <-
  classicalGLM2(
    results = left_files,
    brainstructure = "cortexL",
    gamma = c(0, 0.5, 1),
    correction = "FWER",
    alpha = 0.01
  )
group_right <-
  classicalGLM2(
    results = right_files,
    brainstructure = "cortexR",
    gamma = c(0, 0.5, 1),
    correction = "FWER",
    alpha = 0.01
  )
cifti_classical <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/603_cifti_5k_template_whole.rds")
cifti_classical$data$cortex_left <- group_left$avg_estimate
cifti_classical$data$cortex_right <- group_right$avg_estimate
# >>>>>>  Visual Cue and Tongue (idx = c(1,4)) ----
for(i in c(1,4)) {
  plot(cifti_classical, zlim = c(-1,1), idx = i, hemisphere = "both",
       legend_embed = F,
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
       fname = paste0("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/607_group_classical_",task_file_names[i],"_estimate.png"))
}

# >>>>>> Foot and Hand Tasks (idx = c(2,3)) ----
for(i in c(2,3)) {
  for(hem in c("left","right")) {
    plot(cifti_classical, zlim = c(-1,1), idx = i, hemisphere = hem,
         legend_embed = F,
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         fname = paste0("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/607_group_classical_",hem,"_",task_file_names[i],"_estimate.png"))
  }
}

# >>>> Classical (unsmoothed) ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
group_results_classical <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/group/502_HCP_classical_group_estimates.rds")
cifti_classical <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/603_cifti_5k_template_whole.rds")
cifti_classical$data$cortex_left <- as.matrix(group_results_classical$left)
plot(cifti_classical, hemisphere = "left", view = "lateral", zlim = c(-1,1),
     idx = 4, legend_embed = F,
     fname = file.path(plot_dir, "607_group_classical_unsmoothed_tongue_estimate.png"))

# >> Activations ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
# >>>> Bayesian ----
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples"
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
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
  # cifti_obj <- readRDS("HCP_data/cifti_5k_template_whole.rds")
  cifti_obj <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/603_cifti_5k_template_whole.rds")
  cifti_obj$data$cortex_left <- left_act
  cifti_obj$data$cortex_right <- right_act
  if(task_idx %in% c(1,4)) {
    plot(cifti_obj, color_mode = "qualitative", colors = col_pal, hemisphere = "both",
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         fname = file.path(plot_dir,paste0('607_group_bayes_',task_names[task_idx],"_activations.png")))
  }
  if(task_idx %in% c(2,3)) {
    for(hem in c("left","right")) {
      plot(cifti_obj, color_mode = "qualitative", colors = col_pal, hemisphere = hem,
           surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
           surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
           fname = file.path(plot_dir,paste0('607_group_bayes_',hem,"_",task_names[task_idx],"_activations.png")))
    }
  }
}


# >>>> Classical ----
# As a note, this would not be necessary in the case of the current (2021-02-10)
# state of the BayesfMRI package, but some retesting using the classical GLM
# is necessary in order to obtain the standard errors of the estimates.
# results_dir <- "HCP_results/5k_results/group"
# result_files <- list.files(results_dir,full.names = T)
# FDR_result <- readRDS(grep("classical_activations_PW_FDR", result_files, value= T))
# FWER_result <- readRDS("HCP_results/5k_results/group/502_HCP_classical_activations_PW_visit1_FWER.rds")
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
classical_result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed/"
classical_result_files <- list.files(classical_result_dir, full.names = T) |>
  grep(pattern = "classical", value = T) |>
  grep(pattern = "visit1", value = T) |>
  grep(pattern = "sessionLR", value = T, invert = T) |>
  grep(pattern = "FWHM6", value = T)
left_files <- grep(pattern = "left", classical_result_files, value = T)
right_files <- grep(pattern = "right", classical_result_files, value = T)
library(BayesfMRI)
group_left <-
  classicalGLM2(
    results = left_files,
    brainstructure = "cortexL",
    gamma = c(0, 0.5, 1),
    correction = "FWER",
    alpha = 0.01
  )
group_right <-
  classicalGLM2(
    results = right_files,
    brainstructure = "cortexR",
    gamma = c(0, 0.5, 1),
    correction = "FWER",
    alpha = 0.01
  )
left_active <- sapply(group_left$active_result, function(x) x$active, simplify = F) |>
  Reduce(f = `+`)
left_active[left_active == 0] <- NA
right_active <- sapply(group_right$active_result, function(x) x$active, simplify = F) |>
  Reduce(f = `+`)
right_active[right_active == 0] <- NA
cifti_classical <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/603_cifti_5k_template_whole.rds")
library(BayesfMRI)
library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path',"/Applications/workbench")
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
col_pal <- c(light_orange,"red","purple")
cifti_classical$data$cortex_left <- left_active
cifti_classical$data$cortex_right <- right_active
task_names <- c("visual_cue","foot","hand","tongue")
for(task_idx in c(1,4)) {
  plot(cifti_classical, hemisphere = "both", idx = task_idx,
       color_mode = "qualitative", colors = col_pal,
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
       fname = file.path(plot_dir,paste0('/607_group_classical_',task_names[task_idx],"_activations.png")))
}

for(task_idx in c(2,3)) {
  for(hem in c("left","right")) {
    plot(cifti_classical, hemisphere = hem, idx = task_idx,
         color_mode = "qualitative", colors = col_pal,
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         fname = file.path(plot_dir,paste0('607_group_classical_',hem,"_",task_names[task_idx],"_activations.png")))
  }
}

# >>>> Classical (unsmoothed) ----
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
FWER_result <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/group/502_HCP_classical_activations_FWER.rds")
FWER_result <- FWER_result$left
FWER_result <- Reduce(`+`, FWER_result)
FWER_result[FWER_result == 0] <- NA
library(BayesfMRI)
library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path',"/Applications/workbench")
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
# col_pal <- c("yellow",light_orange,"red","purple")
col_pal <- c(light_orange,"red","purple")
cifti_unsmoothed <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/603_cifti_5k_template_whole.rds")
cifti_unsmoothed <- remove_xifti(cifti_unsmoothed, remove = "cortex_right")
cifti_unsmoothed$data$cortex_left <- FWER_result
plot(cifti_unsmoothed, colors = col_pal, color_mode = "qualitative",
     hemisphere = 'left', view = "lateral", idx = 4,
     fname = file.path(plot_dir,"607_group_classical_unsmoothed_tongue_activations.png"))

# FIGURE 9: Group Test-Retest MSE and Correlation ----
# >> MSE ----
library(tidyverse)
truth <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/502_HCP_classical_group_PW_estimates_visit2.rds")
classical_est <- readRDS("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed/group/03_HCP_classical_group_PW_estimates_visit1.rds")
classical_subsample_est <- readRDS("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed/group/03_HCP_classical_subsample_estimates.rds")
classical_est_unsmoothed <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/502_HCP_classical_group_PW_estimates_visit1.rds")
classical_subsample_est_unsmoothed <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/group/502_HCP_classical_subsample_estimates.rds")
bayes_est <- list(
  left = readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/group/501_HCP_BayesGLM2_45subj_result_left_thresh0_20210317.rds")$estimates,
  right = readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/group/501_HCP_BayesGLM2_45subj_result_right_thresh0_20210318.rds")$estimates
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

mse_classical_unsmoothed <- mapply(function(est_hem,tru_hem){
  mapply(function(est,tru) {
    mean((est - tru)^2)
  }, est = split(est_hem, col(est_hem)), tru = split(tru_hem, col(tru_hem)))
}, est_hem = classical_est_unsmoothed, tru_hem = truth)

rownames(mse_bayes) <- rownames(mse_classical) <- rownames(mse_classical_unsmoothed) <- c("cue","foot","hand","tongue")

mse_bayes
mse_classical
mse_classical_unsmoothed

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

mse_classical_subsample_unsmoothed <- sapply(classical_subsample_est_unsmoothed, function(size_ests) {
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
         model = "Classical (FWHM = 6mm)",
         num_subjects = "45") #%>%
  # filter(task == "tongue")

mse_classical_unsmoothed_df <-
  reshape2::melt(mse_classical_unsmoothed) %>%
  pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
  mutate(tongue = (left_tongue + right_tongue)/2,
         cue = (left_cue + right_cue)/2) %>%
  select(-ends_with("_tongue"), -ends_with("_cue")) %>%
  pivot_longer(cols = everything() ,names_to = "task", values_to = "MSE") %>%
  mutate(task = sub("_"," ", task),
         model = "Classical (unsmoothed)",
         num_subjects = "45") #%>%
# filter(task == "tongue")

mse_bayes_df <-
  reshape2::melt(mse_bayes) %>%
  pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
  mutate(tongue = (left_tongue + right_tongue)/2,
         cue = (left_cue + right_cue)/2) %>%
  select(-ends_with("_tongue"), -ends_with("_cue")) %>%
  pivot_longer(cols = everything() ,names_to = "task", values_to = "MSE") %>%
  mutate(task = sub("_"," ", task),
         model = "Bayesian (unsmoothed)",
         num_subjects = "45") #%>%
  # filter(task == "tongue")

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
         model = "Classical (FWHM = 6mm)")

classical_mse_subsample_unsmoothed_df <- reshape2::melt(mse_classical_subsample_unsmoothed) %>%
  mutate(num_subjects = sub("subjects","", L1),
         sample_num = sub("sample","",L2)) %>%
  select(-L1, -L2) %>%
  pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
  mutate(tongue = (left_tongue + right_tongue)/2,
         cue = (left_cue + right_cue)/2) %>%
  select(-ends_with("_tongue"), -ends_with("_cue")) %>%
  pivot_longer(-c(num_subjects,sample_num),names_to = "task", values_to = "MSE") %>%
  mutate(task = sub("_"," ", task),
         model = "Classical (unsmoothed)")

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
         model = "Bayesian (unsmoothed)")

median_mse_df <-
  full_join(classical_mse_subsample_df,classical_mse_subsample_unsmoothed_df) %>%
  full_join(bayes_mse_subsample_df) %>%
  group_by(model,task,num_subjects) %>%
  summarize(MSE = median(MSE)) %>%
  full_join(mse_bayes_df) %>%
  full_join(mse_classical_df) %>%
  full_join(mse_classical_unsmoothed_df)
  # filter(task == "tongue")

# Used for the OHBM 2021 poster
# color_pal <- c(
#   rgb(39,86,197, maxColorValue = 255),
#   rgb(169,99,36, maxColorValue = 255)
# )

library(RColorBrewer)
my_pal <- brewer.pal(6, "Paired")[c(6,5,2)]

mse_plot <-
  full_join(classical_mse_subsample_df,classical_mse_subsample_unsmoothed_df) %>%
  full_join(bayes_mse_subsample_df) %>%
  mutate(model = factor(model, levels = c("Classical (unsmoothed)","Classical (FWHM = 6mm)", "Bayesian (unsmoothed)"))) %>%
  # filter(task == "tongue") %>%
  ggplot() +
  # geom_boxplot(aes(x = task, y = MSE, fill = num_subjects, color = num_subjects)) +
  geom_jitter(aes(x = num_subjects, y = MSE, color = model), width = 0.1, height = 0, alpha = 0.3) +
  geom_line(aes(x = as.numeric(factor(num_subjects)), y = MSE, color = model), data = median_mse_df) +
  geom_point(aes(x = num_subjects, y = MSE, color = model), data = mse_classical_unsmoothed_df, size = 4) +
  geom_point(aes(x = num_subjects, y = MSE, color = model), data = mse_classical_df, size = 4) +
  geom_point(aes(x = num_subjects, y = MSE, color = model), data = mse_bayes_df, size = 4) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  # geom_segment(aes(x = as.numeric(factor(num_subjects))-0.3,xend = as.numeric(factor(num_subjects))+0.3, y = MSE,yend = MSE, color = model), data = median_mse_df) +
  labs(x = "Subsample Size", y = "Test-Retest MSE") +
  scale_color_discrete("", type = my_pal)+
  # scale_color_manual("",values = color_pal, labels = c("Bayesian","Classical")) +
  facet_wrap(~task, scales = "free_y") +
  # guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_classic() +
  # guides(color = FALSE) +
  theme(legend.position = "top",
        # legend.margin = margin(0,0,-5,-35),
        text = element_text(size = 16)
        )

mse_plot
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
ggsave(file.path(plot_dir,"607_group_mse.png"),plot = mse_plot, width = 8, height = 5.5)
 # ggsave("~/Desktop/SMI 2021 poster images/5_group_mse.png",plot = mse_plot, width = 2.5, height = 3.5)

# >> Correlation ----
library(tidyverse)
# Calculate the correlation for the Group estimates
truth <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/502_HCP_classical_group_PW_estimates_visit2.rds")
classical_est <- readRDS("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed/group/03_HCP_classical_group_PW_estimates_visit1.rds")
classical_subsample_est <- readRDS("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed/group/03_HCP_classical_subsample_estimates.rds")
classical_est_unsmoothed <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/502_HCP_classical_group_PW_estimates_visit1.rds")
classical_subsample_est_unsmoothed <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/group/502_HCP_classical_subsample_estimates.rds")
bayes_est <- list(
  left = readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/group/501_HCP_BayesGLM2_45subj_result_left_thresh0_20210317.rds")$estimates,
  right = readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/group/501_HCP_BayesGLM2_45subj_result_right_thresh0_20210318.rds")$estimates
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

cor_classical_unsmoothed <- mapply(function(est_hem,tru_hem){
  mapply(function(est,tru) {
    cor(est,tru)
  }, est = split(est_hem, col(est_hem)), tru = split(tru_hem, col(tru_hem)))
}, est_hem = classical_est_unsmoothed, tru_hem = truth)

cor_bayes
cor_classical
cor_classical_unsmoothed

rownames(cor_bayes) <- rownames(cor_classical) <- rownames(cor_classical_unsmoothed) <- c("cue","foot","hand","tongue")

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

cor_classical_unsmoothed_subsample <- sapply(classical_subsample_est_unsmoothed, function(size_ests) {
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
         model = "Classical (FWHM = 6mm)",
         num_subjects = "45")

cor_classical_unsmoothed_df <-
  reshape2::melt(cor_classical_unsmoothed) %>%
  pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
  mutate(tongue = (left_tongue + right_tongue)/2,
         cue = (left_cue + right_cue)/2) %>%
  select(-ends_with("_tongue"), -ends_with("_cue")) %>%
  pivot_longer(cols = everything() ,names_to = "task", values_to = "corr") %>%
  mutate(task = sub("_"," ", task),
         model = "Classical (unsmoothed)",
         num_subjects = "45")

cor_bayes_df <-
  reshape2::melt(cor_bayes) %>%
  pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
  mutate(tongue = (left_tongue + right_tongue)/2,
         cue = (left_cue + right_cue)/2) %>%
  select(-ends_with("_tongue"), -ends_with("_cue")) %>%
  pivot_longer(cols = everything() ,names_to = "task", values_to = "corr") %>%
  mutate(task = sub("_"," ", task),
         model = "Bayesian (unsmoothed)",
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
         model = "Classical (FWHM = 6mm)")

classical_cor_subsample_unsmoothed_df <- reshape2::melt(cor_classical_unsmoothed_subsample) %>%
  mutate(num_subjects = sub("subjects","", L1),
         sample_num = sub("sample","",L2)) %>%
  select(-L1, -L2) %>%
  pivot_wider(names_from = c(Var2,Var1), values_from = value) %>%
  mutate(tongue = (left_tongue + right_tongue)/2,
         cue = (left_cue + right_cue)/2) %>%
  select(-ends_with("_tongue"), -ends_with("_cue")) %>%
  pivot_longer(-c(num_subjects,sample_num),names_to = "task", values_to = "corr") %>%
  mutate(task = sub("_"," ", task),
         model = "Classical (unsmoothed)")

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
         model = "Bayesian (unsmoothed)")

median_cor_df <- full_join(classical_cor_subsample_df,classical_cor_subsample_unsmoothed_df) %>%
  full_join(bayes_cor_subsample_df) %>%
  group_by(model,task,num_subjects) %>%
  summarize(corr = median(corr)) %>%
  full_join(cor_bayes_df) %>%
  full_join(cor_classical_df) %>%
  full_join(cor_classical_unsmoothed_df)

library(RColorBrewer)
my_pal <- brewer.pal(6, "Paired")[c(6,5,2)]

cor_plot <- full_join(classical_cor_subsample_df,classical_cor_subsample_unsmoothed_df) %>%
  full_join(bayes_cor_subsample_df) %>%
  mutate(model = factor(model, levels = c("Classical (unsmoothed)","Classical (FWHM = 6mm)","Bayesian (unsmoothed)"))) %>%
  ggplot() +
  # geom_boxplot(aes(x = task, y = MSE, fill = num_subjects, color = num_subjects)) +
  geom_jitter(aes(x = num_subjects, y = corr, color = model), width = 0.1, height = 0, alpha = 0.3) +
  geom_line(aes(x = as.numeric(factor(num_subjects)), y = corr, color = model), data = median_cor_df) +
  geom_point(aes(x = num_subjects, y = corr, color = model), data = cor_classical_df, size = 4) +
  geom_point(aes(x = num_subjects, y = corr, color = model), data = cor_classical_unsmoothed_df, size = 4) +
  geom_point(aes(x = num_subjects, y = corr, color = model), data = cor_bayes_df, size = 4) +
  geom_hline(aes(yintercept = 1), lty = 2) +
  # geom_segment(aes(x = as.numeric(factor(num_subjects))-0.3,xend = as.numeric(factor(num_subjects))+0.3, y = MSE,yend = MSE, color = model), data = median_mse_df) +
  labs(x = "Subsample Size", y = "Correlation") +
  scale_color_discrete("", type = my_pal) +
  # scale_fill_manual(name = "Subsample Size", values = my_pal) +
  # scale_color_manual(name = "Subsample Size", values = my_pal) +
  facet_wrap(~task, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "top",
        text = element_text(size = 16))
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
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
ggsave(file.path(plot_dir,"607_group_cor.png"),plot = cor_plot, width = 8, height = 5.5)

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
# result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/classical"
result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
library(abind)
result_files <- list.files(result_dir, full.names = T) |>
  grep(pattern = "_visit1_", value = T) |>
  grep(pattern = "FWHM6", value = T) |>
  grep(pattern = "classical", value = T) |>
  grep(pattern = "sessionLR", value = T, invert = T)
# result_files <- grep("_visit1_", result_files, value = T)
# load("HCP_data/subjects.Rdata")
load("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/subjects.Rdata")
# classical_estimates <- sapply(c('left','right'), function(hem) {
#   sapply(subjects, function(subject) {
#     filen <- grep(paste0(subject, "_visit1_",hem), result_files, value = T)
#     if(length(filen) > 1) filen = filen[1]
#     result_obj <- readRDS(filen)
#     out <-
#       abind(result_obj$betas_classical$LR$data[[paste0("cortex_", hem)]],
#             result_obj$betas_classical$RL$data[[paste0("cortex_", hem)]],
#             along = 3)
#     return(out)
#   }, simplify = "array")
# }, simplify = F)

subsample_subjects <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples/501_sample_subjects_20210319.rds")
subsample_subjects[[4]] <- as.matrix(subjects)
task_names <- c("cue","foot","hand","tongue")
library(tidyverse)

class_num_act <- sapply(subsample_subjects, function(samp_n) {
  sub_size <- sapply(seq(ncol(samp_n)), function(samp_num) {
    files <- sapply(samp_n[,samp_num], function(x) grep(pattern =x, result_files,value = T))
    group_left <-
      classicalGLM2(
        results = grep(pattern = "left", files, value = T),
        brainstructure = "cortexL",
        gamma = c(0, 0.5),
        correction = "FWER",
        alpha = 0.01
      )
    group_right <-
      classicalGLM2(
        results = grep(pattern = "right", files, value = T),
        brainstructure = "cortexR",
        gamma = c(0, 0.5),
        correction = "FWER",
        alpha = 0.01
      )
    num_active <- sapply(list(group_left, group_right), function(hem) {
      sapply(hem$active_result, function(gam) { as.matrix(apply(gam$active, 2, sum))}, simplify = F)
    }, simplify = F)
    names(num_active) <- c("left","right")
    sum_active <- reshape2::melt(num_active) %>%
      mutate(task = task_names[Var1],
             gamma = as.numeric(sub("gamma = ","",L2)),
             other_hem = ifelse(L1 == "left","right","left"),
             task = ifelse(task == "foot" | task =="hand", paste(other_hem,task),task)) %>%
      group_by(gamma, task) %>%
      summarize(value = sum(value)) %>%
      mutate(samp_num = samp_num,
             num_subjects = nrow(samp_n)) %>%
      ungroup()
    return(sum_active)
  }, simplify = F)
  sub_size_df <- Reduce(full_join, sub_size)
  return(sub_size_df)
}, simplify = F)
class_act_df <- Reduce(full_join,class_num_act)

# load("HCP_data/subjects.Rdata")
# subsample_idx <- sapply(subsample_subjects, function(sam){
#   out <- apply(sam,1:2, function(x) which(subjects == x))
#   return(out)
# }, simplify = F)
# threshs <- c(0,0.5,1)
# class_num_act <- sapply(subsample_idx, function(ssize) {
#   sample_num_active <- sapply(1:10,function(samp) {
#     sapply(c('left','right'), function(hem) {
#       samp_est <- classical_estimates[[hem]][,,,ssize[,samp]]
#       num_locs <- dim(samp_est)[1]
#       bonferroni_cutoff <- 0.01 / num_locs # alpha = 0.01
#       data_df <- reshape2::melt(samp_est)
#       vertex_lists <- split(data_df,data_df$Var1)
#       thresh_active <- sapply(c(0,0.5,1), function(thr) {
#         active_all <- sapply(vertex_lists, function(v_df){
#           vt_df <- split(v_df, v_df$Var2)
#           active_vertex <- sapply(vt_df, function(vt) {
#             ttest_res <- t.test(x = vt$value, mu = thr, alternative = "greater")
#             active_v <- ifelse(ttest_res$p.value < bonferroni_cutoff, 1, 0)
#             return(active_v)
#           }, simplify = "array")
#           return(active_vertex)
#         }, simplify = 'array')
#         out <- apply(active_all,1,sum)
#         return(as.matrix(out))
#       }, simplify = F)
#       names(thresh_active) <- paste0(threshs,"%")
#       return(thresh_active)
#     }, simplify = F)
#   }, simplify = F)
#   names(sample_num_active) <- as.character(1:10)
#   return(sample_num_active)
# }, simplify = F)
# names(class_num_act) <- c("10","20","30")

# class_num_act_45 <- sapply("45", function(ssize) {
#   sample_num_active <- sapply(1,function(samp) {
#     sapply(c('left','right'), function(hem) {
#       samp_est <- classical_estimates[[hem]]
#       num_locs <- dim(samp_est)[1]
#       bonferroni_cutoff <- 0.01 / num_locs # alpha = 0.01
#       data_df <- reshape2::melt(samp_est)
#       vertex_lists <- split(data_df,data_df$Var1)
#       thresh_active <- sapply(c(0,0.5,1), function(thr) {
#         active_all <- sapply(vertex_lists, function(v_df){
#           vt_df <- split(v_df, v_df$Var2)
#           active_vertex <- sapply(vt_df, function(vt) {
#             ttest_res <- t.test(x = vt$value, mu = thr, alternative = "greater")
#             active_v <- ifelse(ttest_res$p.value < bonferroni_cutoff, 1, 0)
#             return(active_v)
#           }, simplify = "array")
#           return(active_vertex)
#         }, simplify = 'array')
#         out <- apply(active_all,1,sum)
#         return(as.matrix(out))
#       }, simplify = F)
#       names(thresh_active) <- paste0(threshs,"%")
#       return(thresh_active)
#     }, simplify = F)
#   }, simplify = F)
#   names(sample_num_active) <- "1"
#   return(sample_num_active)
# }, simplify = F)

library(tidyverse)
classical_num_act <-
  reshape2::melt(class_num_act) %>%
  full_join(reshape2::melt(class_num_act_45)) %>%
  mutate(model = "Classical GLM")

# >> Plot ----
col_pal <- scales::hue_pal()(5)
library(viridis)
col_pal <- viridis(7,option = "C")[1:5]
col_pal <- rgb(c(230, 86, 0, 0,204),c(159,180,158,114,121),c(0,233,115,178,167), maxColorValue = 255)
library(tidyverse)
task_names <- c("cue","foot","hand","tongue")

bayes_num_act[[4]] <- bayes_num_act_45[[1]]
names(bayes_num_act)[4] <- "45"
num_act_df <-
  reshape2::melt(bayes_num_act) %>%
  mutate(task = task_names[Var1],
         task = ifelse(task %in% c("foot","hand"),paste(L3,task),task),
         gamma = as.numeric(sub("gamma = ","", L4)),
         samp_num = as.numeric(L2),
         num_subjects = as.numeric(L1)) %>%
  group_by(task, gamma, num_subjects, samp_num) %>%
  summarize(value = sum(value)) %>%
  ungroup() %>%
  filter(gamma != 1) %>%
  mutate(model = "Bayesian GLM") %>%
  full_join(mutate(class_act_df, model = "Classical GLM")) %>%
  filter(task != "cue")

num_act_plot <- ggplot(mapping = aes(x = num_subjects, y = value, color = task)) +
  geom_jitter(width = 1, height = 1, alpha = 0.3,
              data = filter(num_act_df, num_subjects < 45)) +
  geom_point(data = filter(num_act_df, num_subjects == 45), size = 4) +
  geom_line(data = group_by(num_act_df, num_subjects, gamma, model, task) %>% summarize(value = mean(value)), lwd = 1.6) +
  scale_color_manual("", values = col_pal) +
  facet_grid(gamma~model, scales = "free_y", labeller = label_bquote(rows = gamma == .(gamma)~"%")) +
  labs(x = "Number of Subjects", y = "Number of Active Locations") +
  # theme_classic() +
  theme(legend.position = "right",
        text = element_text(size = 14),
        panel.spacing = unit(0.2, "lines"),
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black",size=1,fill = NA))
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
ggsave(file.path(plot_dir,"607_num_act_plots_thr0_and_thr05.png"), plot = num_act_plot, width = 7, height = 6)


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
    # full_join(classical_num_act) %>%
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

plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
ggsave(file.path(plot_dir,"607_num_act_plots_thr0.png"), plot = num_act_plots[[1]], width = 7, height = 6)
ggsave(file.path(plot_dir,"607_num_act_plots_thr0_5.png"), plot = num_act_plots[[2]], width = 7, height = 6)

# >> Now as one plot ----

num_act_df <-
  reshape2::melt(bayes_num_act) %>%
  full_join(reshape2::melt(bayes_num_act_45)) %>%
  mutate(model = "Bayesian GLM") %>%
  full_join(classical_num_act) %>%
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
  select(-L1,-L4) %>%
  pivot_longer(cols = -c(threshold,num_subjects,model, L2), names_to = "task", values_to = "num_active") %>%
  mutate(both_hems = ifelse(task %in% c("cue","tongue"),"Two Hemispheres","One Hemisphere")) %>%
  filter(#both_hems == hemi,
    task != "cue",
    threshold != 1
  ) %>%
  rename(Sample = L2)
sub_df <- filter(num_act_df, num_subjects < 45)
means_df <- group_by(num_act_df,num_subjects, task, threshold, model) %>%
  summarize(num_active = mean(num_active))
num_act_plot <- ggplot(sub_df, aes(x = num_subjects, y = num_active, color = task)) +
  geom_jitter(width = 1, height = 0, alpha = 0.3) +
  geom_line(data = means_df,lwd = 1.6) +
  geom_point(data = filter(means_df,num_subjects == 45), size = 4) +
  facet_grid(threshold~model, scales = "free", labeller = label_bquote(rows = gamma == .(threshold) ~'%')) +
  # facet_wrap(~threshold + model, scales = "free_y", nrow = 2, ncol = 2, byrow = T) +
  # labs(x = "Number of Subjects", y = "Number of Active Locations", title = expression(paste0("Activation Threshold (",gamma,"=",thr,"%)"))) +
  labs(x = "Number of Subjects", y = "Number of Active Locations") +
  scale_color_manual("", values = col_pal) +
  # geom_hline(yintercept = 0) +
  theme_bw() +
  theme(legend.position = "right",
        text = element_text(size = 14),
        panel.grid = element_blank())
num_act_plot

plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
ggsave(file.path(plot_dir,"607_num_act_plots_thr0_and_thr05.png"), plot = num_act_plot, width = 7, height = 6)

# >> Here's a test to see if the Bayesian 0.5% activations are greater than the classical GLM 0% ----
num_act_df %>%
  filter(threshold != 1,
         num_subjects != 45) %>%
  mutate(model = sub(" GLM", "", model)) %>%
  filter((model == "Bayesian" & threshold == 0.5) | (model == "Classical" & threshold == 0)) %>%
  pivot_wider(names_from = c(model,threshold), values_from = num_active) %>%
  group_by(num_subjects,task) %>%
  summarize(estimate = t.test(x = Bayesian_0.5, y = Classical_0, paired = T)$estimate,
            p_value = t.test(x = Bayesian_0.5, y = Classical_0, paired = T)$p.value)


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
load("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/subjects.Rdata")
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
        print(paste(subject,vis,hem,thr))
        result_file <- grep(
          paste0(subject,"_",vis,"_",hem,"_thr",thr,".rds"),
          bayes_subject_files, value = T
        )
        if(length(result_file) > 1) result_file <- result_file[1]
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
classical_estimates <- readRDS("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed/05_avg_classical_estimates.rds")
# classical_estimates <- readRDS("/Volumes/GoogleDrive/.shortcut-targets-by-id/1UetPPzvfP-rJYdT-9Kp9YBW5a9io2CcA/BayesGLM_Validation_old/5k_results/602_avg_estimates_PW_classical.rds")
threshs <- c(0,0.5,1)
combined_active_FWER <- sapply(classical_estimates, function(hem_res) {
  num_locs <- dim(hem_res[[1]])[1]
  bonferroni_cutoff <- 0.01 / num_locs # alpha = 0.01
  data_df <- reshape2::melt(hem_res)
  vertex_lists <- split(data_df,data_df$Var1)
  thresh_active <- sapply(threshs, function(thr) {
    active_all <- sapply(vertex_lists, function(v_df){
      vt_df <- split(v_df, v_df$Var2)
      active_vertex <- sapply(vt_df, function(vt) {
        ttest_res <- t.test(x = vt$value, mu = thr, alternative = "greater")
        active_v <- ifelse(ttest_res$p.value < bonferroni_cutoff, 1, 0)
        return(active_v)
      }, simplify = "array")
      return(active_vertex)
    }, simplify = 'array')
    return(t(active_all))
  }, simplify = F)
  names(thresh_active) <- paste0(threshs,"%")
  return(thresh_active)
}, simplify = FALSE)
saveRDS(combined_active_FWER, "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/502_HCP_classical_activations_PW_FWER.rds")
# >>>> Group ----
classical_group_num_active <- readRDS("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed/group/502_HCP_classical_activations_PW_visit1_FWER.rds")
# classical_group_num_active <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/502_HCP_classical_activations_PW_FWER.rds")
classical_group_num_active <- sapply(classical_group_num_active, function(hem_act) {
  out <- sapply(hem_act, function(thr_act) {
    as.matrix(apply(thr_act,2,sum)) # Have to do as.matrix for reshape2::melt
  }, simplify = F)
  names(out) <- c(0,0.5,1)
  return(out)
}, simplify = F)

reshape2::melt(classical_group_num_active)

# >>>> Subject ----
load("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/subjects.Rdata")
classical_subject_files <- list.files(
  "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed",
  full.names = T) |>
  grep(pattern = "FWHM6", value = T) |>
  grep(pattern = "classical", value = T) |>
  grep(pattern = "visit1", value = T)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench')
library(BayesfMRI)
subject <- subjects[1]; hem <- 'left'; threshold <- 0
save_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed/activations"
for(subject in subjects) {
  for(hem in c("left","right")) {
    subj_hem_file <- grep(pattern = paste0(subject,"_visit1_",hem),
                         classical_subject_files, value = T)
    if(length(subj_hem_file) > 1) subj_hem_file <- subj_hem_file[1]
    subj_hem_obj <- readRDS(subj_hem_file)
    for(threshold in c(0,0.5,1)) {
      thr_chr <- paste0("thr",sub("\\.","", as.character(threshold)),"_")
      subj_hem_active <-
        id_activations.classical(
          subj_hem_obj$GLMs_classical[[paste0("cortex",toupper(substring(hem,1,1)))]],
          alpha = 0.01,
          threshold = threshold,
          correction = "FWER"
        )
      saveRDS(subj_hem_active,file = file.path(save_dir,paste0("05_visit1_subject_",subject,"_",hem,"_",thr_chr,".rds")))
    }
  }
}

classical_subject_activation_files <- list.files("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed/activations", full.names = T)

classical_subject_num_active <- sapply(subjects, function(subject) {
  sapply(paste0("visit",1), function(vis) {
    sapply(c('left','right'), function(hem) {
      L_or_R <- toupper(substring(hem,1,1))
      all_thr <- sapply(c(0,0.5,1), function(thr) {
        thr_chr <- paste0("thr",sub("\\.","", as.character(thr)),"_")
        print(paste(subject,vis,hem,thr))
        result_file <- grep(
          paste0(vis,"_subject_",subject,"_",hem,"_",thr_chr),
          classical_subject_activation_files, value = T
        )
        if(length(result_file) > 1) result_file <- result_file[1]
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
                         model == "Bayes") #& thresh == 0.5
group_classical_df <- filter(group_activations_df,
                             model == "Classical") #& thresh == 0.5
group_df <- full_join(group_bayes_df, group_classical_df) %>%
  mutate(model = ifelse(model == "Bayes", "Bayesian GLM","Classical GLM"),
         model_thresh = paste0(model," (",thresh,"%)"))
subj_bayes_df <- filter(subject_activations_df,
                        model == "Bayes") #, thresh == 0.5
subj_classical_df <- filter(subject_activations_df,
                            model == "Classical") #, thresh == 0.5
subj_df <- full_join(subj_bayes_df, subj_classical_df)

col_pal <- rgb(c(230, 86, 0, 0,204),c(159,180,158,114,121),c(0,233,115,178,167), maxColorValue = 255)

compare_num_plot <-
  subj_df %>%
  filter(task != "cue") %>%
  mutate(model = ifelse(model == "Bayes", "Bayesian GLM","Classical GLM"),
         model_thresh = paste0(model," (",thresh,"%)")) %>%
  # ggplot() +
  ggplot(aes(x = task, y = activations, color = model, group = model)) +
  # geom_boxplot(aes(x = model, y = activations + 1, color = model)) +
  # geom_boxplot(aes(x = task, y = activations), outlier.alpha = 0, width = 0.4) +
  geom_boxplot(aes(group = interaction(task, model)), color = "black",outlier.alpha = 0, width = 0.4) +
  # geom_jitter(aes(x = task, y = activations, color = task), width = 0.2, height = 0, alpha = 0.3) +
  # geom_jitter(width = 0.2, height = 0, alpha = 0.3) +
  geom_point(alpha = 0.3, position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.4)) +
  # geom_point(aes(x = task, y = activations, color = task), shape = 19, data = group_df, size = 4, show.legend = FALSE) +
  # geom_point(shape = 19, data = group_df, size = 4, show.legend = FALSE) +
  geom_point(shape = 19, data = group_df, size = 4, show.legend = FALSE, position = position_dodge(width = 0.4)) +
  geom_point(color = 'black', shape = 1, data = group_df, size = 4, show.legend = FALSE, position = position_dodge(width = 0.4)) +
  # scale_color_manual("", values = col_pal) +
  # scale_shape_manual("",values = c("B","C")) +
  scale_color_discrete("") +
  facet_wrap(~thresh, scales = "free_y", labeller = label_bquote(cols = gamma == .(thresh)~'%')) +
  # facet_wrap(~task) +
  # facet_grid(~thresh, scales = "free", labeller = label_bquote(cols = gamma == .(thresh)~'%')) +
  # facet_grid(model_thresh~., scales = "fixed") +
  geom_hline(yintercept = 0) +
  labs(x = "",y = "Number of Active Locations") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4))) +
  # scale_y_log10() +
  theme_classic() +
  # theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5), legend.position = "none")
  theme(legend.position = "top",
        text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 0.6, hjust = 0.4),
        axis.ticks.x = element_blank())

compare_num_plot
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
ggsave(file.path(plot_dir,"5_compare_num_plot.png"), width = 8, height = 4.5)

library(gridExtra)
marrangeGrob(grobs = num_activations_plots, nrow = 1, ncol = 3, top ="")

library(ggpubr)
num_activations_plot <- ggarrange(plotlist = num_activations_plots, nrow = 1, ncol = 3)
num_activations_plot
ggsave("plots/5_num_activations_plot.png", width = 10, height = 7)

# <<<<< APPENDIX FIGURES >>>>> ----

# FIGURE D.4: Permutation activations ----
# The classical analysis code
library(ciftiTools)
ciftiTools::ciftiTools.setOption("wb_path", "/Applications/workbench/bin_macosx64/wb_command") # Dan's Macbook Pro
wb_cmd <- "/Applications/workbench/bin_macosx64/wb_command" # Dan's Macbook Pro
library(INLA)
inla.setOption(pardiso.license = "~/licenses/pardiso.lic") # Dan's Macbook Pro
library(BayesfMRI)
main_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/visit1_data"
data_dir <- main_dir
result_dir <- "~/Desktop" # Dan's Macbook Pro
load("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/subjects.Rdata") # Macbook Pro
tasks <- c('cue','lf','lh','rf','rh','t') # Task data frame columns
names_tasks <- c('cue','left_foot','left_hand','right_foot','right_hand','tongue')
colors_tasks <- c('black',RColorBrewer::brewer.pal(5, 'Set2'))
cols_LH <- c(1,4:6) #cue, right foot, right hand, tongue
cols_RH <- c(1:3,6) #cue, left foot, left hand, tongue
cols_list <- list(cols_LH, cols_RH)
TR = 0.72 #temporal resolution of data
thetas <- NULL # No starting values for precision parameters
subject <- subjects[2]; visit <- 1; h <- 1
dir_s <- file.path(data_dir, subject, 'MNINonLinear', 'fsaverage_LR32k')
fname_gifti_left <- file.path(dir_s, paste0(subject,'.L.midthickness.32k_fs_LR.surf.gii'))
fname_gifti_right <- file.path(dir_s, paste0(subject,'.R.midthickness.32k_fs_LR.surf.gii'))
dir1_s <- file.path(data_dir,subject, 'MNINonLinear', 'Results', 'tfMRI_MOTOR_LR')
dir2_s <- file.path(data_dir,subject, 'MNINonLinear', 'Results', 'tfMRI_MOTOR_RL')
fname1_ts <- file.path(dir1_s,'tfMRI_MOTOR_LR_Atlas_FWHM6.dtseries.nii')
fname2_ts <- file.path(dir2_s,'tfMRI_MOTOR_RL_Atlas_FWHM6.dtseries.nii')
hem <- c('left','right')[h]
cols_h <- cols_list[[h]]
tasks_h <- tasks[cols_h]
names_tasks_h <- names_tasks[cols_h]
K <- length(tasks_h)
#Set up list of onsets
onsets1 <- onsets2 <- vector('list', length=K)
names(onsets1) <- names(onsets2) <- names_tasks_h
for(t in tasks_h){
  ind_t <- which(tasks_h==t)
  fname1_t <- file.path(dir1_s,'EVs',paste0(t,'.txt'))
  onsets1[[ind_t]] <- read.table(fname1_t, header=FALSE)
  fname2_t <- file.path(dir2_s,'EVs',paste0(t,'.txt'))
  onsets2[[ind_t]] <- read.table(fname2_t, header=FALSE)
}
#Set up nuisance regressors
motion1 <- as.matrix(read.table(file.path(dir1_s,'Movement_Regressors.txt'), header=FALSE))
motion2 <- as.matrix(read.table(file.path(dir2_s,'Movement_Regressors.txt'), header=FALSE))
start_time <- proc.time()[3]
result_svh <- BayesGLM_cifti(cifti_fname = c(fname1_ts, fname2_ts), # Multi-session
                             # cifti_fname = fname1_ts, # single-session
                             surfL_fname = fname_gifti_left,
                             surfR_fname = fname_gifti_right,
                             brainstructures = hem,
                             design = NULL,
                             onsets = list(onsets1, onsets2), # Multi-session
                             # onsets = onsets1, # single-session
                             TR = TR,
                             nuisance = list(motion1, motion2), # Multi-session
                             # nuisance = motion1, # single-session
                             nuisance_include = c('drift','dHRF'),
                             scale_BOLD = TRUE,
                             scale_design = TRUE,
                             GLM_method = 'classical',
                             ar_order = 6,
                             ar_smooth = 6,
                             session_names = c('LR','RL'), # Multiple sessions
                             # session_names = c('LR'), # single session
                             resamp_res = 5000, # Don't forget to change this
                             num.threads = 6, # Remember the tradeoff here (speed/memory) 4 to 6 threads seems optimal based on testing
                             verbose = TRUE,
                             outfile = NULL,
                             return_INLA_result = T,
                             avg_sessions = T,
                             trim_INLA = T,
                             num_permute = 1000)
total_time <- proc.time()[3] - start_time
result_svh$total_time <- total_time
saveRDS(result_svh, file=file.path(result_dir,paste0("500_",subject,"_visit",visit,"_",hem,"_5k_classical_FWHM6_",format(Sys.Date(),"%Y%m%d"),".rds")))

library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
col_pal <- c(light_orange,"red","purple")
result <- readRDS("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed/permutations/500_114823_visit1_left_5k_classical_FWHM6_20211101.rds")
library(BayesfMRI)
active_perm <- sapply(c(0,0.5,1), function(thr) {
  id_out <- id_activations_cifti(result,alpha = 0.01,method = 'classical',threshold = thr,correction = 'permutation')
  return(id_out$activations_xifti)
},simplify = F)

active_perm <- Reduce(`+`,active_perm)
active_perm$data$cortex_left[active_perm$data$cortex_left == 0] <- NA
plot(active_perm, idx = 4, colors = col_pal, color_mode = "qualitative",
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii", hemisphere = 'left', view = 'lateral',
     fname = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/05_classical_114823_activations_permutation.png"
     )

active_FWER <- sapply(c(0,0.5,1), function(thr) {
  id_out <- id_activations_cifti(result,alpha = 0.01,method = 'classical',threshold = thr,correction = 'FWER')
  return(id_out$activations_xifti)
},simplify = F)
active_FWER <- Reduce(`+`,active_FWER)
active_FWER$data$cortex_left[active_FWER$data$cortex_left == 0] <- NA
plot(active_FWER, idx = 4, colors = col_pal, color_mode = "qualitative",
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii", hemisphere = 'left', view = 'lateral', fname = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/05_classical_activations_FWER.png")

active_FDR <- sapply(c(0,0.5,1), function(thr) {
  id_out <- id_activations_cifti(result,alpha = 0.01,method = 'classical',threshold = thr,correction = 'FDR')
  return(id_out$activations_xifti)
},simplify = F)
active_FDR <- Reduce(`+`,active_FDR)
active_FDR$data$cortex_left[active_FDR$data$cortex_left == 0] <- NA
plot(active_FDR, idx = 4, colors = col_pal, color_mode = "qualitative",
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii", hemisphere = 'left', view = 'lateral', fname = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/05_classical_114823_activations_FDR.png"
     )

active_bayes <- list.files("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/individual/PW/activations", full.names = T) |>
  grep(pattern = "visit1_subject_103818_left", value = T) |>
  grep(pattern = "_classical", value = T, invert = T) |>
  sapply(FUN = function(x) {y <- readRDS(x); return(y$activations_xifti)},simplify = F) |>
  Reduce(f = `+`)
active_bayes$data$cortex_left[active_bayes$data$cortex_left == 0] <- NA
plot(active_bayes, idx = 4, colors = col_pal, color_mode = "qualitative",
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii", hemisphere = 'left', view = 'lateral')# fname = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/05_bayes_activations.png")

# FIGURE F.7: Subject test-retest metrics comparison to full-res classical analyses -----

bayes_result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
# classical_result_dir_unsmoothed <- "HCP_results/32k_results"
classical_result_dir_unsmoothed <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/32k_results"
# classical_result_dir_smoothed <- "HCP_results/32k_results/smoothed"
classical_result_dir_smoothed <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/32k_results/smoothed"
classical_result_dir_5k <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/classical"
classical_result_dir_5k_unsmoothed <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
# sphere_5k_result_dir <- "HCP_results/5k_results/individual/PW/spherical_analyses"
sphere_5k_result_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/individual/PW/spherical_analyses"

# load("HCP_data/subjects.Rdata")
load("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/subjects.Rdata")

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
# saveRDS(classical_visit2_5k_sphere,"HCP_results/5k_results/5_classical_visit2_sphere_est.rds")
# classical_visit2_5k_sphere <- readRDS("HCP_results/5k_results/5_classical_visit2_sphere_est.rds")

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
# saveRDS(bayes_visit1_est, "HCP_results/5k_results/5_bayes_visit1_est.rds")
# bayes_visit1 <- readRDS("HCP_results/5k_results/5_bayes_visit1_est.rds")
bayes_visit1 <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/5_bayes_visit1_est.rds")
bayes_visit1 <- bayes_visit1[keep_subjects]
# HERE ----
# bayes_visit1_sphere <- list()
# for(subject in subjects) {
#   bayes_file_left <- grep(paste0(subject,"_visit1_left"), list.files(sphere_5k_result_dir, full.names = T), value = T)
#   bayes_file_right <- grep(paste0(subject,"_visit1_right"), list.files(sphere_5k_result_dir, full.names = T), value = T)
#   ests <- list(
#     left = readRDS(bayes_file_left)$betas_Bayesian$avg$data$cortex_left,
#     right = readRDS(bayes_file_right)$betas_Bayesian$avg$data$cortex_right
#   )
#   bayes_visit1_sphere[[subject]] <- ests
# }
# saveRDS(bayes_visit1_sphere,"HCP_results/5k_results/5_bayes_visit1_sphere_est.rds")
# bayes_visit1_sphere <- readRDS("HCP_results/5k_results/5_classical_visit2_sphere_est.rds")
bayes_visit1_sphere <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/individual/PW/spherical_analyses")

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

ggsave("plots/5_mse_vs_full_res_classical.png",plot = subj_mse_plot, width = 7, height = 5)

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

ggsave("plots/5_cor_vs_full_res_classical.png",plot = subj_cor_plot, width = 7, height = 4)


# FIGURE H.10: Single run amplitude estimates ----
# Bayesian ----
# Classical ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
library(BayesfMRI)
main_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results"
result_dir <- file.path(main_dir, "smoothed")
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
load(file.path(main_dir,"../subjects.RData"))
subjects <- subjects[c(1,2,4)]
subject_files <- list.files(result_dir,full.names = T) |>
  grep(pattern = "FWHM6", value = T) |>
  grep(pattern ="classical", value = T) |>
  grep(pattern = "visit1", value = T) |>
  sapply(X = subjects,
         FUN = function(m,sf) grep(pattern = m, x = sf, value = T),
         simplify = T) |>
  c()
for(subject in subjects) {
  filens <- grep(pattern = subject, subject_files, value = T)
  result_obj_left <- readRDS(grep(pattern = "left", filens, value = T))
  result_obj_right <- readRDS(grep(pattern = "right", filens, value = T))
  plot_obj <- combine_xifti(result_obj_left$betas_classical$LR,
                            result_obj_right$betas_classical$LR)
  plot(plot_obj, zlim = c(-1,1), idx = 4, legend_embed = F,
       fname = file.path(plot_dir,paste0("607_classical_",subject,"_visit1_sessionLR_tongue_estimate.png")))
}

# FIGURE H.11: Single-subject estimates ----
# >> Bayesian ----
# >> Classical ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
library(BayesfMRI)
main_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results"
result_dir <- file.path(main_dir, "smoothed")
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
subject_files <- list.files(result_dir,full.names = T) |>
  grep(pattern = "FWHM6", value = T) |>
  grep(pattern ="classical", value = T) |>
  grep(pattern = "visit1", value = T) |>
  grep(pattern = "103818", value = T)
task_names <- c("visual_cue","foot","hand","tongue")
subject <- "103818"

result_obj_left <- readRDS(grep(pattern = "left", subject_files, value = T))
result_obj_right <- readRDS(grep(pattern = "right", subject_files, value = T))
plot_obj <- combine_xifti(result_obj_left$betas_classical$LR,
                          result_obj_right$betas_classical$LR)
for(task_idx in c(1,4)) {
  plot(plot_obj, zlim = c(-1,1), idx = task_idx, legend_embed = F,
       fname = file.path(plot_dir,
                         paste0("607_classical_",subject,"_visit1_sessionLR_",
                                task_names[task_idx],"_estimate.png")))
}

for(task_idx in c(2,3)) {
  for(hem in c("left","right")) {
    other_hem <- grep(hem, c("left","right"), value = T, invert = T)
    plot(plot_obj, zlim = c(-1,1), idx = task_idx, legend_embed = F,
         hemisphere = hem,
         fname = file.path(plot_dir,
                           paste0("607_classical_",subject,"_visit1_sessionLR_",
                                  other_hem,"_",
                                  task_names[task_idx],"_estimate.png")))
  }
}

# FIGURE H.12: Average Estimates ----
# >> Bayesian ----
# >> Classical ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
library(BayesfMRI)
main_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results"
result_dir <- file.path(main_dir, "smoothed")
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
subject_files <- list.files(result_dir,full.names = T) |>
  grep(pattern = "FWHM6", value = T) |>
  grep(pattern ="classical", value = T) |>
  grep(pattern = "visit1", value = T) |>
  grep(pattern = "103818", value = T)
task_names <- c("cue","foot","hand","tongue")
subject <- "103818"

result_obj_left <- readRDS(grep(pattern = "left", subject_files, value = T))
result_obj_right <- readRDS(grep(pattern = "right", subject_files, value = T))
plot_obj <- combine_xifti(result_obj_left$betas_classical$avg,
                          result_obj_right$betas_classical$avg)
for(task_idx in c(1,4)) {
  plot(plot_obj, zlim = c(-1,1), idx = task_idx, legend_embed = F,
       fname = file.path(plot_dir,
                         paste0("600_subject_",subject,"_",
                                task_names[task_idx],"_classical_estimates.png")))
}

for(task_idx in c(2,3)) {
  for(hem in c("left","right")) {
    other_hem <- grep(hem, c("left","right"), value = T, invert = T)
    plot(plot_obj, zlim = c(-1,1), idx = task_idx, legend_embed = F,
         hemisphere = hem,
         fname = file.path(plot_dir,
                           paste0("600_subject_",subject,"_",
                                  other_hem,"_",
                                  task_names[task_idx],"_classical_estimates.png")))
  }
}

# FIGURE H.14: ICC values, all tasks ----
library(ciftiTools)
ciftiTools.setOption('wb_path',"/Applications/workbench")
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
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
result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
avg_estimates_classical <- readRDS(file.path(result_dir,"05_avg_classical_estimates.rds"))
library(abind)
ICC_weights <- sapply(avg_estimates_classical, function(hem_est) {
  combined_est <- Reduce(abind,hem_est)
  avg_est <- apply(combined_est,1:2,mean)
  weight_est <- apply(avg_est,2,function(avg_task) abs(avg_task) / sum(abs(avg_task)))
  return(weight_est)
}, simplify = F)
task_file_names <- c("visual_cue","foot","hand","tongue")

# >> Bayesian ----
result_dir <- "HCP_results/5k_results"
avg_estimates <- readRDS(file.path(result_dir,"5_avg_estimates.rds"))
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

bayesian_icc_avg_cifti <- readRDS("HCP_data/cifti_5k_template_whole.rds")
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
    fname = paste0("plots/5_",task_file_names[task_idx],"_icc_bayesian.png")
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
      fname = paste0("plots/5_",hem[h],"_",task_file_names[ti],"_icc_bayesian.png")
    )
  }
}

# >> Classical ----
avg_estimates_classical <- readRDS(file.path(result_dir,"05_avg_classical_estimates.rds"))
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

classical_icc_avg_cifti <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/603_cifti_5k_template_whole.rds")
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
    fname = file.path(plot_dir,paste0("602_",task_file_names[task_idx],"_icc_classical.png"))
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
      fname = file.path(plot_dir,paste0("602_",hem[h],"_",task_file_names[ti],"_icc_classical.png"))
    )
  }
}

# FIGURE E.10: Single-run activations ----
result_dir_classical <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
library(BayesfMRI)
library(ciftiTools)
ciftiTools.setOption('wb_path',"/Applications/workbench")
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
col_pal <- c(light_orange,"red","purple")
task_names <- c("cue","foot","hand","tongue")
subject <- "103818"
# >> Bayesian ----
# >> Classical ----
subject_files <- grep(pattern = "103818",
                      list.files(result_dir_classical, full.names = TRUE),
                      value = TRUE) |>
  grep(pattern = "FWHM6", value = T) |>
  grep(pattern = "visit1", value = T) |>
  grep(pattern = "sessionLR", value = T) |>
  grep(pattern = "classical", value = T)
left_obj <- readRDS(grep("left",subject_files, value = T))
left_cifti <- left_obj$betas_classical$LR
left_active <- sapply(c(0,0.5,1), function(thr) {
  id_activations_cifti(left_obj,alpha = 0.01,method = 'classical',threshold = thr,session_name = "LR")
}, simplify = F)
left_active <- left_active[[1]]$activations$cortexL$active +
  left_active[[2]]$activations$cortexL$active +
  left_active[[3]]$activations$cortexL$active
left_active[left_active == 0] <- NA
left_cifti$data$cortex_left <- left_active
right_obj <- readRDS(grep("right",subject_files, value = T))
right_cifti <- right_obj$betas_classical$LR
right_active <- sapply(c(0,0.5,1), function(thr) {
  id_activations_cifti(right_obj,alpha = 0.01,method = 'classical',threshold = thr, session_name = "LR")
}, simplify = F)
right_active <- right_active[[1]]$activations$cortexR$active +
  right_active[[2]]$activations$cortexR$active +
  right_active[[3]]$activations$cortexR$active
right_active[right_active == 0] <- NA
right_cifti$data$cortex_right <- right_active
plot_obj <- combine_xifti(left_cifti, right_cifti)

for(task_idx in c(1,4)) {
  plot(plot_obj, idx = task_idx, legend_embed = F, color_mode = "qualitative",
       colors = col_pal,
       fname = file.path(plot_dir,
                         paste0("607_classical_",subject,"_visit1_sessionLR_",
                                task_names[task_idx],"_activations.png")))
}

for(task_idx in c(2,3)) {
  for(hem in c("left","right")) {
    other_hem <- grep(hem, c("left","right"), value = T, invert = T)
    plot(plot_obj, idx = task_idx, legend_embed = F, colors = col_pal,
         hemisphere = hem, color_mode = "qualitative",
         fname = file.path(plot_dir,
                           paste0("607_classical_",subject,"_visit1_sessionLR_",
                                  other_hem,"_",
                                  task_names[task_idx],"_activations.png")))
  }
}

# FIGURE E.11: Two-run activations ----
result_dir_classical <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
library(BayesfMRI)
library(ciftiTools)
ciftiTools.setOption('wb_path',"/Applications/workbench")
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
col_pal <- c(light_orange,"red","purple")
task_names <- c("cue","foot","hand","tongue")
subject <- "103818"
# >> Bayesian ----
# >> Classical ----
subject_files <- grep(pattern = "103818",
                      list.files(result_dir_classical, full.names = TRUE),
                      value = TRUE) |>
  grep(pattern = "FWHM6", value = T) |>
  grep(pattern = "visit1", value = T)
left_obj <- readRDS(grep("left",subject_files, value = T))
left_cifti <- left_obj$betas_classical$avg
left_active <- sapply(c(0,0.5,1), function(thr) {
  id_activations_cifti(left_obj,alpha = 0.01,method = 'classical',threshold = thr)
}, simplify = F)
left_active <- left_active[[1]]$activations$cortexL$active +
  left_active[[2]]$activations$cortexL$active +
  left_active[[3]]$activations$cortexL$active
left_active[left_active == 0] <- NA
left_cifti$data$cortex_left <- left_active
right_obj <- readRDS(grep("right",subject_files, value = T))
right_cifti <- right_obj$betas_classical$avg
right_active <- sapply(c(0,0.5,1), function(thr) {
  id_activations_cifti(right_obj,alpha = 0.01,method = 'classical',threshold = thr)
}, simplify = F)
right_active <- right_active[[1]]$activations$cortexR$active +
  right_active[[2]]$activations$cortexR$active +
  right_active[[3]]$activations$cortexR$active
right_active[right_active == 0] <- NA
right_cifti$data$cortex_right <- right_active
plot_obj <- combine_xifti(left_cifti, right_cifti)

for(task_idx in c(1,4)) {
  plot(plot_obj, idx = task_idx, legend_embed = F, color_mode = "qualitative",
       colors = col_pal,
       fname = file.path(plot_dir,
                         paste0("607_classical_",subject,"_visit1_",
                                task_names[task_idx],"_activations.png")))
}

for(task_idx in c(2,3)) {
  for(hem in c("left","right")) {
    other_hem <- grep(hem, c("left","right"), value = T, invert = T)
    plot(plot_obj, idx = task_idx, legend_embed = F, colors = col_pal,
         hemisphere = hem, color_mode = "qualitative",
         fname = file.path(plot_dir,
                           paste0("607_classical_",subject,"_visit1_",
                                  other_hem,"_",
                                  task_names[task_idx],"_activations.png")))
  }
}


# FIGURE: Correlations boxplot for smoothed classical results ----
classical_result_dir_5k <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/classical"

classical_smoothed_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
smoothed_files <- list.files(classical_smoothed_dir, full.names = T) |>
  grep(pattern = "_5k_",value = T)
# How far are the analyses right now?
smoothed_files_8mm <- grep(pattern = "_FWHM8_",smoothed_files, value = TRUE) |>
  grep(pattern = "visit1", value = TRUE) |>
  grep(pattern = "right", value = TRUE)
tail(smoothed_files_8mm,1) # Done up to subject 917255 (That's all of them!)
load("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/subjects.Rdata")
subjects_smooth <- subjects[seq(1,which(subjects == "917255"))]

classical_visit1_5k <- sapply(subjects_smooth, function(subject) {
  classical_file_left <- grep(paste0(subject,"_visit1_left"), list.files(classical_result_dir_5k, full.names = T), value = T)
  classical_file_right <- grep(paste0(subject,"_visit1_right"), list.files(classical_result_dir_5k, full.names = T), value = T)
  class_visit1_est <- list(
    left = readRDS(classical_file_left)$betas_classical$avg$data$cortex_left,
    right = readRDS(classical_file_right)$betas_classical$avg$data$cortex_right
  )
  return(class_visit1_est)
}, simplify = FALSE)

classical_visit1_5k_smoothed <- sapply(seq(3,10), function(fwhm) {
  sapply(subjects_smooth, function(subject) {
    print(paste(subject,fwhm))
    classical_file_left <- grep(paste0(subject,"_visit1_left"), smoothed_files, value = T) |> grep(pattern = paste0("FWHM",fwhm), value = T)
    classical_file_right <- grep(paste0(subject,"_visit1_right"), smoothed_files, value = T) |> grep(pattern = paste0("FWHM",fwhm), value = T)
    if(length(classical_file_left) > 1) classical_file_left <- classical_file_left[1]
    if(length(classical_file_right) > 1) classical_file_right <- classical_file_right[1]
    class_visit1_est <- list(
      left = readRDS(classical_file_left)$betas_classical$avg$data$cortex_left,
      right = readRDS(classical_file_right)$betas_classical$avg$data$cortex_right
    )
    return(class_visit1_est)
  }, simplify = FALSE)
}, simplify = F)
names(classical_visit1_5k_smoothed) <- paste0("FWHM",seq(3,10))

classical_visit1_5k_smoothed$FWHM0 <- classical_visit1_5k

classical_visit2_5k <- sapply(subjects, function(subject) {
  classical_file_left <- grep(paste0(subject,"_visit2_left"), list.files(classical_result_dir_5k, full.names = T), value = T)
  classical_file_right <- grep(paste0(subject,"_visit2_right"), list.files(classical_result_dir_5k, full.names = T), value = T)
  class_visit2_est <- list(
    left = readRDS(classical_file_left)$betas_classical$avg$data$cortex_left,
    right = readRDS(classical_file_right)$betas_classical$avg$data$cortex_right
  )
  return(class_visit2_est)
}, simplify = FALSE)

smooth_dims <- sapply(classical_visit1_5k_smoothed, function(x) {
  sapply(x, function(xx) {
    sapply(xx, function(xxx) as.matrix(dim(xxx)), simplify = F)
  }, simplify = F)
}, simplify = F)
library(tidyverse)
reshape2::melt(smooth_dims, value.name = "dimension") %>%
  group_by(L3, Var1) %>%
  summarize(num_diff = length(unique(dimension)))

classical_cor_5k <- sapply(classical_visit1_5k_smoothed, function(cv1_5ks) {
  mapply(function(est,tru) {
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
  }, est = cv1_5ks, tru = classical_visit2_5k, SIMPLIFY = F)
}, simplify = FALSE)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

color_pal <- c(gg_color_hue(2)[1], "red", gg_color_hue(2)[2], "turquoise2","blue3")

smooth_cor_plot <-
  reshape2::melt(classical_cor_5k, value.name = "Correlation") %>%
  mutate(variable = sub("_"," ", variable),
         FWHM = sub("FWHM","",L1),
         FWHM = factor(FWHM, levels = as.character(0:10))) %>%
  ggplot() +
  geom_boxplot(aes(x = FWHM, y = Correlation, color = FWHM)) +
  labs(
    y = "Correlation",
    x = ""
  ) +
  scale_color_discrete("Smoothing FWHM") +
  facet_wrap(~variable, scales = "fixed",nrow = 1) +
  guides(color = guide_legend(nrow = 1)) +
  theme_classic() +
  theme(text=element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.box.spacing = unit(-15,"points"),
        # legend.title = element_blank(),
        plot.margin = margin(0,1,1,1,"pt"))
# theme(plot.title = element_text(hjust = 0.5,size = 7),
#       axis.title = element_text(size = 9))

smooth_cor_plot

ggsave(filename = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/05_smooth_cor_plot.png",plot = smooth_cor_plot, width = 7, height = 4)

# <<<<< REVIEWER RESPONSE FIGURES >>>>> ----

# FIGURE: Comparison of Bayesian Estimates for smoothed and unsmoothed data ----
smoothed_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
smoothed_files <- list.files(smoothed_dir, full.names = T) |>
  grep(pattern = "Bayes", value = TRUE) |>
  grep(pattern = "visit1", value = TRUE) |>
  grep(pattern = "103818", value = TRUE)
unsmoothed_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/PW"
unsmoothed_files <- list.files(unsmoothed_dir, full.names = T) |>
  grep(pattern = "visit1", value = TRUE) |>
  grep(pattern = "103818", value = TRUE)

unsmoothed_cifti <- list(
  left = list(),
  right = list()
)
smoothed_cifti <- list(
  left = list(),
  right = list()
)
for(hem in c('left','right')) {
  unsmoothed_result_obj <- readRDS(grep(hem,unsmoothed_files, value = T))
  unsmoothed_cifti[[hem]] <- unsmoothed_result_obj$betas_Bayesian$avg
  smoothed_result_obj <- readRDS(grep(hem,smoothed_files, value = T))
  smoothed_cifti[[hem]] <- smoothed_result_obj$betas_Bayesian$avg
}

unsmoothed_cifti <- combine_xifti(unsmoothed_cifti$left,unsmoothed_cifti$right)
smoothed_cifti <- combine_xifti(smoothed_cifti$left,smoothed_cifti$right)

plot(unsmoothed_cifti, idx = 4, zlim = c(-1,1), title = "Unsmoothed")
plot(smoothed_cifti, idx = 4, zlim = c(-1,1), title = "Smoothed")

# FIGURE: Plots of hyperparameter densities for lateral tasks ----
result_obj_left <- readRDS("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/all_tasks/500_103818_visit1_left_5k_6tasks_20210917.rds")
result_obj_right <- readRDS("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/all_tasks/500_103818_visit1_right_5k_6tasks_20210917.rds")
library(tidyverse)
hyperpar_df <-
  mutate(result_obj_left$GLMs_Bayesian$cortexL$theta_posteriors,
         hemisphere = 'left') %>%
  full_join(
    mutate(result_obj_right$GLMs_Bayesian$cortexR$theta_posteriors,
           hemisphere = 'right')
  ) %>%
  filter(grepl(pattern = "hand",beta) | grepl("foot",beta))

my_labels <- c(
  `kappa left foot` = expression(kappa["left foot"]),
  `kappa right foot` = expression(kappa["right foot"]),
  `kappa left hand` = expression(kappa["left hand"]),
  `kappa right hand` = expression(kappa["right hand"]),
  `tau left foot` = expression(tau["left foot"]),
  `tau right foot` = expression(tau["right foot"]),
  `tau left hand` = expression(tau["left hand"]),
  `tau right hand` = expression(tau["right hand"])
)

kappa_tau_plot <-
  hyperpar_df %>%
  mutate(beta = sub("_"," ", beta),
         param = sub("log_","",param),
         param_beta = paste(param,beta),
         param_beta = factor(param_beta,
                             levels = c("kappa left foot", "kappa right foot",
                                        "kappa left hand", "kappa right hand",
                                        "tau left foot", "tau right foot",
                                        "tau left hand", "tau right hand"),
                             labels = c('kappa[plain(left~foot)]',
                                        'kappa[plain(right~foot)]',
                                        'kappa[plain(left~hand)]',
                                        'kappa[plain(right~hand)]',
                                        'tau[plain(left~foot)]',
                                        'tau[plain(right~foot)]',
                                        'tau[plain(left~hand)]',
                                        'tau[plain(right~hand)]')),
         hemisphere = paste(hemisphere, "hemisphere")) %>%
  ggplot(aes(x = exp(value), y = density, fill = hemisphere, group = hemisphere)) +
  geom_line() +
  geom_ribbon(aes(ymax = density, ymin = 0), alpha = 0.3) +
  labs(x = "value", y = "posterior density") +
  scale_fill_discrete("") +
  facet_wrap(~ param_beta, scales = "free", nrow = 2, labeller = label_parsed) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text = element_text(size = 14))

kappa_tau_plot

ggsave(filename = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/05_lateral_task_posterior_density_kappa_tau.png",plot = kappa_tau_plot, width = 9, height = 5)

# FIGURE: Errors across hemisphere ----
library(Matrix)
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
classical_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
result_files <- list.files(classical_dir, full.names = TRUE) |>
  grep(pattern = "FWHM6", value = TRUE) |>
  grep(pattern = "Bayes", value = TRUE) |>
  grep(pattern = "visit1", value = TRUE) |>
  grep(pattern = "103818", value = TRUE)

result_files_classical <- list.files(classical_dir, full.names = TRUE) |>
  grep(pattern = "FWHM6", value = TRUE) |>
  grep(pattern = "classical", value = TRUE) |>
  grep(pattern = "visit1", value = TRUE) |>
  grep(pattern = "103818", value = TRUE)

# hem <- 'left'
err_list <- list(
  left = list(
    session1 = matrix(NA,4443, 284),
    session2 = matrix(NA,4443, 284)
  ),
  right = list(
    session1 = matrix(NA,4444, 284),
    session2 = matrix(NA,4444, 284)
  )
)
for(hem in c('left','right')) {
  L_or_R <- toupper(substring(hem,1,1))
  hem_file <- grep(hem,result_files, value = T)
  hem_obj <- readRDS(hem_file)
  hem_file_classical <- grep(hem, result_files_classical, value = T)
  hem_obj_classical <- readRDS(hem_file_classical)
  ntime <- nrow(hem_obj$design[[1]])
  nvox <- sum(hem_obj$GLMs_Bayesian[[paste0("cortex",L_or_R)]]$mask)
  nsess <- length(hem_obj$session_names)
  pw_y <- hem_obj$GLMs_Bayesian[[paste0("cortex",L_or_R)]]$y
  pw_X <- hem_obj$GLMs_Bayesian[[paste0("cortex",L_or_R)]]$X
  # Bayesian
  # betas <- hem_obj$GLMs_Bayesian[[paste0("cortex",L_or_R)]]$beta_estimates
  # betas <- sapply(betas,function(b) c(b[hem_obj$GLMs_Bayesian[[paste0("cortex",L_or_R)]]$mask,]), simplify = F)
  # Classical
  betas <- sapply(hem_obj_classical$GLMs_classical[[paste0("cortex",L_or_R)]][1:2],
                  function(x) c(x$estimates[x$mask,]), simplify = F)
  fit_pwy <- mapply(`%*%`, x = pw_X, y = betas)
  fit_pwy <- Reduce(function(x,y) c(x@x,y@x),fit_pwy)
  err <- pw_y - fit_pwy
  err_array <- array(err, dim = c(nvox,ntime,nsess))
  err_sess <- asplit(err_array,3)
  err_list[[hem]]$session1 <- err_sess[[1]]
  err_list[[hem]]$session2 <- err_sess[[2]]
}

err_list <- sapply(err_list, Reduce, f = cbind, simplify = F)
err_mat <- Reduce(rbind,err_list)
err_corr <- cor(t(err_mat))
corr_upper_tri <- err_corr[upper.tri(err_corr)]

# random_vert <- sample(1:4443, size = 3)
# These were found randomly earlier, but I didn't set a seed to make the plots.
# Fortunately, the file names include the vertex numbers, so that is why they
# are being specified in this way here.
random_vert <- c(1973,2775, 3423)
cifti_corr <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_data/603_cifti_5k_template_whole.rds")
cifti_corr$data$cortex_left <- t(err_corr[random_vert, 1:4443])
cifti_corr$data$cortex_right <- t(err_corr[random_vert, 4444:8887])
vert_idx <- 1
vert_mat <- matrix(0,8887,3)
for(i in 1:3) {
  vert_mat[random_vert[i],i] <- 1
}
cifti_vert <- newdata_xifti(cifti_corr, vert_mat)
# cifti_vert <- convert_xifti(cifti_vert, to = "dlabel")
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
for(vert_idx in 1:3) {
  plot(cifti_corr, zlim = c(-0.5,0.5), idx = vert_idx, legend_embed = F,
       title = paste("Correlations: Vertex",random_vert[vert_idx]),
       fname = file.path(plot_dir,paste0("05_correlations_vertex",
                                         random_vert[vert_idx],
                                         "_classical",
                                         ".png")))

  # plot(cifti_vert, color_mode = "qualitative", idx = vert_idx,
  #      title = paste("Location, Vertex",random_vert[vert_idx]),
  #      # borders = TRUE, edge_color = 'black',
  #      colors = c("yellow","blue"),
  #      fname = file.path(plot_dir,paste0("05_vertex",random_vert[vert_idx],".png")))
}

# FIGURE: Classical subject estimates FWHM 4,8,10----
library(ciftiTools)
ciftiTools.setOption("wb_path","/Applications/workbench")
result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
result_files <- list.files(result_dir, full.names = T)
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
load("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
for(subject in subjects) {
  for(fwhm in c(4,8,10)) {
    result_file <- grep(subject, result_files, value = T) |>
      grep(pattern = "visit1", value = T) |>
      grep(pattern = "left", value = T) |>
      grep(pattern = paste0("FWHM",fwhm), value = T)
    if(length(result_file) > 1) result_file <- result_file[1]
    result <- readRDS(result_file)
    plot(result$betas_classical$avg, hemisphere = "left", view = "lateral",
         idx = 4, zlim = c(-1,1), legend_embed = F,
         fname = file.path(plot_dir,
                           paste0("600_subject_",subject,
                                  "_tongue_classical_estimates_FWHM",
                                  fwhm,".png")))
  }
}

# FIGURE: Classical subject activations FWHM 4,8,10----
library(ciftiTools)
ciftiTools.setOption("wb_path","/Applications/workbench")
library(INLA)
inla.setOption(pardiso.license = "~/licenses/pardiso.lic")
library(BayesfMRI)
result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
result_files <- list.files(result_dir, full.names = T)
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
col_pal <- c(light_orange,"red","purple")
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
load("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
for(subject in subjects) {
  for(fwhm in c(4,8,10)) {
    result_file <- grep(subject, result_files, value = T) |>
      grep(pattern = "visit1", value = T) |>
      grep(pattern = "left", value = T) |>
      grep(pattern = paste0("FWHM",fwhm), value = T)
    if(length(result_file) > 1) result_file <- result_file[1]
    result <- readRDS(result_file)
    subj_act <- sapply(c(0,0.5,1), function(g) {
      out <- id_activations_cifti(result,alpha = 0.01,method = "classical",threshold = g,correction = "FWER")
      return(out$activations_xifti)
    }, simplify = F)
    subj_act <- Reduce(`+`,subj_act)
    subj_act$data$cortex_left[subj_act$data$cortex_left == 0] <- NA
    plot(subj_act, hemisphere = "left", view = "lateral",
         idx = 4, color_mode = "qualitative", colors = col_pal, legend_embed = F,
         fname = file.path(plot_dir,
                           paste0("600_subject_",subject,
                                  "_tongue_classical_activations_FWHM",
                                  fwhm,".png"))
         )
  }
}

# FIGURE: Classical group estimates and activations FWHM ----
library(ciftiTools)
ciftiTools.setOption("wb_path","/Applications/workbench")
library(INLA)
inla.setOption(pardiso.license = "~/licenses/pardiso.lic")
library(BayesfMRI)
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
col_pal <- c(light_orange,"red","purple")
result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed"
for(fwhm in paste0("FWHM",c(4,6,8,10))) {
  result_files <- list.files(result_dir, full.names = T) |>
    grep(pattern = "classical", value = T) |>
    grep(pattern = fwhm, value = T) |>
    grep(pattern = "visit1", value = T) |>
    grep(pattern = "left", value = T)
  classical_group <-
    classicalGLM2(
      results = result_files,
      brainstructure = "cortexL",
      session_name = "avg",
      gamma = c(0, 0.5, 1),
      correction = "FWER",
      alpha = 0.01
    )
  cifti_obj <- readRDS(result_files[1])$betas_classical$avg
  cifti_est <- cifti_obj
  cifti_est$data$cortex_left <- classical_group$avg_estimate
  cifti_act <- cifti_obj
  all_act <- sapply(classical_group$active_result, function(x) x$active, simplify = F)
  all_act <- Reduce(`+`,all_act)
  all_act[all_act == 0] <- NA
  cifti_act$data$cortex_left <- all_act
  plot(cifti_est, view = "lateral", hemisphere = "left", zlim = c(-1,1), idx = 4,
       legend_embed = F,
       fname = file.path(plot_dir,paste0("607_group_classical_tongue_estimates_",fwhm,".png")))
  plot(cifti_act, view = "lateral", hemisphere = "left", idx = 4,
       legend_embed = F, color_mode = "qualitative", colors = col_pal,
       fname = file.path(plot_dir,paste0("607_group_classical_tongue_activations_",fwhm,".png")))
}

# FIGURE: Group-level permutation testing ----
result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed/"
result_files <- list.files(result_dir, full.names = T) |>
  grep(pattern = "FWHM6", value = T) |>
  grep(pattern = "classical", value = T) |>
  grep(pattern = "visit1", value = T) |>
  grep(pattern = "left", value = T)

subject_estimates <- sapply(result_files, function(x) readRDS(x)$betas_classical$avg$data$cortex_left)
# # 5000 is the number of permutations we will perform
# set.seed(47408)
# reorderings <- sapply(seq(5000), function(s) sample(c(1,-1),size = 45, replace = T), simplify = F)
# library(parallel)
# cl <- makeCluster(6)
# # Number of tasks is 4, there are 4443 locations, and 45 subjects
# null_dist <- parSapplyLB(cl,reorderings, function(ro,se){
#   perm_est <- t(t(se) * ro)
#   perm_est <- array(perm_est, dim = c(4443,4,45))
#   max_t <- apply(perm_est,1:2,function(x) t.test(x)$statistic)
#   return(apply(max_t,2,max))
# }, se = subject_estimates)
# stopCluster(cl)
# saveRDS(null_dist, "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed/03_max_tvalues.rds")
null_dist <- readRDS("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed/03_max_tvalues.rds")
alpha <- 0.01
thresholds <- apply(null_dist, 1, function(x) sort(x, decreasing = T)[ceiling(alpha*length(x))])

actual_t <- sapply(c(0,0.5,1), function(gam){
  out <- apply(array(subject_estimates, dim = c(4443,4,45)),1:2, function(x) t.test(x,mu = gam)$statistic)
  return(out)
}, simplify = F)
active_perm <- sapply(actual_t, function(act_t){
  mapply(function(t_k,thr) as.numeric(t_k > thr), t_k = split(act_t,col(act_t)), thr = thresholds)
}, simplify = F)
active_perm <- Reduce(`+`, active_perm)
active_perm[active_perm == 0] <- NA

# actual_t <- apply(array(subject_estimates, dim = c(4443,4,45)),1:2, function(x) t.test(x)$statistic)
# active_perm <- mapply(function(t_k,thr) as.numeric(t_k > thr), t_k = split(actual_t,col(actual_t)), thr = thresholds)
active_cifti <- readRDS(result_files[1])$betas_classical$avg
active_cifti$data$cortex_left <- active_perm
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
light_orange <- colorRampPalette(c("white","orange"))(3)[2]
my_pal <- c(light_orange,"red","purple")
plot(
  active_cifti,
  hemisphere = "left",
  idx = 4,
  view = 'lateral',
  color_mode = "qualitative",
  colors = my_pal,
  fname = file.path(
    plot_dir,
    "607_group_classical_tongue_activations_permutation.png"
  )
)

# Here's a plot of the FDR activations with group testing
library(BayesfMRI)
fdr_act <- classicalGLM2(results = result_files,
                brainstructure = "cortexL",
                gamma = c(0,0.5,1),
                correction = "FDR",
                alpha = 0.01)
fdr_act <- sapply(fdr_act$active_result, function(x) x$active, simplify = F)
fdr_act <- Reduce(`+`, fdr_act)
fdr_act[fdr_act == 0] <- NA
active_cifti <- readRDS(result_files[1])$betas_classical$avg
active_cifti$data$cortex_left <- fdr_act
plot_dir <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots"
light_orange <- colorRampPalette(c("white","orange"))(3)[2]
my_pal <- c(light_orange,"red","purple")
plot(
  active_cifti,
  hemisphere = "left",
  idx = 4,
  view = 'lateral',
  color_mode = "qualitative",
  colors = my_pal,
  fname = file.path(
    plot_dir,
    "607_group_classical_tongue_activations_fdr.png"
  )
)

# <<<<< NOT USED IN MANUSCRIPT >>>>>----
# FIGURE: First few frames from subject 562345 ----
library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench')
cifti_fname <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/visit2_data/562345/MNINonLinear/Results/tfMRI_MOTOR_LR/tfMRI_MOTOR_LR_Atlas.dtseries.nii"
surfL <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/visit1_data/562345/MNINonLinear/fsaverage_LR32k/562345.L.midthickness.32k_fs_LR.surf.gii"
surfR <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/visit1_data/562345/MNINonLinear/fsaverage_LR32k/562345.R.midthickness.32k_fs_LR.surf.gii"

cifti_obj <- read_cifti(cifti_fname = cifti_fname,surfL_fname = surfL,surfR_fname = surfR)

ntime <- ncol(cifti_obj$data$cortex_left)
for(tt in seq(ntime)) {
  plot(cifti_obj, idx = tt, legend_embed = F,
       fname = paste0("~/Desktop/subject562345_visit2LR/05_frame_",sprintf('%03d',tt),".png"))
}

# FIGURE: Plot of hyperparameter densities across hemispheres ----
result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/PW"
result_files <- list.files(result_dir, full.names = TRUE) |>
  grep(pattern = "103818",value = TRUE) |>
  grep(pattern = "visit1", value = TRUE)
left_result <- readRDS(result_files[1])
right_result <- readRDS(result_files[2])
library(tidyverse)
mutate(left_result$GLMs_Bayesian$cortexL$theta_posteriors, hemisphere = 'left') %>%
  full_join(
    mutate(right_result$GLMs_Bayesian$cortexR$theta_posteriors, hemisphere = 'right')
  )

# FIGURE: Comparison of cutoffs from permutation testing and permutation testing ----
# >> 5k data ----
model_obj <- readRDS("~/Desktop/500_103818_visit1_left_5k_classical_1000permutations_20210908.rds")
model_obj <- model_obj$GLMs_classical$cortexL
sess_ind <- 3
beta_est <- model_obj[[sess_ind]]$estimates
se_beta <- model_obj[[sess_ind]]$SE_estimates
DOF <- model_obj[[sess_ind]]$DOF

alpha <- 0.01
t_stats <- (model_obj[[sess_ind]]$null_estimates) / model_obj[[sess_ind]]$null_SE_estimates
max_tstats <- apply(t_stats,2:3,max, na.rm = TRUE)
null_thresholds <- apply(max_tstats, 1, quantile, probs = (1 - alpha))

n_comps <- sum(!is.na(beta_est[,1]))
bonf_cutoff <- qt(p = 1 - alpha/n_comps,df = DOF,lower.tail = T)

cutoffs_df <- data.frame(
  Task = c("cue","right foot","right hand","tongue"),
  cutoff = c(null_thresholds)
)

col_pal <- rgb(c(230, 86, 0, 0,204),c(159,180,158,114,121),c(0,233,115,178,167), maxColorValue = 255)

library(tidyverse)
g <- reshape2::melt(max_tstats, value.name = "Tstats") %>%
  mutate(Task = c("cue","right foot","right hand","tongue")[Var1],
         Task = factor(Task, levels = c("cue","right foot","right hand","tongue","Bonferroni"))) %>%
  ggplot() +
  geom_density(aes(x = Tstats, fill = Task, color = Task), alpha = 0.3) +
  geom_vline(aes(xintercept = cutoff, color = Task), data = cutoffs_df) +
  geom_vline(aes(xintercept = bonf_cutoff)) +
  labs(x = "Maximum null distribution test statistics", y = "Density", title = "Permutation Cutoffs - 5K data") +
  lims(y = c(0,1.4), x = c(2,7)) +
  geom_text(aes(x = bonf_cutoff, y = 1, label = "Bonferroni"), angle = 90, vjust = 0) +
  geom_hline(aes(yintercept = 0)) +
  scale_fill_discrete("Task") +
  scale_color_discrete("Task") +
  theme_classic()
# g

# >> 5k smoothed data ----
# model_obj <- readRDS("~/Desktop/500_103818_visit1_left_5k_classical_FWHM6_20210928.rds")
model_obj <- readRDS("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/5k_results/smoothed/permutations/500_103818_visit1_left_5k_classical_FWHM6_20211101.rds")
model_obj <- model_obj$GLMs_classical$cortexL
sess_ind <- 3
beta_est <- model_obj[[sess_ind]]$estimates
se_beta <- model_obj[[sess_ind]]$SE_estimates
DOF <- model_obj[[sess_ind]]$DOF

alpha <- 0.01
t_stats <- (model_obj[[sess_ind]]$null_estimates) / model_obj[[sess_ind]]$null_SE_estimates
max_tstats <- apply(t_stats,2:3,max, na.rm = TRUE)
null_thresholds <- apply(max_tstats, 1, quantile, probs = (1 - alpha))

n_comps <- sum(!is.na(beta_est[,1]))
bonf_cutoff <- qt(p = 1 - alpha/n_comps,df = DOF,lower.tail = T)

cutoffs_df <- data.frame(
  Task = c("cue","right foot","right hand","tongue"),
  cutoff = c(null_thresholds)
)

col_pal <- rgb(c(230, 86, 0, 0,204),c(159,180,158,114,121),c(0,233,115,178,167), maxColorValue = 255)

library(tidyverse)
g4 <- reshape2::melt(max_tstats, value.name = "Tstats") %>%
  mutate(Task = c("cue","right foot","right hand","tongue")[Var1],
         Task = factor(Task, levels = c("cue","right foot","right hand","tongue","Bonferroni"))) %>%
  ggplot() +
  geom_density(aes(x = Tstats, fill = Task, color = Task), alpha = 0.3) +
  geom_vline(aes(xintercept = cutoff, color = Task), data = cutoffs_df) +
  geom_vline(aes(xintercept = bonf_cutoff)) +
  labs(x = "Maximum null distribution test statistics", y = "Density", title = "Permutation Cutoffs - 5K smoothed data") +
  lims(y = c(0,1.4), x = c(2,7)) +
  geom_text(aes(x = bonf_cutoff, y = 1, label = "Bonferroni"), angle = 90, vjust = 0) +
  geom_hline(aes(yintercept = 0)) +
  scale_fill_discrete("Task") +
  scale_color_discrete("Task") +
  theme_classic()
# g4

# >> 32k data ----
model_obj <- readRDS("~/Desktop/500_103818_visit1_left_32k_classical_1000permutations_20210916.rds")
model_obj <- model_obj$GLMs_classical$cortexL
sess_ind <- 3
beta_est <- model_obj[[sess_ind]]$estimates
se_beta <- model_obj[[sess_ind]]$SE_estimates
DOF <- model_obj[[sess_ind]]$DOF

alpha <- 0.01
t_stats <- (model_obj[[sess_ind]]$null_estimates) / model_obj[[sess_ind]]$null_SE_estimates
max_tstats <- apply(t_stats,2:3,max, na.rm = TRUE)
null_thresholds <- apply(max_tstats, 1, quantile, probs = (1 - alpha))

n_comps <- sum(!is.na(beta_est[,1]))
bonf_cutoff_32k <- qt(p = 1 - alpha/n_comps,df = DOF,lower.tail = T)

cutoffs_df <- data.frame(
  Task = c("cue","right foot","right hand","tongue"),
  cutoff = c(null_thresholds)
)

col_pal <- rgb(c(230, 86, 0, 0,204),c(159,180,158,114,121),c(0,233,115,178,167), maxColorValue = 255)

library(tidyverse)
g2 <- reshape2::melt(max_tstats, value.name = "Tstats") %>%
  mutate(Task = c("cue","right foot","right hand","tongue")[Var1],
         Task = factor(Task, levels = c("cue","right foot","right hand","tongue","Bonferroni"))) %>%
  ggplot() +
  geom_density(aes(x = Tstats, fill = Task, color = Task), alpha = 0.3) +
  geom_vline(aes(xintercept = cutoff, color = Task), data = cutoffs_df) +
  geom_vline(aes(xintercept = bonf_cutoff_32k)) +
  labs(x = "Maximum null distribution test statistics", y = "Density", title = "Permutation Cutoffs - 32K data") +
  lims(y = c(0,1.4), x = c(2,7)) +
  geom_text(aes(x = bonf_cutoff_32k, y = 1, label = "Bonferroni"), angle = 90, vjust = 0) +
  geom_hline(aes(yintercept = 0)) +
  scale_fill_discrete("Task") +
  scale_color_discrete("Task") +
  theme_classic()
# g2

# >> 32k smoothed data ----
model_obj <- readRDS("~/Desktop/500_103818_visit1_left_32k_smoothed_FWHM5_classical_1000permutations_20210916.rds")
model_obj <- model_obj$GLMs_classical$cortexL
sess_ind <- 3
beta_est <- model_obj[[sess_ind]]$estimates
se_beta <- model_obj[[sess_ind]]$SE_estimates
DOF <- model_obj[[sess_ind]]$DOF

alpha <- 0.01
t_stats <- (model_obj[[sess_ind]]$null_estimates) / model_obj[[sess_ind]]$null_SE_estimates
max_tstats <- apply(t_stats,2:3,max, na.rm = TRUE)
null_thresholds <- apply(max_tstats, 1, quantile, probs = (1 - alpha))

n_comps <- sum(!is.na(beta_est[,1]))
bonf_cutoff_32k <- qt(p = 1 - alpha/n_comps,df = DOF,lower.tail = T)

cutoffs_df <- data.frame(
  Task = c("cue","right foot","right hand","tongue"),
  cutoff = c(null_thresholds)
)

col_pal <- rgb(c(230, 86, 0, 0,204),c(159,180,158,114,121),c(0,233,115,178,167), maxColorValue = 255)

library(tidyverse)
g3 <- reshape2::melt(max_tstats, value.name = "Tstats") %>%
  mutate(Task = c("cue","right foot","right hand","tongue")[Var1],
         Task = factor(Task, levels = c("cue","right foot","right hand","tongue","Bonferroni"))) %>%
  ggplot() +
  geom_density(aes(x = Tstats, fill = Task, color = Task), alpha = 0.3) +
  geom_vline(aes(xintercept = cutoff, color = Task), data = cutoffs_df) +
  geom_vline(aes(xintercept = bonf_cutoff_32k)) +
  labs(x = "Maximum null distribution test statistics", y = "Density", title = "Permutation Cutoffs - 32K smoothed data") +
  lims(y = c(0,1.4),x = c(2,7)) +
  geom_text(aes(x = bonf_cutoff_32k, y = 1, label = "Bonferroni"), angle = 90, vjust = 0) +
  geom_hline(aes(yintercept = 0)) +
  scale_fill_discrete("Task") +
  scale_color_discrete("Task") +
  theme_classic()
# g3

# >> Put the three together ----
library(gridExtra)
g_all <- grid.arrange(g + theme(legend.position = 'none'),g2 + theme(legend.position = 'none'),g4 + theme(legend.position = 'none'),g3 + theme(legend.position = 'none'), nrow = 2, ncol = 2)
g_all

ggsave(filename = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/05_permutation_cutoffs.png",plot = g_all, width = 10, height = 8)

# FIGURE: Residual plot ----
result <- readRDS("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/HCP_results/5k_results/individual/PW/500_103818_visit1_left_5k_20210125.rds")
y <- result$GLMs_Bayesian$cortexL$y
X <- Reduce(rbind,result$GLMs_Bayesian$cortexL$X)
beta <- c(result$betas_Bayesian$avg$data$cortex_left)
pred <- (X %*% beta)@x
resid <- y - pred
resid_mat <- matrix(resid, 4443, 568)

cifti_obj <- result$betas_Bayesian$avg
cifti_obj$data$cortex_left[,1] <- resid_mat[,1]

library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
plot(cifti_obj, hemisphere = 'left',idx = 1, view = 'lateral')


# FIGURE: Comparison of smoothed/unsmoothed 5k/32k results ----
result_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Gambling_Task_Dan"
subject_files <- list.files(result_dir, full.names = TRUE) |>
  grep(pattern = "103818", value = T)

cifti_list <- list(
  bayes_5k = readRDS(grep("runLR",grep(
    "5k_Bayes_classical_gambling", subject_files, value = T
  ), value = T))$betas_Bayesian$LR,
  classical_5k = readRDS(grep("runLR",grep(
    "5k_Bayes_classical_gambling", subject_files, value = T
  ), value = T))$betas_classical$LR,
  classical_5k_smoothed = readRDS(
    grep("5k_smoothed_classical_gambling", subject_files, value = T)
  )$betas_classical$LR,
  bayes_32k = readRDS(grep(
    "32k_Bayes_classical_gambling", subject_files, value = T
  ))$betas_Bayesian$LR,
  classical_32k = readRDS(grep(
    "32k_classical_gambling", subject_files, value = T
  ))$betas_classical$LR,
  classical_32k_smoothed = readRDS(
    grep("32k_smoothed_classical_gambling", subject_files, value = T)
  )$betas_classical$LR
)

library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')

mapply(function(cif,nm) {
  plot(
    cif,
    hemisphere = 'left',
    view = 'lateral',
    idx = 1,
    zlim = c(-1, 1),
    legend_embed = F,
    surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
    surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
    fname = paste0("/Volumes/GoogleDrive/My Drive/MEJIA_LAB_Dan/BayesGLM_Validation/plots/05_",nm,"_103818_gambling_estimate.png")
  )
}, cif = cifti_list, nm = names(cifti_list))



# FIGURE: Full-res classical ----
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
         fname = paste0("plots/5_classical_",
                        subject,"_",task_names[task_idx],
                        "_32k_unsmoothed_estimate.png"))
  }
  for(task_idx in 2:3) {
    for(hem in c('left','right')) {
      plot(cifti_obj, idx = task_idx, zlim = c(-1,1), legend_embed = F,
           hemisphere = hem,
           surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
           surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
           fname = paste0("plots/5_classical_",
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
    cifti_target_fname = paste0("HCP_results/32k_results/5_",
                                subject,"_visit1_32k_smoothed_estimates.dscalar.nii"),
    surf_FWHM = 6,
    surfL_fname = paste0("HCP_data/subject_surfaces/",subject,
                         ".L.midthickness.32k_fs_LR.surf.gii"),
    surfR_fname = paste0("HCP_data/subject_surfaces/",subject,
                         ".R.midthickness.32k_fs_LR.surf.gii")
  )
  cifti_obj <- read_cifti(paste0("HCP_results/32k_results/5_", subject,
                                 "_visit1_32k_smoothed_estimates.dscalar.nii"))
  for(task_idx in c(1,4)) {
    plot(cifti_obj, idx = task_idx, zlim = c(-1,1), legend_embed = F,
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         fname = paste0("plots/5_classical_",
                        subject,"_",task_names[task_idx],
                        "_32k_smoothed_estimate.png"))
  }
  for(task_idx in 2:3) {
    for(hem in c('left','right')) {
      plot(cifti_obj, idx = task_idx, zlim = c(-1,1), legend_embed = F,
           hemisphere = hem,
           surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
           surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
           fname = paste0("plots/5_classical_",
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
         fname = paste0("plots/5_smoothed_classical_",
                        subject,"_",task_names[task_idx],
                        "_32k_estimate.png"))
  }
  for(task_idx in 2:3) {
    for(hem in c('left','right')) {
      plot(cifti_obj, idx = task_idx, zlim = c(-1,1), legend_embed = F,
           hemisphere = hem,
           surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
           surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
           fname = paste0("plots/5_smoothed_classical_",
                          subject,"_",hem,"_",task_names[task_idx],
                          "_32k_estimate.png"))
    }
  }
}

# FIGURE: Spherical vs midthickness surfaces----
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
             fname = paste0("plots/5_bayes_",c('midthickness','spherical','diff')[surf],"_",subject,"_",visit,"_",task_names[task_idx],"_estimate.png"))
      }
      for(task_idx in 2:3) {
        for(hem in c('left','right')) {
          plot(surf_ciftis[[surf]], idx = task_idx, zlim = zlims[[surf]], hemisphere = hem,
               legend_embed = F,
               surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
               surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
               fname = paste0("plots/5_bayes_",c('midthickness','spherical','diff')[surf],"_",subject,"_",visit,"_",hem,"_",task_names[task_idx],"_estimate.png"))
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
             fname = paste0("plots/5_bayes_",c('midthickness','spherical','diff')[surf],"_",subject,"_",visit,"_",task_names[task_idx],"_activation.png"))
      }
      for(task_idx in 2:3) {
        for(hem in c('left','right')) {
          use_colors <-
            sort(unique(surf_ciftis[[surf]]$data[[paste0("cortex_", hem)]][, task_idx]))
          plot(surf_ciftis[[surf]], idx = task_idx, color_mode = "qualitative",
               hemisphere = hem, legend_embed = F, colors = col_pal[use_colors],
               surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
               surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
               fname = paste0("plots/5_bayes_",c('midthickness','spherical','diff')[surf],"_",subject,"_",visit,"_",hem,"_",task_names[task_idx],"_activation.png"))
        }
      }
    }
  }
}

# FIGURE: ICC per task ----
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
avg_estimates_classical <- readRDS(file.path(result_dir,"5_avg_estimates_classical.rds"))
ICC_weights <- sapply(avg_estimates_classical, function(hem_est) {
  combined_est <- Reduce(abind,hem_est)
  avg_est <- apply(combined_est,1:2,mean)
  weight_est <- apply(avg_est,2,function(avg_task) abs(avg_task) / sum(abs(avg_task)))
  return(weight_est)
}, simplify = F)

# >> Bayesian ----
result_dir <- "HCP_results/5k_results"
avg_estimates <- readRDS(file.path(result_dir,"5_avg_estimates.rds"))
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
bayesian_icc_avg_cifti <- readRDS("HCP_data/cifti_5k_template_whole.rds")
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
    fname = paste0("plots/5_",task_file_names[task_idx],"_icc_bayesian.png")
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
      fname = paste0("plots/5_",hem[h],"_",task_file_names[ti],"_icc_bayesian.png")
    )
  }
}

# >> Classical ----
avg_estimates_classical <- readRDS(file.path(result_dir,"5_avg_estimates_classical.rds"))
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

classical_icc_avg_cifti <- readRDS("HCP_data/cifti_5k_template_whole.rds")
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
    fname = paste0("plots/5_",task_file_names[task_idx],"_icc_classical")
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
      fname = paste0("plots/5_",hem[h],"_",task_file_names[ti],"_icc_classical.png")
    )
  }
}

# FIGURE: weighted ICC by task ----

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
avg_estimates_classical <- readRDS(file.path(result_dir,"5_avg_estimates_classical.rds"))
ICC_weights <- sapply(avg_estimates_classical, function(hem_est) {
  combined_est <- Reduce(abind,hem_est)
  avg_est <- apply(combined_est,1:2,mean)
  weight_est <- apply(avg_est,2,function(avg_task) abs(avg_task) / sum(abs(avg_task)))
  return(weight_est)
}, simplify = F)


# >> Bayesian ----
result_dir <- "HCP_results/5k_results"
session_estimates <- readRDS(file.path(result_dir,"5_session_estimates.rds"))
avg_estimates <- readRDS(file.path(result_dir,"5_avg_estimates.rds"))
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
avg_estimates_classical <- readRDS(file.path(result_dir,"5_avg_estimates_classical.rds"))
session_estimates_classical <- readRDS(file.path(result_dir,"5_session_estimates_classical.rds"))
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

ggsave("plots/5_wicc_task_scatterplot.png", plot = wicc_scatter, width = 5, height = 4)


# FIGURE: Group Estimates and Activations----
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

bayes_estimates <- readRDS("HCP_data/cifti_5k_template_whole.rds")
bayes_estimates$data$cortex_left <- left_result$estimates
bayes_estimates$data$cortex_right <- right_result$estimates

plot(bayes_estimates, idx = 4, zlim = c(-1,1),
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
     fname = "plots/5_group_tongue_estimates.png")

# >>>> Classical ----
class_group <- readRDS("HCP_results/5k_results/502_HCP_classical_group_PW_estimates.rds")
class_estimates <- readRDS("HCP_data/cifti_5k_template_whole.rds")
class_estimates$data$cortex_left <- class_group$left
class_estimates$data$cortex_right <- class_group$right
plot(class_estimates, zlim = c(-1,1), idx = 4,
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
     fname = "plots/5_group_tongue_estimates_classical.png")

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
cifti_obj <- readRDS("HCP_data/cifti_5k_template_whole.rds")
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
         "~/Desktop/SMI 2021 poster images/5_",
         model_names[i],
         "_",task_names[task_idx],
         ".png"
         ),
       view = "lateral", legend_embed = F, zlim = c(-1,1),
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")
}

# FIGURE: ICC Thresholded Images ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench/')


task_file_names <- c("visual_cue","foot","hand","tongue")
bayesian_icc_avg_cifti <- readRDS("HCP_data/cifti_5k_template_whole.rds")
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
    fname = paste0("plots/5_",task_file_names[task_idx],"_icc_quality_bayesian.png")
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
      fname = paste0("plots/5_",hem[h],"_",task_file_names[ti],"_icc_quality_bayesian.png")
    )
  }
}

# >> Classical ----
result_dir <- "HCP_results/5k_results"
avg_estimates <- readRDS(file.path(result_dir,"5_average_estimates_classical.rds"))
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
classical_icc_avg_cifti <- readRDS("HCP_data/cifti_5k_template_whole.rds")
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
    fname = paste0("plots/5_",task_file_names[task_idx],"_icc_quality_classical.png")
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
      fname = paste0("plots/5_",hem[h],"_",task_file_names[ti],"_icc_quality_classical.png")
    )
  }
}





# FIGURE: Group areas of activation ----


# FIGURE: Single-session estimates ----

# FIGURE: Single-session activations ----
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
      cifti_obj <- readRDS("HCP_data/cifti_5k_template_whole.rds")
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
              "plots/5_bayes_",
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
                "plots/5_bayes_",
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
        cifti_obj <- readRDS("HCP_data/cifti_5k_template_whole.rds")
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
              "plots/5_classical_",
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
                "plots/5_classical_",
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

# FIGURE: Single- vs. Multi-session activation ----
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

# FIGURE: Dice x area for single vs. multisession data ----
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

# FIGURE: Motor Cortex Parcellation ----
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
  cifti_out <- readRDS("HCP_data/cifti_5k_template_whole.rds")
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

# FIGURE: AR coefficient map ----
left_obj <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session/500_103818_visit1_left_5k_sessionLR_20210322.rds")
right_obj <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session/500_103818_visit1_right_5k_sessionLR_20210322.rds")
cifti_obj <- readRDS("HCP_data/cifti_5k_template_whole.rds")
cifti_obj$data$cortex_left <- left_obj$prewhitening_info$left$AR_coeffs
cifti_obj$data$cortex_right <- right_obj$prewhitening_info$right$AR_coeffs
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
plot(cifti_obj, hemisphere = 'left',idx = 1,
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")

# FIGURE: Mesh structure ----
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
  fname = "plots/5_mesh_cifti.png")

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
  fname = "plots/5_spherical_mesh_cifti.png")

# FIGURE: Full-resolution classical results ----
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
               fname = paste0("plots/5_15k_classical_",
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
                 fname = paste0("plots/5_15k_classical_",
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
         fname = paste0("plots/5_classical_",subject,"_",visit,"_15k_estimate.png"))
    # }
  }
}

#  FIGURES: Number of runs x Number of activations (Subject) ----
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

ggsave("plots/5_num_act_plots_thr0.png", plot = num_act_plots[[1]], width = 7, height = 6)
ggsave("plots/5_num_act_plots_thr0.5.png", plot = num_act_plots[[2]], width = 7, height = 6)
