# This is a script to calculate the intra-class correlation coefficient (ICC)
# for the HCP data. This will be done using averages to calculate the between-
# subject variance and the individual sessions to calculate the within-subject
# variation.

result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
result_files <- list.files(result_dir,full.names = TRUE)
result_files <- grep(".rds", result_files, value = T)
# Grab estimates ----
# >>  Bayesian ----
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
# saveRDS(session_estimates, file.path(result_dir,"602_session_estimates.rds"))
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
saveRDS(avg_estimates, file.path(result_dir,"602_avg_PW_estimates.rds"))

# >>  Classical ----
# Individual sessions
session_estimates_classical <- sapply(c('left','right'), function(hem) {
  out_list_l1 <- sapply(c('visit1','visit2'), function(visit) {
    file_locs <- grep(hem, grep(visit, result_files,value = TRUE), value = TRUE)
    session_estimates <-
      sapply(file_locs,
             function(file_n){
               result_obj <- readRDS(file_n)
               which_cortex <- paste0('cortex_',hem)
               dim_est <- dim(result_obj$betas_classical$LR$data[[which_cortex]])
               out_obj <- array(NA,dim = c(dim_est,2)) # Third index corresponds to session (1-LR, 2-RL)
               out_obj[,,1] <- result_obj$betas_classical$LR$data[[which_cortex]]
               out_obj[,,2] <- result_obj$betas_classical$RL$data[[which_cortex]]
               return(out_obj)
             }, simplify = 'array')
  }, simplify = F)
  return(out_list_l1)
}, simplify = F)
# saveRDS(session_estimates_classical, file.path(result_dir,"602_session_estimates_classical.rds"))
# Averages
avg_estimates_classical <- sapply(c('left','right'), function(hem) {
  out_list_l1 <- sapply(c('visit1','visit2'), function(visit) {
    file_locs <- grep(hem, grep(visit, result_files,value = TRUE), value = TRUE)
    session_estimates <-
      sapply(file_locs,
             function(file_n){
               result_obj <- readRDS(file_n)
               which_cortex <- paste0('cortex_',hem)
               out_obj <- result_obj$betas_classical$avg$data[[which_cortex]]
               return(out_obj)
             }, simplify = 'array')
  }, simplify = F)
  return(out_list_l1)
}, simplify = F)
saveRDS(avg_estimates_classical, file.path(result_dir,"602_avg_estimates_PW_classical.rds"))

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

# ICC Calculation ----
# >> Bayesian ----
session_estimates <- readRDS(file.path(result_dir,"602_session_estimates.rds"))
avg_estimates <- readRDS(file.path(result_dir,"602_avg_estimates.rds"))
# >>>> ICC using session estimates ----
library(abind)
abind_f <- function(...) abind(..., along = 3)
ICC_values_separate <- sapply(session_estimates, function(hem_est) {
  combined_visits <- Reduce(abind_f, hem_est)
  vertex_ICC <- apply(combined_visits, 1:2, function(tX) {
    return(ICC(t(tX)))
  })
  return(vertex_ICC)
}, simplify = F)
sapply(ICC_values_separate, function(x) summary(c(x)))
par(mfrow=c(2,4))
for(h in c('left','right')) {
  for(ta in 1:4) {
    ta_name <- c('visual','foot','hand','tongue')
    hist(ICC_values_separate[[h]][,ta], main = paste(h,ta_name[ta]), xlab = 'ICC')
  }
}

# >>>> ICC using average estimates ----
library(abind)
abind_g <- function(...) abind(...,along = 4)
ICC_values_average <- sapply(avg_estimates, function(hem_est) {
  combined_visits <- Reduce(abind_g,hem_est)
  vertex_ICC <- apply(combined_visits, 1:2, ICC)
  return(vertex_ICC)
}, simplify = F)
sapply(ICC_values_average, function(x) summary(c(x)))
par(mfrow=c(2,4))
for(h in c('left','right')) {
  for(ta in 1:4) {
    ta_name <- c('visual','foot','hand','tongue')
    hist(ICC_values_average[[h]][,ta], main = paste(h,ta_name[ta]), xlab = 'ICC')
  }
}

# >> Classical ----
session_estimates_classical <- readRDS(file.path(result_dir,"602_session_estimates_classical.rds"))
avg_estimates_classical <- readRDS(file.path(result_dir,"602_avg_estimates_classical.rds"))
# >>>> ICC using session estimates ----
library(abind)
abind_f <- function(...) abind(..., along = 3)
ICC_values_separate_classical <- sapply(session_estimates_classical, function(hem_est) {
  combined_visits <- Reduce(abind_f, hem_est)
  vertex_ICC <- apply(combined_visits, 1:2, function(tX) {
    return(ICC(t(tX)))
  })
  return(vertex_ICC)
}, simplify = F)
sapply(ICC_values_separate_classical, function(x) summary(c(x)))
par(mfrow=c(2,4))
for(h in c('left','right')) {
  for(ta in 1:4) {
    ta_name <- c('visual','foot','hand','tongue')
    hist(ICC_values_separate_classical[[h]][,ta], main = paste(h,ta_name[ta]), xlab = 'ICC')
  }
}

# >>>> ICC using average estimates ----
library(abind)
abind_g <- function(...) abind(...,along = 4)
ICC_values_average_classical <- sapply(avg_estimates_classical, function(hem_est) {
  combined_visits <- Reduce(abind_g,hem_est)
  vertex_ICC <- apply(combined_visits, 1:2, ICC)
  return(vertex_ICC)
}, simplify = F)
sapply(ICC_values_average_classical, function(x) summary(c(x)))
par(mfrow=c(2,4))
for(h in c('left','right')) {
  for(ta in 1:4) {
    ta_name <- c('visual','foot','hand','tongue')
    hist(ICC_values_average_classical[[h]][,ta], main = paste(h,ta_name[ta]), xlab = 'ICC')
  }
}
# Weighted ICC Calculation ----
# The weights here all come from the averaged classical estimates, as they are unbiased
ICC_weights <- sapply(avg_estimates_classical, function(hem_est) {
  combined_est <- Reduce(abind,hem_est)
  avg_est <- apply(combined_est,1:2,mean)
  weight_est <- apply(avg_est,2,function(avg_task) abs(avg_task) / sum(abs(avg_task)))
  return(weight_est)
}, simplify = F)
# >> Bayesian ----
bayes_wICC_separate <- round(apply(mapply(function(wt,icc) {
  return(apply(wt*icc,2,sum))
}, wt = ICC_weights, icc = ICC_values_separate),1,sum),3)
bayes_wICC_avg <- round(apply(mapply(function(wt,icc) {
  return(apply(wt*icc,2,sum))
}, wt = ICC_weights, icc = ICC_values_average),1,sum),3)
# >> Classical ----
classical_wICC_separate <- round(apply(mapply(function(wt,icc) {
  return(apply(wt*icc,2,sum))
}, wt = ICC_weights, icc = ICC_values_separate_classical),1,sum),3)
classical_wICC_avg <- round(apply(mapply(function(wt,icc) {
  return(apply(wt*icc,2,sum))
}, wt = ICC_weights, icc = ICC_values_average_classical),1,sum),3)

# Plotting ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
plot_dir <- "~/github/BayesGLM_Validation/HCP_results/plots"
# >> Bayesian ----
# >>>> Separate Estimates ----
bayesian_icc_sep_cifti <- readRDS(result_files[1])$betas_Bayesian[[1]]
bayesian_icc_sep_cifti$data$cortex_left <- ICC_values_separate$left
bayesian_icc_sep_cifti$data$cortex_right <- ICC_values_separate$right
bayesian_icc_sep_cifti$meta$cortex$medial_wall_mask$right <-
  readRDS(result_files[2])$betas_Bayesian[[1]]$meta$cortex$medial_wall_mask$right
plot(
  bayesian_icc_sep_cifti,
  idx = 1,
  title = paste("Bayesian Visual Cue ICC, wICC =", bayes_wICC_separate[1]),
  color_mode = 'sequential',
  colors = c('black','red','yellow','white'),
  zlim = c(0,0.3,0.6,0.9),
  save = TRUE,
  surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
  surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
  fname = "plots/602_visual_cue_icc_bayesian_separate"
)

# >>>> Average Estimates ----
task_file_names <- c("visual_cue","foot","hand","tongue")
bayesian_icc_avg_cifti <- readRDS(result_files[1])$betas_Bayesian[[1]]
bayesian_icc_avg_cifti$data$cortex_left <- ICC_values_average$left
bayesian_icc_avg_cifti$data$cortex_right <- ICC_values_average$right
bayesian_icc_avg_cifti$meta$cortex$medial_wall_mask$right <-
  readRDS(result_files[2])$betas_Bayesian[[1]]$meta$cortex$medial_wall_mask$right
for(i in 1:4) {
  plot(
    bayesian_icc_avg_cifti,
    idx = i,
    hemisphere = 'left',
    title = paste("wICC =", bayes_wICC_avg[i]),
    color_mode = 'sequential',
    colors = c('black','red','yellow','white'),
    zlim = c(0,0.3,0.6,0.9),
    save = TRUE,
    surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
    surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
    fname = paste0("plots/602_",task_file_names[i],"_icc_bayesian_avg_left")
  )
}


# >> Classical ----
# >>>> Separate Estimates ----
classical_icc_sep_cifti <- readRDS(result_files[1])$betas_classical[[1]]
classical_icc_sep_cifti$data$cortex_left <- ICC_values_separate_classical$left
classical_icc_sep_cifti$data$cortex_right <- ICC_values_separate_classical$right
classical_icc_sep_cifti$meta$cortex$medial_wall_mask$right <-
  readRDS(result_files[2])$betas_classical[[1]]$meta$cortex$medial_wall_mask$right
plot(
  classical_icc_sep_cifti,
  idx = 1,
  title = paste("Classical Visual Cue ICC, wICC =", classical_wICC_separate[1]),
  color_mode = 'sequential',
  colors = c('black','red','yellow','white'),
  zlim = c(0,0.3,0.6,0.9),
  save = TRUE,
  surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
  surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
  fname = "plots/602_visual_cue_icc_classical_separate"
)

# >>>> Average Estimates ----
classical_icc_avg_cifti <- readRDS(result_files[1])$betas_classical[[1]]
classical_icc_avg_cifti$data$cortex_left <- ICC_values_average_classical$left
classical_icc_avg_cifti$data$cortex_right <- ICC_values_average_classical$right
classical_icc_avg_cifti$meta$cortex$medial_wall_mask$right <-
  readRDS(result_files[2])$betas_classical[[1]]$meta$cortex$medial_wall_mask$right
for(i in 1:4) {
  plot(
    classical_icc_avg_cifti,
    idx = i,
    hemisphere = 'left',
    title = paste("wICC =", classical_wICC_avg[i]),
    color_mode = 'sequential',
    colors = c('black','red','yellow','white'),
    zlim = c(0,0.3,0.6,0.9),
    save = TRUE,
    surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
    surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
    fname = paste0("plots/602_",task_file_names[i],"_icc_classical_avg_left")
  )
}
