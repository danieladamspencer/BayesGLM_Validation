# This is a script meant to provide reproducible code for each figure in the
# manuscript

# FIGURE 1: Single-subject estimates ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench/')
library(INLA)
inla.setOption(pardiso.lic = "~/pardiso.lic")
library(BayesfMRI)
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
load("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
result_files <- list.files(result_dir, full.names = T)
result_files <- c(sapply(subjects,grep, x = result_files, value = T))
vis <- "visit1"
task_idx <- 4
task_name <- "Tongue"
result_files <- grep(vis, result_files, value = T)
# >> Bayesian ----
for(s in subjects) {
  cifti_template <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
  left_result <- readRDS(grep('left',grep(s,result_files, value = T), value = T))
  cifti_template$data$cortex_left <- left_result$betas_Bayesian$avg$data$cortex_left
  rm(left_result)
  right_result <- readRDS(grep('right',grep(s,result_files, value = T), value = T))
  cifti_template$data$cortex_right <- right_result$betas_Bayesian$avg$data$cortex_right
  rm(right_result)
  plot(cifti_template, idx = task_idx, #title = paste("Subject",s,task_name,"Task"),
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
       fname = paste0("plots/600_subject_",s,"_",tolower(task_name),"_estimates.png"),
       zlim = c(-1,1))
}
# >> Classical ----
for(s in subjects) {
  cifti_template <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
  left_result <- readRDS(grep('left',grep(s,result_files, value = T), value = T))
  cifti_template$data$cortex_left <- left_result$betas_classical$avg$data$cortex_left
  rm(left_result)
  right_result <- readRDS(grep('right',grep(s,result_files, value = T), value = T))
  cifti_template$data$cortex_right <- right_result$betas_classical$avg$data$cortex_right
  rm(right_result)
  plot(cifti_template, idx = task_idx, #title = paste("Subject",s,task_name,"Task"),
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
       fname = paste0("plots/600_subject_",s,"_",tolower(task_name),"_classical_estimates.png"),
       zlim = c(-1,1))
}

# FIGURE 2: Single-subject areas of activation ----
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

# FIGURE 3: Single-subject activation overlap ----
# >> Bayesian ----
result_dir <- "HCP_results/5k_results/individual/PW/activations"
load("HCP_data/subjects.Rdata")
subjects <- subjects[c(1,2)]
task_idx <- 4 # This is for the tongue
thr <- c(0,0.5)
thr_chr <- paste0("_thr",sub("\\.","",thr),"_")
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
subject_files <- unlist(sapply(thr_chr, grep, x = subject_files,  value=T))
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

# >> Classical ----
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

# FIGURE 6: Dice Bootstrap CIs and Area x Dice ----
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
dice_result_files <- grep("604_", dice_result_files, value = T)
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

library(tidyverse)
dice_ci_plot <- reshape2::melt(bstrap_ci_results) %>%
  mutate(task = task_names[Var2],
         threshold = factor(L2, levels = c("0%","0.5%","1%"))) %>%
  pivot_wider(names_from = Var1, values_from = value) %>%
  ggplot() +
  geom_errorbar(aes(ymin = Low, ymax = High, x = threshold, color = L1), position = position_dodge(width = 0.2), width = 0.2) +
  geom_point(aes(y = Mean, x = threshold, color = L1), position = position_dodge(width = 0.2)) +
  facet_wrap(vars(task), scales = "free", nrow = 3, ncol = 2) +
  labs(x = "Activation Threshold", y = "Dice Coefficient") +
  scale_color_discrete("") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "top")

dice_ci_plot

ggsave("plots/607_dice_ci_plot.png", width = 3, height = 7, plot = dice_ci_plot)


area_result_files <- list.files("HCP_results/5k_results", full.names = T)
area_result_files <- grep("604_", area_result_files, value = T)
area_result_files <- grep("_PW_", area_result_files, value = T)
area_result_files <- grep("_overlap_", area_result_files, value = T)
area_files <- list(Bayesian = grep("classical",area_result_files, invert = T, value = T)[c(2,1,3)],
                   Classical = grep("FWER",area_result_files, value = T)[c(2,1,3)])
area_results <- sapply(area_files, function(df) sapply(df,readRDS, simplify = F), simplify = F)
names(area_results$Bayesian) <- c("0%","0.5%","1%")
names(area_results$Classical) <- c("0%","0.5%","1%")

area_by_dice <- reshape2::melt(dice_results, value.name = "dice") %>%
  cbind(overlap = reshape2::melt(area_results, value.name = "overlap")$overlap) %>%
  mutate(Var2 = task_names[Var2],
         Var2 = as.factor(Var2),
         L2 = factor(L2, levels = c("0%","0.5%","1%"))) %>%
  ggplot() +
  geom_point(aes(x = overlap/100, y = dice, color = L1)) +
  facet_grid(Var2 ~ L2, scales = "free_x") +
  scale_color_discrete("") +
  labs(y = "Dice Coefficient", x = "Size of Overlap (100s of vertices)") +
  theme_bw() +
  theme(legend.position = "top", panel.grid = element_blank())

area_by_dice

ggsave("plots/607_area_by_dice.png", width = 4.5, height = 7, plot = area_by_dice)

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

# FIGURE 8: ICC Thresholded Images ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench/')
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

# FIGURE 9: ICC quality barplots ----
task_file_names <- c("visual_cue","foot","hand","tongue")
task_names <- c("Visual Cue","Foot","Hand","Tongue")

library(tidyverse)
ICC_df <- reshape2::melt(ICC_quality) %>%
  mutate(Model = "Bayesian") %>%
  full_join(
    reshape2::melt(ICC_quality_classical) %>%
      mutate(Model = "Classical")
  ) %>%
  mutate(task = task_names[Var2],
         quality = factor(value, labels = c("Poor","Fair","Good","Excellent")))

ICC_quality_df <- ICC_df %>%
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
  scale_fill_discrete("ICC Quality") +
  labs(x = "", y = "Proportion of Data Locations") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_flip()

icc_quality_plot

ggsave("plots/607_icc_quality_plot.png",plot = icc_quality_plot, height = 2.5, width = 10)

# FIGURE 10: Group Estimates ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench/')
task_file_names <- c("visual_cue","foot","hand","tongue")
# >> Bayesian ----
group_results_left <- readRDS("HCP_results/5k_results/group/501_HCP_BayesGLM2_45subj_result_left_thresh0_20210317.rds")
group_results_right <- readRDS("HCP_results/5k_results/group/501_HCP_BayesGLM2_45subj_result_right_thresh0_20210318.rds")
cifti_both <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
cifti_both$data$cortex_left <- as.matrix(group_results_left$estimates)
cifti_both$data$cortex_right <- as.matrix(group_results_right$estimates)
# >>>>  Visual Cue and Tongue (idx = c(1,4)) ----
for(i in c(1,4)) {
  plot(cifti_both, zlim = c(-1,1), idx = i, hemisphere = "both",
       legend_embed = F,
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
       fname = paste0("plots/607_group_bayes_",task_file_names[i],"_estimate.png"))
}

# >>>> Foot and Hand Tasks (idx = c(2,3)) ----
for(i in c(2,3)) {
  for(hem in c("left","right")) {
    plot(cifti_both, zlim = c(-1,1), idx = i, hemisphere = hem,
         legend_embed = F,
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         fname = paste0("plots/607_group_bayes_",hem,"_",task_file_names[i],"_estimate.png"))
  }
}

# >> Classical ----
group_results_classical <- readRDS("HCP_results/5k_results/group/502_HCP_classical_group_estimates.rds")
cifti_classical <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
cifti_classical$data$cortex_left <- as.matrix(group_results_classical$left)
cifti_classical$data$cortex_right <- as.matrix(group_results_classical$right)
# >>>>  Visual Cue and Tongue (idx = c(1,4)) ----
for(i in c(1,4)) {
  plot(cifti_classical, zlim = c(-1,1), idx = i, hemisphere = "both",
       legend_embed = F,
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
       fname = paste0("plots/607_group_classical_",task_file_names[i],"_estimate.png"))
}

# >>>> Foot and Hand Tasks (idx = c(2,3)) ----
for(i in c(2,3)) {
  for(hem in c("left","right")) {
    plot(cifti_classical, zlim = c(-1,1), idx = i, hemisphere = hem,
         legend_embed = F,
         surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
         surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
         fname = paste0("plots/607_group_classical_",hem,"_",task_file_names[i],"_estimate.png"))
  }
}

# FIGURE 11: Group areas of activation ----
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
# >> Bayesian ----
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW"
# result_dir <- "HCP_results/5k_results/individual/PW"
result_files <- list.files(result_dir, full.names = T)
result_files <- grep("501_HCP_BayesGLM2_45subj", result_files, value = T)
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

# >> Classical ----
# As a note, this would not be necessary in the case of the current (2021-02-10)
# state of the BayesfMRI package, but some retesting using the classical GLM
# is necessary in order to obtain the standard errors of the estimates.
results_dir <- "HCP_results/5k_results/group"
result_files <- list.files(results_dir,full.names = T)
FDR_result <- readRDS(grep("classical_activations_PW_FDR", result_files, value= T))
FWER_result <- readRDS(grep("classical_activations_PW_FWER", result_files, value = T))
library(BayesfMRI)
library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path',"/Applications/workbench")
light_orange <- grDevices::colorRampPalette(c("orange","white"))(3)[2]
col_pal <- c("yellow",light_orange,"red","purple")

left_active <- FDR_result$left$`0%`
left_active <- left_active + Reduce(`+`, FWER_result$left)
left_active[left_active == 0] <- NA
right_active <- FDR_result$right$`0%`
right_active <- right_active + Reduce(`+`, FWER_result$right)
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

