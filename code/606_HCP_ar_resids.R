# This is a script to examine the results from the classical GLM
result_dir_classical <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/classical"
result_dir_Bayesian <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
result_files_classical <- list.files(result_dir_classical, full.names = T)
result_files_classical <- grep("500_", result_files_classical, value = T)
result_files_Bayesian <- list.files(result_dir_Bayesian, full.names = T)
result_files_Bayesian <- grep("500_", result_files_Bayesian, value = T)
load("HCP_data/subjects.Rdata")
library(Matrix)
for(subj in subjects){
  for(visit in 1:2) {
    subj_cifti <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
    for(hem in c("left","right")) {
      # Need to read in the Bayesian results first because they contain the
      # prewhitened data, which is much faster than prewhitening all over again
      bayes_file <- grep(paste0("500_",subj,"_visit",visit,"_",hem,"_5k"),result_files_Bayesian, value = T)
      bayes_obj <- readRDS(bayes_file)
      pw_y <- bayes_obj$GLMs_Bayesian[[paste0("cortex",toupper(substring(hem,1,1)))]]$y
      pw_X <- bayes_obj$GLMs_Bayesian[[paste0("cortex",toupper(substring(hem,1,1)))]]$X
      rm(bayes_obj)
      classical_file <- grep(paste0("500_",subj,"_visit",visit,"_",hem,"_5k"),result_files_classical, value = T)
      result_obj <- readRDS(classical_file)
      K <- length(result_obj$beta_names)
      for(j in seq(length(result_obj$session_names))) {
        sess_X <- pw_X[[j]]
        NT <- nrow(sess_X)
        sess_y <- pw_y[(seq(NT) + (j-1)*NT)]
        sess_beta <- result_obj$betas_classical[[result_obj$session_names[j]]]$data[[paste0("cortex_",hem)]]
        N <- nrow(sess_beta)
        ntime <- NT / N
        y_hat <- as.numeric(sess_X %*% c(sess_beta))
        resids <- matrix(sess_y - y_hat,ntime, N)
        sig_ar <- apply(resids, 2, function(x) {
          ar_fit <- ar.yw(x, aic = T)
          clim <- qnorm((1+0.95)/2)/sqrt(ar_fit$n.used)
          ar_sig <- which(abs(ar_fit$ar) > clim)
          return(ar_fit$ar[ar_sig])
        })
        resid_ar <- t(apply(resids, 2, function(x) ar.yw(x, order.max = 6, aic = F)$ar))
        if(j == 1) {
          subj_cifti$data[[paste0("cortex_",hem)]] <- resid_ar
        } else {
          subj_cifti$data[[paste0("cortex_",hem)]] <- cbind(subj_cifti$data[[paste0("cortex_",hem)]],resid_ar)
        }
      }
    }
    saveRDS(subj_cifti, paste0("HCP_results/5k_results/ar_resids/606_subject_", subj, "_visit",visit,"_ar_resids.rds"))
  }
}

# View the ar coefficients for the residuals ----
library(ciftiTools)
ciftiTools.setOption('wb_path', "/Applications/workbench")
# >> Individual subject ----
ar_obj <- readRDS("HCP_results/5k_results/ar_resids/606_subject_103818_visit1_ar_resids.rds")
plot(ar_obj,
     idx = 7,
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")
# >> Averages ----
ar_dir <- "HCP_results/5k_results/ar_resids"
ar_files <- list.files(ar_dir, full.names = T)
ar_obj <- readRDS(ar_files[1])
for(fn in 2:length(ar_files)) {
  ar_obj <- ar_obj + readRDS(ar_files[fn])
}
ar_sum_1 <- select_xifti(ar_obj, 1:6)
ar_sum_2 <- select_xifti(ar_obj, 7:12)
ar_avg <- (ar_sum_1 + ar_sum_2) / (length(ar_files)*2)
plot(ar_avg,
     idx = 1:6,
     widget = T,
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")
# These look pretty good when I look at the overall averages. The residual AR
# coefficients have patterns of higher values along the ridges of the cifti
# surfaces, and the AR6 coefficient has some consistent higher values in the
# occipital regions, but the values are almost all between -0.01 and 0.01
