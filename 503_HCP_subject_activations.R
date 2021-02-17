# This is a script to run the id_activations on the subject_level data


# BAYESIAN ----
# Load libraries and result directory ----
library(BayesfMRI)
library(excursions)
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
save_dir <- "HCP_results/5k_results/individual/PW/activations"
load("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/subjects.Rdata")
# subjects <- subjects[seq(5)]
subject_files <- unlist(sapply(subjects, grep, x = list.files(result_dir, full.names = T), value = T))
subject_files <- grep("visit2", subject_files, value=T)

# Run the id_activations ----
threshs <- c(0.5)
alpha <- 0.01
hem <- c("left","right")

for(subject in subjects){
  for(h in hem) {
    L_or_R <- toupper(substring(h,1,1))
    result_obj <- readRDS(grep(h, grep(subject, subject_files, value = T), value = T))
    V <- nrow(result_obj$betas_Bayesian$avg$data[[paste0("cortex_",h)]])
    K <- ncol(result_obj$betas_Bayesian$avg$data[[paste0("cortex_",h)]])
    pred_idx <- result_obj$GLMs_Bayesian[[paste0("cortex",L_or_R)]]$INLA_result$misc$configs$config[[1]]$pred_idx
    mu_pred <- result_obj$GLMs_Bayesian[[paste0("cortex",L_or_R)]]$INLA_result$misc$configs$config[[1]]$mean[pred_idx]
    Q_pred <- result_obj$GLMs_Bayesian[[paste0("cortex",L_or_R)]]$INLA_result$misc$configs$config[[1]]$Q[pred_idx,pred_idx]
    for(thr in threshs) {
      for(a in alpha) {
        exc_obj <- excursions(alpha = a, u = thr, mu = mu_pred, Q = Q_pred, type = ">", method = "EB")
        act_obj <- list(
          active = matrix(exc_obj$E, V,K),
          excur_result = exc_obj
        )
        saveRDS(act_obj,
                file.path(save_dir,
                          paste0("503_HCP_visit2_subject_",
                                 subject,
                                 "_",h,
                                 "_thr",
                                 sub("\\.","",as.character(thr)),
                                 "_alpha",
                                 sub("\\.","",as.character(a)),
                                 "_activations.rds")))
      }
    }
  }
}


# CLASSICAL ----
# library(BayesfMRI)
# subject_files <- c(
#   "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_103818_visit1_left_5k_20210125.rds",
#   "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/500_103818_visit1_right_5k_20210125.rds"
# )
# classical_activations <- sapply(c("left","right"), function(hem) {
#   L_or_R <- toupper(substring(hem,1,1))
#   fn <- grep(hem,subject_files,value =T)
#   result <- readRDS(fn)
#   c_act <- id_activations.classical(result, alpha = 0.01)
#   tongue_act <- cbind(c_act[[paste0("cortex",L_or_R)]]$FDR_act[,4],c_act[[paste0("cortex",L_or_R)]]$FWER_act[,4])
#   tongue_act <- apply(tongue_act,1,sum)
#   tongue_act <- ifelse(tongue_act == 0, NA, tongue_act)
#   return(as.matrix(tongue_act))
# })
#
# saveRDS(classical_activations, "HCP_results/5k_results/503_classical_tongue_activations.rds")
