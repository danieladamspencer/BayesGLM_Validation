# This is a script for calculating the Dice coefficient for each subject's
# areas of activation from the average estimates for each visit.

# Grab contextual information ----
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
result_files <- list.files(result_dir, full.names = TRUE)
load("HCP_data/subjects.Rdata") # The subjects object (character vector)
result_files <- grep(".rds", result_files, value = T)
hems <- c('left','right')
threshold <- 0
subject_dice <- vector('numeric', length = length(subjects))
library(BayesfMRI)
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
library(excursions)

tasks <- c("visual cue","tongue","right foot","right hand","left foot","left hand")

# Calculate subject Dice coefficients ----
# BAYESIAN ----
# subjects_dice <- list()
# for(subj in subjects) {
#   visit_exc <- list()
#   for(v in 1:2) {
#     hems_exc <- list()
#     for(h in hems) {
#       L_or_R <- toupper(substring(h,1,1))
#       filen <- grep(paste0(subj,"_visit",v,"_",h), result_files, value = T)
#       result_subj <- readRDS(filen)
#       model_obj <- result_subj$GLMs_Bayesian[[paste0("cortex",L_or_R)]]
#       nvox <- model_obj$mesh$n
#       K <- length(model_obj$beta_names)
#       config <- excursions:::private.get.config(model_obj$INLA_result,1)
#       pred_inds <- model_obj$INLA_result$misc$configs$config[[1]]$pred_idx
#       # Remember that u is the threshold here
#       avg_excursions <- excursions(alpha = 0.01, u = threshold, mu =
#                            model_obj$INLA_result$misc$configs$config[[1]]$mean[pred_inds],
#                          Q = model_obj$INLA_result$misc$configs$config[[1]]$Q[pred_inds,pred_inds],
#                          type = ">")
#       mat_act <- matrix(avg_excursions$E,nvox,K)
#       mat_out <- matrix(0, nvox, (K+2))
#       if(h == 'left') mat_out[,1:4] <- mat_act[,c(1,4,2,3)]
#       if(h == 'right') mat_out[,c(1,2,5,6)] <- mat_act[,c(1,4,2,3)]
#       hems_exc[[h]] <- mat_out
#     }
#     visit_exc[[v]] <- Reduce(rbind,hems_exc)
#   }
#   subjects_dice[[subj]] <- mapply(function(v1,v2) {sum(v1 * v2) / mean(c(sum(v1),sum(v2))) },
#          v1 = split(visit_exc[[1]],col(visit_exc[[1]])),
#          v2 = split(visit_exc[[2]],col(visit_exc[[2]])))
# }
# all_dice <- Reduce(rbind, subjects_dice)
# saveRDS(all_dice, paste0("HCP_results/5k_results/604_Dice_coefficient_PW_threshold",threshold,".rds"))
# CLASSICAL ----
subjects_dice <- list()
for(subj in subjects) {
  visit_exc <- list()
  for(v in 1:2) {
    hems_exc <- list()
    for(h in hems) {
      L_or_R <- toupper(substring(h,1,1))
      filen <- grep(paste0(subj,"_visit",v,"_",h), result_files, value = T)
      result_subj <- readRDS(filen)
      model_obj <- result_subj$GLMs_Bayesian[[paste0("cortex",L_or_R)]]
      nvox <- model_obj$mesh$n
      K <- length(model_obj$beta_names)
      class_act <- id_activations.classical(result_subj,alpha = 0.01)
      # hems_exc[[h]] <- class_act[[paste0("cortex",L_or_R)]]$FWER_act
      mat_out <- matrix(0, nvox, (K+2))
      if(h == 'left') mat_out[,1:4] <- class_act[[paste0("cortex",L_or_R)]]$FDR_act[,c(1,4,2,3)]
      if(h == 'right') mat_out[,c(1,2,5,6)] <- class_act[[paste0("cortex",L_or_R)]]$FDR_act[,c(1,4,2,3)]
      hems_exc[[h]] <- mat_out
    }
    visit_exc[[v]] <- Reduce(rbind,hems_exc)
  }
  subjects_dice[[subj]] <- mapply(function(v1,v2) {sum(v1 * v2) / mean(c(sum(v1),sum(v2))) },
                                  v1 = split(visit_exc[[1]],col(visit_exc[[1]])),
                                  v2 = split(visit_exc[[2]],col(visit_exc[[2]])))
}
all_dice <- Reduce(rbind, subjects_dice)
saveRDS(all_dice, paste0("HCP_results/5k_results/604_classical_FDR_Dice_coefficient_PW_threshold",threshold,".rds"))
