# This is a script for calculating the Dice coefficient for each subject's
# areas of activation from the average estimates for each visit.

# Grab contextual information ----
result_dir <- "HCP_results/5k_results/individual"
result_files <- list.files(result_dir, full.names = TRUE)
load("HCP_data/subjects.Rdata") # The subjects object (character vector)
hems <- c('left','right')
visits <- 1:2
thresholds <- c(0,0.005)
subject_dice <- vector('numeric', length = length(subjects))
library(BayesfMRI)
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')

# Calculate subject Dice coefficients ----
for(subj in subjects) {
  for(v in visits) {
    for(h in hems) {
      L_or_R <- toupper(substring(h,1,1))
      filen <- grep(paste0(subj,"_visit",v,"_",h), result_files, value = T)
      config <- excursions:::private.get.config(model_obj$INLA_result,1)
      result_subj <- readRDS(filen)
      for(thr in thresholds) {
        excursions:::inla.output.indices()
        pred_inds <- model_obj$INLA_result$misc$configs$config[[1]]$pred_idx
        test <- excursions(alpha = 0.01, u = 0, mu =
                             model_obj$INLA_result$misc$configs$config[[1]]$mean,
                           Q = model_obj$INLA_result$misc$configs$config[[1]]$Q,
                           type = ">", ind = pred_inds)
        result_active <- id_activations(result_subj$GLMs_Bayesian[[paste0('cortex',L_or_R)]],
                                        threshold = thr, alpha = 0.01, type = NULL)
      }
    }
  }
}
