# This is a script to run the id_activations on the subject_level data

# Load libraries and result directory ----
library(BayesfMRI)
result_dir <- "HCP_results/5k_results/individual/PW"
glm_files <- grep("500_", list.files(result_dir,full.names = T), value = T)
subject <- "103818"
subject_files <- grep(subject,glm_files, value = T)

# Run the id_activations ----
threshs <- c(0,0.5,1)
alpha <- 0.01
hem <- c("left","right")
for(h in hem) {
  result_obj <- readRDS(grep(h, subject_files, value = T))
  L_or_R <- toupper(substring(h,1,1))
  for(thr in threshs) {
    for(a in alpha) {
      act_obj <-
        id_activations.posterior(result_obj$GLMs_Bayesian[[paste0("cortex", L_or_R)]],
                                 threshold = thr,
                                 alpha = a)
      saveRDS(act_obj,
              file.path(result_dir,
                        paste0("503_HCP_subject_",
                               subject,
                               "_thr",
                               sub("\\.","",as.character(thr)),
                               "_alpha",
                               sub("\\.","",as.character(a)),
                               "_activations.rds")))
    }
  }
}

