# This is a script to run the id_activations on the subject_level data

# # BAYESIAN ----
# # Load libraries and result directory ----
# library(BayesfMRI)
# library(excursions)
# result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW"
# save_dir <- "HCP_results/5k_results/individual/PW/activations"
# load("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/subjects.Rdata")
# subject_files <- unlist(sapply(subjects, grep, x = list.files(result_dir, full.names = T), value = T))
# # subject_files <- grep("visit2", subject_files, value=T)
#
# # Run the id_activations ----
# threshs <- c(0,0.5,1)
# alpha <- 0.01
# hem <- c("left","right")
#
# for(subject in subjects){
#   for(v in paste0("visit",1:2)) {
#     for(h in hem) {
#       L_or_R <- toupper(substring(h,1,1))
#       result_obj <- readRDS(grep(paste0(subject,"_",v,"_",h), subject_files, value = T))
#       for(thr in threshs) {
#         for(a in alpha) {
#           active_out <-
#             id_activations_cifti(result_obj,
#                                  alpha = a,
#                                  method = "Bayesian",
#                                  threshold = thr)
#           saveRDS(active_out,
#                   file.path(save_dir,
#                             paste0("503_HCP_",v,"_subject_",
#                                    subject,
#                                    "_",h,
#                                    "_thr",
#                                    sub("\\.","",as.character(thr)),
#                                    "_alpha",
#                                    sub("\\.","",as.character(a)),
#                                    "_activations.rds")))
#         }
#       }
#     }
#   }
# }


# CLASSICAL ----
# Load libraries and result directory ----
library(BayesfMRI)
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/classical"
save_dir <- "HCP_results/5k_results/individual/PW/activations"
load("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/subjects.Rdata")
# subjects <- subjects[4]
subject_files <- unlist(sapply(subjects, grep, x = list.files(result_dir, full.names = T), value = T))
# subject_files <- grep("visit1", subject_files, value=T)

# Run the id_activations ----
threshs <- c(0,0.5,1)
alpha <- 0.01
hem <- c("left","right")

for(subject in subjects){
  for(h in hem) {
    for(v in 1:2) {
      L_or_R <- toupper(substring(h,1,1))
      result_obj <-
        readRDS(grep(
          paste0("500_", subject, "_visit", v, "_", h),
          subject_files,
          value = T
        ))
      for(thr in threshs) {
        for(a in alpha) {
          act_obj <-
            id_activations.classical(
              result_obj$GLMs_classical[[paste0("cortex", L_or_R)]],
              alpha = alpha,
              threshold = thr,
              correction = "FWER"
            )
          saveRDS(act_obj,
                  file.path(save_dir,
                            paste0("503_HCP_visit",v,
                                   "_subject_",
                                   subject,
                                   "_",h,
                                   "_thr",
                                   sub("\\.","",as.character(thr)),
                                   "_alpha",
                                   sub("\\.","",as.character(a)),
                                   "_activations_classical.rds")))
        }
      }
    }
  }
}
