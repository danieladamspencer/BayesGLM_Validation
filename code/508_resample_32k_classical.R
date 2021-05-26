# Resample the full-res classical results ----
load("HCP_data/subjects.Rdata")
library(parallel)
cl <- makeCluster(4)
parSapply(cl, X = subjects, function(subject) {
  classical_dir <- "HCP_results/32k_results/smoothed"
  surface_dir <- "HCP_data/subject_surfaces"
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
  library(ciftiTools)
  ciftiTools.setOption('wb_path','/Applications/workbench')
  subject_surfL <- file.path(surface_dir,paste0(subject,".L.midthickness.32k_fs_LR.surf.gii"))
  subject_surfR <- file.path(surface_dir,paste0(subject,".R.midthickness.32k_fs_LR.surf.gii"))
  for(visit in paste0("visit",1:2)) {
    subject_left <- readRDS(grep(paste0("500_",subject,"_",visit,"_left"), list.files(classical_dir,full.names = T), value = T))$betas_classical$avg
    subject_right <- readRDS(grep(paste0("500_",subject,"_",visit,"_right"), list.files(classical_dir,full.names = T), value = T))$betas_classical$avg
    subject_cifti <- xifti_combine(subject_left,subject_right)
    resample_cifti <- resample_xifti(subject_cifti, resamp_res = 5000)
    saveRDS(resample_cifti, file.path(classical_dir,paste0("500_",subject,"_",visit,"_both_resampled_classical.rds")))
  }
})
stopCluster(cl)
