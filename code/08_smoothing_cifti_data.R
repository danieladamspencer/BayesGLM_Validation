# This is a more correct method for smoothing the data
# First, a test about the smoothing methods used
# library(ciftiTools)
# ciftiTools.setOption('wb_path', "/Applications/workbench")
# surfL_fname <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/visit1_data/103818/MNINonLinear/fsaverage_LR32k/103818.L.midthickness.32k_fs_LR.surf.gii"
# surfR_fname <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/visit1_data/103818/MNINonLinear/fsaverage_LR32k/103818.R.midthickness.32k_fs_LR.surf.gii"
# unsmoothed_cifti <- read_cifti(cifti_fname = "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/visit1_data/103818/MNINonLinear/Results/tfMRI_MOTOR_LR/tfMRI_MOTOR_LR_Atlas.dtseries.nii", surfL_fname = surfL_fname, surfR_fname = surfR_fname)
# for(fwhm in 3:6) {
#   smoothed_cifti <- smooth_cifti(unsmoothed_cifti, cifti_target_fname = "~/Desktop/cifti_smoothing_test.dtseries.nii", surf_FWHM = fwhm, surfL_fname = surfL_fname, surfR_fname = surfR_fname)
#   unsmoothed_cifti <- smoothed_cifti
# }
# plot(smoothed_cifti, fname = "~/Desktop/smoothed_3_4_5_6mmFHWM_subcortical_103818_visit1_idx1.png", surfL= surfL_fname, surfR= surfR_fname)
#
# old_smoothed_cifti <- read_cifti(cifti_fname = "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/visit1_data/103818/MNINonLinear/Results/tfMRI_MOTOR_LR/tfMRI_MOTOR_LR_Atlas_FWHM6.dtseries.nii", surfL_fname = surfL_fname, surfR_fname = surfR_fname)
#
# plot(old_smoothed_cifti, fname = "~/Desktop/old_smoothed_6mmFWHM_103818_visit1_idx1.png")
#
# unsmoothed_cifti <- read_cifti(cifti_fname = "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/visit1_data/103818/MNINonLinear/Results/tfMRI_MOTOR_LR/tfMRI_MOTOR_LR_Atlas.dtseries.nii", surfL_fname = surfL_fname, surfR_fname = surfR_fname)
#
# correct_6mm_smooth <- smooth_cifti(unsmoothed_cifti, surf_FWHM = 6, surfL_fname = surfL_fname, surfR_fname = surfR_fname)
# plot(correct_6mm_smooth, fname = "~/Desktop/correct_6mm_smoothing_103818_visit1_LR_idx1.png", surfL= surfL_fname, surfR= surfR_fname)

# NOW TO CORRECTLY RUN THE SMOOTHING ----
load("/Volumes/Lab_Data_Drive/users/danspen/HCP_Motor_Task_Dan/subjects.Rdata")
library(parallel)
cl <- makeCluster(10)
parSapplyLB(cl, subjects, function(subject) {
  library(ciftiTools)
  ciftiTools.setOption('wb_path', "/Applications/workbench")
  visit1_dir <- "/Volumes/Lab_Data_Drive/users/danspen/HCP_Motor_Task_Dan/visit1_data"
  surfL <- file.path(visit1_dir, subject,"MNINonLinear","fsaverage_LR32K",
                     paste0(subject,".L.midthickness.32k_fs_LR.surf.gii"))
  surfR <- file.path(visit1_dir, subject,"MNINonLinear","fsaverage_LR32K",
                     paste0(subject,".R.midthickness.32k_fs_LR.surf.gii"))
  unsmoothed_fname_LR <- file.path(visit1_dir, subject,"MNINonLinear","Results",
                                   "tfMRI_MOTOR_LR","tfMRI_MOTOR_LR_Atlas.dtseries.nii")
  unsmoothed_fname_RL <- file.path(visit1_dir, subject,"MNINonLinear","Results",
                                   "tfMRI_MOTOR_RL","tfMRI_MOTOR_RL_Atlas.dtseries.nii")
  for(fwhm in 2:10) {
    smoothed_cifti_LR <- smooth_cifti(x = unsmoothed_fname_LR,
                                      cifti_target_fname =
                                        file.path(visit1_dir, subject,
                                                  "MNINonLinear","Results",
                                                  "tfMRI_MOTOR_LR",
                                                  paste0("tfMRI_MOTOR_LR_Atlas_FWHM",fwhm,".dtseries.nii")),
                                      surf_FWHM = fwhm, vol_FWHM = fwhm,
                                      surfL_fname = surfL, surfR_fname = surfR)
    smoothed_cifti_RL <- smooth_cifti(x = unsmoothed_fname_RL,
                                      cifti_target_fname =
                                        file.path(visit1_dir, subject,
                                                  "MNINonLinear","Results",
                                                  "tfMRI_MOTOR_RL",
                                                  paste0("tfMRI_MOTOR_RL_Atlas_FWHM",fwhm,".dtseries.nii")),
                                      surf_FWHM = fwhm, vol_FWHM = fwhm,
                                      surfL_fname = surfL, surfR_fname = surfR)
  }
  return(NULL)
})
