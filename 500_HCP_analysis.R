# This is a script for running the HCP data analysis
library(ciftiTools)
# wb_cmd <- "/Applications/workbench" # Mac Pro
# ciftiTools::ciftiTools.setOption("wb_path", wb_cmd) # Mac Pro
ciftiTools::ciftiTools.setOption("wb_path", "/Applications/workbench/bin_macosx64/wb_command") # Dan's Macbook Pro
wb_cmd <- "/Applications/workbench/bin_macosx64/wb_command" # Dan's Macbook Pro
library(INLA)
inla.setOption(pardiso.license = "~/licenses/pardiso.lic") # Dan's Macbook Pro
library(BayesfMRI)
main_dir <- "~/github/BayesGLM_Validation" # Dan's Macbook Pro
# main_dir <- "/Volumes/Macintosh HD/Users/Shared/HCP" # Mac Pro
data_dir <- "~/github/BayesGLM_Validation/HCP_data" # Dan's Macbook Pro
# data_dir <- "/Volumes/Macintosh HD/Users/Shared/HCP/visit1_data/" # Mac Pro
result_dir <- "~/github/BayesGLM_Validation/HCP_results/1k_results/PW" # Dan's Macbook Pro
# result_dir <- "/Volumes/Macintosh HD/Users/Shared/HCP/5k_results/PW" # Mac Pro
load(file.path(data_dir,"subjects.Rdata")) # Macbook Pro
# load(file.path(main_dir,"subjects.Rdata")) # Mac Pro
tasks <- c('cue','lf','lh','rf','rh','t') # Task data frame columns
names_tasks <- c('cue','left foot','left hand','right foot','right hand','tongue')
colors_tasks <- c('black',RColorBrewer::brewer.pal(5, 'Set2'))
cols_LH <- c(1,4:6) #cue, right foot, right hand, tongue
cols_RH <- c(1:3,6) #cue, left foot, left hand, tongue
cols_list <- list(cols_LH, cols_RH)
TR = 0.72 #temporal resolution of data
thetas <- NULL # No starting values for precision parameters
subjects <- subjects[1:2]
for(subject in subjects) {
  dir_s <- file.path(data_dir, subject, 'MNINonLinear', 'fsaverage_LR32k')
  fname_gifti_left <- file.path(dir_s, paste0(subject,'.L.midthickness.32k_fs_LR.surf.gii'))
  fname_gifti_right <- file.path(dir_s, paste0(subject,'.R.midthickness.32k_fs_LR.surf.gii'))
  for(visit in 1:2) {
    if(visit == 1) {
      dir1_s <- file.path(data_dir,subject, 'MNINonLinear', 'Results', 'tfMRI_MOTOR_LR')
      dir2_s <- file.path(data_dir,subject, 'MNINonLinear', 'Results', 'tfMRI_MOTOR_RL')
      fname1_ts <- file.path(dir1_s,'tfMRI_MOTOR_LR_Atlas.dtseries.nii')
      fname2_ts <- file.path(dir2_s,'tfMRI_MOTOR_RL_Atlas.dtseries.nii')
    }
    if(visit == 2) {
      # dir1_s <- file.path(data_dir,"retest",subject, 'MNINonLinear', 'Results', 'tfMRI_MOTOR_LR')
      dir1_s <- file.path(main_dir,"visit2_data",subject, 'MNINonLinear', 'Results', 'tfMRI_MOTOR_LR')
      # dir2_s <- file.path(data_dir,"retest",subject, 'MNINonLinear', 'Results', 'tfMRI_MOTOR_RL')
      dir2_s <- file.path(main_dir,"visit2_data",subject, 'MNINonLinear', 'Results', 'tfMRI_MOTOR_RL')
      fname1_ts <- file.path(dir1_s,'tfMRI_MOTOR_LR_Atlas.dtseries.nii')
      fname2_ts <- file.path(dir2_s,'tfMRI_MOTOR_RL_Atlas.dtseries.nii')
    }
    #analyze hemispheres separately due to different in set of tasks
    for(h in 1:2){

      #h=1 -- left hemisphere
      #h=2 -- right hemisphere
      hem <- c('left','right')[h]

      cols_h <- cols_list[[h]]
      tasks_h <- tasks[cols_h]
      names_tasks_h <- names_tasks[cols_h]
      K <- length(tasks_h)

      #Set up list of onsets
      onsets1 <- onsets2 <- vector('list', length=K)
      names(onsets1) <- names(onsets2) <- names_tasks_h
      for(t in tasks_h){
        ind_t <- which(tasks_h==t)
        fname1_t <- file.path(dir1_s,'EVs',paste0(t,'.txt'))
        onsets1[[ind_t]] <- read.table(fname1_t, header=FALSE)
        fname2_t <- file.path(dir2_s,'EVs',paste0(t,'.txt'))
        onsets2[[ind_t]] <- read.table(fname2_t, header=FALSE)
      }

      #Set up nuisance regressors
      motion1 <- as.matrix(read.table(file.path(dir1_s,'Movement_Regressors.txt'), header=FALSE))
      motion2 <- as.matrix(read.table(file.path(dir2_s,'Movement_Regressors.txt'), header=FALSE))

    start_time <- proc.time()[3]
    result_svh <- BayesGLM_cifti(cifti_fname = c(fname1_ts, fname2_ts), # Multi-session
                                 surfL_fname = fname_gifti_left,
                                 surfR_fname = fname_gifti_right,
                                 brainstructures = hem,
                                 # wb_path = wb_cmd,
			                           design = NULL,
                                 onsets = list(onsets1, onsets2), # Multi-session
                                 TR = TR,
                                 nuisance = list(motion1, motion2), # Multi-session
                                 nuisance_include = c('drift','dHRF'),
			                           scale_BOLD = TRUE,
			                           scale_design = TRUE,
                                 GLM_method = 'both',
			                           prewhiten = TRUE,
			                           ar_order = 6,
			                           surface_sigma = 6 / (2*sqrt(2*log(2))),
                                 session_names = c('LR','RL'), # Multiple sessions
                                 resamp_res = 1000, # Don't forget to change this
                                 num.threads = 4, # Remember the tradeoff here (speed/memory) 4 to 6 threads seems optimal based on testing
                                 verbose = TRUE,
                                 outfile = NULL,
                                 return_INLA_result = T,
                                 avg_betas_over_sessions = T,
			                           trim_INLA = T)
    total_time <- proc.time()[3] - start_time
    result_svh$total_time <- total_time
    saveRDS(result_svh, file=file.path(result_dir,paste0("500_",subject,"_visit",visit,"_",hem,"_1k_",format(Sys.Date(),"%Y%m%d"),".rds")))
    }
  }
}
