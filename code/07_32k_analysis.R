# # This is a script for running 32k analysis using the INLA implementation of
# # the Bayesian GLM
#
# # Simulated data ----
# library(ciftiTools)
# ciftiTools.setOption('wb_path',"/Applications/workbench")
# library(INLA)
# inla.setOption(pardiso.license = "~/licenses/pardiso.lic")
# library(BayesfMRI)
# library(brainSim)
# # Simulate the data at 32k
# set.seed(47401)
# # Testing at a lower resolution to see how computation time scales
# resolution <- 32000
# sim_data <-
#   simulate_cifti_multiple(
#     wb_path = "/Applications/workbench",
#     n_subjects = 1,
#     n_sessions = 1,
#     n_runs = 1,
#     ntasks = 1,
#     ntime = 300,
#     resamp_res = NULL,
#     max_amplitude = 2,
#     onsets = NULL,
#     durations = NULL,
#     TR = 1,
#     subject_var = NULL,
#     session_var = NULL,
#     run_var = NULL,
#     ar_error = NULL,
#     surfL = NULL,
#     surfR = NULL
#   )
# # Write the data to a directory for analysis
# cifti_fname <- paste0("~/Desktop/simulated_data/07_",resolution,"resolution_1subject_1session_1run_1task_T300.dtseries.nii")
# write_cifti(
#   xifti = sim_data$simulated_cifti$`Subject 1 Session 1 Run 1`,
#   cifti_fname = cifti_fname,
#   surfL_fname = NULL,
#   surfR_fname = NULL,
#   verbose = FALSE
# )
# # Now try running the analyses
# inla_result <-
#   BayesGLM_cifti(
#     cifti_fname = cifti_fname,
#     surfL_fname = sim_data$simulated_cifti$`Subject 1 Session 1 Run 1`$surf$cortex_left,
#     surfR_fname = sim_data$simulated_cifti$`Subject 1 Session 1 Run 1`$surf$cortex_right,
#     brainstructures = 'left',
#     design = sim_data$design,
#     onsets = NULL,
#     TR = 1,
#     nuisance = NULL,
#     nuisance_include = NULL,
#     scale_BOLD = TRUE,
#     scale_design = TRUE,
#     GLM_method = "Bayesian",
#     ar_order = 0,
#     ar_smooth = 0,
#     session_names = NULL,
#     resamp_res = NULL,
#     num.threads = 4,
#     verbose = TRUE,
#     outfile = NULL,
#     return_INLA_result = TRUE,
#     avg_sessions = FALSE,
#     trim_INLA = TRUE,
#     tol = 1e-3
#   )
# # Save the analysis
# save_fname <- paste0("~/Desktop/simulated_data/07_",resolution,"resolution_1subject_1session_1run_1task_T300_inla.rds")
# saveRDS(inla_result, file = save_fname)


# Gambling data ----
library(ciftiTools)
ciftiTools.setOption('wb_path',"/Applications/workbench")
library(INLA)
inla.setOption(pardiso.license = "~/licenses/pardiso.lic")
library(BayesfMRI)
data_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Gambling_Task_Dan/visit1_data"
surfaces_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/visit1_data"
load("/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan/subjects.Rdata")
TR <- 0.72
do_smooth <- TRUE
subject <- subjects[1]; visit <- 1; hem <- 'left'; run <- "LR"
# subjects <- subjects[c(1,2,4)]
for(subject in subjects) {
  subject_surface_files <-
    list.files(file.path(surfaces_dir, subject, "MNINonLinear/fsaverage_LR32k"),
               full.names = T)
  subject_surface_left <- grep("\\.L\\.", subject_surface_files, value = T)
  subject_surface_right <- grep("\\.R\\.", subject_surface_files, value = T)
  for(visit in 1) {
    data_dir_s <- file.path(data_dir, subject)
    if(visit == 2) data_dir_s <- file.path(data_dir,"retest",subject)
    # Run hemispheres separately
    for(hem in c("left")) {
      cifti_fnames <- character()
      onsets <- list()
      mvmnt <- list()
      for(run in c("LR","RL")) {
        run_dir <- file.path(
          data_dir_s,
          paste0("MNINonLinear/Results/tfMRI_GAMBLING_", run))
        # Set up cifti
        run_cifti_fname <- list.files(run_dir, full.names = T)
        run_cifti_fname <- grep(".dtseries.nii", run_cifti_fname, value = T)
        # Set up win - loss EV
        ev_files <- list.files(file.path(run_dir,"EVs"), full.names = TRUE)
        ev_files <- grep("_event.txt",ev_files, value = T)
        ev_names <- sub(paste0(file.path(run_dir,"EVs"),"/"),"",ev_files)
        ev_names <- sub(".txt","",ev_names)
        run_evs <- sapply(ev_files, read.table, header = F, simplify = F)
        run_evs <- sapply(run_evs, function(x) x[,-3], simplify = F)
        names(run_evs) <- ev_names
        run_onsets <- Reduce(rbind,run_evs)
        run_onsets <- list(as.matrix(run_onsets[order(run_onsets$V1),-3]))
        # Set up movement
        run_movement <-
          read.table(file.path(run_dir, "Movement_Regressors.txt"), header = F)
        # Smooth data, if desired
        if(do_smooth) {
          smoothed_fname <- file.path(run_dir, paste0("tfMRI_GAMBLING_",run,"_Atlas_smoothed.dtseries.nii"))
          smooth_cifti(
            x = run_cifti_fname,
            cifti_target_fname = smoothed_fname,
            surf_FWHM = 5,
            vol_FWHM = 5,
            surfL_fname = subject_surface_left,
            surfR_fname = subject_surface_right,
            cerebellum_fname = NULL,
            subcortical_zeroes_as_NA = FALSE,
            cortical_zeroes_as_NA = FALSE,
            subcortical_merged = FALSE
          )
          run_cifti_fname <- smoothed_fname
        }
        # Add run data to lists
        cifti_fnames <- c(cifti_fnames, run_cifti_fname)
        onsets[[run]] <- run_onsets
        mvmnt[[run]] <- as.matrix(run_movement)
      }
      # Set up analysis with BayesGLM_cifti
      result_svh <- BayesGLM_cifti(cifti_fname = cifti_fnames,
                                   surfL_fname = subject_surface_left,
                                   surfR_fname = subject_surface_right,
                                   brainstructures = hem,
                                   design = NULL,
                                   onsets = onsets,
                                   TR = TR,
                                   nuisance = mvmnt,
                                   nuisance_include = c('drift','dHRF'),
                                   scale_BOLD = TRUE,
                                   scale_design = TRUE,
                                   GLM_method = 'classical',
                                   ar_order = 6,
                                   ar_smooth = 6,
                                   session_names = c('LR','RL'), # Multiple sessions
                                   resamp_res = 5000, # Don't forget to change this
                                   num.threads = 6, # Remember the tradeoff here (speed/memory) 4 to 6 threads seems optimal based on testing
                                   verbose = TRUE,
                                   outfile = NULL,
                                   return_INLA_result = T,
                                   avg_sessions = T,
                                   trim_INLA = T)
      saveRDS(result_svh, file = file.path("/Volumes/Lab_Data_Drive/users/danspen/HCP_Gambling_Task_Dan", paste0("07_32k_smoothed_classical_gambling_subject",subject,"_visit",visit,"_",hem,"_",format(Sys.Date(),"%Y%m%d"),".rds")))
    }
  }
}
