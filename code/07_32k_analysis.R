# This is a script for running 32k analysis using the INLA implementation of
# the Bayesian GLM

# Simulated data ----
library(ciftiTools)
ciftiTools.setOption('wb_path',"/Applications/workbench")
library(INLA)
inla.setOption(pardiso.license = "~/licenses/pardiso.lic")
library(BayesfMRI)
library(brainSim)
# Simulate the data at 32k
set.seed(47401)
# Testing at a lower resolution to see how computation time scales
resolution <- 10000
sim_data <-
  simulate_cifti_multiple(
    wb_path = "/Applications/workbench",
    n_subjects = 1,
    n_sessions = 1,
    n_runs = 1,
    ntasks = 1,
    ntime = 300,
    resamp_res = 10000,
    max_amplitude = 2,
    onsets = NULL,
    durations = NULL,
    TR = 1,
    subject_var = NULL,
    session_var = NULL,
    run_var = NULL,
    ar_error = NULL,
    surfL = NULL,
    surfR = NULL
  )
# Write the data to a directory for analysis
cifti_fname <- paste0("~/Desktop/simulated_data/07_",resolution,"resolution_1subject_1session_1run_1task_T300.dtseries.nii")
write_cifti(
  xifti = sim_data$simulated_cifti$`Subject 1 Session 1 Run 1`,
  cifti_fname = cifti_fname,
  surfL_fname = NULL,
  surfR_fname = NULL,
  verbose = FALSE
)
# Now try running the analyses
inla_result <-
  BayesGLM_cifti(
    cifti_fname = cifti_fname,
    surfL_fname = sim_data$simulated_cifti$`Subject 1 Session 1 Run 1`$surf$cortex_left,
    surfR_fname = sim_data$simulated_cifti$`Subject 1 Session 1 Run 1`$surf$cortex_right,
    brainstructures = 'left',
    design = sim_data$design,
    onsets = NULL,
    TR = 1,
    nuisance = NULL,
    nuisance_include = NULL,
    scale_BOLD = TRUE,
    scale_design = TRUE,
    GLM_method = "Bayesian",
    ar_order = 0,
    ar_smooth = 0,
    session_names = NULL,
    resamp_res = NULL,
    num.threads = 4,
    verbose = TRUE,
    outfile = NULL,
    return_INLA_result = TRUE,
    avg_sessions = FALSE,
    trim_INLA = TRUE,
    tol = 1e-3
  )
# Save the analysis
save_fname <- paste0("~/Desktop/simulated_data/07_",resolution,"resolution_1subject_1session_1run_1task_T300_inla.rds")
saveRDS(inla_result, file = save_fname)
