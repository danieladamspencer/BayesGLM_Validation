# This is a script to plot the activation maps for the classical analysis of the
# HCP data.
library(ciftiTools)
ciftiTools.setOption('wb_path',"/Applications/workbench")

# Plot each task's FDR & FWER activation ----
FDR_activations <- readRDS("HCP_results/5k_results/group/502_HCP_classical_activations_FDR.rds")
FWER_activations <- readRDS("HCP_results/5k_results/group/502_HCP_classical_activations_FWER.rds")
thresh <- 0
hem <- c('left','right')
task_names <- c('visual cue','foot','hand','tongue')
activation_maps <- vector('list',length(task_names))
names(activation_maps) <- task_names

# Changing this to plot both hemispheres of the tongue task ONLY
for(task in c(4)) {
  # Load a template for the task
  activation_maps[[task]] <- readRDS("HCP_data/603_cifti_5k_template_whole.rds")
  for(h in hem) {
    # Set the data to match the dimensions needed for each task
    activation_maps[[task]]$data[[paste0('cortex_',h)]] <-
      matrix(NA,
             nrow = nrow(activation_maps[[task]]$data[[paste0('cortex_',h)]]),
             ncol = 2)
    activation_maps[[task]]$data[[paste0('cortex_',h)]][,1] <-
      ifelse(FDR_activations[[h]]$`0%`[,task] == 0, NA, FDR_activations[[h]]$`0%`[,task])
    activation_maps[[task]]$data[[paste0('cortex_',h)]][,2] <-
      ifelse(FWER_activations[[h]]$`0%`[,task] == 0, NA, 2*FWER_activations[[h]]$`0%`[,task])
    activation_maps[[task]]$data[[paste0('cortex_',h)]] <-
      as.matrix(
        apply(
          activation_maps[[task]]$data[[paste0('cortex_', h)]],
          1,
          function(x) {
            max_x <- max(x,na.rm = TRUE)
            max_x <- ifelse(max_x == -Inf, NA, max_x)
            return(max_x)
          }))
  }
}

for(task in task_names[4]) {
  plot(activation_maps[[task]], color_mode = 'qualitative', hemisphere = 'both', colors = c('darkblue','forestgreen'),
       surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
       surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii",
       save = TRUE, fname = paste0("plots/604_",sub(" ","_", task),"_classical_activation_map"))
}
