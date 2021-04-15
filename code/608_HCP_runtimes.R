# Find the runtimes for the BayesGLM2 functions
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples"
result_files <- list.files(result_dir, full.names = T)
runtimes <- sapply(paste0("_",c(10,20,30,45),"subj_"), function(ss) {
  size_files <- grep(ss,result_files, value = T)
  out <- sapply(size_files, function(sf) readRDS(sf)$total_time)
}, simplify = F)
saveRDS(runtimes, "HCP_results/5k_results/group/608_BayesGLM2_runtimes.rds")
runtimes <- readRDS("HCP_results/5k_results/group/608_BayesGLM2_runtimes.rds")
# Quality check
hist(runtimes[[1]]) # Weird how the subsamples from earlier runs take longer. Since the figures are only run using the first 60 entries, I'll just use those
hist(runtimes[[2]]) # Same here
hist(runtimes[[3]])
hist(runtimes[[4]])

c(median(runtimes[[1]][seq(60)]), sd(runtimes[[1]][seq(60)])) /60
c(median(runtimes[[2]][seq(60)]), sd(runtimes[[2]][seq(60)])) /60
c(median(runtimes[[3]][seq(60)]), sd(runtimes[[3]][seq(60)])) /60
c(median(runtimes[[4]]), sd(runtimes[[4]])) /60

library(tidyverse)
subj_files_45 <- list.files("HCP_results/5k_results/group", full.names = T) %>%
  grep(pattern = "BayesGLM2_45subj",x = .,value = T)
runtimes45 <- sapply(subj_files_45, function(x) {readRDS(x)$total_time})
c(median(runtimes45), sd(runtimes45)) /60
