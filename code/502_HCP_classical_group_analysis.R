# This is a script to perform the t-tests necessary to get the p-values in
# order to determine areas of activation for the classical estimates
classical_estimates <-
  readRDS("HCP_results/5k_results/602_avg_estimates_PW_classical.rds")

# Calculate the average estimates ----
library(abind)
average_estimates <- sapply(classical_estimates, function(hem_est) {
  combined_estimates <- abind(hem_est[[1]], hem_est[[2]], rev.along = 0)
  avg_est <- apply(combined_estimates,1:2,mean)
  return(avg_est)
}, simplify = F)
saveRDS(average_estimates, "HCP_results/5k_results/502_HCP_classical_group_PW_estimates.rds")

# Average estimates by visit
for(v in 1:2) {
  visit_estimates <- sapply(classical_estimates, function(ce) {
    visit_est <- ce[[paste0("visit",v)]]
    avg_visit_est <- apply(visit_est,1:2,mean)
    return(avg_visit_est)
  })
  saveRDS(visit_estimates,
          paste0("HCP_results/5k_results/502_HCP_classical_group_PW_estimates_visit",v,".rds"))
}

# Average estimates by subsample
sample_subjects <- readRDS("/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/group/PW/subsamples/501_sample_subjects_20210319.rds")
classical_estimates <-
  readRDS("HCP_results/5k_results/602_avg_estimates_PW_classical.rds")
load("HCP_data/subjects.Rdata")
subsamp_ests <- sapply(sample_subjects, function(num_sub) {
  num_sampl <- apply(num_sub, 2, function(sample_sub) {
    use_idx <- sapply(sample_sub, function(subj) which(subjects == subj))
    out <- sapply(classical_estimates, function(hem_results) {
      apply(hem_results$visit1[,,use_idx],1:2, mean)
    })
    return(out)
  })
  names(num_sampl) <- paste0("sample",1:10)
  return(num_sampl)
}, simplify = F)
names(subsamp_ests) <- paste0(c(10,20,30),"subjects")
saveRDS(subsamp_ests, "HCP_results/5k_results/group/502_HCP_classical_subsample_estimates.rds")

# FWER ----
# bonferroni_cutoff <- 0.01 / (4443+4444) # alpha = 0.01
threshs <- c(0,0.5,1)
combined_active_FWER <- sapply(classical_estimates, function(hem_res) {
  num_locs <- dim(hem_res[[1]])[1]
  bonferroni_cutoff <- 0.01 / num_locs # alpha = 0.01
  data_df <- reshape2::melt(hem_res)
  vertex_lists <- split(data_df,data_df$Var1)
  thresh_active <- sapply(threshs, function(thr) {
    active_all <- sapply(vertex_lists, function(v_df){
      vt_df <- split(v_df, v_df$Var2)
      active_vertex <- sapply(vt_df, function(vt) {
        ttest_res <- t.test(x = vt$value, mu = thr, alternative = "greater")
        active_v <- ifelse(ttest_res$p.value < bonferroni_cutoff, 1, 0)
        return(active_v)
      }, simplify = "array")
      return(active_vertex)
    }, simplify = 'array')
    return(t(active_all))
  }, simplify = F)
  names(thresh_active) <- paste0(threshs,"%")
  return(thresh_active)
}, simplify = FALSE)
saveRDS(combined_active_FWER, "HCP_results/5k_results/group/502_HCP_classical_activations_PW_FWER.rds")

# FDR ----
# BH_FDR <- function(p, FDR = 0.05) {
#   p_rank <- rank(-p, na.last=T)
#   BH_cutoff <- FDR * p_rank / length(p_rank)
#   out <- ifelse(is.na(p), 1, p)
#   out <- ifelse(out < BH_cutoff, 1,0)
#   return(out)
# }

threshs <- c(0)
combined_active_FDR <- sapply(classical_estimates, function(hem_res) {
  data_df <- reshape2::melt(hem_res)
  vertex_lists <- split(data_df,data_df$Var1)
  thresh_active <- sapply(threshs, function(thr) {
    pval_all <- sapply(vertex_lists, function(v_df){
      vt_df <- split(v_df, v_df$Var2)
      pval_vertex <- sapply(vt_df, function(vt) {
        ttest_res <- t.test(x = vt$value, mu = thr, alternative = "greater")
        return(ttest_res$p.value)
      }, simplify = "array")
      return(pval_vertex)
    }, simplify = 'array')
    pval_adj <- apply(pval_all,1,p.adjust, method = "BH")
    active_all <- apply(pval_adj,2,function(x) as.numeric(x < 0.01))
    # active_all2 <- apply(pval_all,1,BH_FDR, FDR = 0.01) # This might be wrong
    return(active_all)
  }, simplify = F)
  names(thresh_active) <- paste0(threshs,"%")
  return(thresh_active)
}, simplify = FALSE)
saveRDS(combined_active_FDR, "HCP_results/5k_results/group/502_HCP_classical_activations_PW_FDR.rds")
