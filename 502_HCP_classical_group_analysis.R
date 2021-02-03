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

# FWER ----
bonferroni_cutoff <- 0.01 / (4443+4444) # alpha = 0.01
threshs <- c(0)
combined_active_FWER <- sapply(classical_estimates, function(hem_res) {
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
saveRDS(combined_active_FWER, "HCP_results/5k_results/502_HCP_classical_activations_PW_FWER.rds")

# FDR ----
BH_FDR <- function(p, FDR = 0.05) {
  p_rank <- rank(-p, na.last=T)
  BH_cutoff <- FDR * p_rank / length(p_rank)
  out <- ifelse(is.na(p), 1, p)
  out <- ifelse(out < BH_cutoff, 1,0)
  return(out)
}

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
    active_all <- apply(pval_all,1,BH_FDR, FDR = 0.01)
    return(active_all)
  }, simplify = F)
  names(thresh_active) <- paste0(threshs,"%")
  return(thresh_active)
}, simplify = FALSE)
saveRDS(combined_active_FDR, "HCP_results/5k_results/502_HCP_classical_activations_PW_FDR.rds")
