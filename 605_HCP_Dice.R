# # This is a script for calculating the Dice coefficient for each subject's
# # areas of activation from the average estimates for each visit.
#
# # Calculate subject Dice coefficients ----
# # BAYESIAN ----
# result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session"
# result_files <- list.files(result_dir, full.names = TRUE)
# load("HCP_data/subjects.Rdata") # The subjects object (character vector)
# result_files <- grep(".rds", result_files, value = T)
# hems <- c('left','right')
# threshold <- 1
# subject_dice <- vector('numeric', length = length(subjects))
# library(BayesfMRI)
# library(ciftiTools)
# ciftiTools.setOption('wb_path','/Applications/workbench')
# library(excursions)
#
# tasks <- c("visual cue","tongue","right foot","right hand","left foot","left hand")
# subjects_dice <- list()
# for(subj in subjects) {
#   visit_exc <- list()
#   for(v in 1:2) {
#     hems_exc <- list()
#     for(h in hems) {
#       L_or_R <- toupper(substring(h,1,1))
#       filen <- grep(paste0(subj,"_visit",v,"_",h), result_files, value = T)
#       result_subj <- readRDS(filen)
#       model_obj <- result_subj$GLMs_Bayesian[[paste0("cortex",L_or_R)]]
#       nvox <- model_obj$mesh$n
#       K <- length(model_obj$beta_names)
#       config <- excursions:::private.get.config(model_obj$INLA_result,1)
#       pred_inds <- model_obj$INLA_result$misc$configs$config[[1]]$pred_idx
#       # Remember that u is the threshold here
#       avg_excursions <- excursions(alpha = 0.01, u = threshold, mu =
#                            model_obj$INLA_result$misc$configs$config[[1]]$mean[pred_inds],
#                          Q = model_obj$INLA_result$misc$configs$config[[1]]$Q[pred_inds,pred_inds],
#                          type = ">")
#       mat_act <- matrix(avg_excursions$E,nvox,K)
#       mat_out <- matrix(0, nvox, (K+2))
#       if(h == 'left') mat_out[,1:4] <- mat_act[,c(1,4,2,3)]
#       if(h == 'right') mat_out[,c(1,2,5,6)] <- mat_act[,c(1,4,2,3)]
#       hems_exc[[h]] <- mat_out
#     }
#     visit_exc[[v]] <- Reduce(rbind,hems_exc)
#   }
#   # subjects_dice[[subj]] <- mapply(function(v1,v2) {sum(v1 * v2) / mean(c(sum(v1),sum(v2))) },
#   #        v1 = split(visit_exc[[1]],col(visit_exc[[1]])),
#   #        v2 = split(visit_exc[[2]],col(visit_exc[[2]])))
#   subjects_dice[[subj]] <- mapply(function(v1,v2) {sum(v1 * v2)},
#                                   v1 = split(visit_exc[[1]],col(visit_exc[[1]])),
#                                   v2 = split(visit_exc[[2]],col(visit_exc[[2]])))
# }
# all_dice <- Reduce(rbind, subjects_dice)
# # saveRDS(all_dice, paste0("HCP_results/5k_results/604_Dice_coefficient_PW_threshold",threshold,".rds"))
# saveRDS(all_dice, paste0("HCP_results/5k_results/604_overlap_size_PW_threshold",threshold,".rds"))
# # CLASSICAL ----
# result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/classical"
# result_files <- list.files(result_dir, full.names = TRUE)
# load("HCP_data/subjects.Rdata") # The subjects object (character vector)
# result_files <- grep(".rds", result_files, value = T)
# result_files <- grep("500_", result_files, value = T)
# hems <- c('left','right')
# threshold <- 1
# subject_dice <- vector('numeric', length = length(subjects))
# library(BayesfMRI)
# library(ciftiTools)
# ciftiTools.setOption('wb_path','/Applications/workbench')
# subjects_dice <- list()
# for(subj in subjects) {
#   visit_exc <- list()
#   for(v in 1:2) {
#     hems_exc <- list()
#     for(h in hems) {
#       L_or_R <- toupper(substring(h,1,1))
#       filen <- grep(paste0(subj,"_visit",v,"_",h), result_files, value = T)
#       result_subj <- readRDS(filen)
#       model_obj <- result_subj$GLMs_classical[[paste0("cortex",L_or_R)]]
#       nvox <- nrow(result_subj$betas_classical[[1]]$data[[paste0("cortex_",h)]])
#       K <- length(result_subj$beta_names)
#       class_act <- id_activations.classical(model_obj,alpha = 0.01,correction = "FWER",threshold = threshold)
#       in_mask <- model_obj[[1]]$mask
#       mat_out <- matrix(0, nvox, (K+2))
#       if(h == 'left') mat_out[,1:4] <- class_act$active[in_mask,c(1,4,2,3)]
#       if(h == 'right') mat_out[,c(1,2,5,6)] <- class_act$active[in_mask,c(1,4,2,3)]
#       hems_exc[[h]] <- mat_out
#     }
#     visit_exc[[v]] <- Reduce(rbind,hems_exc)
#   }
#   # subjects_dice[[subj]] <- mapply(function(v1,v2) {sum(v1 * v2) / mean(c(sum(v1),sum(v2))) },
#   #                                 v1 = split(visit_exc[[1]],col(visit_exc[[1]])),
#   #                                 v2 = split(visit_exc[[2]],col(visit_exc[[2]])))
#   subjects_dice[[subj]] <- mapply(function(v1,v2) {sum(v1 * v2)},
#                                   v1 = split(visit_exc[[1]],col(visit_exc[[1]])),
#                                   v2 = split(visit_exc[[2]],col(visit_exc[[2]])))
# }
# all_dice <- Reduce(rbind, subjects_dice)
# # saveRDS(all_dice, paste0("HCP_results/5k_results/604_classical_FWER_Dice_coefficient_PW_threshold",threshold,".rds"))
# saveRDS(all_dice, paste0("HCP_results/5k_results/604_classical_FWER_overlap_size_PW_threshold",threshold,".rds"))
#
# # # Plot the Dice Coefficients ----
# dice_result_files <- list.files("HCP_results/5k_results", full.names = T)
# dice_result_files <- grep("604_", dice_result_files, value = T)
# dice_result_files <- grep("_PW_", dice_result_files, value = T)
# dice_result_files <- grep("_Dice_", dice_result_files, value = T)
# dice_files <- list(Bayesian = grep("classical",dice_result_files, invert = T, value = T)[c(2,1,3)],
#                    Classical = grep("FWER",dice_result_files, value = T)[c(2,1,3)])
# dice_results <- sapply(dice_files, function(df) sapply(df,readRDS, simplify = F), simplify = F)
# task_names <- c("visual cue","tongue","right foot","right hand","left foot","left hand")
# names(dice_results$Bayesian) <- c("0%","0.5%","1%")
# names(dice_results$Classical) <- c("0%","0.5%","1%")
# library(tidyverse)
# #
# # # >> Boxplots ----
# reshape2::melt(dice_results) %>%
#   mutate(Var2 = task_names[Var2],
#          Var2 = as.factor(Var2),
#          L2 = factor(L2, levels = c("0%","0.5%","1%"))) %>%
#   ggplot() +
#   geom_boxplot(aes(x = value, y = L1, fill = L1)) +
#   labs(y = "Dice Coefficients", x = "") +
#   scale_fill_discrete("") +
#   facet_grid(Var2~L2, scales = "free_y") +
#   theme_bw() +
#   theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())
#
# reshape2::melt(dice_results) %>%
#   filter(L1 == "Bayesian") %>%
#   mutate(Var2 = task_names[Var2],
#          Var2 = as.factor(Var2),
#          L2 = factor(L2, levels = c("0%","0.5%","1%"))) %>%
#   ggplot() +
#   geom_boxplot(aes(x = L2, y = value)) +
#   labs(y = "Dice Coefficients", x = "Activation Threshold") +
#   facet_grid(Var2~., scales = "free_y") +
#   theme_bw() #+
#   theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())
#
# # >>>> Paired Hypothesis Testing ----
# bstrap_test <- function(x, mu = 0, num_resamp = 5000, alternative = "greater") {
#   bsamp_means <- sapply(seq(num_resamp), function(s){
#     out <- mean(sample(x, length(x), replace = T), na.rm = T)
#     return(out)
#   })
#   leq_geq <- NULL
#   if(alternative == "greater") leq_geq <- `>=`
#   if(alternative == "less") leq_geq <- `<=`
#   p_out <- sum(leq_geq(bsamp_means,mu), na.rm = T) / num_resamp
#   return(1 - p_out)
# }
# task_names <- c("visual cue","tongue","right foot","right hand","left foot","left hand")
# sig_df <- sapply(dice_results, function(method_results) {
#   sapply(method_results, function(thr_results) {
#     rownames(thr_results) <- NULL
#     return(thr_results)
#   }, simplify = F)
# }, simplify = F) %>%
#   reshape2::melt() %>%
#   pivot_wider(names_from = L1) %>%
#   mutate(diff = Bayesian - Classical) %>%
#   group_by(Var2, L2) %>%
#   summarize(
#     pvalue = bstrap_test(diff)
#   ) %>%
#   ungroup() %>%
#   mutate(aster = ifelse(pvalue < 0.01, "***", ifelse(
#     pvalue < 0.05, "**", ifelse(pvalue < 0.1, "*", "")
#   )),
#   L1 = "Bayesian", Var2 = task_names[Var2],
#   L2 = factor(L2, levels = c("0%","0.5%","1%")),
#   pval = paste("p =",format(round(pvalue,4), nsmall = 4, scientific = F)),
#   pval = ifelse(pval == "p = 0.0000","p < 0.0001",pval))
#
# # # >> Violin plots ----
# GeomSplitViolin <- ggproto(
#   "GeomSplitViolin",
#   GeomViolin,
#   draw_group = function(self, data, ..., draw_quantiles = NULL) {
#     data <-
#       transform(
#         data,
#         xminv = x - violinwidth * (x - xmin),
#         xmaxv = x + violinwidth * (xmax - x)
#       )
#     grp <- data[1, "group"]
#     newdata <-
#       plyr::arrange(transform(data, x = if (grp %% 2 == 1)
#         xminv
#         else
#           xmaxv), if (grp %% 2 == 1)
#             y
#         else-y)
#     newdata <-
#       rbind(newdata[1,], newdata, newdata[nrow(newdata),], newdata[1,])
#     newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <-
#       round(newdata[1, "x"])
#
#     if (length(draw_quantiles) > 0 &
#         !scales::zero_range(range(data$y))) {
#       stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
#                                                 1))
#       quantiles <-
#         ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
#       aesthetics <-
#         data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
#       aesthetics$alpha <-
#         rep(1, nrow(quantiles))
#       both <- cbind(quantiles, aesthetics)
#       quantile_grob <-
#         GeomPath$draw_panel(both, ...)
#       ggplot2:::ggname("geom_split_violin",
#                        grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
#     }
#     else {
#       ggplot2:::ggname("geom_split_violin",
#                        GeomPolygon$draw_panel(newdata, ...))
#     }
#   }
# )
#
# geom_split_violin <-
#   function(mapping = NULL,
#            data = NULL,
#            stat = "ydensity",
#            position = "identity",
#            ...,
#            draw_quantiles = NULL,
#            trim = TRUE,
#            scale = "area",
#            na.rm = FALSE,
#            show.legend = NA,
#            inherit.aes = TRUE) {
#     layer(
#       data = data,
#       mapping = mapping,
#       stat = stat,
#       geom = GeomSplitViolin,
#       position = position,
#       show.legend = show.legend,
#       inherit.aes = inherit.aes,
#       params = list(
#         trim = trim,
#         scale = scale,
#         draw_quantiles = draw_quantiles,
#         na.rm = na.rm,
#         ...
#       )
#     )
#   }
# vplot <- reshape2::melt(dice_results) %>%
#   mutate(Var2 = task_names[Var2],
#          Var2 = as.factor(Var2),
#          L2 = factor(L2, levels = c("0%","0.5%","1%"))) %>%
#   ggplot(aes(y = value, x = 1, fill = L1)) +
#   geom_split_violin() +
#   labs(y = "Dice Coefficients", x = "") +
#   scale_fill_discrete("") +
#   facet_grid(Var2~L2) +
#   theme_bw() +
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
#   coord_flip()
#
# final_vplot <- vplot + geom_text(aes(label = aster, y = 0.85, x = 1.4),data = sig_df) +
#   geom_text(aes(label = pval, y = 0.85, x = 1.33),data = sig_df, size = 2.2)
# ggsave("plots/604_Dice_violin_plots.png", width = 7, height = 8, plot = final_vplot)
#
# # Plot Overlap Size x Dice Coefficient ----
# dice_result_files <- list.files("HCP_results/5k_results", full.names = T)
# dice_result_files <- grep("604_", dice_result_files, value = T)
# dice_result_files <- grep("_PW_", dice_result_files, value = T)
# dice_result_files <- grep("_Dice_", dice_result_files, value = T)
# dice_files <- list(Bayesian = grep("classical",dice_result_files, invert = T, value = T)[c(2,1,3)],
#                    Classical = grep("FWER",dice_result_files, value = T)[c(2,1,3)])
# dice_results <- sapply(dice_files, function(df) sapply(df,readRDS, simplify = F), simplify = F)
# task_names <- c("visual cue","tongue","right foot","right hand","left foot","left hand")
# names(dice_results$Bayesian) <- c("0%","0.5%","1%")
# names(dice_results$Classical) <- c("0%","0.5%","1%")
# area_result_files <- list.files("HCP_results/5k_results", full.names = T)
# area_result_files <- grep("604_", area_result_files, value = T)
# area_result_files <- grep("_PW_", area_result_files, value = T)
# area_result_files <- grep("_overlap_", area_result_files, value = T)
# area_files <- list(Bayesian = grep("classical",area_result_files, invert = T, value = T)[c(2,1,3)],
#                    Classical = grep("FWER",area_result_files, value = T)[c(2,1,3)])
# area_results <- sapply(area_files, function(df) sapply(df,readRDS, simplify = F), simplify = F)
# names(area_results$Bayesian) <- c("0%","0.5%","1%")
# names(area_results$Classical) <- c("0%","0.5%","1%")
# library(tidyverse)
#
# reshape2::melt(dice_results, value.name = "dice") %>%
#   cbind(overlap = reshape2::melt(area_results, value.name = "overlap")$overlap) %>%
#   mutate(Var2 = task_names[Var2],
#          Var2 = as.factor(Var2),
#          L2 = factor(L2, levels = c("0%","0.5%","1%"))) %>%
#   ggplot() +
#   geom_point(aes(x = overlap, y = dice, color = L1)) +
#   facet_grid(Var2 ~ L2, scales = "free_x") +
#   scale_color_discrete("") +
#   labs(y = "Dice Coefficient", x = "Size of Overlapping Area\n(# data locations)") +
#   theme_bw()

# Calculate all activations ----
library(BayesfMRI)
library(parallel)
cl <- makeCluster(3)
result_dir <- "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session"
# result_dir <- "~/Desktop"
result_files <- list.files(result_dir, full.names = T)
load("HCP_data/subjects.Rdata")
# load("/Users/Shared/Lab_Data/HCP_Motor_Task_Dan/subjects.Rdata")
subjects <- subjects[c(1,2,4)]
all_activations <-
  # sapply(subjects, function(subject) {
  parSapply(cl, subjects, function(subject, result_files){
    library(BayesfMRI)
    subject_activations <- sapply(c("LR","RL"), function(session) {
      session_activations <- sapply(paste0("visit",1:2), function(visit) {
        hem_activations <- sapply(c("left","right"), function(hem) {
          L_or_R <- toupper(substring(hem,1,1))
          result_obj <-
            readRDS(grep(
              paste0(subject, "_", visit, "_", hem, "_5k_session", session),
              result_files,
              value = T
            ))
          threshold_activations <- sapply(c(0,0.5,1), function(threshold) {
            bayes_activations <-
              id_activations_cifti(
                model_obj = result_obj,
                alpha = 0.01,
                method = "Bayesian",
                threshold = threshold
              )
            classical_activations <-
              id_activations_cifti(
                model_obj = result_obj,
                alpha = 0.01,
                method = "classical",
                threshold = threshold
              )
            class_act <- classical_activations$activations[[paste0("cortex",L_or_R)]]$active
            class_act <- matrix(as.numeric(class_act), nrow = nrow(class_act), ncol = ncol(class_act))
            colnames(class_act) <- colnames(bayes_activations$activations[[paste0("cortex",L_or_R)]]$active)
            bayes_act <- bayes_activations$activations[[paste0("cortex",L_or_R)]]$active
            return(
              list(
                Bayes = bayes_act,
                Classical = class_act
              )
            )
          }, simplify = FALSE)
          return(threshold_activations)
        }, simplify = FALSE)
        return(hem_activations)
      }, simplify = FALSE)
      return(session_activations)
    }, simplify = FALSE)
    return(subject_activations)
  }, result_files = result_files, simplify = FALSE)

  # })

saveRDS(all_activations, "/Volumes/GoogleDrive/My Drive/BayesGLM_Validation/5k_results/individual/PW/single_session/605_all_activations.rds")
# saveRDS(all_activations, file = "~/Desktop/605_all_activations.rds")
