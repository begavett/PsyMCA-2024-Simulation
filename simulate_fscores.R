library(pacman)
p_load(dplyr, magrittr, ggplot2, psych, data.table, tidyr, janitor, hablar, mirt)
source("src/simulate_fscores.R")
source("src/estimate_models.R")
source("src/run_simulations.R")

# Setup item parameters and factors

factor_names <- c("Factor1", "Factor2", "Factor3", "Factor4")
intercept_list = vector(mode = "list", length = length(factor_names))
names(intercept_list) <- factor_names

fac_item_names <- list(Factor1 = paste0("F1_item_", 1:10),
                       Factor2 = paste0("F2_item_", 1:10),
                       Factor3 = paste0("F3_item_", 1:10),
                       Factor4 = paste0("F4_item_", 1:10))

fac_item_thresh <- list(Factor1 = c(3, 2, 1, 4, 2, 5, 1, 3, 4, 5),
                        Factor2 = c(3, 2, 1, 4, 2, 5, 1, 3, 4, 5),
                        Factor3 = c(3, 2, 1, 4, 2, 5, 1, 3, 4, 5),
                        Factor4 = c(3, 2, 1, 4, 2, 5, 1, 3, 4, 5))

for(i in factor_names) {
  intercept_list[[i]] <- data.matrix(matrix(NA, nrow = lengths(fac_item_names)[i], ncol = max(fac_item_thresh[[i]])))
  
  for(j in 1:lengths(fac_item_names)[i]) {
    intercept_list[[i]][j, 1:fac_item_thresh[[i]][j]] <- sort(rnorm(fac_item_thresh[[i]][j], mean = 0, sd = 1), decreasing = TRUE)
  }
}


simulate_fscores <- function(n = 1000,
                             fac_mus = c(0, 0, 0, 0),
                             fac_sigmas = c(1, 1, 1, 1),
                             fac_names = factor_names,
                             fac_rhos = c(r12 = .5, r13 = .5, r14 = .5, r23 = .5, r24 = .5, r34 = .5), # set rhos with the upper right triangle
                             item_names = fac_item_names,
                             n_thresholds = fac_item_thresh,
                             intercepts = intercept_list,
                             slopes = list(Factor1 = data.matrix(rnorm(length(item_names[["Factor1"]]), mean = 2, sd = 1)),
                                           Factor2 = data.matrix(rnorm(length(item_names[["Factor2"]]), mean = 2, sd = 1)),
                                           Factor3 = data.matrix(rnorm(length(item_names[["Factor3"]]), mean = 2, sd = 1)),
                                           Factor4 = data.matrix(rnorm(length(item_names[["Factor4"]]), mean = 2, sd = 1))),
                             emp = FALSE,
                             quiet = TRUE) {
  
  if(any(duplicated(unlist(fac_item_names)))) {
    stop ("Duplicate item_names detected!")
  }
  
  sim_thetas <- faux::rnorm_multi(n = n,
                                  vars = length(fac_names),
                                  mu = fac_mus,
                                  sd = fac_sigmas,
                                  r = fac_rhos,
                                  varnames = fac_names,
                                  empirical = emp,
                                  as.matrix = FALSE) %>%
    janitor::clean_names(case = "none")
  
  sim_item_resp_list <- mirt_model_list <- sim_fs_list <- vector(mode = "list", length = length(fac_names))
  names(sim_item_resp_list) <- names(mirt_model_list) <- names(sim_fs_list) <- fac_names
  
  for(f in fac_names){
    
    sim_item_resp_list[[f]] <- mirt::simdata(a = slopes[[f]],
                                             d = intercepts[[f]],
                                             itemtype = ifelse(n_thresholds[[f]] == 1, "dich", "graded"),
                                             Theta = sim_thetas[, f]) %>%
      as_tibble() %>%
      setnames(item_names[[f]])
    
    mirt_model_list[[f]] <- mirt::mirt(data = sim_item_resp_list[[f]],
                                       model = 1,
                                       itemtype = "graded",
                                       verbose = !quiet) %>%
      { if(quiet) suppressMessages(.) else .}
    
    sim_fs_list[[f]] <- mirt::fscores(object = mirt_model_list[[f]], 
                                      full.scores.SE = TRUE) %>%
      data.frame() %>%
      set_names(paste0(f, c("_FS", "_FS_SE"))) %>%
      bind_cols(mirt::fscores(object = mirt_model_list[[f]], method = "plausible") %>%
                  data.frame() %>%
                  set_names(paste0(f, "_PV"))) %>%
      bind_cols(sim_item_resp_list[[f]])
  }
  return(sim_fs_list %>% bind_cols())
}
