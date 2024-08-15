library(pacman)
p_load(dplyr, magrittr, ggplot2, psych, data.table, tidyr, janitor, hablar, mirt, faux)

# Setup item parameters and factors

factor_names <- c("Factor1", "Factor2")
intercept_list = vector(mode = "list", length = length(factor_names))
names(intercept_list) <- factor_names

fac_item_names <- list(Factor1 = paste0("F1_item_", 1:10),
                       Factor2 = paste0("F2_item_", 1:10))

fac_item_thresh <- list(Factor1 = c(3, 2, 1, 4, 2, 5, 1, 3, 4, 5),
                        Factor2 = c(3, 2, 1, 4, 2, 5, 1, 3, 4, 5))

for(i in factor_names) {
  intercept_list[[i]] <- data.matrix(matrix(NA, nrow = lengths(fac_item_names)[i], ncol = max(fac_item_thresh[[i]])))
  
  for(j in 1:lengths(fac_item_names)[i]) {
    intercept_list[[i]][j, 1:fac_item_thresh[[i]][j]] <- sort(rnorm(fac_item_thresh[[i]][j], mean = 0, sd = 1), decreasing = TRUE)
  }
}


simulate_fscores <- function(n = 50000,
                             fac_mus = c(0, -0.03),
                             fac_sigmas = c(1, 0.15),
                             fac_names = factor_names,
                             fac_rhos = c(r12 = .3), # set rhos with the upper right triangle
                             item_names = fac_item_names,
                             n_thresholds = fac_item_thresh,
                           #  intercepts = intercept_list,
                            # slopes = list(Factor1 = data.matrix(rnorm(length(item_names[["Factor1"]]), mean = 2, sd = 1)),
                                           # Factor2 = data.matrix(rnorm(length(item_names[["Factor2"]]), mean = 2, sd = 1)),
                                           # Factor3 = data.matrix(rnorm(length(item_names[["Factor3"]]), mean = 2, sd = 1)),
                                           # Factor4 = data.matrix(rnorm(length(item_names[["Factor4"]]), mean = 2, sd = 1))),
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
  return(sim_thetas)
 
  
}


#Generate a simulation dataset based on N

N <- 50000
num_assess <- 10
sim_ids <- sort(rep(1:N, num_assess))
assess_times <- rep(0:9, N)

sim_data <- as.data.frame(bind_cols(sim_ids, assess_times))
names(sim_data) <- c("sim_id", "assess_time")

foo <- foo %>% mutate(sim_id = row_number())

sim_data <- inner_join(sim_data, foo, by="sim_id") %>%
  mutate(cog = Factor1 + assess_time*Factor2)


