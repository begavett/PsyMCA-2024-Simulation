library(pacman)
p_load(dplyr, magrittr, ggplot2, psych, data.table, tidyr, janitor, hablar,
       mirt, faux, Hmisc, arules, discretization, synthpop, lavaan, stringr, 
       sqldf, semTools)

source(here::here("recode.R"))

#---- Load and preprocess the data ----
# Use the dialog box to open ~/Dropbox/Projects/PsyMCA-2024-Simulation/data/HCAP_harmonized_data_PsyMCA-v6.rdata
#load(file.choose())

original_data <- `HCAP_harmonized_data_PsyMCA-v6`
# synthetic_data <- syn(original_data %>% select(-c(id_hrs, id_elsa, id_mexcog,
#                                                   id_lasidad, id_haalsi, 
#                                                   id_charls,
#                                                   demographics,
#                                                   COGNITIVE_FACTORS,
#                                                   economics,
#                                                   HEALTH_AND_FUNCTIONING,
#                                                   CESD_ITEMS,
#                                                   behaviors,
#                                                   COGNITIVE_TEST_INDICATORS,
#                                                   admin)))

hcap_haalsi_memory <- original_data %>%
  select(id, study, age:educattain_resp,
         fmem_bayes3, # Factor score for mem (Bayes with 3 posterior draws)
         fmem, # Factor score for mem
         fmem_se, # SE of Factor score for mem
         cerad_imm_1 = u201, 
         cerad_imm_2 = u202, 
         cerad_del1 = u210, 
         cerad_del2 = u211, 
         cerad_recog1 = u213, 
         cerad_recog2 = u214,
         cerad_cpd = u212, 
         lmi = u203, 
         lmd = u207, 
         lm_recog = u215, 
         mmse_3wi = u258,
         mmse_3wd = u209,
         brave_imm = u232, 
         brave_del = u230,
         hearing_aid1,
         hearing,
         hyper,
         smokestat,
         mweight,
         mheight,
         cesd1:cesd8, # Threshold = 4
         t2diab,
         drinksperwk,
         vision_dis,
         vision_near,
         fqmdactx,
         fqvgactx,
         mingrp) %>%
  filter(study %in% c(1, 5)) %>%
  mutate(cohort = ifelse(study == 1, "HRS-HCAP", "HAALSI-HCAP"),
         less_edu = ifelse(educattain_resp >= 3, 0, 1), # using upper secondary as threshold
         bmi = mweight/mheight^2, 
         obese = ifelse(bmi >= 30, 1, 0),
         excessive_drink = ifelse(drinksperwk > 12, 1, 0),
         eversmoker = ifelse(smokestat == 3, 0, 1),
         depression = case_when(is.na(cesd1) & is.na(cesd2) & is.na(cesd3) &
                                  is.na(cesd4) & is.na(cesd5) & is.na(cesd6) &
                                  is.na(cesd7) & is.na(cesd8) ~ NA_real_,
                                rowSums(across(cesd1:cesd8), na.rm = T) >= 4 ~ 1,
                                rowSums(across(cesd1:cesd8), na.rm = T) < 4 ~ 0),
         vision_loss = case_when(
           vision_dis == 5 | vision_near == 5 ~ 1,
           is.na(vision_dis) & is.na(vision_near) ~ NA_real_,
           TRUE ~ 0))

#---- Set population level parameters ----
pop_hyppar <- hcap_haalsi_memory %>%
  group_by(study) %>%
  summarise(n = n(),
            m_mem = mean_(fmem),
            sd_mem = sd_(fmem)) %>%
  ungroup() %>%
  mutate(m_mem = m_mem - m_mem[study == 1],
         sd_mem = sd_mem / sd_mem[study == 1],
         var_mem = sd_mem^2)


pop_hyppar
cormat_hcap <- hcap_haalsi_memory %>%
  filter(study == 1) %>%
  select(fmem, lmi, lmd, lm_recog, cerad_cpd, hyper, depression) %>%
  cor(use = "pairwise")


hcap_haalsi_memory %>%
  filter(study == 5) %>%
  select(fmem, lmi, lmd, lm_recog, cerad_cpd, hyper, depression) %>%
  cor(use = "pairwise")

mu_haalsi <- hcap_haalsi_memory %>%
  filter(study == 5) %>%
  select(fmem, lmi, lmd, lm_recog, cerad_cpd, hyper, depression) %>%
  colMeans(na.rm = TRUE)

mu_haalsi["fmem"] <- mu_haalsi["fmem"] - {
  hcap_haalsi_memory %>%
    filter(study == 1) %>%
    pull(fmem) %>%
    mean_()
}

sd_haalsi <- hcap_haalsi_memory %>%
  filter(study == 5) %>%
  select(fmem, lmi, lmd, lm_recog, cerad_cpd, hyper, depression) %>%
  apply(., 2, sd, na.rm = TRUE)

sd_haalsi["fmem"] <- sd_haalsi["fmem"] / {
  hcap_haalsi_memory %>%
    filter(study == 1) %>%
    pull(fmem) %>%
    sd_()
}

pop_mem_m <- '

mem =~ lmi + lmd + lm_recog + cerad_cpd

mem ~ hyper + depression

'

pop_mem_1g_f <- sem(pop_mem_m, 
                    data = hcap_haalsi_memory %>% filter(cohort == "HRS-HCAP"),
                    estimator = "mlr", 
                    std.lv = TRUE)

pop_mem_2g_f <- sem(pop_mem_m, 
                    data = hcap_haalsi_memory,
                    estimator = "mlr", 
                    std.lv = TRUE,
                    group = "cohort")

summary(pop_mem_1g_f, fit.measures = TRUE, standardized = TRUE)
summary(pop_mem_2g_f)
partable(pop_mem_1g_f)
pop_mem_pars <- partable(pop_mem_2g_f) %>%
  mutate(formula = paste(lhs, op, rhs))
partable_1g <- partable(pop_mem_1g_f) %>%
  mutate(formula = paste(lhs, op, rhs))
partable_1g <- partable_1g %>%
  filter(!str_detect(formula, "depression"),
         !str_detect(formula, "hyper"))

for(i in 1:nrow(partable_1g)) {
  pop_mem_pars$est[pop_mem_pars$formula == partable_1g$formula[i]] <- partable_1g$est[i]
}

pop_mem_pars$est[pop_mem_pars$formula == "mem ~1 " & pop_mem_pars$group == 2] <- mu_haalsi["fmem"]
pop_mem_pars$est[pop_mem_pars$formula == "mem ~~ mem" & pop_mem_pars$group == 2] <- sd_haalsi["fmem"]


pop_mem_pars$free <- 0
pop_mem_pars$formula <- NULL


#---- Simulation and discretizatin function ----
process_data <- function(suffix, iter) {
  cat(paste("Simulating data:", suffix))
  set.seed(iter)
  # Simulate the data based on population parameters
  sim_data <- simulateData(model = pop_mem_pars, 
                           model.type = "sem", 
                           std.lv = TRUE,
                           estimator = "mlr",
                           sample.nobs = c(20*pop_hyppar$n[pop_hyppar$study == 1], 
                                           20*pop_hyppar$n[pop_hyppar$study == 5])) %>%
    mutate(across(lmi:cerad_cpd, round),
           cohort = ifelse(group == 1, "HRS", "HAALSI"))
  
  
  mem_cfa <- '
  mem =~ lmi + lmd + lm_recog + cerad_cpd
  '
  
  test.seq <- c("loadings","intercepts")
  meq.list <- list()
  for (i in 0:length(test.seq)) {
    if (i == 0L) {
      meq.label <- "configural"
      group.equal <- ""
      long.equal <- ""
    } else {
      meq.label <- test.seq[i]
      group.equal <- test.seq[1:i]
      long.equal <- test.seq[1:i]
    }
    meq.list[[meq.label]] <- measEq.syntax(configural.model = mem_cfa,
                                           data = sim_data,
                                           ID.fac = "std.lv",
                                           group = "cohort",
                                           group.equal = group.equal,
                                           return.fit = TRUE)
  }
  
  mi_out <- compareFit(meq.list)
  summary(mi_out)
  
  # Discretization
  for (n_obs in 10:500){
    if (suffix == "ef"){
      test <- sim_data %>%
        mutate(across(c(lmi, lmd, lm_recog, cerad_cpd),
                      \(x) as.numeric(cut2(x, g = 10, m = n_obs)), .names = "{.col}_ef")) %>%
        select(group, ends_with("ef"))
    } else if (suffix == "ei"){
      test <- recodeOrdinal(sim_data,
                            varlist_orig = c("lmi", "lmd", "lm_recog", "cerad_cpd"),
                            varlist_tr = c("lmi_ei", "lmd_ei", "lm_recog_ei", "cerad_cpd_ei"),
                            ncat = 10,
                            nobs = n_obs) %>%
        select(group, ends_with("ei"))
    }
    
    count_tib <- test %>%
      filter(group == 2) %>%
      select(ends_with(suffix)) %>%
      pivot_longer(cols = ends_with(suffix)) %>%
      group_by(name, value) %>%
      dplyr::count()
    
    if(all(count_tib$n >= 10)){
      print(n_obs)
      break
    } else if (n_obs == 64){
      stop("No solution found n obs for bin")
    }
  }
  
  
  #n_obs <- 2
  
  if (suffix == "ef"){
    sim_data_r <- sim_data %>%
      mutate(across(c(lmi, lmd, lm_recog, cerad_cpd), 
                    \(x) as.numeric(cut2(x, g = 10, m = n_obs)), 
                    .names = "{.col}_ef")) 
  } else if (suffix == "ei"){
    sim_data_r <- recodeOrdinal(
      sim_data, 
      varlist_orig = c("lmi", "lmd", "lm_recog", "cerad_cpd"),
      varlist_tr = c("lmi_ei", "lmd_ei", "lm_recog_ei", "cerad_cpd_ei"),
      ncat = 10, nobs = n_obs)
  }
  
  # Define column names based on the suffix
  item_names <- c(paste0("lmi_", suffix), paste0("lmd_", suffix), 
                  paste0("lm_recog_", suffix), paste0("cerad_cpd_", suffix))
  
  # Define model formula based on the suffix
  model_formula <- sprintf('
  mem =~ %s
  
  mem ~ hyper + depression
  ', paste(item_names, collapse = " + "))
  
  # Process the data
  cat(paste("Processing suffix:", suffix))
  
  # Fit SEM model
  sem_model <- sem(model_formula, data = sim_data_r, ordered = TRUE, 
                   estimator = "wlsmv", std.lv = TRUE)
  sem_summary <- summary(sem_model, fit.measures = TRUE, standardize = TRUE)
  
  # Fit MIRT model
  mirt_model <- mirt(sim_data_r %>% select(ends_with(suffix)), 
                     model = 1, itemtype = "graded")
  mirt_coef <- coef(mirt_model, simplify = TRUE)
  
  # Randomly select a linking item
  linking_item_number <- sample(1:4, 1)
  linking_item <- item_names[linking_item_number]
  
  # Fit HRS model
  fit_hrs <- mirt(sim_data_r %>% filter(cohort == "HRS") %>% select(ends_with(suffix)),
                  model = 1, itemtype = "graded", SE = FALSE, 
                  technical = list(NCYCLES = 50000))
  bank_hrs <- mod2values(fit_hrs)
  
  # Update study column
  sim_data_r <- sim_data_r %>%
    mutate(study = ifelse(group == 1, "A_HRS", "B_HAALSI"))
  
  # Fit HAALSI model
  fit_haalsi <- mirt(sim_data_r %>% filter(cohort == "HAALSI") %>% select(ends_with(suffix)),
                     model = 1, itemtype = "graded", SE = FALSE, 
                     technical = list(NCYCLES = 50000))
  bank_haalsi <- mod2values(fit_haalsi)
  
  # Fit multiple group model
  fit_mg <- multipleGroup(sim_data_r %>% select(ends_with(suffix)), 
                          model = 1, itemtype = "graded", 
                          group = sim_data_r$study, SE = FALSE, 
                          invariance = c(linking_item, "free_means", "free_var"),
                          technical = list(NCYCLES = 50000))
  mg_coef <- coef(fit_mg, simplify = TRUE)
  
  # Perform DIF analysis
  dif_res <- DIF(fit_mg, which.par = c("a1", paste0("d", 1:9)), 
                 scheme = "add_sequential", 
                 technical = list(NCYCLES = 50000))
  dif_items <- rownames(dif_res)
  
  anchor_items <- setdiff(item_names, dif_items)
  
  # Get parameter table from bifactor model
  hz_mg_pars <- mod2values(fit_mg)
  
  # Constrain all HRS item parameters to their banked values
  for(i in item_names) {
    hz_mg_pars[hz_mg_pars$item == i, "value" ] <- rep(bank_hrs[bank_hrs$item == i, "value"], 2)
  }
  
  if(length(dif_items) > 0) {
    # Set DIF items
    for(i in setdiff(item_names, anchor_items)) {
      hz_mg_pars[hz_mg_pars$item == i & hz_mg_pars$group == "A_HRS", "value"] <- bank_hrs[bank_hrs$item == i, "value"]
      
      if (length(bank_haalsi[bank_haalsi$item == i, "value"]) == length(hz_mg_pars[hz_mg_pars$item == i & hz_mg_pars$group == "B_HAALSI", "value"])) {
        hz_mg_pars[hz_mg_pars$item == i & hz_mg_pars$group == "B_HAALSI", "value"] <- bank_haalsi[bank_haalsi$item == i, "value"]
      } else {
        use_pars <- paste("B_HAALSI", bank_haalsi$item, bank_haalsi$name, sep = "_")
        hz_mg_pars <- hz_mg_pars %>%
          mutate(label = paste(group, item, name, sep = "_")) %>%
          filter(label %in% use_pars) %>%
          select(-label)
        hz_mg_pars[hz_mg_pars$item == i & hz_mg_pars$group == "B_HAALSI", "value"] <- bank_haalsi[bank_haalsi$item == i, "value"]
      }
    }
  }
  
  hz_mg_pars$est <- FALSE
  hz_mg_pars <- hz_mg_pars %>%
    mutate(est = ifelse(item %in% setdiff(item_names, anchor_items) & value != 0, TRUE, est)) %>%
    mutate(est = ifelse(group == "B_HAALSI" & name %in% c("MEAN_1", "COV_11"), TRUE, est))
  
  # Check item parameters
  hz_mg_pars_check <- hz_mg_pars %>%
    filter(item %in% item_names, name %in% c("a1", paste0("d", 1:9))) %>%
    arrange(item, group, name)
  
  # Fit final multiple group model
  hz_model <- multipleGroup(sim_data_r %>% select(ends_with(suffix)),
                            model = 1, group = sim_data_r$study,
                            itemtype = "graded", pars = hz_mg_pars,
                            technical = list(NCYCLES = 50000), SE = FALSE)
  hz_model_coef <- coef(hz_model, simplify = TRUE)
  
  # Fix all parameters to their estimates and generate factor scores
  hz_model_pars <- mod2values(hz_model)
  hz_model_pars$est <- FALSE
  hz_model1 <- multipleGroup(sim_data_r %>% select(ends_with(suffix)),
                             model = 1, group = sim_data_r$study,
                             itemtype = "graded", pars = hz_model_pars,
                             technical = list(NCYCLES = 50000))
  hz_model_fscores <- fscores(hz_model1, full.scores.SE = FALSE, QMC = TRUE) %>%
    as_tibble()
  

  sim_data_r_fs <- sim_data_r %>%
    bind_cols(hz_model_fscores)
  

  lm_model <- lm(F1 ~ depression + hyper, data = sim_data_r_fs) %>% 
    summary()
  
  lm_model_hrs <- lm(F1 ~ depression + hyper, data = sim_data_r_fs %>% 
                       filter(study == "A_HRS")) %>% 
    summary()
  
  lm_model_haalsi <- lm(F1 ~ depression + hyper, data = sim_data_r_fs %>% 
                          filter(study == "B_HAALSI")) %>% 
    summary()
  
  
  return(list(sim_data = sim_data, sem_summary = sem_summary, mirt_coef = mirt_coef, mg_coef = mg_coef, 
              dif_items = dif_items, hz_mg_pars_check = hz_mg_pars_check, 
              hz_model_coef = hz_model_coef, lm_model = lm_model,
              lm_model_hrs = lm_model_hrs, lm_model_haalsi = lm_model_haalsi))
}

lm_full_coef <- tibble()
lm_haalsi_coef <- tibble()
lm_hrs_coef <- tibble()
dif_items <- tibble()
for (s in 1:12){
  for (type in c("ef", "ei")){
    lm_full_coef %<>%
      rbind(broom::tidy(results[[paste0("iter_", s, "_", type)]]$lm_model) %>%
              mutate(num = s, type = type))
    lm_haalsi_coef %<>%
      rbind(broom::tidy(results[[paste0("iter_", s, "_", type)]]$lm_model_haalsi) %>%
              mutate(num = s, type = type))
    lm_hrs_coef %<>%
      rbind(broom::tidy(results[[paste0("iter_", s, "_", type)]]$lm_model_hrs) %>%
              mutate(num = s, type = type))
    dif_items %<>%
      rbind(tibble(num = s, type = type,
                   dif_items = results[[paste0("iter_", s, "_", type)]]$dif_items))
    
  }
}

lm_haalsi_coef %>%
  group_by(term, type) %>%
  dplyr::summarize(mean_est = mean(estimate),
                   mean_se = mean(std.error))

lm_hrs_coef %>%
  group_by(term, type) %>%
  dplyr::summarize(mean_est = mean(estimate),
                   mean_se = mean(std.error))

# results_ef <- process_data("ef", iter = 1)
# results_ei <- process_data("ei", iter = 1)
# 
# # check if the same simulate data is used
# diffdf::diffdf(results_ef$sim_data, results_ei$sim_data)


results <- list()
nsims <- 20
successes <- 0
tries <- 10000
max_tries <- tries+100

while (successes < nsims & tries < max_tries) {
  
  tries <- tries + 1
  
  results_ef <- tryCatch(process_data("ef", iter = tries), error = function(e) e, finally = print("Done"))
  
  if(!"error" %in% class(results_ef)) {
    results_ei <- tryCatch(process_data("ei", iter = tries), error = function(e) e, finally = print("Done"))
    
    if(!"error" %in% class(results_ei)) {
      successes <- successes + 1
      results[[paste0("iter_", successes, "_ef")]] <- results_ef
      results[[paste0("iter_", successes, "_ei")]] <- results_ei
    }
  }
  cat(paste0("Tries = ", tries, "; Successes = ", successes, "\n"))
}

# # Parallel 
# # This wouldn't do parallel without errors yet
# p_load(doParallel, parallel)
# reqpackage <- c("dplyr", "ggplot2", "psych",
#                 "data.table", "tidyr", "janitor", "hablar", "mirt", "faux",
#                 "Hmisc", "arules", "discretization", "synthpop", "lavaan",
# "stringr", "sqldf")
# # Process data for both _ef and _ei suffixes
# cl <- makeCluster(4)
# registerDoParallel(cl)
# 
# foreach (i = 1:nsims, .packages = reqpackage) %dopar% {
#   # Run for both _ef and _ei suffixes
#   results_ef <- process_data("ef", iter = i)
#   results_ei <- process_data("ei", iter = i)
# 
#   # Store the results
#   results[[paste0("iter_", i, "_ef")]] <- results_ef
#   results[[paste0("iter_", i, "_ei")]] <- results_ei
# 
#   # # Collect DIF items for both suffixes
#   # dif_items_df <- rbind(dif_items_df,
#   #                       data.frame(Iteration = i, Suffix = "ef",
#   #                                  DIF_Item = results_ef$dif_items,
#   #                                  stringsAsFactors = FALSE))
#   # dif_items_df <- rbind(dif_items_df,
#   #                       data.frame(Iteration = i, Suffix = "ei",
#   #                                  DIF_Item = results_ei$dif_items,
#   #                                  stringsAsFactors = FALSE))
#   saveRDS(results, here::here("output", paste0("results_", iter, ".rds")))
# }
# stopCluster(cl)



# You can access the results using results_ef and results_ei
###loop through sapply to get a table of dif_items
