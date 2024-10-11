library(pacman)
p_load(dplyr, magrittr, ggplot2, psych, data.table, tidyr, janitor, hablar,
       mirt, faux, Hmisc, arules, discretization, synthpop, lavaan, stringr, 
       sqldf, semTools, sjmisc, future.apply, readr, progressr, semhelpinghands,
       furrr, filelock)

write_to_log <- function(iteration, log_file) {
  lock <- filelock::lock(log_file)
  cat(paste(Sys.time(), "- Iteration:", iteration, "completed\n"), file = log_file, append = TRUE)
  filelock::unlock(lock)
}

#---- Load and preprocess the data ----
# Use the dialog box to open ~/Dropbox/Projects/PsyMCA-2024-Simulation/data/HCAP_harmonized_data_PsyMCA-v6.rdata
#load(file.choose())

load("~/Dropbox/PsyMCA2024_Data_Sharing_HCAP/HCAP_harmonized_data_PsyMCA-v6.rdata")

original_data <- `HCAP_harmonized_data_PsyMCA-v6`

# Continuous defined as items with > 10 response categories

hrs_haalsi_cog <- original_data %>%
  select(id, study, age:educattain_resp,
         fmem_bayes3, # Factor score for mem (Bayes with 3 posterior draws)
         fmem, # Factor score for mem
         fmem_se, # SE of Factor score for mem
         dom_cat = u105,
         mo_cat = u106,
         yr_cat = u107,
         dow_cat = u108,
         time_cat = u109,
         country_cat = u112,
         state_cat = u113,
         county_cat = u114,
         city_cat = u115,
         season_cat = u116,
         floor_cat = u117,
         address_cat = u118,
         bldg_cat = u140,
         cerad_imm1_cont = u201,
         cerad_imm2_cont = u202, 
         cerad_del1_cont = u210, 
         cerad_del2_cont = u211, 
         cerad_recog1_cont = u213,
         cerad_recog2_cont = u214,
         cerad_cpd_cont = u212,
         lmi_cont = u203,
         lmd_cont = u207,
         lm_recog_cont = u215,
         mmse_3wi_cat = u258,
         mmse_3wd_cat = u209,
         brave_imm_cat = u232,
         brave_del_cat = u230,
         ravens_cont = u301,
         num_series_cont = u302,
         traila_cont = u401,
         trailb_cont = u303,
         sim_cat = u304,
         token_cat = u306,
         spa_span_f_cat = u711,
         spa_span_b_cat = u712,
         gonogo_cat = u308,
         mot_prog_cat = u713,
         mmse_bac_sp_cat = u409,
         back_count_cont = u403,
         sdmt_cont = u402,
         cancel_sym_cont = u405,
         cancel_let_cont = u414,
         ser7s_cat = u412,
         bdaynam_cat = u404,
         fdaynam_cat = u725,
         cents_cat = u709,
         animals_cont = u501,
         cactus_cat = u502,
         scissors_cat = u504,
         watch_cat = u505,
         pencil_cat = u506,
         elbow_cat = u508,
         write_sentence_cat = u509,
         close_eyes_cat = u511,
         repetition_cat = u513,
         hammer_cat = u515,
         point2_cat = u517,
         market_cat = u518,
         threestep_cat = u523,
         pres_pm_cat = u701,
         phon_flu_cont = u714,
         bnt_cont = u715,
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
           TRUE ~ 0),
         cerad_imm_cont = ifelse(cohort == "HRS-HCAP", cerad_imm1_cont, cerad_imm2_cont),
         cerad_del_cont = ifelse(cohort == "HRS-HCAP", cerad_del1_cont, cerad_del2_cont),
         cerad_recog_cont = ifelse(cohort == "HRS-HCAP", cerad_recog1_cont, cerad_recog2_cont)
  )


hrs_haalsi_cog %>%
  select(ends_with(c("_cont"))) %>%
  lapply(table)

hrs_haalsi_cog %>%
  select(ends_with(c("_cont"))) %>%
  psych::describe()

sk <- hrs_haalsi_cog %>%
  select(ends_with(c("_cont"))) %>%
  psych::describe() %>%
  data.frame() %>%
  select(skew, kurtosis)

sk %>%
  cor()


# Population Parameters ---------------------------------------------------


# Establish the population
# Mu & Sigma of reference group is 0(1)
# Mu & Sigma of focal group is from HAALSI memory factor

mu_haalsi <- hrs_haalsi_cog %>%
  filter(study == 5) %>%
  select(fmem, cerad_imm_cont, cerad_del_cont, cerad_recog_cont, cerad_cpd_cont, 
         lmi_cont, lmd_cont, lm_recog_cont, hyper, depression) %>%
  colMeans(na.rm = TRUE)

mu_haalsi["fmem"] <- mu_haalsi["fmem"] - {
  hrs_haalsi_cog %>%
    filter(study == 1) %>%
    pull(fmem) %>%
    mean_()
}

sd_haalsi <- hrs_haalsi_cog %>%
  filter(study == 5) %>%
  select(fmem, cerad_imm_cont, cerad_del_cont, cerad_recog_cont, cerad_cpd_cont, 
         lmi_cont, lmd_cont, lm_recog_cont, hyper, depression) %>%
  apply(., 2, sd, na.rm = TRUE)

sd_haalsi["fmem"] <- sd_haalsi["fmem"] / {
  hrs_haalsi_cog %>%
    filter(study == 1) %>%
    pull(fmem) %>%
    sd_()
}


## Test and population characteristics -------------------------------------

n_items <- 10
max_categories <- 10
min_cell <- 50
pop_n_sim_RG <- 5e05 # 500k
pop_n_sim_FG <- 5e05 # 500k
set.seed(1234)
loadings <- rnorm(n = n_items) %>%
  abs() %>%
  fisherz2r() %>%
  round(3)
set.seed(245)
skew_kurt <- MASS::mvrnorm(10, mu = colMeans(sk), Sigma = data.matrix(sk %>% cov())) %>%
  data.frame()


# Specify population model in reference group
population_model_rg <- paste0('cog =~ ', paste0(paste0(loadings, "*x", 1:n_items), collapse = " + "), 
                              "\n cog ~0*1
                              \n cog ~~ 1*cog")

cat(population_model_rg)

# Specify population model in focal group
population_model_fg <- paste0('cog =~ ', paste0(paste0(loadings, "*x", 1:n_items), collapse = " + "), 
                              "\n cog ~-.90*1
                              \n cog ~~ 0.75*cog")

cat(population_model_fg)

# Simulate population data for the reference group
pop_data_rg <- simulateData(model = population_model_rg, 
                            sample.nobs = pop_n_sim_RG,
                            model.type = "cfa",
                            estimator = "mlr", 
                            skewness = skew_kurt$skew,
                            kurtosis = skew_kurt$kurtosis,
                            seed = 48291) %>%
  mutate(group = "Reference") %>%
  as.data.table()

# Simulate population data for the focal group
pop_data_fg <- simulateData(model = population_model_fg, 
                            sample.nobs = pop_n_sim_FG,
                            model.type = "cfa",
                            estimator = "mlr",
                            skewness = skew_kurt$skew,
                            kurtosis = skew_kurt$kurtosis,
                            seed = 47819) %>%
  mutate(group = "Focal") %>%
  as.data.table()


pop_data_all <- rbind(pop_data_rg, pop_data_fg)
rm(pop_data_rg, pop_data_fg)

pop_data_all %>%
  sample_n(10000) %>%
  pivot_longer(cols = starts_with("x"),
               names_to = "item",
               values_to = "score") %>%
  mutate(item = factor(item, levels = paste0("x", 1:max(parse_number(item))))) %>%
  ggplot(aes(x = score, fill = group)) +
  geom_histogram(position = "dodge") +
  facet_wrap(~item, ncol = 5) +
  xlim(-5, 5)

psych::describeBy(pop_data_all, pop_data_all$group)
skew_kurt

cfa_1f <- paste0("cog =~ NA*", paste0("x", 1:10, collapse = " + "))

# assess fit of the model to the population data using multiple groups analysis
cfa(cfa_1f, data = pop_data_all, estimator = "mlr", group = "group",
    std.lv = TRUE,
    group.equal = c("loadings", "intercepts")) %>%
  summary(fit.measures = TRUE, standardized = TRUE)

pop_data_npsy <- copy(pop_data_all)

# Convert population data to a format that is more typical of actual neuropsych data - integers with lower and upper bounds of 0 and 20, respectively, some with non-normal distributions
# First 5 items have positive skew
# pop_data_npsy[, (paste0("x", 1:(n_items/2))) := lapply(.SD, function(x) replace(round(2 * (x + 0.5*abs(min(x))), 0), which(round(2 * (x + 0.5*abs(min(x))), 0) < 0), 0)), .SDcols = paste0("x", 1:(n_items/2))]
# 
# # Second 5 items have negative skew
# pop_data_npsy[, (paste0("x", (1 + (n_items/2)):n_items)) := lapply(.SD, function(x) replace(round(2 * (x + 1.5*abs(min(x))), 0), which(round(2 * (x + 1.5*abs(min(x))), 0) > 20), 20)), .SDcols = paste0("x", (1 + (n_items/2)):n_items)]
# 
pop_data_npsy[, paste0("x", 1:n_items) := lapply(.SD, function(x) round(3*(x - min(x)), 0)),  .SDcols = paste0("x", 1:n_items)]

psych::describeBy(pop_data_npsy, pop_data_npsy$group)
lapply(pop_data_npsy, table)

set.seed(444)
pop_data_npsy %>%
  sample_n(10000) %>%
  pivot_longer(cols = starts_with("x"),
               names_to = "item",
               values_to = "score") %>%
  mutate(item = factor(item, levels = paste0("x", 1:max(parse_number(item))))) %>%
  ggplot(aes(x = score, fill = group)) +
  geom_histogram(position = "dodge", binwidth = 1) +
  facet_wrap(~item, ncol = 5, scales = "free")

pop_cfa_npsy <- cfa(cfa_1f, data = pop_data_npsy, estimator = "mlr", group = "group",
                    std.lv = TRUE,
                    group.equal = c("loadings", "intercepts"))

summary(pop_cfa_npsy, fit.measures = TRUE, standardized = TRUE)

pop_cfa_model <- partable(pop_cfa_npsy)



# Recoding ----------------------------------------------------------------

# Recoding

# Equal Frequency

equal_frequency_recoding <- pop_data_npsy[, paste0("x", 1:n_items, "_ef") := lapply(.SD, function(x) {
  as.numeric(Hmisc::cut2(x, g = max_categories, m = min_cell))
}), .SDcols = paste0("x", 1:n_items)]

lapply(select(equal_frequency_recoding, ends_with("_ef")), table)

eq_freq_scheme <- data.frame(item = paste0("x", 1:n_items),
                             scheme = NA_character_)

for(i in 1:nrow(eq_freq_scheme)) {
  
  eq_freq_scheme$scheme[eq_freq_scheme$item == paste0("x", i)] <- table(equal_frequency_recoding[, get(paste0("x", i))], equal_frequency_recoding[, get(paste0("x", i, "_ef"))]) %>%
    data.frame() %>%
    filter(Freq > 0) %>%
    mutate(across(Var1:Var2, \(x) as.numeric(as.character(x)))) %>%
    dplyr::group_by(Var2) %>%
    dplyr::summarise(scheme = paste0(first(Var1), ":", last(Var1), "=", first(Var2), ";")) %>%
    mutate(scheme = if_else(row_number() == n(), gsub(":[0-9]+=", ":hi=", scheme), scheme)) %>%
    pull(scheme) %>%
    paste(collapse = "\n") %>%
    gsub("^[0-9]+:", "lo:", .)
  
}

eq_freq_scheme

# Equal Interval

equal_interval_recoding <- pop_data_npsy[, paste0("x", 1:n_items, "_ei") := lapply(.SD, function(x) {
  cut(x, breaks = max_categories, labels = FALSE, include.lowest = TRUE)
}), .SDcols = paste0("x", 1:n_items)]

lapply(select(equal_interval_recoding, ends_with("ei")), table)


eq_int_scheme <- data.frame(item = paste0("x", 1:n_items),
                            scheme = NA_character_)

for(i in 1:nrow(eq_int_scheme)) {
  
  eq_int_scheme$scheme[eq_int_scheme$item == paste0("x", i)] <- table(equal_interval_recoding[, get(paste0("x", i))], equal_interval_recoding[, get(paste0("x", i, "_ei"))]) %>%
    data.frame() %>%
    filter(Freq > 0) %>%
    mutate(across(Var1:Var2, \(x) as.numeric(as.character(x)))) %>%
    dplyr::group_by(Var2) %>%
    dplyr::summarise(scheme = paste0(first(Var1), ":", last(Var1), "=", first(Var2), ";")) %>%
    mutate(scheme = if_else(row_number() == n(), gsub(":[0-9]+=", ":hi=", scheme), scheme)) %>%
    pull(scheme) %>%
    paste(collapse = "\n") %>%
    gsub("^[0-9]+:", "lo:", .)
  
}

eq_int_scheme

# k-means clustering

kmeans_recoding <- pop_data_npsy[, paste0("x", 1:n_items, "_km") := lapply(.SD, function(x) {
  kmeans(x, centers = max_categories)$cluster
}), .SDcols = paste0("x", 1:n_items)]

lapply(select(kmeans_recoding, ends_with("km")), table)

kmeans_scheme <- data.frame(item = paste0("x", 1:n_items),
                            scheme = NA_character_)

for(i in 1:nrow(kmeans_scheme)) {
  
  kmeans_scheme$scheme[kmeans_scheme$item == paste0("x", i)] <- table(kmeans_recoding[, get(paste0("x", i))], kmeans_recoding[, get(paste0("x", i, "_km"))]) %>%
    data.frame() %>%
    arrange(Var1) %>%
    filter(Freq > 0) %>%
    mutate(Var2 = cumsum(!duplicated(Var2))) %>%
    mutate(across(Var1:Var2, \(x) as.numeric(as.character(x)))) %>%
    dplyr::group_by(Var2) %>%
    dplyr::summarise(scheme = paste0(first(Var1), ":", last(Var1), "=", first(Var2), ";")) %>%
    mutate(scheme = if_else(row_number() == n(), gsub(":[0-9]+=", ":hi=", scheme), scheme)) %>%
    pull(scheme) %>%
    paste(collapse = "\n") %>%
    gsub("^[0-9]+:", "lo:", .)
  
  
}

kmeans_scheme


# Manual scheme





# Combine recoding schemes

all_recoding_schemes <- eq_freq_scheme %>%
  mutate(type = "ef") %>%
  bind_rows(eq_int_scheme %>%
              mutate(type = "ei")) %>%
  bind_rows(kmeans_scheme %>%
              mutate(type = "km"))

rm(equal_frequency_recoding, equal_interval_recoding, kmeans_recoding)



# Functions ---------------------------------------------------------------


## Automated MIMIC modeling for DIF detection ------------------------------


mimic_dif <- function(config_model = cfa_1f, data = samp_data_npsy, grp = "group", dif_alpha = .01) {
  
  mimic0_syntax <- paste0(config_model, "\n cog ~ ", grp, "\n ", paste0("x", 1:(ncol(data) - 1), " ~ 0*", grp, collapse = "\n"))
  
  m0 <- cfa(mimic0_syntax, data = data, estimator = "mlr", std.lv = TRUE)
  
  i <- 0
  pval <- 0
  
  while(pval < dif_alpha) {
    ref <- get(paste0("m", i))
    mi_i <- modindices(ref, sort. = TRUE) %>%
      filter(rhs == grp) %>%
      slice(1) %>%
      pull(lhs)
    
    assign(paste0("m", i, "mi"), modindices(ref, sort. = TRUE) %>%
             filter(rhs == grp) %>%
             slice(1) %>%
             pull(lhs))
    
    assign(paste0("m", i + 1),
           cfa(gsub(paste0(get(paste0("m", i, "mi")), " ~ 0\\*group"), 
                    paste0(get(paste0("m", i, "mi")), " ~ NA\\*group"), 
                    ptable_to_syntax(partable(ref))), 
               data = data, 
               estimator = "mlr", 
               std.lv = TRUE)
    )
    
    assign(paste0("pval", i + 1), anova(ref, get(paste0("m", i + 1)))$`Pr(>Chisq)`[2])
    
    temp <- data.frame(candidate_item = get(paste0("m", i, "mi")),
                       p = round(get(paste0("pval", i + 1)), 4))
    if(i == 0) {
      results_df <- temp
    } else {
      results_df <- bind_rows(results_df, temp)
    }
    pval <- get(paste0("pval", i + 1))
    i <- i + 1
    if(i == ncol(data)) pval <- 1
  }
  return(results_df)
}

# mimic_dif()

## Automated Multiple groups analysis for DIF detection ------------------------------

mg_dif <- function(config_model = cfa_1f, data = samp_data_npsy, grp = "group", dif_alpha = .01) {
  
  mg0_syntax <- paste0(config_model, "\n ", paste0("x", 1:(ncol(data) - 1), " ~ c(ig1.", 1:(ncol(data) - 1), ", ig1.", 1:(ncol(data) - 1), ")*1", collapse = "\n"))
  
  mg0_syntax <- gsub("=~ NA\\*", "=~ c(lg1.1, lg1.1)*", mg0_syntax)
  
  for(itno in 2:(str_count(config_model, "\\+") + 1)) {
    mg0_syntax <- gsub(paste0("\\+ x", itno), paste0("\\+ c(lg1.", itno, ", lg1.", itno, ")*x", itno), mg0_syntax)
  }
  
  mg0_syntax <- paste0(mg0_syntax, "\n cog ~~ c(1, NA)*cog \n cog ~ c(0, NA)*1", collapse = "\n")
  
  m0 <- cfa(mg0_syntax, data = data, estimator = "mlr", std.lv = TRUE, group = "group")
  unadj_factor_scores <- lavPredict(m0, se = "standard", assemble = TRUE)
  unadj_factor_scores$se <- as.numeric(attr(lavPredict(m0, se = "standard"), "se")[[1]])
  names(unadj_factor_scores) <- c("unadj_cog_fs", "group", "unadj_cog_fs_se")
  i <- 0
  pval <- 0
  
  while(pval < dif_alpha) {
    ref <- get(paste0("m", i))
    ref_syntax <- get(paste0("mg", i, "_syntax"))
    
    # ts_i <- lavTestScore(ref)$uni %>%
    #   data.frame() %>%
    #   arrange(p.value) %>%
    #   mutate(adj.p = p.adjust(p.value, method = p.method)) %>%
    #   slice(1) %>%
    #   pull(lhs)
    
    modindices(ref)
    
    assign(paste0("m", i, "ts"), lavTestScore(ref)$uni %>%
             data.frame() %>%
             arrange(p.value) %>%
             slice(1) %>%
             select(lhs, p.value))
    
    assign(paste0("m", i, "free"), partable(ref) %>%
             filter(plabel == get(paste0("m", i, "ts"))$lhs) %>%
             mutate(code = paste(lhs, op, rhs)) %>%
             pull(code))
    
    # Old code, which only freed the intercept or the loading, but not both
    # assign(paste0("mg", i + 1, "_syntax"), ifelse(str_detect(get(paste0("m", i, "free")), " =~ "), gsub(paste0("1.", parse_number(get(paste0("m", i, "free"))), ")\\*x", parse_number(get(paste0("m", i, "free")))),
    #                                                                                                   paste0("2.", parse_number(get(paste0("m", i, "free"))), ")\\*x", parse_number(get(paste0("m", i, "free")))),
    #                                                                                                   ref_syntax),
    #                                             ifelse(str_detect(get(paste0("m", i, "free")), " ~"), gsub(paste0("1.", parse_number(get(paste0("m", i, "free"))), ")\\*1"),
    #                                                                                    paste0("2.", parse_number(get(paste0("m", i, "free"))), ")\\*1"),
    #                                                                                    ref_syntax), warning(paste0("syntax error found with", get(paste0("m", i, "free")))))))
    
    # Updated code, which frees both the intercept and loading of an item simultaneously if either its intercept or loading has evidence of DIF.
    assign(paste0("mg", i + 1, "_syntax"), 
           gsub(paste0("1.", parse_number(get(paste0("m", i, "free"))), ")\\*x", parse_number(get(paste0("m", i, "free")))),
                paste0("2.", parse_number(get(paste0("m", i, "free"))), ")\\*x", parse_number(get(paste0("m", i, "free")))), gsub(paste0("1.", parse_number(get(paste0("m", i, "free"))), ")\\*1"),
                                                                                                                                  paste0("2.", parse_number(get(paste0("m", i, "free"))), ")\\*1"),
                                                                                                                                  ref_syntax)))
    
    assign(paste0("m", i + 1),
           cfa(get(paste0("mg", i + 1, "_syntax")), 
               data = data, 
               estimator = "mlr", 
               std.lv = TRUE,
               group = "group")
    )
    
    assign(paste0("pval", i + 1), get(paste0("m", i, "ts"))$p.value)
    
    adj_factor_scores <- lavPredict(get(paste0("m", i + 1)), se = "standard", assemble = TRUE)
    adj_factor_scores$se <- as.numeric(attr(lavPredict(get(paste0("m", i + 1)), se = "standard"), "se")[[1]])
    names(adj_factor_scores) <- c("adj_cog_fs", "group", "adj_cog_fs_se")
    cn_factor_scores <- bind_cols(unadj_factor_scores, adj_factor_scores %>% select(-group)) %>%
      mutate(fs_diff = unadj_cog_fs - adj_cog_fs,
             pooled_se = sqrt(unadj_cog_fs_se * adj_cog_fs_se),
             scaled_diff = fs_diff / pooled_se,
             diff_gt_03 = abs(fs_diff) > 0.3,
             diff_gt_1se = abs(scaled_diff) > 1)
    
    
    temp <- data.frame(candidate_item = paste0("x", parse_number(get(paste0("m", i, "free")))),
                       # par_type = ifelse(str_detect(get(paste0("m", i, "free")), " =~ "), "Loading",
                       #                               ifelse(str_detect(get(paste0("m", i, "free")), " ~"), "Intercept", 
                       #                                      warning(paste0("syntax error found with", get(paste0("m", i, "free")))))),
                       p = round(get(paste0("pval", i + 1)), 4),
                       salient_diff_pct_03 = sum_(cn_factor_scores$diff_gt_03) / length(cn_factor_scores$diff_gt_03),
                       salient_diff_pct_1se = sum_(cn_factor_scores$diff_gt_1se) / length(cn_factor_scores$diff_gt_1se),
                       salient_diff_pct_03_rg = sum_(cn_factor_scores$diff_gt_03[cn_factor_scores$group == "Reference"]) / length(cn_factor_scores$diff_gt_03[cn_factor_scores$group == "Reference"]),
                       salient_diff_pct_1se_rg = sum_(cn_factor_scores$diff_gt_1se[cn_factor_scores$group == "Reference"]) / length(cn_factor_scores$diff_gt_1se[cn_factor_scores$group == "Reference"]),
                       salient_diff_pct_03_fg = sum_(cn_factor_scores$diff_gt_03[cn_factor_scores$group == "Focal"]) / length(cn_factor_scores$diff_gt_03[cn_factor_scores$group == "Focal"]),
                       salient_diff_pct_1se_fg = sum_(cn_factor_scores$diff_gt_1se[cn_factor_scores$group == "Focal"]) / length(cn_factor_scores$diff_gt_1se[cn_factor_scores$group == "Focal"]))
    
    
    
    if(i == 0) {
      results_df <- temp
    } else {
      results_df <- bind_rows(results_df, temp)
    }
    pval <- get(paste0("pval", i + 1))
    i <- i + 1
    if(i == ncol(data)) pval <- 1
  }
  
  
  return(results_df)
}

#mg_dif()

# DIF detection ---------------------------------------------------------

sim_dif_detection <- function(iter, 
                              suffixes = c("ef", "ei", "km"), 
                              n_RG = 3496, 
                              n_FG = 631, 
                              model_syntax = cfa_1f, 
                              model_partable = pop_cfa_model, 
                              recoding_scheme = all_recoding_schemes, 
                              min_cell = 10, 
                              dif_alpha = .01) {
  
  require(lavaan)
  require(dplyr)
  require(mirt)
  require(data.table)
  
  cat(paste(iter, ". Simulating data \n\n"))
  
  # Sample from population
  samp_data <- simulateData(model = model_partable, 
                            sample.nobs = c(n_RG, n_FG),
                            model.type = "cfa",
                            estimator = "mlr",
                            seed = iter) %>%
    mutate(group = case_when(group == 1 ~ "Reference",
                             group == 2 ~ "Focal"))
  
  samp_data_npsy <- samp_data %>%
    mutate(across(c(everything(), -group), \(x) round(x, 0)))
  
  n_items <- ncol(samp_data_npsy) - 1
  
  # cfa(model_syntax, data = samp_data_npsy, estimator = "mlr", group = "group",
  #     std.lv = TRUE,
  #     group.equal = c("loadings", "intercepts")) %>%
  #   summary(fit.measures = TRUE, standardized = TRUE)
  
  cat("Continuous model \n\n")
  
  # dif_continuous <- mimic_dif(config_model = model_syntax,
  #                             data = samp_data_npsy,
  #                             grp = "group",
  #                             dif_alpha = .05) %>%
  #   bind_rows(data.frame(candidate_item = paste0("x", 1:n_items),
  #                        p = NA)) %>%
  #   distinct(candidate_item, .keep_all = TRUE) %>%
  #   mutate(DIF = case_when(p < dif_alpha ~ 1,
  #                          p >= dif_alpha ~ 0,
  #                          is.na(p) ~ 0),
  #          type = "cn")
  
  dif_continuous <- mg_dif(config_model = model_syntax,
                           data = samp_data_npsy,
                           grp = "group",
                           dif_alpha = .01) %>%
    bind_rows(data.frame(candidate_item = paste0("x", 1:n_items),
                         p = NA)) %>%
    distinct(candidate_item, .keep_all = TRUE) %>%
    mutate(DIF = case_when(p < dif_alpha ~ 1,
                           p >= dif_alpha ~ 0,
                           is.na(p) ~ 0),
           type = "cn")
  
  
  
  for(suffix in suffixes) {
    
    samp_data_recoded <- samp_data_npsy %>%
      mutate(group = factor(group) %>%
               relevel(ref = "Reference"))
    
    # Recode
    for(i in 1:(ncol(samp_data_npsy) - 1)) {
      
      # Apply recoding scheme
      samp_data_recoded[paste0("x", i, "_", suffix)] <- car::recode(samp_data_npsy[,paste0("x", i)], recodes = all_recoding_schemes %>%
                                                                      filter(type == suffix,
                                                                             item == paste0("x", i)) %>%
                                                                      pull(scheme))
      
      # Collapse sparse cells (minimum n per cell required for both groups)
      samp_data_recoded[paste0("x", i, "_", suffix)] <- lordif::collapse(samp_data_recoded[, paste0("x", i, "_", suffix)], 
                                                                         group = as.numeric(factor(samp_data_recoded$group)),
                                                                         minCell = min_cell)
    }
    
    # Multiple group IRT model (GRM) -  all items as anchors
    mg_mod <- mirt::multipleGroup(data = samp_data_recoded %>%
                                    select(ends_with(suffix)),
                                  group = samp_data_recoded$group,
                                  model = 1,
                                  itemtype = "graded",
                                  #SE = TRUE,
                                  invariance = c(paste0("x", 1:n_items, "_", suffix), "free_mean", "free_var"),
                                  verbose = FALSE)
    
    assign(paste0(suffix, "_unadj_factor_scores"), fscores(mg_mod, full.scores.SE = TRUE) %>%
             as_tibble() %>%
             bind_cols(samp_data_recoded %>% select(group)) %>%
             set_names(c("unadj_cog_fs", "unadj_cog_fs_se", "group")))
    
    cat("DIF with MIRT - ", suffix, "\n\n")
    
    # DIF testing (constrained baseline model, iteratively remove anchor items)
    dif_mod <- mirt::DIF(mg_mod, 
                         which.par = coef(mg_mod, simplify = TRUE)$Reference$items %>% attr(., "dimnames") %>% purrr::pluck(2),
                         items2test = paste0("x", 1:n_items, "_", suffix),
                         scheme = "drop_sequential",
                         seq_stat = dif_alpha,
                         technical = list(NCYCLES = 50000),
                         verbose = FALSE,
                         max_run = n_items)
    
    if(nrow(dif_mod) > 0) {
      dif_mod_ret <- mirt::DIF(mg_mod, 
                               which.par = coef(mg_mod, simplify = TRUE)$Reference$items %>% attr(., "dimnames") %>% purrr::pluck(2),
                               items2test = paste0("x", 1:n_items, "_", suffix),
                               scheme = "drop_sequential",
                               seq_stat = dif_alpha,
                               technical = list(NCYCLES = 50000),
                               verbose = FALSE,
                               max_run = n_items,
                               return_seq_model = TRUE)
      
      assign(paste0(suffix, "_adj_factor_scores"), fscores(dif_mod_ret, full.scores.SE = TRUE) %>%
               as_tibble() %>%
               bind_cols(samp_data_recoded %>% select(group)) %>%
               set_names(c("adj_cog_fs", "adj_cog_fs_se", "group")))
      
      assign(paste0(suffix, "_factor_scores"),  bind_cols(get(paste0(suffix, "_unadj_factor_scores")), get(paste0(suffix, "_adj_factor_scores")) %>% select(-group)) %>%
               mutate(fs_diff = unadj_cog_fs - adj_cog_fs,
                      pooled_se = sqrt(unadj_cog_fs_se * adj_cog_fs_se),
                      scaled_diff = fs_diff / pooled_se,
                      diff_gt_03 = abs(fs_diff) > 0.3,
                      diff_gt_1se = abs(scaled_diff) > 1))
      
      assign(paste0("dif_results_", suffix), dif_mod %>%
               mutate(candidate_item = rownames(.)) %>%
               bind_rows(data.frame(candidate_item = paste0("x", 1:n_items, "_", suffix),
                                    p = NA)) %>%
               distinct(candidate_item, .keep_all = TRUE) %>%
               mutate(DIF = case_when(p < dif_alpha ~ 1,
                                      p >= dif_alpha ~ 0,
                                      is.na(p) ~ 0),
                      type = suffix, 
                      salient_diff_pct_03 = sum_(get(paste0(suffix, "_factor_scores"))$diff_gt_03) / length(get(paste0(suffix, "_factor_scores"))$diff_gt_03),
                      salient_diff_pct_1se = sum_(get(paste0(suffix, "_factor_scores"))$diff_gt_1se) / length(get(paste0(suffix, "_factor_scores"))$diff_gt_1se),
                      salient_diff_pct_03_rg = sum_(get(paste0(suffix, "_factor_scores"))$diff_gt_03[get(paste0(suffix, "_factor_scores"))$group == "Reference"]) / length(get(paste0(suffix, "_factor_scores"))$diff_gt_03[get(paste0(suffix, "_factor_scores"))$group == "Reference"]),
                      salient_diff_pct_1se_rg = sum_(get(paste0(suffix, "_factor_scores"))$diff_gt_1se[get(paste0(suffix, "_factor_scores"))$group == "Reference"]) / length(get(paste0(suffix, "_factor_scores"))$diff_gt_1se[get(paste0(suffix, "_factor_scores"))$group == "Reference"]),
                      salient_diff_pct_03_fg = sum_(get(paste0(suffix, "_factor_scores"))$diff_gt_03[get(paste0(suffix, "_factor_scores"))$group == "Focal"]) / length(get(paste0(suffix, "_factor_scores"))$diff_gt_03[get(paste0(suffix, "_factor_scores"))$group == "Focal"]),
                      salient_diff_pct_1se_fg = sum_(get(paste0(suffix, "_factor_scores"))$diff_gt_1se[get(paste0(suffix, "_factor_scores"))$group == "Focal"]) / length(get(paste0(suffix, "_factor_scores"))$diff_gt_1se[get(paste0(suffix, "_factor_scores"))$group == "Focal"])))
    } else {
      
      assign(paste0("dif_results_", suffix), dif_mod %>%
               mutate(candidate_item = rownames(.)) %>%
               bind_rows(data.frame(candidate_item = paste0("x", 1:n_items, "_", suffix),
                                    p = NA)) %>%
               distinct(candidate_item, .keep_all = TRUE) %>%
               mutate(DIF = case_when(p < dif_alpha ~ 1,
                                      p >= dif_alpha ~ 0,
                                      is.na(p) ~ 0),
                      type = suffix))
    }
  }
  
  dif_results <- rbindlist(mget(paste0("dif_results_", suffixes)), fill = TRUE) %>%
    bind_rows(dif_continuous) %>%
    mutate(candidate_item = str_split_i(candidate_item, "_", 1),
           sim = iter)
  
  
  write_to_log(iter, log_file)
  
  return(dif_results)
  
}



# Run Simulations ---------------------------------------------------------

log_file <- "sim_dif_log.txt"

test1 <- sim_dif_detection(iter = 31)

iter_start <- 1
nsims <- 20
iter_end <- iter_start + nsims - 1

test2 <- Map(sim_dif_detection, iter = iter_start:iter_end)


plan(multisession, workers = 4)

# Enable progress reporting
handlers(global = TRUE)



with_progress({
  p <- progressor(along = iter_start:iter_end)
  
  test2 <- future_Map(function(i) {
    p()  # Report progress
    sim_dif_detection(iter = i)
  }, 
  i = iter_start:iter_end,
  future.seed = TRUE)
})

plan(sequential)


# Compile Results ---------------------------------------------------------



res2 <- test2 %>% 
  rbindlist(fill = TRUE)
#rm(test2)
#saveRDS(res2, "Output/res2.Rds")
#res2 <- readRDS("Output/res2.Rds")

res2

res2 %>%
  filter(DIF == 1)

res2 %>%
  group_by(candidate_item, type) %>%
  summarise(sims = length(DIF),
            DIF_n = sum(DIF),
            .groups = "drop") %>%
  mutate(DIF_pct = round(100*DIF_n/sims, 3)) %>%
  pivot_wider(id_cols = c(candidate_item, sims), 
              names_from = type,
              values_from = c(DIF_n, DIF_pct)) %>%
  mutate(EF_vs_cont = DIF_pct_ef - DIF_pct_cn,
         EI_vs_cont = DIF_pct_ei - DIF_pct_cn,
         KM_vs_cont = DIF_pct_km - DIF_pct_cn) %>%
  arrange(parse_number(candidate_item)) %>%
  data.frame() 


res2 %>%
  group_by(candidate_item, type) %>%
  summarise(sims = length(DIF),
            DIF_n = sum(DIF),
            .groups = "drop") %>%
  mutate(DIF_pct = round(100*DIF_n/sims, 3),
         DIF_pct_se = sqrt((DIF_pct/100 * (1-DIF_pct/100))/sims),
         DIF_pct_lo = DIF_pct - 100*DIF_pct_se,
         DIF_pct_hi = DIF_pct + 100*DIF_pct_se,
         candidate_item = factor(candidate_item, levels = paste0("x", 1:max(parse_number(candidate_item))))) %>%
  ggplot(aes(x = type, y = DIF_pct, fill = type)) +
  geom_col() +
  geom_errorbar(aes(ymin = DIF_pct_lo, ymax = DIF_pct_hi)) +
  facet_wrap(~candidate_item, ncol = 5) +
  ylab("Percent of Simulations with DIF") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_discrete(name = "Type") +
  ylim(0, 100)


res2 %>%
  group_by(type, sim) %>%
  mutate(anyDIF = ifelse(any(DIF == 1), 1, 0)) %>%
  ungroup() %>%
  distinct(type, sim, .keep_all = TRUE) %>%
  group_by(type) %>%
  summarise(sims = n(),
            sims_w_DIF = sum(anyDIF == 1),
            DIF_pct = sims_w_DIF / sims,
            DIF_pct_se = sqrt((DIF_pct * (1-DIF_pct))/sims),
            DIF_pct = DIF_pct*100,
            DIF_pct_lo = DIF_pct - 100*DIF_pct_se,
            DIF_pct_hi = DIF_pct + 100*DIF_pct_se,
            .groups = "drop") %>%
  ggplot(aes(x = type, y = DIF_pct, fill = type)) +
  geom_col() +
  geom_errorbar(aes(ymin = DIF_pct_lo, ymax = DIF_pct_hi)) +
  ylab("Percent of Simulations with DIF") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_discrete(name = "Type") +
  ylim(0, 100)


res2 %>%
  group_by(type, sim) %>%
  mutate(anyDIF = ifelse(any(DIF == 1), 1, 0)) %>%
  ungroup() %>%
  distinct(type, sim, .keep_all = TRUE) %>%
  group_by(type) %>%
  summarise(sims = n(),
            sims_w_DIF = sum(anyDIF == 1),
            .groups = "drop") %>%
  mutate(DIF_pct = round(100*sims_w_DIF/sims, 3))

table(res2$candidate_item, res2$type, res2$DIF)

table(res2$candidate_item, res2$type, res2$DIF) %>%
  proportions(margin = c(1, 2))
