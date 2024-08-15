library(pacman)
p_load(dplyr, magrittr, ggplot2, psych, data.table, tidyr, janitor, hablar,
       mirt, faux, Hmisc, arules, discretization, synthpop, lavaan, stringr)
source("~/Dropbox/Research/Code/recode.R")

# Use the dialog box to open ~/Dropbox/Projects/PsyMCA-2024-Simulation/data/HCAP_harmonized_data_PsyMCA-v4.rdata
#load(file.choose())

use_items <- c("lmi_ef", "lmd_ef", "lm_recog_ef", "cerad_cpd_ef")

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
         currsmoke,
         mweight,
         mheight,
         cesd1:cesd8, # Threshold = 4
         t2diab,
         drinksperwk,
         vision_dis,
         vision_near,
         fqmdactx,
         fqvgactx) %>%
  mutate(bmi = mweight/mheight^2) %>%
  filter(study %in% c(1, 5)) %>%
  mutate(cohort = ifelse(study == 1, "HRS-HCAP", "HAALSI-HCAP"),
         depression = case_when(is.na(cesd1) & is.na(cesd2) & is.na(cesd3) &
                                  is.na(cesd4) & is.na(cesd5) & is.na(cesd6) &
                                  is.na(cesd7) & is.na(cesd8) ~ NA_real_,
                                rowSums(across(cesd1:cesd8), na.rm = T) >= 4 ~ 1,
                                rowSums(across(cesd1:cesd8), na.rm = T) < 4 ~ 0))

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

sim_data <- simulateData(model = pop_mem_pars, 
                         model.type = "sem", 
                         std.lv = TRUE,
                         estimator = "mlr",
                         sample.nobs = c(pop_hyppar$n[pop_hyppar$study == 1], 
                                         pop_hyppar$n[pop_hyppar$study == 5])) %>%
  mutate(across(lmi:cerad_cpd, round),
         cohort = ifelse(group == 1, "HRS", "HAALSI"))


sim_data_long <- sim_data %>%
  pivot_longer(cols = lmi:depression)


ggplot(sim_data_long, aes(x = value, fill = cohort)) +
  geom_density(position = "dodge", alpha = .5) +
  facet_wrap(~name, scales = "free")

hcap_haalsi_memory_long <- hcap_haalsi_memory %>%
  select(cohort, lmi, lmd, lm_recog, cerad_cpd, hyper, depression) %>%
  pivot_longer(cols = lmi:depression)


psych::describe(sim_data)

hcap_haalsi_memory %>%
  select(cohort, lmi, lmd, lm_recog, cerad_cpd, hyper, depression) %>%
  psych::describe()


sim_data_r <- sim_data %>%
  mutate(across(c(lmi, lmd, lm_recog, cerad_cpd), \(x) as.numeric(cut2(x, g = 10, m = 5)), .names = "{.col}_ef"))

sim_data_r <- recodeOrdinal(sim_data_r, 
                            varlist_orig = c("lmi", "lmd", "lm_recog", "cerad_cpd"),
                            varlist_tr = c("lmi_ei", "lmd_ei", "lm_recog_ei", "cerad_cpd_ei"))

lapply(sim_data_r, table)

ggplot(sim_data_r, aes(x = lmi, y = lmi_ei, colour = cohort)) +
  geom_jitter()

ggplot(sim_data_r, aes(x = lmd, y = lmd_ei, colour = cohort)) +
  geom_jitter()

ggplot(sim_data_r, aes(x = lm_recog, y = lm_recog_ei, colour = cohort)) +
  geom_jitter()

ggplot(sim_data_r, aes(x = cerad_cpd, y = cerad_cpd_ei, colour = cohort)) +
  geom_jitter()



ggplot(sim_data_r, aes(x = lmi, y = lmi_ef, colour = cohort)) +
  geom_jitter()

ggplot(sim_data_r, aes(x = lmd, y = lmd_ef, colour = cohort)) +
  geom_jitter()

ggplot(sim_data_r, aes(x = lm_recog, y = lm_recog_ef, colour = cohort)) +
  geom_jitter()

ggplot(sim_data_r, aes(x = cerad_cpd, y = cerad_cpd_ef, colour = cohort)) +
  geom_jitter()

psych::describe(sim_data_r %>% select(ends_with("ef")))

lapply(sim_data_r %>% select(ends_with("ef")), table)

mem_ef_m <- '
mem =~ lmi_ef + lmd_ef + lm_recog_ef + cerad_cpd_ef

mem ~ hyper + depression
'


sem(mem_ef_m, data = sim_data_r, ordered = TRUE, estimator = "wlsmv", std.lv = TRUE) %>%
  summary(fit.measures = TRUE, standardize = TRUE)


mirt(sim_data_r %>% select(ends_with("ef")), 
     model = 1, 
     itemtype = "graded") %>%
  coef(simplify = TRUE)

linking_item_number <- sample(1:4, 1)
linking_item <- c("lmi_ef", "lmd_ef", "lm_recog_ef", "cerad_cpd_ef")[linking_item_number]

fit_hrs <- mirt(sim_data_r %>% filter(cohort == "HRS") %>% select(ends_with("ef")),
                model = 1,
                itemtype = "graded",
                SE = TRUE,
                technical = list(NCYCLES = 50000))

bank_hrs <- mod2values(fit_hrs)

sim_data_r <- sim_data_r %>%
  mutate(study = ifelse(group == 1, "A_HRS", "B_HAALSI"))

fit_haalsi <- mirt(sim_data_r %>% filter(cohort == "HAALSI") %>% select(ends_with("ef")),
                model = 1,
                itemtype = "graded",
                SE = TRUE,
                technical = list(NCYCLES = 50000))

bank_haalsi <- mod2values(fit_haalsi)

fit_mg <- multipleGroup(sim_data_r %>% select(ends_with("ef")), 
                        model = 1,
                        itemtype = "graded",
                        group = sim_data_r$study,
                        SE = TRUE,
                        invariance = c(linking_item, "free_means", "free_var"),
                        technical = list(NCYCLES = 50000))

coef(fit_mg, simplify = TRUE)

dif_res <- DIF(fit_mg, 
               which.par = c("a1", paste0("d", 1:9)),
               scheme = "add_sequential",
               technical = list(NCYCLES = 50000))

dif_items <- rownames(dif_res)


anchor_items <- setdiff(c("lmi_ef", "lmd_ef", "lm_recog_ef", "cerad_cpd_ef"), dif_items)

## Fix banked item parameters except those with DIF (freely estimated separately for each cohort)

# Get parameter table from bifactor model
hz_mg_pars <- mod2values(fit_mg)

# First, constrain all HRS item parameters to their banked values
for(i in use_items) {
  hz_mg_pars[hz_mg_pars$item == i, "value" ] <- rep(bank_hrs[bank_hrs$item == i, "value"], 2)
}

# Next, for the linking items with DIF, give them starting values from their respective banks
# This will overwrite above making sure that linking items with DIF get different starting values by cohort
for(i in setdiff(use_items, anchor_items)) {
  print(i)
  hz_mg_pars[hz_mg_pars$item == i & hz_mg_pars$group == "A_HRS", "value"] <- bank_hrs[bank_hrs$item == i, "value"]
  hz_mg_pars[hz_mg_pars$item == i & hz_mg_pars$group == "B_HAALSI", "value"] <- bank_haalsi[bank_haalsi$item == i, "value"]
}

hz_mg_pars$est <- FALSE

# Freely estimate parameters for items with DIF
hz_mg_pars <- hz_mg_pars %>%
  mutate(est = ifelse(item %in% setdiff(use_items, anchor_items) & value != 0, TRUE, est))

# Finally, make sure that the mean and variance of the global factor in HAALSI are free 
hz_mg_pars <- hz_mg_pars %>%
  mutate(est = ifelse(group == "B_HAALSI" & name %in% c("MEAN_1", "COV_11"), TRUE, est))

# Check item parameters to see that linking items with DIF are free
# and DIF-free linking items are fixed
hz_mg_pars %>%
  filter(item %in% use_items,
         name %in% c("a1", paste0("d", 1:9))) %>%
  arrange(item, group, name)

# This should estimate item parameters for linking items with DIF (in both groups separately) and VIP hyperparameters
hz_model <- multipleGroup(sim_data_r %>% select(ends_with("ef")),
                 model = 1,
                 group = sim_data_r$study,
                 itemtype = "graded",
                 pars = hz_mg_pars,
                 technical = list(NCYCLES = 50000),
                 SE = TRUE)

coef(hz_model, simplify = TRUE)

## Fix all parameters to their estimates from above and generate factor scores

hz_model_pars <- mod2values(hz_model)
hz_model_pars$est <- FALSE
hz_model1 <- multipleGroup(sim_data_r %>% select(ends_with("ef")),
                           model = 1,
                           group = sim_data_r$study,
                           itemtype = "graded",
                           pars = hz_model_pars,
                           technical = list(NCYCLES = 50000))

hz_model_fscores <- fscores(hz_model1, full.scores.SE = TRUE, QMC = TRUE) %>%
  as_tibble()

hcap_haalsi_memory_fs <- hcap_haalsi_memory %>%
  bind_cols(hz_model_fscores)


ggplot(hcap_haalsi_memory_fs, aes(x = fmem, y = F1, colour = cohort)) +
  geom_point() +
  geom_smooth(method = "lm")


lm(F1 ~ depression + hyper, data = hcap_haalsi_memory_fs) %>% 
  summary()
