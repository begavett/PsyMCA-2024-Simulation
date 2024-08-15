library(pacman)
p_load(dplyr, magrittr, ggplot2, psych, data.table, tidyr, janitor, hablar,
       mirt, faux, Hmisc, arules, discretization, synthpop)
source("~/Dropbox/Research/Code/recode.R")

# Use the dialog box to open ~/Dropbox/Projects/PsyMCA-2024-Simulation/data/HCAP_harmonized_data_PsyMCA-v4.rdata
load(file.choose())

original_data <- `HCAP_harmonized_data_PsyMCA-v5`
synthetic_data <- syn(original_data %>% select(-c(id_hrs, id_elsa, id_mexcog,
                                                  id_lasidad, id_haalsi, 
                                                  id_charls,
                                                  demographics,
                                                  COGNITIVE_FACTORS,
                                                  economics,
                                                  HEALTH_AND_FUNCTIONING,
                                                  CESD_ITEMS,
                                                  behaviors,
                                                  COGNITIVE_TEST_INDICATORS,
                                                  admin)))

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
  mutate(cohort = ifelse(study == 1, "HRS-HCAP", "HAALSI-HCAP"))

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

hcap_haalsi_memory %>%
  summarise(n = n(),
            m_mem = mean_(fmem),
            sd_mem = sd_(fmem))

hcap_haalsi_memory %>%
  select(cohort, lmi, lmd, lm_recog, cerad_cpd) %>%
  group_by(cohort) %>%
  summarise(across(everything(), .fns = list(sd = sd_, var = var_)))

hcap_haalsi_memory <- hcap_haalsi_memory %>%
  mutate(across(c(lmi, lmd, lm_recog, cerad_cpd), \(x) as.numeric(scale(x)), .names = "{.col}_z"))

mem_m <- '
mem =~ lmi_z + lmd_z + lm_recog_z + cerad_cpd_z
'

mem_f0 <- cfa(mem_m, 
             data = hcap_haalsi_memory, 
             #ordered = c("mmse_3wi", "mmse_3wd"),
             std.lv = FALSE, estimator = "mlr")

mem_f <- cfa(mem_m, 
             data = hcap_haalsi_memory, 
             #ordered = c("mmse_3wi", "mmse_3wd"),
             group = "cohort",
             std.lv = TRUE, estimator = "mlr")

summary(mem_f, fit.measures = TRUE, standardize = TRUE)
pop_partable <- partable(mem_f)

pop_partable <- pop_partable %>%
  mutate(est = case_when(
    # loadings
    # lhs == "mem" & rhs == "mmse_3wi" ~ 0.51,
    # lhs == "mem" & rhs == "mmse_3wd" ~ 0.76,
    lhs == "mem" & rhs == "lmi_z" ~ 0.71,
    lhs == "mem" & rhs == "lmd_z" ~ 0.74,
    lhs == "mem" & rhs == "lm_recog_z" ~ 0.62,
    lhs == "mem" & rhs == "cerad_cpd_z" ~ 0.70,
# 
#     # intercepts
#     lhs == "lmi" & op == "~1" ~ 9.83,
#     lhs == "lmd" & op == "~1" ~ 7.39,
#     lhs == "lm_recog" & op == "~1" ~ 10.28,
#     lhs == "cerad_cpd" & op == "~1" ~ 5.81,
#     
#     # variances
#     lhs == "lmi" & op == "~~" & rhs == "lmi" & group == 2 ~ 17.5,
#     lhs == "lmd" & op == "~~" & rhs == "lmi" & group == 2 ~ 17.7,
#     lhs == "lm_recog" & op == "~~" & rhs == "lmi" & group == 2 ~ 3.68,
#     lhs == "cerad_cpd" & op == "~~" & rhs == "lmi" & group == 2 ~ 7.23,
#     
#     lhs == "lmi" & op == "~~" & rhs == "lmi" & group == 1 ~ 26.1,
#     lhs == "lmd" & op == "~~" & rhs == "lmi" & group == 1 ~ 29.6,
#     lhs == "lm_recog" & op == "~~" & rhs == "lmi" & group == 1 ~ 7.97,
#     lhs == "cerad_cpd" & op == "~~" & rhs == "lmi" & group == 1 ~ 10.6,
#     
    # group 2 latent variable means & sds
    lhs == "mem" & op == "~1" & group == "2" ~ -0.899,
    lhs == "mem" & op == "~~" & rhs == "mem" ~ 0.569,
    TRUE ~ est),
    free = 0)

cfa(model = pop_partable,
    data = hcap_haalsi_memory, 
    #ordered = c("mmse_3wi", "mmse_3wd"),
    group = "cohort",
    std.lv = TRUE, estimator = "mlr",
    optim.force.converged = TRUE) %>%
  summary(standardize = TRUE)

sim_data <- simulateData(model = pop_partable,
             model.type = "cfa",
             estimator = "mlr",
             sample.nobs = c(pop_hyppar$n[pop_hyppar$study == 1], pop_hyppar$n[pop_hyppar$study == 5])) %>%
  mutate(cohort = case_when(group == 1 ~ "HRS_HCAP",
                            group == 2 ~ "HAALSI_HCAP"),
         lmi = ifelse(cohort == "HRS_HCAP",
                      lmi_z * sd_(hcap_haalsi_memory$lmi[hcap_haalsi_memory$study == 1]) + mean_(hcap_haalsi_memory$lmi[hcap_haalsi_memory$study == 1]),
                      lmi_z * sd_(hcap_haalsi_memory$lmi[hcap_haalsi_memory$study == 5]) + mean_(hcap_haalsi_memory$lmi[hcap_haalsi_memory$study == 5])) %>%
           round(),
         lmd = ifelse(cohort == "HRS_HCAP",
                      lmd_z * sd_(hcap_haalsi_memory$lmd[hcap_haalsi_memory$study == 1]) + mean_(hcap_haalsi_memory$lmd[hcap_haalsi_memory$study == 1]),
                      lmd_z * sd_(hcap_haalsi_memory$lmd[hcap_haalsi_memory$study == 5]) + mean_(hcap_haalsi_memory$lmd[hcap_haalsi_memory$study == 5])) %>%
           round(),
         lm_recog = ifelse(cohort == "HRS_HCAP",
                           lm_recog_z * sd_(hcap_haalsi_memory$lm_recog[hcap_haalsi_memory$study == 1]) + mean_(hcap_haalsi_memory$lm_recog[hcap_haalsi_memory$study == 1]),
                           lm_recog_z * sd_(hcap_haalsi_memory$lm_recog[hcap_haalsi_memory$study == 5]) + mean_(hcap_haalsi_memory$lm_recog[hcap_haalsi_memory$study == 5])) %>%
           round(),
         cerad_cpd = ifelse(cohort == "HRS_HCAP",
                            cerad_cpd_z * sd_(hcap_haalsi_memory$cerad_cpd[hcap_haalsi_memory$study == 1] + mean_(hcap_haalsi_memory$cerad_cpd[hcap_haalsi_memory$study == 1])),
                            cerad_cpd_z * sd_(hcap_haalsi_memory$cerad_cpd[hcap_haalsi_memory$study == 5]) + mean_(hcap_haalsi_memory$cerad_cpd[hcap_haalsi_memory$study == 5])) %>%
           round())
         


sim_data_long <- sim_data %>%
  pivot_longer(cols = lmi:cerad_cpd)


ggplot(sim_data_long, aes(x = value, fill = cohort)) +
  geom_density(position = "dodge", alpha = .5) +
  facet_wrap(~name, scales = "free")

hcap_haalsi_memory_long <- hcap_haalsi_memory %>%
  select(cohort, lmi, lmd, lm_recog, cerad_cpd) %>%
  pivot_longer(cols = lmi:cerad_cpd)


psych::describe(sim_data)

hcap_haalsi_memory %>%
  select(cohort, lmi, lmd, lm_recog, cerad_cpd) %>%
  psych::describe()


sim_data_r <- sim_data %>%
  mutate(across(c(lmi, lmd, lm_recog, cerad_cpd), \(x) as.numeric(cut2(x, g = 10, m = 10)), .names = "{.col}_ef"))

sim_data_r <- recodeOrdinal(sim_data_r, 
                            varlist_orig = c("lmi", "lmd", "lm_recog", "cerad_cpd"),
                            varlist_tr = c("lmi_ei", "lmd_ei", "lm_recog_ei", "cerad_cpd_ei"))

lapply(sim_data_r, table)

ggplot(sim_data_r, aes(x = lmi, y = lmi_ei)) +
  geom_point()

ggplot(sim_data_r, aes(x = lmd, y = lmd_ei)) +
  geom_point()

ggplot(sim_data_r, aes(x = lm_recog, y = lm_recog_ei)) +
  geom_point()

ggplot(sim_data_r, aes(x = cerad_cpd, y = cerad_cpd_ei)) +
  geom_point()



ggplot(sim_data_r, aes(x = lmi, y = lmi_ef)) +
  geom_point()

ggplot(sim_data_r, aes(x = lmd, y = lmd_ef)) +
  geom_point()

ggplot(sim_data_r, aes(x = lm_recog, y = lm_recog_ef)) +
  geom_point()

ggplot(sim_data_r, aes(x = cerad_cpd, y = cerad_cpd_ef)) +
  geom_point()

psych::describe(sim_data_r %>% select(ends_with("ef")))

lapply(sim_data_r %>% select(ends_with("ef")), table)

mem_ef_m <- '
mem =~ lmi_ef + lmd_ef + lm_recog_ef + cerad_cpd_ef
'


cfa(mem_ef_m, data = sim_data_r, ordered = TRUE, estimator = "wlsmv", std.lv = TRUE) %>%
  summary(fit.measures = TRUE, standardize = TRUE)


mirt(sim_data_r %>% select(ends_with("ef")), 
     model = 1, 
     itemtype = "graded") %>%
  coef(simplify = TRUE)

fit_mg <- multipleGroup(sim_data_r %>% select(ends_with("ef")), 
     model = 1,
     itemtype = "graded",
     group = sim_data_r$cohort,
     SE = TRUE,
     invariance = c("lmd_ef", "free_means", "free_var"),
     technical = list(NCYCLES = 50000))

coef(fit_mg, simplify = TRUE)

DIF(fit_mg, 
    which.par = c("a1", paste0("d", 1:9)),
    scheme = "add_sequential")
