library(pacman)
p_load(dplyr, magrittr, ggplot2, psych, data.table, tidyr, janitor, hablar,
       mirt, faux)

# Use the dialog box to open ~/Dropbox/Projects/PsyMCA-2024-Simulation/data/HCAP_harmonized_data_PsyMCA-v4.rdata
load(file.choose())

original_data <- `HCAP_harmonized_data_PsyMCA-v4`

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
  filter(study %in% c(1, 5))

