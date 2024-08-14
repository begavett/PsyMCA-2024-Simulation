#this is a test

library(haven)
library(tidyverse)

load("/Users/ryanandrews/BOSTON UNIVERSITY Dropbox/Ryan Andrews/PsyMCA2024_Data_Sharing_HCAP/HCAP_harmonized_data_PsyMCA-v6.rdata")
hrs_core <- haven::read_dta("/Users/ryanandrews/Dropbox (BOSTON UNIVERSITY)/PsyMCA2024_Data_Sharing_HCAP/interim-HcapHarmoniz-308-ScoresIncludingHRS2016Core.dta")


## IN 2016 when HCAP started, create a dataset based on that year that has HRS TICS variables
#for all subjects and HCAP variables for HCAP subjects.

## Then, create a factor model (1-factor model) that would be based upon what they did with
#HCAP, but would include the HRS variables. 

## You could filter out people who don't have HCAP

analytic <- hrs_core %>% filter(study ==1 | study ==10) %>%
  mutate(orientation_hrshcap = u105 + u106 + u107 + u108,
         orientation_hcap = u113 + u114 + u115 + u116 + u117 + u118,
         #502 and 504 TICS separately (available in both)
         naming_tics = u502 + u504,
         naming_tics_hcap = u505 + u506 + u508,
         naming_tics_hcap = case_when(
           naming_tics_hcap <=2 ~ 1,
           TRUE ~ 2
         )
  )

varlist_orig <- c("u201", "u203", "u207", "u213", "u215", "u301", 
                  "u302", "u303", "u401", "u402", "u403", "u414", 
                  "u501", "u601", "u216")
varlist_tr <- c("u201r", "u203r", "u207r", "u213r", "u215r", "u301r", 
                "u302r", "u303r", "u401r", "u402r", "u403r", "u414r", 
                "u501r", "u601r", "u216r")
#source()
analytic2 <- as.data.frame(analytic)
recoded_analytic = recodeOrdinal(df = analytic2,
              varlist_orig = varlist_orig,
              varlist_tr = varlist_tr,
              type = "interval",
              ncat = 15)


item_names <- c(varlist_tr,
                "naming_tics_hcap",
                "naming_tics",
                "orientation_hcap",
                "orientation_hrshcap",
                "u209",
                "u258",
                "u409",
                "u515",
                "u517",
                "u518",
                "u602",
                "u701",
                "u705",
                "u410",
                "u412")

library(mirt)

model1f <- "cog=1-30"
fit1f <- mirt(data = recoded_analytic %>% select(all_of(item_names)),
              model = model1f)

summary(fit1f)
fit_ind1f <- M2(fit1f,calcNull = TRUE, CI = 0.90, type = "C2", 
                QMC=TRUE, na.rm = TRUE)

lavmod1f <- "cog=~u201r + u203r + u207r + u209 + u213r + u215r + u258 + u301r + u302r + u303r + u401r + u402r + u403r + u409 + u414r + u501r + orientation_hrshcap + orientation_hcap + naming_tics + naming_tics_hcap + u515 + u517 + u518 + u601r + u602 + u701 + u705 + u410 + u216r + u412"
fit1f_lav <- lavaan::cfa(model = lavmod1f,
                         data = recoded_analytic)

hist(foo$u201)
hist(foo$u201r)
hist(foo$u203)
hist(foo$u203r)
#u105 - u108 are orientation items




# hrshcap <- `HCAP_harmonized_data_PsyMCA-v6` %>% filter(study==1)
# hrs_core2016 <- hrs_core %>% filter(study==1)
