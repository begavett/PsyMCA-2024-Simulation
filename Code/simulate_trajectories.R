require(lme4)
require(mirt)
require(tidyverse)
require(ggplot2)

# show size of environment objects
# sort( sapply(ls(),function(x){object.size(get(x))})) 

# rm(list = ls()[grepl("^fit",ls())])



#   extractMIRTParm is a function that extracts item parameters from a mirt
#       estimate multidimensional (uncorrelated dimensions) model and converts
#       discrimination (a) parameters to a vector of item parameters and 
#       threshhold (d) parameters to a matrix of d values.
#   Parameters
#       model - estimated mirt model object
#       max_thresh - maximum number of threshold paramters across items
#   Value list with 2 elements
#       a - vector of a parameters, 1 for each item
#       d - matrix of d parameters, rows correspons to items, columns to threshholds

extractMIRTParm <- function(model,max_thresh = 9) {
    require(mirt)
    
    parm <- mod2values(model)
    
    parmw <- parm %>% dplyr::select(item, name, value) %>% 
        pivot_wider(id_cols = item, names_from = name, values_from = value) %>% 
        filter(!item == "GROUP") %>% mutate(
            d1 = case_when(
                is.na(d1) ~ d,
                TRUE ~ d1
            )
        )
    
    a <- as.matrix(parmw %>% dplyr::select(a1))
    d <- as.matrix(parmw %>% dplyr::select(paste0("d",1:max_thresh)))
    return(list(a,d))
}


#   simTraj is a function that 1) estimates a mirt model using a simulated
#       item response dataset, 2) calculates factor scores based on the estimated 
#       model, 3) estimates linear mixed effects models with true cognition and  
#       simulated cognition (factor scores from simulated item responses) as dependent 
#       variables, time in study as a fixed effect variable, and person and 
#       person-by-time random effects, and 4) 4) returns person specific estimates 
#       (random effects) of cognitive components slopes and intercepts.
#   Parameters
#       data - data frame that contains true cognition value, simulated item responses, 
#           id variable, and time in study variable
#       mirt_mod - mirt model specification
#   Value - list with 5 elements:
#       data frame with estimated random effects for mirt_mod applied to data. 
#           Includes intercept and slope random effects for true cognition (_true) 
#           and simulated cognition factor scores (_sim)

# simTraj <- function(data = df, item_labels = item_labels, mirt_mod = mirt_all,
#     iteration = iteration) {
simTraj <- function(data = df, sim_cog_var = "sim_cog_all", frml = frml,
    iteration = iteration) {
    require(mirt)
    require(tidyverse)
    require(lme4)
    
    # # estimate mirt model using simulated responses and generate factor scores
    re_empty <- data.frame(id = NA, int_true = NA, slope_true = NA,
        int_true_se = NA,slope_true_se = NA, int_sim = NA, slope_sim = NA,
        int_sim_se = NA, slope_sim_se = NA)
    # 
    # tryCatch({
    #     mod <- mirt(data[,item_labels], mirt_mod)
    #     data$sim_cog <- fscores(mod)
    #     data <- data %>% relocate(sim_cog, .after = true_cog)
    # }, error=function(e){return(re_empty)})
    # data$itertion <- iteration
    
    # frml <- "true_cog ~ time + agebl_75 + slope_true_qrtl +
    # time:agebl_75 + time:slope_true_qrtl + (1 + time | id)"
    # frml <- "true_cog ~ time + agebl_75 + 
    # time:agebl_75 + (1 | id)"
    # 
    tryCatch({
        # mixed effects longitudinal model - true cognition
        # res_long_true <- lmer(true_cog ~ time + (1 + time | id), data = data)
        res_long_true <- lmer(formula = frml, data = data)
        # summary(res_long_true)
        data$cog_true_pred <- predict(res_long_true, newdata = data, type = "response")
        re_true <- data.frame(ranef(res_long_true))
        re_true <- re_true %>% 
            pivot_wider(id_cols = grp, names_from = term,values_from = c(condval,condsd))
        names(re_true) <- c("id","int_true","slope_true","int_true_se","slope_true_se")
        re_true$id <- as.integer(levels(re_true$id))[re_true$id]
        fe_true <- data.frame(t(fixef(res_long_true)))
        names(fe_true) <- paste0(names(fe_true),"_true")
        
        
        data$sim_cog <- data[,sim_cog_var]
        #res_long_sim <- lmer(sim_cog ~ time + (1 + time | id), data = data)
        frml <- sub("true_cog","sim_cog",frml)
        res_long_sim <- lmer(formula = frml, data = data)
        # summary(res_long_sim)
        
        data$cog_sim_pred <- predict(res_long_sim, newdata = data, type = "response")
        data <- data %>% dplyr::select(-sim_cog)
        re_sim <- data.frame(ranef(res_long_sim))
        re_sim <- re_sim %>% 
            pivot_wider(id_cols = grp, names_from = term,values_from = c(condval,condsd))
        names(re_sim) <- c("id","int_sim","slope_sim","int_sim_se","slope_sim_se")
        re_sim$id <- as.integer(levels(re_sim$id))[re_sim$id] 
        fe_sim <- data.frame(t(fixef(res_long_sim)))
        names(fe_sim) <- paste0(names(fe_sim),"_sim")
        
        re <- re_true %>% left_join(re_sim, by="id")
        fe <- bind_cols(fe_true,fe_sim)
        
    }, error=function(e){re <- re_empty})
    re$iteration <- iteration
    
    return(list(data,re,fe))
    
}

#   simulateTrajectories is a function that 1) simulates item responses for 
#   true cognition vectors for different components of the HRS-HCAP cognitive test
#   battery (MMSE, HRS TICS, HCAP, MMSE+TICS+HCAP), and 2) generates simulated observed (factor)
#   scores for each of the cognitive components (simTraj()), 3) estimates linear mixed effects models
#   with true cognition and simulated cognition (factor scores from simulated item responses)
#   as dependent variables, time in study as a fixed effect variable, and person and 
#   person-by-time random effects (simTraj()), 4) collates person specific estimates (random effects)
#   of cognitive components slopes and intercepts.
# 
#   Parameters
#       true_sim - data frame with true cognition values, rows correspond to assessments,
#           columns correspond to simulation iterations.
#       a_par - vector/data frame containing a (discrimination) item parameters 
#           from mirt co-calibration model
#       d_par - vector/data frame containing d (threshhold) item parameters 
#           from mirt co-calibration model
#       item_labels - item labels
#       niters - number of simulation iterations
#
#   Value - list with 5 elements:
#       ds - list of datasets for each niter simulations. Each dataset contains
#           the true cognition value (true_cog) and the simulated item responses conditional
#           on true_cog
#       re_all - estimated random effects for the full HRS-HCAP cognitive test battery
#           (MMSE+TICS+HCAP). Includes intercept and slope random effects for true
#           cognition (_true) and simulated cognition factor scores (_sim)
#       re_mmse - estimated random effects for the MMSE component of the HRS-HCAP 
#           cognitive test battery. Includes intercept and slope random effects for true
#           cognition (_true) and simulated cognition factor scores (_sim)
#       re_tics - estimated random effects for the TICS component of the HRS-HCAP 
#           cognitive test battery. Includes intercept and slope random effects for true
#           cognition (_true) and simulated cognition factor scores (_sim)
#       re_hcap - estimated random effects for the HCAP component of the HRS-HCAP 
#           cognitive test battery. Includes intercept and slope random effects for true
#           cognition (_true) and simulated cognition factor scores (_sim)

simulateTrajectories <- function(true_sim = s4, a_par = a, d_par = d, item_labels, 
    niter = 100, iter_group = 10, out_dir = "Analysis/Simulation Results/") {
    require(mirt)
    require(tidyverse)
    for (itgrp in 6:iter_group) {
        ds <- list()
        re_all <- list()
        re_mmse <- list()
        re_tics <- list()
        re_hcap <- list()
        fe_all <- list()
        fe_mmse <- list()
        fe_tics <- list()
        fe_hcap <- list()
        set.seed(092724)
        for (iter in 1:niter) {
            cat(paste0("Group - ",itgrp, ", Iteration - ",iter,"\n"))
            
            # simulate item responses
            iteration <- ((itgrp - 1) * niter) + iter
            df <- data.frame(simdata(a = a_par, d = d_par, Theta = true_sim[,paste0("vm",iteration)], 
                itemtype = 'graded'))
            names(df) <- item_labels
            df$id <- true_sim$id
            df$time <- true_sim$time
            df$agebl_75 <- true_sim$agebl_75
            df$true_cog <- true_sim[,paste0("vm",iteration)]
            df$iteration <- iteration
            df <- df %>% relocate(c(id,iteration,time,true_cog))
            
            true_re <- readRDS("~/Psychometrics Conference/2024/Simulation WG/PsyMCA-2024-Simulation/Data/true_random_effects.rds")
            df <- df %>% left_join((true_re %>% dplyr::select(id,int_true_qrtl,
                slope_true_qrtl)),by="id")
            
            # ds[[iteration]] <- df
            
            flag_low_response <- FALSE
            for (itnm in item_labels) {
                if (length(table(df[,itnm])) < 2) {
                    flag_low_response <- TRUE
                }
            }
            if (flag_low_response == TRUE) {
                iter = iter - 1
                next
            }
            
            # mirt_models
            
            mirt_mmse <- "mmse = 1-11"
            mirt_tics <- "tics = 1-8"
            mirt_hcap <- "hcap = 1-18"
            mirt_all <- "cog = 1-37"
            
            ### estimate mirt model using simulated responses and generate factor scores
            re_empty <- data.frame(id = NA, int_true = NA, slope_true = NA, 
                int_true_se = NA,slope_true_se = NA, int_sim = NA, slope_sim = NA, 
                int_sim_se = NA, slope_sim_se = NA)
            
            tryCatch({
                mod <- mirt(df[,item_labels], mirt_all)
                df$sim_cog_all <- fscores(mod)
                mod <- mirt(df[,item_labels[1:11]], mirt_mmse)
                df$sim_cog_mmse <- fscores(mod)
                mod <- mirt(df[,item_labels[12:19]], mirt_tics)
                df$sim_cog_tics <- fscores(mod)
                mod <- mirt(df[,item_labels[20:37]], mirt_hcap)
                df$sim_cog_hcap <- fscores(mod)
                df <- df %>% 
                    relocate(c(sim_cog_all,sim_cog_mmse,sim_cog_tics,sim_cog_hcap), 
                        .after = true_cog)
            }, error=function(e){return(re_empty)})
            df$iteration <- iteration
            
            # rescale factor scores to equate metric to true cognition
            df <- df %>% mutate(
                sim_cog_mmse_rs = (((sim_cog_mmse - mean(sim_cog_mmse)) / 
                    sd(sim_cog_mmse)) * sd(true_cog)) + mean(true_cog),
                sim_cog_tics_rs = (((sim_cog_tics - mean(sim_cog_tics)) / 
                    sd(sim_cog_tics)) * sd(true_cog)) + mean(true_cog),
                sim_cog_hcap_rs = (((sim_cog_hcap - mean(sim_cog_hcap)) / 
                    sd(sim_cog_hcap)) * sd(true_cog)) + mean(true_cog),
                sim_cog_all_rs = (((sim_cog_all - mean(sim_cog_all)) / 
                    sd(sim_cog_all)) * sd(true_cog)) + mean(true_cog),
            )
            
            # res <- lmer(formula = frml,data = df)
            # summary(res)
            # fixef(res)
            # 
            
            
            ### calculate simulated random effects
            
            # frml <- as.formula("true_cog ~ time + (1 | id)")
            frml <- "true_cog ~ time + agebl_75 + slope_true_qrtl +
            time:agebl_75 + time:slope_true_qrtl + (1 + time | id)"
            # re <- simTraj(data = df, item_labels = item_labels, mirt_mod = mirt_all,
            #     iteration = iteration)
            res <- simTraj(data = df, sim_cog_var = "sim_cog_all_rs",frml = frml,
                iteration = iteration)
            df <- res[[1]]
            df <- df %>% rename(cog_true_pred_all = cog_true_pred,
                cog_sim_pred_all = cog_sim_pred)
            re <- res[[2]]
            re$label <- "MMSE+TICS+HCAP"
            re_all[[paste0("iteration-",iteration)]] <- re
            fe <- res[[3]]
            fe$label <- "MMSE+TICS+HCAP"
            fe_all[[paste0("iteration-",iteration)]] <- fe
            
            # re <- simTraj(data = df, item_labels = item_labels[1:11], mirt_mod = mirt_mmse,
            #     iteration = iteration)
            res <- simTraj(data = df, sim_cog_var = "sim_cog_mmse_rs",frml = frml,
                iteration = iteration)
            df <- res[[1]]
            df <- df %>% rename(cog_true_pred_mmse = cog_true_pred,
                 cog_sim_pred_mmse = cog_sim_pred)
            re <- res[[2]]
            re$label <- "MMSE"
            re_mmse[[paste0("iteration-",iteration)]] <- re
            fe <- res[[3]]
            fe$label <- "MMSE"
            fe_mmse[[paste0("iteration-",iteration)]] <- fe
            
            # re <- simTraj(data = df, item_labels = item_labels[12:19], mirt_mod = mirt_tics,
            #     iteration = iteration)
            res <- simTraj(data = df, sim_cog_var = "sim_cog_tics_rs",frml = frml,
                iteration = iteration)
            df <- res[[1]]
            df <- df %>% rename(cog_true_pred_tics = cog_true_pred,
                cog_sim_pred_tics = cog_sim_pred)
            re <- res[[2]]
            re$label <- "TICS"
            re_tics[[paste0("iteration-",iteration)]] <- re
            fe <- res[[3]]
            fe$label <- "TICS"
            fe_tics[[paste0("iteration-",iteration)]] <- fe
            
            # re <- simTraj(data = df, item_labels = item_labels[20:37], mirt_mod = mirt_hcap,
            #     iteration = iteration)
            res <- simTraj(data = df, sim_cog_var = "sim_cog_hcap_rs",frml = frml,
                iteration = iteration)
            df <- res[[1]]
            df <- df %>% rename(cog_true_pred_hcap = cog_true_pred,
                cog_sim_pred_hcap = cog_sim_pred)
            re <- res[[2]]
            re$label <- "HCAP"
            re_hcap[[paste0("iteration-",iteration)]] <- re
            fe <- res[[3]]
            fe$label <- "HCAP"
            fe_hcap[[paste0("iteration-",iteration)]] <- fe
            
            ds[[paste0("iteration-",iteration)]] <- df
            
            
        } # end for iter
        
        saveRDS(ds, file = paste0(out_dir,"ds_iteration_group_", itgrp, ".rds"))
        saveRDS(re_mmse, file = paste0(out_dir,"re_mmse_iteration_group_", itgrp, ".rds"))
        saveRDS(re_tics, file = paste0(out_dir,"re_tics_iteration_group_", itgrp, ".rds"))
        saveRDS(re_hcap, file = paste0(out_dir,"re_hcap_iteration_group_", itgrp, ".rds"))
        saveRDS(re_all, file = paste0(out_dir,"re_all_iteration_group_", itgrp, ".rds"))
        saveRDS(fe_mmse, file = paste0(out_dir,"fe_mmse_iteration_group_", itgrp, ".rds"))
        saveRDS(fe_tics, file = paste0(out_dir,"fe_tics_iteration_group_", itgrp, ".rds"))
        saveRDS(fe_hcap, file = paste0(out_dir,"fe_hcap_iteration_group_", itgrp, ".rds"))
        saveRDS(fe_all, file = paste0(out_dir,"fe_all_iteration_group_", itgrp, ".rds"))
        
    } # end for itgrp
    
    for (itgrp in 1:iter_group) {
        if (itgrp == 1) {
            re_mmse <- readRDS(paste0(out_dir,"re_mmse_iteration_group_", itgrp, ".rds"))
            re_mmse <- re_mmse %>% bind_rows()
            re_tics <- readRDS(paste0(out_dir,"re_tics_iteration_group_", itgrp, ".rds"))
            re_tics <- re_tics %>% bind_rows()
            re_hcap <- readRDS(paste0(out_dir,"re_hcap_iteration_group_", itgrp, ".rds"))
            re_hcap <- re_hcap %>% bind_rows()
            re_all <- readRDS(paste0(out_dir,"re_all_iteration_group_", itgrp, ".rds"))
            re_all <- re_all %>% bind_rows()
            
            fe_mmse <- readRDS(paste0(out_dir,"fe_mmse_iteration_group_", itgrp, ".rds"))
            fe_mmse <- fe_mmse %>% bind_rows()
            fe_tics <- readRDS(paste0(out_dir,"fe_tics_iteration_group_", itgrp, ".rds"))
            fe_mmse <- fe_mmse %>% bind_rows()
            fe_hcap <- readRDS(paste0(out_dir,"fe_hcap_iteration_group_", itgrp, ".rds"))
            fe_mmse <- fe_mmse %>% bind_rows()
            fe_all <- readRDS(paste0(out_dir,"fe_all_iteration_group_", itgrp, ".rds"))
            fe_mmse <- fe_mmse %>% bind_rows()
            
        } else {
            re_mmse_temp <- readRDS(paste0(out_dir,"re_mmse_iteration_group_", itgrp, ".rds"))
            re_mmse_temp <- re_mmse_temp %>% bind_rows()
            re_mmse <- bind_rows(re_mmse, re_mmse_temp)
            re_tics_temp <- readRDS(paste0(out_dir,"re_tics_iteration_group_", itgrp, ".rds"))
            re_tics_temp <- re_tics_temp %>% bind_rows()
            re_tics <- bind_rows(re_tics, re_tics_temp)
            re_hcap_temp <- readRDS(paste0(out_dir,"re_hcap_iteration_group_", itgrp, ".rds"))
            re_hcap_temp <- re_hcap_temp %>% bind_rows()
            re_hcap <- bind_rows(re_hcap, re_hcap_temp)
            re_all_temp <- readRDS(paste0(out_dir,"re_all_iteration_group_", itgrp, ".rds"))
            re_all_temp <- re_all_temp %>% bind_rows()
            re_all <- bind_rows(re_all, re_all_temp)
            
            fe_mmse_temp <- readRDS(paste0(out_dir,"fe_mmse_iteration_group_", itgrp, ".rds"))
            fe_mmse_temp <- fe_mmse_temp %>% bind_rows()
            fe_mmse <- bind_rows(fe_mmse, fe_mmse_temp)
            fe_tics_temp <- readRDS(paste0(out_dir,"fe_tics_iteration_group_", itgrp, ".rds"))
            fe_tics_temp <- fe_tics_temp %>% bind_rows()
            fe_tics <- bind_rows(fe_tics, fe_tics_temp)
            fe_hcap_temp <- readRDS(paste0(out_dir,"fe_hcap_iteration_group_", itgrp, ".rds"))
            fe_hcap_temp <- fe_hcap_temp %>% bind_rows()
            fe_hcap <- bind_rows(fe_hcap, fe_hcap_temp)
            fe_all_temp <- readRDS(paste0(out_dir,"fe_all_iteration_group_", itgrp, ".rds"))
            fe_all_temp <- fe_all_temp %>% bind_rows()
            fe_all <- bind_rows(fe_all, fe_all_temp)
            
        }
    }
    saveRDS(re_mmse, file = paste0(out_dir,"re_mmse", ".rds"))
    saveRDS(re_tics, file = paste0(out_dir,"re_tics", ".rds"))
    saveRDS(re_hcap, file = paste0(out_dir,"re_hcap", ".rds"))
    saveRDS(re_all, file = paste0(out_dir,"re_all", ".rds"))

    saveRDS(fe_mmse, file = paste0(out_dir,"fe_mmse", ".rds"))
    saveRDS(fe_tics, file = paste0(out_dir,"fe_tics", ".rds"))
    saveRDS(fe_hcap, file = paste0(out_dir,"fe_hcap", ".rds"))
    saveRDS(fe_all, file = paste0(out_dir,"fe_all", ".rds"))
    
    
    # return(list("ds" = ds, "re_all" = re_all, "re_mmse" = re_mmse,
    #     "re_tics" = re_tics, "re_hcap" = re_hcap))
}


#   ************************* Run/Test Simulation Code *************************

#   Code in this segment can be uncommented to run and test the functions in this file.
#   This code extracts longitudinal "true" cognition values from
#   "simulated_longitudinal_true_cognition.rds". This file is a synthetic dataset
#   with 792 individuals uniformly distributed in age for 60-95, each with baseline 
#   assessment and 10 follow-up assessments at 1 year intervals. Parameters to build 
#   this dataset came from an empirical analysis of a large dataset representative 
#   of older adults. This dataset is a long-form dataset that includes a person id
#   variable (id) and time-in-study variable (time) that are used in longitudinal
#   mixed effects model generated in the simTraj() function. It contains 4000
#   simulations of true cognition for the longitudinal data srtructure (labelled
#   vm1-vm4000). The simulation of true cognition is described in
#   longitudinal_cognition_simulation.<html/Rmd>.

# # load simulated "true" cognition
# true_sim <- readRDS("~/Psychometrics Conference/2024/Simulation WG/PsyMCA-2024-Simulation/Data/simulated_longitudinal_true_cognition.rds")
# 
# # Item names for analysis
# mmse_nmsr <- c("ornttime_mmser","orntplace_mmser","imrec_mmser","delrec_mmser",
#                "spellbck_mmser","naming_mmser","phrase_mmser","commfoll_mmser","writesent_mmser",
#                "draw_mmser","comm3step_mmser")
# 
# tics_nmsr <- c("wlimrc_ticsr","wldlrc_ticsr","count_ticsr","animfl_ticsr",
#                "numser_ticsr","serial7_ticsr","ornttime_ticsr","naming_ticsr")
# 
# hcap_nmsr <- c("naming_hcapr","ceradtot_hcapr","ceraddel_hcapr","ceradrcg_hcapr",
#                "animfl_hcapr","bmim_hcapr","bmdel_hcapr","lmim_hcapr","lmdel_hcapr","lmrcg_hcapr",
#                "conprxim_hcapr","conprxdel_hcapr","dsmt_hcapr","numser_hcapr","ravens_hcapr",
#                "trailsa_hcapir","trailsb_hcapir","spatial_hcapr")
# item_labels <- c(mmse_nmsr,tics_nmsr,hcap_nmsr)

# load("~/Psychometrics Conference/2024/Simulation WG/PsyMCA-2024-Simulation/Analysis/co_calibration_results.Rdata")
# 
# a <- extractMIRTParm((res_md_1a))[[1]]
# d <- extractMIRTParm((res_md_1a))[[2]]
# 
# sim20 <- simulateTrajectories(true_sim = true_sim, a_par = a, d_par = d,
#     item_labels=item_labels, iter_group = 5, niter = 4)
# 
# simulateTrajectories(true_sim = true_sim, a_par = a, d_par = d,
#     item_labels=item_labels, iter_group = 5, niter = 3)
# code to read and merge iteration group files
# for (itgrp in 1:iter_group) {
#     if (itgrp == 1) {
#         re_mmse <- readRDS(paste0(out_dir,"re_mmse_iteration_group_", itgrp, ".rds"))
#         re_tics <- readRDS(paste0(out_dir,"re_tics_iteration_group_", itgrp, ".rds"))
#         re_hcap <- readRDS(paste0(out_dir,"re_hcap_iteration_group_", itgrp, ".rds"))
#         re_all <- readRDS(paste0(out_dir,"re_all_iteration_group_", itgrp, ".rds"))
#     } else {
#         re_mmse_temp <- readRDS(paste0(out_dir,"re_mmse_iteration_group_", itgrp, ".rds"))
#         re_mmse <- bind_rows(re_mmse, re_mmse_temp)
#         re_tics_temp <- readRDS(paste0(out_dir,"re_tics_iteration_group_", itgrp, ".rds"))
#         re_tics <- bind_rows(re_tics, re_tics_temp)
#         re_hcap_temp <- readRDS(paste0(out_dir,"re_hcap_iteration_group_", itgrp, ".rds"))
#         re_hcap <- bind_rows(re_hcap, re_hcap_temp)
#         re_all_temp <- readRDS(paste0(out_dir,"re_all_iteration_group_", itgrp, ".rds"))
#         re_all <- bind_rows(re_all, re_all_temp)
#     }
# }
# saveRDS(re_all, file = paste0(out_dir,"re_all", ".rds"))

# re_mmse <- readRDS(paste0(out_dir,"re_mmse", ".rds"))
# re_tics <- readRDS(paste0(out_dir,"re_tics", ".rds"))
# re_hcap <- readRDS(paste0(out_dir,"re_hcap", ".rds"))
# re_all <- readRDS(paste0(out_dir,"re_all", ".rds"))

# code to examine correlations of true cognition across iterations
# dsg1 <- readRDS(paste0(out_dir,"ds_iteration_group_1.rds"))
# tc <- dsg1[[1]] %>% dplyr::select(id,time,true_cog) %>% 
#     left_join((dsg1[[5]] %>% dplyr::select(id,time,true_cog)) %>% 
#         rename(true_cog2 = true_cog), by = c("id","time"))
# 
# cor(tc[,c("true_cog","true_cog2")])
# cor(s4[,c("vm1","vm5")])
# 
# 
# # saveRDS(sim20,file="~/Psychometrics Conference/2024/Simulation WG/PsyMCA-2024-Simulation/Analysis/sim20.rds")
# 
# ds <- sim20[["ds"]]
# re_all <- sim20[["re_all"]]
# re_mmse <- sim20[["re_mmse"]]
# re_tics <- sim20[["re_tics"]]
# re_hcap <- sim20[["re_hcap"]]
# 
# # iteration 1 chosen arbitrarily
# summary(lm(slope_sim ~ slope_true,data = re_mmse[[1]]))
# summary(lm(slope_sim ~ slope_true,data = re_tics[[1]]))
# summary(lm(slope_sim ~ slope_true,data = re_hcap[[1]]))
# summary(lm(slope_sim ~ slope_true,data = re_all[[1]]))
# 
# # iteration 2 chosen arbitrarily
# ggplot(data=re_all[[2]],aes(slope_true,slope_sim)) + geom_point() +
#     geom_smooth()
# ggplot(data=re_mmse[[2]],aes(slope_true,slope_sim)) + geom_point() +
#     geom_smooth()
# ggplot(data=re_tics[[2]],aes(slope_true,slope_sim)) + geom_point() +
#     geom_smooth()
# ggplot(data=re_hcap[[2]],aes(slope_true,slope_sim)) + geom_point() +
#     geom_smooth()
# 
# # merge re_ datasets across iterations
# re_all_summ <- bind_rows(re_all)
# re_mmse_summ <- bind_rows(re_mmse)
# re_tics_summ <- bind_rows(re_tics)
# re_hcap_summ <- bind_rows(re_hcap)
# ggplot(data=re_all_summ,aes(slope_true,slope_sim)) + geom_point() + geom_smooth()
# ggplot(data=re_mmse_summ,aes(slope_true,slope_sim)) + geom_point() + geom_smooth()
# ggplot(data=re_tics_summ,aes(slope_true,slope_sim)) + geom_point() + geom_smooth()
# ggplot(data=re_hcap_summ,aes(slope_true,slope_sim)) + geom_point() + geom_smooth()
# 
# df <- readRDS(paste0(out_dir,"ds_iteration_group_2.rds"))[[1]]
# plot(df$true_cog,df$sim_cog_mmse)
# plot(df$true_cog,df$sim_cog_tics)
# plot(df$true_cog,df$sim_cog_hcap)
# plot(df$true_cog,df$sim_cog_all)
# summary(lm(sim_cog_all ~ true_cog, data = df))
# summary(df[,c("true_cog","sim_cog_all","sim_cog_mmse","sim_cog_tics","sim_cog_hcap")])
# mod <- mirt(df[,hcap_nmsr], mirt_hcap)
# df$sim_cog <- fscores(mod)
# df$age <- true_sim$age
# 
# 
# set.seed(08222024)
# rids <- sample(unique(df$id),size = 100, replace = FALSE)
# df5 <- df %>%  filter(id %in% rids)
# 
# # vm_pred - mean of 4000 simulations
# ggplot(df5, aes(x = age, y = sim_cog)) + 
#     geom_line(aes(group=id, color=id), show.legend = FALSE) + 
#     geom_smooth(linetype=1,linewidth=1.5,show.legend = FALSE)  +
#     # scale_x_continuous(name = "Time from Baseline (Years)", limits = c(0,10)) + 
#     scale_x_continuous(name = "Age (Years)", limits = c(60,100)) + 
#     scale_y_continuous(name = "Cognitive Score",limits=c(-3,3))
#   

#   ----------------------------- End Run/Test ---------------------------------



