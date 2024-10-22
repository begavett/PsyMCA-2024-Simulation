require(lme4)
require(mirt)
require(tidyverse)
require(ggplot2)

# show size of environment objects
# sort( sapply(ls(),function(x){object.size(get(x))})) 

# rm(list = ls()[grepl("^re",ls())])



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


#   simTraj is a function that inputs a simulated dataset and 1) estimates linear 
#       mixed effects models with true cognition and simulated cognition 
#       (factor scores from simulated item responses) as dependent 
#       variables, time in study as a fixed effect variable, and person and 
#       person-by-time random effects, and 2) returns model predicted estimates of the 
#       dependent variable in the simulated dataset, person specific estimates 
#       (random effects) of the cognitive measure slopes and intercepts, and
#       fixed effects estimates from the model.
#   Parameters
#       data - data frame that contains true cognition value, simulated item responses, 
#           id variable, and time in study variable
#       sim_cog_var - label for variable within dataframe (data) for which 
#           trajectories are being simulated
#       frml - formula specification for (lmer) longitudinal mixed effects model
#           passed from simulateTrajectories
#       iteration = iteration number passed from simulateTrajectories
#   Value - list with 3 elements:
#       data - analytic dataframe (data) with longitudinal model predicted cognition 
#           values for each assessment
#       re - data frame with estimated random effects from longitudinal model. 
#           Includes intercept and slope random effects for true cognition (_true) 
#           and simulated cognition factor scores (_sim)
#       fe - data frame with estimated fixed effects from longitudinal model. 

# simTraj <- function(data = df, item_labels = item_labels, mirt_mod = mirt_all,
#     iteration = iteration) {
simTraj <- function(data = df, sim_cog_var = "sim_cog_all", frml = frml,
    iteration = iteration) {
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
        # str(coef(res_long_true))
        data$cog_true_pred <- predict(res_long_true, newdata = data, type = "response")
        re_true <- data.frame(ranef(res_long_true))
        re_true <- re_true %>% 
            pivot_wider(id_cols = grp, names_from = term,values_from = c(condval,condsd))
        names(re_true) <- c("id","int_true","slope_true","int_true_se","slope_true_se")
        re_true$id <- as.integer(levels(re_true$id))[re_true$id]
        fe_true <- data.frame(t(lme4::fixef(res_long_true)))
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
        fe_sim <- data.frame(t(lme4::fixef(res_long_sim)))
        names(fe_sim) <- paste0(names(fe_sim),"_sim")
        
        re <- re_true %>% left_join(re_sim, by="id")
        fe <- bind_cols(fe_true,fe_sim)
        fe$iteration <- iteration
        
    }, error=function(e){re <- re_empty})
    re$iteration <- iteration
    
    return(list(data,re,fe))
    
}

#   simulateTrajectories is a function that 1) simulates item responses for 
#   true cognition vectors for different measures derived from the HRS-HCAP cognitive test
#   battery (MMSE, HRS TICS, HCAP, MMSE+TICS+HCAP), 2) generates simulated observed (factor)
#   scores for each of the cognitive measures, 3) estimates linear mixed effects models
#   with true cognition and simulated cognition (factor scores from simulated item responses)
#   as dependent variables, time in study as a fixed effect variable, and person and 
#   person-by-time random effects (simTraj()), 4) collates and outputs simulated datasets
#   that include simulated item level responses and factor scores, random effects,
#   and fixed effects from mixed effects models.
#   Up to 10 cognitive measures can be processed in each run.
#
#   Parameters
#       iter_group - number of groups of (niter) iterations
#       niter - number of simulation iterations within iteration groups
#       item_labels - vector that contains variable names of superset of items that occur 
#           in any of the simulated measures (union of items across simulated measures)
#       * vectors with mirt design parameters (must have the same number of elements)
#           models - vector, each element contains mirt model syntax for one simulated measure
#               e.g. "dv1 = 1:5"
#           dv_labels - vector, each element contains factor label for one simulated measure
#               e.g. "dv1"
#           items - vector, each element is a vector that contains item numbers
#               (corresponding to item_labels) for one simulated measure 
#               e.g. "c(1:3,4,7)"
#       in_dir - path to directory for input of previously generated simulated
#           item response datasets. if is.null(in_dir) item responses will be simulated
#           using true_sim, a_par, and d_par
#       * required if is.null(in_dir):
#           true_sim - data frame with true cognition values, rows correspond to assessments,
#               columns correspond to simulation iterations.
#           a_par - vector/data frame containing a (discrimination) item parameters 
#               from mirt co-calibration model
#           d_par - vector/data frame containing d (threshhold) item parameters 
#               from mirt co-calibration model
#       out_dir - path to directory for output of simulation results
#       frml - formula specification for (lmer) longitudinal mixed effects model
#           e.g. "true_cog ~ time + (1 + time | id)"
#
#   Value - Saves lists of results to out_dir - there is a file for each iter_group 
#       that contains a list with niter elements. There are also files that merge
#       random effects and fixed effects across iterations. These files contain granular
#       results from the simulations and can be individually and cumulatively large.
#
#       ds_iteration_group_<iter_group number> - list of datasets for each of the 
#           niter simulations. Each dataset contains the true cognition value (true_cog), 
#           the simulated item responses conditional on true_cog, cognition factor scores, 
#           and longitudinal mixed effect model predicted cognition values for each simulated
#           measure specified in dv_labels. Can be large.
#       re_<dv>_iteration_group_<iter_group number> - list of estimated random effects 
#           for each of the niter simulations  for each simulated measure specified 
#           in dv_labels. Includes intercept and slope random effects for 
#           cognition (_true) and simulated cognition factor scores (_sim)
#       fe_<dv>_iteration_group_<iter_group number> - list of estimated fixed effects 
#           for each of the niter simulations for each simulated measure specified in dv_labels.
#       re_<dv> - dataframe that contains random effects from all of the iterations
#           (iter_group * niter) for specified simulated measure (dv). 
#           Includes intercept and slope random effects for 
#           cognition (_true) and simulated cognition factor scores (_sim).
#           Can be large.
#       fe_<dv> - dataframe that contains fixed effects from all of the iterations
#           (iter_group * niter) for specified simulated measure (dv).

simulateTrajectories <- function(iter_group = 5, niter = 20, item_labels, 
                    true_sim, a_par, d_par, 
                    models, dv_labels, items,
                    in_dir = NULL,
                    out_dir = "Analysis/Simulation Results/",
                    frml = "true_cog ~ time + (1 | id)") {
    require(mirt)
    require(tidyverse)
    
    # design <- data.frame(dv_labels,models,items)
    varnms <- paste0("sim_cog_",dv_labels)
    varnms_rs <- paste0(varnms,"_rs")
    for (itgrp in 1:iter_group) {
        ds <- list()
        for (i in 1:length(dv_labels)) {
            assign(paste0("re_",i),list())
            assign(paste0("fe_",i),list())
        }

        set.seed(092724)
        for (iter in 1:niter) {
            cat(paste0("Group - ",itgrp, ", Iteration - ",iter,"\n"))
            
            # simulate item responses
            iteration <- ((itgrp - 1) * niter) + iter
            
            if (is.null(in_dir)){
                df <- data.frame(simdata(a = a_par, d = d_par, Theta = true_sim[,paste0("vm",iteration)], 
                    itemtype = 'graded'))
                names(df) <- item_labels
                df$id <- true_sim$id
                df$time <- true_sim$time
                df$agebl_75 <- true_sim$agebl_75
                df$true_cog <- true_sim[,paste0("vm",iteration)]
                df$iteration <- iteration
                df <- df %>% relocate(c(id,iteration,time,true_cog))
                
                # true random effect intercept and slope - from calculate_re_true.R 
                true_re <- readRDS("~/Psychometrics Conference/2024/Simulation WG/PsyMCA-2024-Simulation/Data/true_random_effects.rds")
                df <- df %>% left_join((true_re %>% dplyr::select(id,int_true_qrtl,
                        slope_true_qrtl)),by="id")
                
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
            } else {
                df <- readRDS(paste0(in_dir,"ds_iteration_group_",itgrp,".rds"))[[iter]]
                df <- df %>% dplyr::select(c(id:true_cog,agebl_75:slope_true_qrtl,any_of(item_labels)))
            }
            
            for (i in 1:length(dv_labels)) {
                assign(paste0("mirt_dv",i),models[i]) 
            }
            
            ### estimate mirt model using simulated responses and generate factor scores
            # create empty re dataframe to return if error in generation
            re_empty <- data.frame(id = NA, int_true = NA, slope_true = NA, 
                                   int_true_se = NA,slope_true_se = NA, int_sim = NA, slope_sim = NA, 
                                   int_sim_se = NA, slope_sim_se = NA)
            
            tryCatch({
                
                for (i in 1:length(dv_labels)) {
                    itms <- eval(parse(text = items[i]))
                    mod <- mirt(df[,item_labels[itms]], models[i])
                    df[,paste0("sim_cog_",dv_labels[i])] <- fscores(mod)
                }
                # varnms <- paste0("sim_cog_",dv_labels)
                df <- df %>% 
                    relocate(all_of(varnms), 
                             .after = true_cog)
            }, error=function(e){return(re_empty)})

            # rescale factor scores to equate metric to true cognition
            
            # varnms_rs <- paste0(varnms,"_rs")
            
            for (i in 1:length(dv_labels)) {
                df[,varnms_rs[i]] <-
                    (((df[,varnms[i]] - mean(df[,varnms[i]])) / 
                          sd(df[,varnms[i]])) * sd(df$true_cog)) + mean(df$true_cog)
            }
            
            # mixed effects model of simulated longitudinal data
            for (i in 1:length(dv_labels)) {
                res <- simTraj(data = df, sim_cog_var = varnms_rs[i],frml = frml,
                    iteration = iteration)
                
                df <- res[[1]]
                names(df) <- sub("cog_true_pred$",paste0("cog_true_pred_",
                    dv_labels[i]),names(df))
                names(df) <- sub("cog_sim_pred$",paste0("cog_sim_pred_",
                    dv_labels[i]),names(df))
                
                re <- res[[2]]
                re$label <- toupper(dv_labels[i])
                
                fe <- res[[3]]
                fe$label <- toupper(dv_labels[i])
                
                if (i == 1) {
                    re_1[[paste0("iteration-",iteration)]] = re
                    fe_1[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 2) {
                    re_2[[paste0("iteration-",iteration)]] = re
                    fe_2[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 3) {
                    re_3[[paste0("iteration-",iteration)]] = re
                    fe_3[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 4) {
                    re_4[[paste0("iteration-",iteration)]] = re
                    fe_4[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 5) {
                    re_5[[paste0("iteration-",iteration)]] = re
                    fe_5[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 6) {
                    re_6[[paste0("iteration-",iteration)]] = re
                    fe_6[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 7) {
                    re_7[[paste0("iteration-",iteration)]] = re
                    fe_7[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 8) {
                    re_8[[paste0("iteration-",iteration)]] = re
                    fe_8[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 9) {
                    re_9[[paste0("iteration-",iteration)]] = re
                    fe_9[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 10) {
                    re_10[[paste0("iteration-",iteration)]] = re
                    fe_10[[paste0("iteration-",iteration)]] = fe
                }
            }
            
            ds[[paste0("iteration-",iteration)]] <- df
            
            
        } # end for iter
        
        saveRDS(ds, file = paste0(out_dir,"ds_iteration_group_", itgrp, ".rds"))

        for (i in 1:length(dv_labels)) {
            saveRDS(eval(parse(text=paste0("re_",i))), 
                    file = paste0(out_dir,"re_",dv_labels[i],
                                  "_iteration_group_", itgrp, ".rds"))
            saveRDS(eval(parse(text=paste0("fe_",i))), 
                    file = paste0(out_dir,"fe_",dv_labels[i],
                                  "_iteration_group_", itgrp, ".rds"))
        }

    } # end for itgrp
    
    for (itgrp in 1:iter_group) {
        if (itgrp == 1) {
            for (i in 1:length(dv_labels)) {
                assign(paste0("re_",i),
                       readRDS(paste0(out_dir,"re_",dv_labels[i],
                                      "_iteration_group_", itgrp, ".rds")) %>% bind_rows())
                assign(paste0("fe_",i),
                       readRDS(paste0(out_dir,"fe_",dv_labels[i],
                                      "_iteration_group_", itgrp, ".rds")) %>% bind_rows())
            }
            
        } else {
            for (i in 1:length(dv_labels)) {
                re_temp <- readRDS(paste0(out_dir,"re_",dv_labels[i],
                                          "_iteration_group_", itgrp, ".rds")) %>% bind_rows()
                assign(paste0("re_",i),
                       eval(parse(text=paste0("re_",i))) %>% bind_rows(re_temp))
                fe_temp <- readRDS(paste0(out_dir,"fe_",dv_labels[i],
                                          "_iteration_group_", itgrp, ".rds")) %>% bind_rows()
                assign(paste0("fe_",i),
                       eval(parse(text=paste0("fe_",i))) %>% bind_rows(fe_temp))
            }
            
        }
    }
    for (i in 1:length(dv_labels)) {
        saveRDS(eval(parse(text=paste0("re_",i))), 
                file = paste0(out_dir,"re_",dv_labels[i],".rds"))
        saveRDS(eval(parse(text=paste0("fe_",i))), 
                file = paste0(out_dir,"fe_",dv_labels[i],".rds"))
    }
    
}

#   simulateBlendedTrajectories is a function that inputs previously generated cognitive
#       test scores, "blends" results from two measures into a new variable
#       (the first measure will comprise the baseline assessment through blend_time, 
#       the second will include assessments subsequent to blend_time). It estimates
#       linear mixed effect (lmer) models via simTraj() and outputs results that are 
#       comparable to those from simulateTrajectories(). 
#   Up to 10 cognitive measures can be processed in each run.
#
#   Parameters
#       iter_group - number of groups of (niter) iterations
#       niter - number of simulation iterations within iteration groups
#       blend_labels - vector that contains label for measures that will be blended.
#           each element is a vector with two elements that are labels for the
#           first and second measures to be blended (e.g. "tics3-tics13"). A simple,
#           non-blended measure can be specified by repeating the measure 
#           (e.g. "tics3-tics3").
#       blend_time - time-in-study cut-point for assigning first or second measure
#           to blended measure. First is assigned to assessments with time ≤ blend_time,
#           second to time > blend_time.
#       in_dir - path to directory for input of previously generated simulated
#           item response datasets. if is.null(in_dir) item responses will be simulated
#           using true_sim, a_par, and d_par
#       out_dir - path to directory for output of simulation results
#       frml - formula specification for (lmer) longitudinal mixed effects model
#           e.g. "true_cog ~ time + (1 + time | id)"
#
#   Value - Saves lists of results to out_dir - there is a file for each iter_group 
#       that contains a list with niter elements. There are also files that merge
#       random effects and fixed effects across iterations. These files contain granular
#       results from the simulations and can be individually and cumulatively large.
#
#       ds_iteration_group_<iter_group number> - list of datasets for each of the 
#           niter simulations. Each dataset contains the true cognition value (true_cog), 
#           the simulated item responses conditional on true_cog, cognition factor scores, 
#           and longitudinal mixed effect model predicted cognition values for each simulated
#           measure specified in dv_labels. Can be large.
#       re_<dv>_iteration_group_<iter_group number> - list of estimated random effects 
#           for each of the niter simulations  for each simulated measure specified 
#           in dv_labels. Includes intercept and slope random effects for 
#           cognition (_true) and simulated cognition factor scores (_sim)
#       fe_<dv>_iteration_group_<iter_group number> - list of estimated fixed effects 
#           for each of the niter simulations for each simulated measure specified in dv_labels.
#       re_<dv> - dataframe that contains random effects from all of the iterations
#           (iter_group * niter) for specified simulated measure (dv). 
#           Includes intercept and slope random effects for 
#           cognition (_true) and simulated cognition factor scores (_sim).
#           Can be large.
#       fe_<dv> - dataframe that contains fixed effects from all of the iterations
#           (iter_group * niter) for specified simulated measure (dv).




simulateBlendedTrajectories <- function(niter = 20, iter_group = 5,
                    blend_labels = blend_labels,
                    blend_time = blend_time,
                    in_dir = "Analysis/Simulation Results/",
                    out_dir = "Analysis/Simulation Results/",
                    frml = "true_cog ~ time + (1 | id)") {
    require(mirt)
    require(tidyverse)
    require(stringr)
    # source("~/Psychometrics Conference/2024/Simulation WG/PsyMCA-2024-Simulation/Code/simulate_trajectories.R")
    
    dv_labels <- blend_labels
    design <- data.frame(dv_labels)
    design <- design %>% mutate(
        dv1 = str_extract(dv_labels,"^.*(?=-)"),
        dv2 = str_extract(dv_labels,"(?<=-).*")
    )
    
    for (itgrp in 1:iter_group) {
        ds <- list()
        for (i in 1:length(design$dv_labels)) {
            assign(paste0("re_",i),list())
            assign(paste0("fe_",i),list())
        }
        set.seed(092724)
        for (iter in 1:niter) {
            cat(paste0("Group - ",itgrp, ", Iteration - ",iter,"\n"))
            
            #  load simulated item responses
            iteration <- ((itgrp - 1) * niter) + iter
            df <- readRDS(paste0(in_dir,"ds_iteration_group_",itgrp,".rds"))[[iter]]
            df <- df %>% dplyr::select(c(id:true_cog,agebl_75:slope_true_qrtl,any_of(item_labels),
                                         ends_with("_rs")))
            
            for (i in 1:length(design$dv_labels)) {
                lbl <- paste0("sim_cog_",design$dv_labels[i],"_rs")
                lbl1 <- paste0("sim_cog_",design$dv1[i],"_rs")
                lbl2 <- paste0("sim_cog_",design$dv2[i],"_rs")
                df[,lbl] <- case_when(
                    df$time <= blend_time ~ df[,lbl1],
                    df$time > blend_time ~ df[,lbl2]
                )
                
            }
            
            
            # mixed effects model of simulated longitudinal data
            
            varnms_rs <- paste0("sim_cog_",design$dv_labels,"_rs")
            
            for (i in 1:length(design$dv_labels)) {
                res <- simTraj(data = df, sim_cog_var = varnms_rs[i],frml = frml,
                               iteration = iteration)
                
                df <- res[[1]]
                names(df) <- sub("cog_true_pred$",paste0("cog_true_pred_",
                                                         design$dv_labels[i]),names(df))
                names(df) <- sub("cog_sim_pred$",paste0("cog_sim_pred_",
                                                        design$dv_labels[i]),names(df))
                
                re <- res[[2]]
                re$label <- toupper(design$dv_labels[i])
                
                fe <- res[[3]]
                fe$label <- toupper(design$dv_labels[i])
                
                if (i == 1) {
                    re_1[[paste0("iteration-",iteration)]] = re
                    fe_1[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 2) {
                    re_2[[paste0("iteration-",iteration)]] = re
                    fe_2[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 3) {
                    re_3[[paste0("iteration-",iteration)]] = re
                    fe_3[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 4) {
                    re_4[[paste0("iteration-",iteration)]] = re
                    fe_4[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 5) {
                    re_5[[paste0("iteration-",iteration)]] = re
                    fe_5[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 6) {
                    re_6[[paste0("iteration-",iteration)]] = re
                    fe_6[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 7) {
                    re_7[[paste0("iteration-",iteration)]] = re
                    fe_7[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 8) {
                    re_8[[paste0("iteration-",iteration)]] = re
                    fe_8[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 9) {
                    re_9[[paste0("iteration-",iteration)]] = re
                    fe_9[[paste0("iteration-",iteration)]] = fe
                }
                if (i == 10) {
                    re_10[[paste0("iteration-",iteration)]] = re
                    fe_10[[paste0("iteration-",iteration)]] = fe
                }
            }
            
            ds[[paste0("iteration-",iteration)]] <- df
            
            
        } # end for iter
        
        saveRDS(ds, file = paste0(out_dir,"ds_iteration_group_", itgrp, ".rds"))
        
        for (i in 1:length(design$dv_labels)) {
            saveRDS(eval(parse(text=paste0("re_",i))), 
                    file = paste0(out_dir,"re_",design$dv_labels[i],
                                  "_iteration_group_", itgrp, ".rds"))
            saveRDS(eval(parse(text=paste0("fe_",i))), 
                    file = paste0(out_dir,"fe_",design$dv_labels[i],
                                  "_iteration_group_", itgrp, ".rds"))
        }
        
    } # end for itgrp
    
    for (itgrp in 1:iter_group) {
        if (itgrp == 1) {
            for (i in 1:length(design$dv_labels)) {
                assign(paste0("re_",i),
                       readRDS(paste0(out_dir,"re_",design$dv_labels[i],
                                      "_iteration_group_", itgrp, ".rds")) %>% bind_rows())
                assign(paste0("fe_",i),
                       readRDS(paste0(out_dir,"fe_",design$dv_labels[i],
                                      "_iteration_group_", itgrp, ".rds")) %>% bind_rows())
            }
            
            
        } else {
            for (i in 1:length(design$dv_labels)) {
                re_temp <- readRDS(paste0(out_dir,"re_",design$dv_labels[i],
                                          "_iteration_group_", itgrp, ".rds")) %>% bind_rows()
                assign(paste0("re_",i),
                       eval(parse(text=paste0("re_",i))) %>% bind_rows(re_temp))
                fe_temp <- readRDS(paste0(out_dir,"fe_",design$dv_labels[i],
                                          "_iteration_group_", itgrp, ".rds")) %>% bind_rows()
                assign(paste0("fe_",i),
                       eval(parse(text=paste0("fe_",i))) %>% bind_rows(fe_temp))
            }
            
            
            
        }
    }
    
    for (i in 1:length(design$dv_labels)) {
        saveRDS(eval(parse(text=paste0("re_",i))), 
                file = paste0(out_dir,"re_",design$dv_labels[i],".rds"))
        saveRDS(eval(parse(text=paste0("fe_",i))), 
                file = paste0(out_dir,"fe_",design$dv_labels[i],".rds"))
    }
    
}


# dat_cols <- c("id","time","age","agebl_75","slope_true_qrtl","int_true_qrtl")
# res_cols <- c("cog_true_pred_all","cog_sim_pred_all","cog_sim_pred_mmse",
#               "cog_sim_pred_tics","cog_sim_pred_hcap")

#   mergeSimResults is a function that merges model predicted outcomes from lists of 
#       simulated datasets generated by simulateTrajectories() and simulateBlendedTrajectories().
#       It returns a dataframe with number of rows equal to the number of rows in the 
#       simulated dataset and columns for each iteration and each simulated measure.
#
#   Parameters:
#       file_name - template for names of files to be inputted and merged. 
#           "_1" will be programmatically replaced by increasing group numbers.
#       out_dir - path to directory for output of merged file
#       dat_cols - vector containing names of columns (other than cognitive measures) to
#           be included in merged data file
#       res_cols - vector containing names of cognitive measure columns to
#           be included in merged data file
#       iter_group - number of groups of (niter) iterations
#       niter - number of simulation iterations within iteration groups
#
#   Value - Returns dataframe that contains dat_cols variables common to the simulated 
#       datasets and columns for each iteration and each cognitive measure specified
#       in res_cols.


mergeSimResults <-function(file_name = "ds_iteration_group_1.rds",
                           out_dir = out_dir,
                           dat_cols = dat_cols, 
                           res_cols = res_cols,
                           iter_group = 50,
                           niter = 20) {
    res <- (readRDS(paste0(out_dir,file_name))[[1]] %>%
                dplyr::select(any_of(dat_cols)))
    for (i in 1:iter_group) {
        for (n in 1:niter) {
            iteration <- ((i-1) * niter) + n
            fl_nm <- sub("_1",paste0("_",i),file_name)
            ds <- readRDS(paste0(out_dir,fl_nm))[[n]] %>% 
                dplyr::select(all_of(res_cols)) 
            names(ds) <- paste0(res_cols,".",iteration)
            res <- res %>% bind_cols(ds)
        }
    }   
        return(res)
}
    


#   ************************* Run/Test Simulation Code *************************

#   Code in this segment can be uncommented to run and test the functions in this file.
#   This code extracts longitudinal "true" cognition values from
#   "simulated_longitudinal_true_cognition.rds". This file is a synthetic dataset
#   with 792 individuals uniformly distributed in age from 60-95, each with baseline 
#   assessment and 10 follow-up assessments at 1 year intervals. Parameters to build 
#   this dataset came from an empirical analysis of a large dataset representative 
#   of older adults. This dataset is a long-form dataset that includes a person id
#   variable (id) and time-in-study variable (time) that are used in longitudinal
#   mixed effects model generated in the simTraj() function. It contains 4000
#   simulations of true cognition for the longitudinal data structure (labelled
#   vm1-vm4000). The simulation of true cognition is described in
#   longitudinal_cognition_simulation.<html/Rmd>.


# # code to use simulateTrajectories() to generate simulated item response datasets
# #   from simulated true cognition (true_sim) and item parameters (a_par and d_par)
# #   and run mixed effects longitudinal models in which mirt estimated factor scores 
# #   for each simulated cognitive measures are dependent variables
# #   Note that in_dir must be NULL
# 
# source("~/Psychometrics Conference/2024/Simulation WG/PsyMCA-2024-Simulation/Code/simulate_trajectories.R")
# 
# item_labels <- c("wlimrc_ticsr","wldlrc_ticsr","cntbck_ticsr","animfl_ticsr",
#                  "numser_ticsr","serial7_ticsr","ornttime_ticsr","naming_ticsr",
#                  "trailsa_hcapir","trailsb_hcapir")
# dv_labels <- c("tics3","tics10","tics13","tics16")
# models <- c(
#     "tics3 = 1-6",
#     "tics10 = 1-7",
#     "tics13 = 1-8",
#     "tics16 = 1-10"
# )
# items <- c(
#     "c(1:3,6:8)",
#     "c(1:3,5:8)",
#     "c(1:8)",
#     "c(1:10)"
# )
# out_dir = "Analysis/Simulation Results/Temp2/"
# 
# frml <- "true_cog ~ time + agebl_75 + slope_true_qrtl +
#             time:agebl_75 + time:slope_true_qrtl + (1 + time | id)"
# 
# # load simulated "true" cognition
# true_sim <- readRDS("~/Psychometrics Conference/2024/Simulation WG/PsyMCA-2024-Simulation/Data/simulated_longitudinal_true_cognition.rds")
# 
# # load item parameters
# load("~/Psychometrics Conference/2024/Simulation WG/PsyMCA-2024-Simulation/Analysis/co_calibration_results.Rdata")
# 
# a <- extractMIRTParm((res_md_1a))[[1]][12:21]
# d <- extractMIRTParm((res_md_1a))[[2]][12:21,]
# 
# simulateTrajectories(iter_group=3,niter=3,item_labels=item_labels,
#                      true_sim=true_sim,a_par=a,d_par=d,
#                      models=models,dv_labels=dv_labels,items=items,
#                      in_dir=NULL,out_dir=out_dir,frml=frml)


# # code to use simulateTrajectories() to input previously created simulated datasets
# #   and run mixed effects longitudinal models in which mirt estimated factor scores 
# #   for each simulated cognitive measures are dependent variables
# #   Note that in_dir is the folder with the previously created simulated datasets
# 
# source("~/Psychometrics Conference/2024/Simulation WG/PsyMCA-2024-Simulation/Code/simulate_trajectories.R")
# 
# item_labels <- c("wlimrc_ticsr","wldlrc_ticsr","cntbck_ticsr","animfl_ticsr",
#                  "numser_ticsr","serial7_ticsr","ornttime_ticsr","naming_ticsr",
#                  "trailsa_hcapir","trailsb_hcapir")
# dv_labels <- c("tics3","tics10","tics13","tics16")
# models <- c(
#     "tics3 = 1-6",
#     "tics10 = 1-7",
#     "tics13 = 1-8",
#     "tics16 = 1-10"
# )
# items <- c(
#     "c(1:3,6:8)",
#     "c(1:3,5:8)",
#     "c(1:8)",
#     "c(1:10)"
# )
# 
# out_dir = "Analysis/Simulation Results/Temp3/"
# in_dir = "Analysis/Simulation Results/Temp2/"
# 
# frml <- "true_cog ~ time + agebl_75 + slope_true_qrtl +
#             time:agebl_75 + time:slope_true_qrtl + (1 + time | id)"
# 
# 
# simulateTrajectories(iter_group=3,niter=3,item_labels=item_labels,
#                      true_sim=NULL,a_par=NULL,d_par=NULL,
#                      models=models,dv_labels=dv_labels,items=items,
#                      in_dir=in_dir,out_dir=out_dir,frml=frml)


# # code to use simulateBlended Trajectories() to created blended cognitive measures
# #   from previously created simulated datasets
# #   and run mixed effects longitudinal models in which mirt estimated factor scores 
# #   for each of the blended cognitive measures are dependent variables.
# #   Note that in_dir is the folder with the previously created simulated datasets
# 
# blend_labels <- c(
#     c("tics3-tics13"),
#     c("tics3-tics3"),
#     c("tics10-tics13"),
#     c("tics10-tics10"),
#     c("tics13-tics16"),
#     c("tics13-tics13"),
#     c("tics3-tics16"),
#     c("tics16-tics16")
# )
# 
# blend_time <- 5
# 
# in_dir = "Analysis/Simulation Results/Temp2/"
# out_dir = "Analysis/Simulation Results/Temp4/"
# 
# frml <- "true_cog ~ time + agebl_75 + slope_true_qrtl +
#     time:agebl_75 + time:slope_true_qrtl + (1 + time | id)"
# 
# 
# simulateBlendedTrajectories(niter=3,iter_group=3,blend_labels = blend_labels,
#     blend_time = blend_time,in_dir=in_dir,out_dir=out_dir,frml=frml)

# #   code to merge model predicted outcomes from lists of simulated datasets 
# #       generated by simulateTrajectories() and simulateBlendedTrajectories().
# 
# dat_cols <- c("id","time","agebl_75","slope_true_qrtl","int_true_qrtl")
# res_cols <- c("cog_true_pred_tics3-tics13","cog_sim_pred_tics3-tics13",
#               "cog_sim_pred_tics3-tics3","cog_sim_pred_tics10-tics13",
#               "cog_sim_pred_tics10-tics10","cog_sim_pred_tics13-tics16",
#               "cog_sim_pred_tics13-tics13","cog_sim_pred_tics3-tics16",
#               "cog_sim_pred_tics16-tics16")
# out_dir <- "Analysis/Simulation Results/Temp4/"
# 
# r1 <- mergeSimResults(file_name = "ds_iteration_group_1.rds", out_dir = out_dir, 
#                       dat_cols = dat_cols, res_cols = res_cols, iter_group = 3, niter = 3)
# saveRDS(r1,file=paste0(out_dir,"model_predicted_trajectory_simulations.rds"))
# 
# r1 <- readRDS(paste0(out_dir,"model_predicted_trajectory_simulations.rds"))
# 
# rm(r1)



#   ----------------------------- End Run/Test ---------------------------------



