
#   rValeMaurelliModified() is a modification of SimDesign::rValeMaurelli(). It
#   accommodates a multiple group structure. The modification section generates 
#   random normal data for both groups and merges it so that non-normal perturbations
#   are subsequently applied to the merged dataset - that is, skewness and
#   kurtosis parameters are applied to the merged dataset.

rValeMaurelliModified <- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
    skew = rep(0, nrow(sigma)), kurt = rep(0, nrow(sigma)),
    rmean = 0, rsd = 1, focal_mn = NULL, focal_sd = NULL,
    focal_mn_z = NULL, focal_sd_z = NULL) 
{
    stopifnot(!missing(n))
    if (length(sigma) == 1L) 
        sigma <- as.matrix(sigma)
    stopifnot(is.matrix(sigma))
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                     check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != nrow(sigma)) 
        stop("mean and sigma have non-conforming size")
    stopifnot(n > 0)
    stopifnot(ncol(sigma) == length(skew))
    stopifnot(ncol(sigma) == length(kurt))
    sds <- sqrt(diag(sigma))
    cor <- cov2cor(sigma)
    chol_corr <- t(chol(cor))
    k <- ncol(cor)
    for (i in 1:k) {
        if (kurt[i] <= skew[i]^2 - 2) {
            stop("Error: the ", i, " th component of kurtosis is not bigger than skewness squared minus 2.\n")
        }
    }
    constant <- function(sk, ku, start) {
        F <- function(x) {
            F <- numeric(3)
            b <- x[1]
            c <- x[2]
            d <- x[3]
            F[1] <- b^2 + 6 * b * d + 2 * c^2 + 15 * d^2 - 1
            F[2] <- 2 * c * (b^2 + 24 * b * d + 105 * d^2 + 
                                 2) - sk
            F[3] <- 24 * (b * d + c^2 * (1 + b^2 + 28 * b * 
                d) + d^2 * (12 + 48 * b * d + 141 * c^2 + 225 * 
                d^2)) - ku
            F
        }
        obj.fun <- function(par) {
            sum(F(par)^2)
        }
        opt <- nlminb(start = start, objective = obj.fun, control = list(abs.tol = 1e-10, 
            rel.tol = 1e-10, eval.max = 1e+06, iter.max = 1e+06))
        if (opt$converge != 0 || opt$objective > 1e-05) 
            stop("optimizer could not find suitable solution for c0, c1, and c2")
        x0 <- opt$par
        x0
    }
    constants <- matrix(nrow = k, ncol = 4)
    for (i in 1:k) {
        constants[i, 2:4] <- t(constant(skew[i], kurt[i], start = c(1, 0, 0)))
        constants[i, 1] <- -(constants[i, 3])
    }
    solve.p12 <- function(r12, a1, a2, b1, b2, c1, c2, d1, d2) {
        ftn <- function(p12) {
            ((b1 * b2 + 3 * b1 * d2 + 3 * d1 * b2 + 9 * d1 * 
                 d2) * p12) + ((2 * c1 * c2) * p12^2) + ((6 * 
                d1 * d2) * p12^3) - r12
        }
        root <- uniroot(ftn, c(-1, 1), check.conv = TRUE)
        p12 <- root$root
        p12
    }
    inter <- matrix(0, k, k)
    for (i in seq_len(k)) {
        for (j in i:k) {
            if (i == j) 
                next
            inter[i, j] <- solve.p12(cor[i, j], constants[i, 
                1], constants[j, 1], constants[i, 2], constants[j, 
                2], constants[i, 3], constants[j, 3], constants[i, 
                4], constants[j, 4])
            inter[j, i] <- inter[i, j]
        }
    }
    
    ### modification
    if (!is.null(focal_mn_z) & !is.null(focal_sd_z)) {
        rdat <- rnorm(n * (k-1), mean = rmean, sd = rsd) %>% matrix(nrow=n,ncol=(k-1))
        rdat1 <- rnorm(n, mean = rmean, sd = rsd) %>% matrix(nrow=n,ncol=1)
        rdat <- cbind(rdat,rdat1)
        rm(rdat1)
    } else {
        rdat <- rnorm(n * k, mean = rmean, sd = rsd) %>% matrix(nrow=n,ncol=k)
    }
    ### end modification
    
    # rdat <- rnorm(n * k, mean = rmean, sd = rsd) %>% matrix(nrow=n,ncol=k)
    
    Z <- t(chol_corr %*% t(matrix(rdat, ncol=k)))
    Z2 <- Z^2
    Z3 <- Z^3
    ### modification
    Y <- matrix(0, nrow = n, ncol = k)
    for (i in 1:k) {
        if (i == k & !is.null(focal_mn_z) & !is.null(focal_sd_z)) {
            Y[, i] <- ((constants[i, 1] + constants[i, 2] * Z[, i] + 
                constants[i, 3] * Z2[, i] + constants[i, 4] * Z3[,i]) * focal_sd_z) + focal_mn_z
        } else {
            Y[, i] <- ((constants[i, 1] + constants[i, 2] * Z[, i] + 
                constants[i, 3] * Z2[, i] + constants[i, 4] * Z3[,i]) * focal_sd) + focal_mn
        }
    }
    ### end modification
    ### original code
    # Y <- matrix(0, nrow = n, ncol = k)
    # for (i in 1:k) {
    #   Y[, i] <- constants[i, 1] + constants[i, 2] * Z[, i] + 
    #       constants[i, 3] * Z2[, i] + constants[i, 4] * Z3[,i]
    # }
    ### end original
    
    Y <- t(t(Y) * sds + mean)
    
    return(Y)
    
}

#   *********************** Dataset Generation Functions ***********************


#   simNonNormal() is a function that generates random multivariate data that is
#   based on random draws of distributions of factor loadings, intercepts, and from 
#   skewness values derived from HCAP "continuous" tests
#   if focal_mn_z andf focal_sd_z are NULL external z variable will not be generated
simNonNormal <- function(
        N = N, # sample size
        nitems = nitems, # number of items
        seed = seed, # random seed
        rmean = 0, # "true" cognition mean in reference group
        rsd = 1, # "true" cognition sd in reference group
        focal_mn = NULL, # "true" cognition mean in focal group
        focal_sd = NULL,  # "true" cognition sd in focal group
        focal_mn_z = NULL,  # external z variable mean in focal group
        focal_sd_z = NULL) {  # external z variable sd in focal group
    require(ggdist) # rstudent_t
    require(dplyr)
    # source("Code/vale_maurelli_modified.R") # replaced by rValeMaurelliModified
    load("Data/fit_kurtosis.RData") # fit_kurt, change directory as needed
    seed2 <- seed
    i <- 1
    while (i <= nitems) {
        cat(paste0("item-",i,"\n"))
        tryCatch({
            if (i == 1){
                # seed2 <- seed2 + 1
                set.seed(seed2)
                # distribution parameters for loadings, intercepts, and skewness
                # derived from separate analyses of HCAP data
                loadings <- rbeta(1, 30.66796, 13.49961)
                intercepts <- rstudent_t(1, 3.6321458, 0.1240406, 0.0294505)
                skewness <- rstudent_t(10, 1.6562347, -0.3542897, 0.3883936)
                kurtosis <- predict(fit_kurt, newdata = data.frame(skewness))
                # kurtosis <- predict(fit_kurt, newdata = data.frame(poly(skewness,2)))
                skew <- skewness[1]
                kurt <- as.numeric(kurtosis[1])
                
                # predicted skewness in HAALSI derived from regression of HAALSI skewness on
                #   HCAP skewness - continuous and categorical items (fit_skew_hcap_haals)
                skewnessf <- predict(fit_skew_hcap_haals,
                    newdata = data.frame(skewness) %>% rename(skew = skewness)) %>% 
                    as.vector()
                kurtosisf <- predict(fit_kurt, newdata = data.frame(skewnessf))
                skewf <- skewnessf[1]
                kurtf <- as.numeric(kurtosisf[1])
                
                Sigma <- loadings %*% t(loadings)
                diag(Sigma) <- 1
                
                set.seed(seed)
                dfoutr <- data.frame(rValeMaurelliModified(n = N, mean = intercepts, 
                    sigma = Sigma, skew = skew, kurt = kurt, rmean = rmean, 
                    rsd = rsd, focal_mn = rmean, focal_sd = rsd))
                names(dfoutr) <- paste0("x",i)
                set.seed(seed-1)
                dfoutf <- data.frame(rValeMaurelliModified(n = N, mean = intercepts, 
                    sigma = Sigma, skew = skewf, kurt = kurtf, rmean = rmean, 
                    rsd = rsd, focal_mn = focal_mn, focal_sd = focal_sd))
                names(dfoutf) <- paste0("x",i) 

                i <- i + 1
            } else {
                
                seed2 <- seed2 + 1
                set.seed(seed2)
                # cat(paste0("seed=",seed2,", skew=",skew2,"\n"))
                loadings2 <- c(loadings,rbeta(1, 30.66796, 13.49961))
                intercepts2 <- c(intercepts,rstudent_t(1, 3.6321458, 0.1240406, 0.0294505))
                skewness <- rstudent_t(10, 1.6562347, -0.3542897, 0.3883936)
                # fit_kurt derived from 2nd degree polynomial regression of observed
                # kurtosis on observed skewness
                kurtosis <- predict(fit_kurt, newdata = data.frame(skewness))
                skew2 <- c(skew,skewness[1])
                kurt2 <- c(kurt,as.numeric(kurtosis[1]))
                
                skewnessf <- predict(fit_skew_hcap_haals,
                    newdata = data.frame(skewness) %>% rename(skew = skewness)) %>% 
                    as.vector()
                kurtosisf <- predict(fit_kurt, newdata = data.frame(skewnessf) %>% 
                    rename(skewness = skewnessf))
                skew2f <- c(skewf,skewnessf[1])
                kurt2f <- c(kurtf,as.numeric(kurtosisf[1]))
                
                Sigma2 <- loadings2 %*% t(loadings2)
                diag(Sigma2) <- 1
                
                set.seed(seed)
                dfoutr <- data.frame(rValeMaurelliModified(n = N, mean = intercepts2, 
                    sigma = Sigma2, skew = skew2, kurt = kurt2, rmean = rmean, 
                    rsd = rsd, focal_mn = rmean, focal_sd = rsd))
                names(dfoutr) <- paste0("x",1:i)
                set.seed(seed-1)
                dfoutf <- data.frame(rValeMaurelliModified(n = N, mean = intercepts2, 
                    sigma = Sigma2, skew = skew2f, kurt = kurt2f, rmean = rmean, 
                    rsd = rsd, focal_mn = focal_mn, focal_sd = focal_sd))
                names(dfoutf) <- paste0("x",1:i)
                 
                loadings <- loadings2
                intercepts <- intercepts2
                skew <- skew2
                kurt <- kurt2
                skewf <- skew2f
                kurtf <- kurt2f
                i <- i+1
            }
        }, error=function(e){i <- i})
    }
    ### add external variable with correlation with true cognitive ability of 0.30 
    if (!is.null(focal_mn_z) & !is.null(focal_sd_z)) {
        error_flag <- TRUE
        while (error_flag == TRUE) {
            tryCatch({
                loadings2 <- c(loadings,0.30)
                intercepts2 <- c(intercepts,0)
                skew2 <- c(skew,0)
                kurt2 <- c(kurt,0)
                skew2f <- c(skewf,0)
                kurt2f <- c(kurtf,0)
                Sigma2 <- loadings2 %*% t(loadings2)
                diag(Sigma2) <- 1
                # nout <- nrow(dfout)
                set.seed(seed)
                dfoutr <- data.frame(rValeMaurelliModified(n = N, mean = intercepts2,
                    sigma = Sigma2, skew = skew2, kurt = kurt2, rmean = rmean,
                    rsd = rsd, focal_mn = rmean, focal_sd = rsd))
                names(dfoutr) <- c(paste0("x",1:nitems),"z1")
                set.seed(seed-1)
                dfoutf <- data.frame(rValeMaurelliModified(n = N, mean = intercepts2,
                    sigma = Sigma2, skew = skew2f, kurt = kurt2f, rmean = rmean,
                    rsd = rsd, focal_mn = focal_mn, focal_sd = focal_sd, focal_mn_z = focal_mn_z, 
                        focal_sd_z = focal_sd_z))
                names(dfoutf) <- c(paste0("x",1:nitems),"z1")
                 error_flag <- FALSE
            }, error=function(e){error_flag <- TRUE})
        }
        
    }
    ### end external variable
    
    dfout <- rbind(dfoutr,dfoutf)

    return(list("data" = dfout, "loadings" = loadings2, "intercepts" = intercepts2,
        "skewness" = skew2, "kurtosis" = kurt2, "skewnessf" = skew2f, "kurtosisf" = kurt2f))
}

#   ------------------------ End Dataset Generation ----------------------------

# dfout$group <- ifelse(as.numeric(rownames(dfout)) <= 5000,"Reference","Focal")
# summary(dfout[dfout$group == "Reference",])
# summary(dfout[dfout$group == "Focal",])

#   ********** Functions for Test Information from Simulated Datasets **********


infoMG <- function (mg_mod_obj = mg_mod_obj, group_list) {
    require(mirt)
    require(dplyr)
    Theta <- matrix(seq(-6,6,.01))
    for (grp in group_list) {
        if (grp == group_list[1]) {
            info <- data.frame("Theta" = Theta,
                    "info" = testinfo(mg_mod_obj, Theta, group = grp)) %>%
                mutate(label = grp)
        } else {
            info <- bind_rows(info, (data.frame("Theta" = Theta,
                    "info" = testinfo(mg_mod_obj, Theta, group = grp)) %>%
                mutate(label = grp)))
        }
    }
    return(info)
}


#   simulateInfoMG() does one random draw from a 2-group population dataset, does
#   equal frequency and equal interval recoding, uses mirt to estimate a 1-factor
#   graded response model for each recoding type, and generates test information 
#   values

simulateInfoMG <- function(
        pop_data_all = pop_data_all,
        itnms = NULL,
        n_items = n_items,
        nsamp = nsamp,
        seed = 01012025,
        ncat = 10,
        ncell = 10,
        append_collapse = TRUE,
        group_labels = c("Reference","Focal")) {
    require(mirt)
    source("~/Research/Code/recode/recode.R") # change directory as needed
    
    if (is.null(itnms)) {
        varlist_orig <- paste0("x",1:n_items)
    } else {
        varlist_orig <- itnms
    }
    
    varlist_ef <- paste0(varlist_orig,"_ef")
    varlist_ei <- paste0(varlist_orig,"_ei")
    varlist_efc <- paste0(varlist_orig,"_efc")
    varlist_eic <- paste0(varlist_orig,"_eic")
    
    model_mirt <- paste0("cog = 1-",n_items)
    
    set.seed(seed)
    samp_data_rg <- pop_data_all %>%
        filter(group == group_labels[1]) %>%
        slice_sample(n = nsamp, replace = FALSE)
    
    samp_data_fg <- pop_data_all %>%
        filter(group == group_labels[2]) %>%
        slice_sample(n = nsamp, replace = FALSE)
    
    samp_data <- bind_rows(samp_data_rg,samp_data_fg)
    
    # info equal frequency
    
    if(append_collapse){
        varlist_trc <- varlist_efc
        varlist <- varlist_efc
    } else {
        varlist_trc <- NULL
        varlist <- varlist_ef
    }
    samp_data_ef <- recodeOrdinal(df = samp_data,varlist_orig = varlist_orig,
        varlist_tr = varlist_ef,type="quantile",ncat = ncat, nobs=ncell)
    samp_data_ef <- collapseMG(samp_data_ef,varlist_tr = varlist_ef,
        varlist_trc=varlist_trc,group="group",nmin = ncell)
    
    fit_mg_ef <- multipleGroup(data = samp_data_ef[,varlist], model = model_mirt,
        group = samp_data_ef[,"group"],
        invariance = c(varlist[c(1)],"free_means","free_vars"))
    # summary(fit_mg_ef)
    # coef(fit_mg_ef)
    # 
    # plot(fit_mg_ef, type = "info")
    
    info_ef <- infoMG(mg_mod_obj = fit_mg_ef,group_list = group_labels)
    
    
    # info equal interval
    
    if(append_collapse){
        varlist_trc <- varlist_efc
        varlist <- varlist_efc
    } else {
        varlist_trc <- NULL
        varlist <- varlist_ef
    }
    samp_data_ei <- recodeOrdinal(df = samp_data,varlist_orig = varlist_orig,
        varlist_tr = varlist_ei,type="interval",ncat = ncat, nobs=ncell)
    samp_data_ei <- collapseMG(samp_data_ei,varlist_tr = varlist_ei,
        varlist_trc=varlist_trc,group="group",nmin = ncell)
    
    fit_mg_ei <- multipleGroup(data = samp_data_ei[,varlist], model = model_mirt,
        group = samp_data_ei[,"group"],
        invariance = c(varlist[c(1)],"free_means","free_vars"))
    # summary(fit_mg_ei)
    # coef(fit_mg_ei)
    # 
    # plot(fit_mg_ei, type = "info")
    
    info_ei <- infoMG(mg_mod_obj = fit_mg_ei,group_list = group_labels)
    
    info_mg <- bind_rows(
        (info_ef %>% mutate(
            label = paste0(label,"_EF")
        )),
        (info_ei %>% mutate(
            label = paste0(label,"_EI")
        ))
    )
    
    return(info_mg)
}


#   simInfoMGMult() is a function to control multiple draws from a population 
#   data using simulateInfoMG() and compiles results from all draws

simInfoMGMult <- function(
        pop_data_all = pop_data_all,
        nsim = 100,
        n_items = n_items,
        nsamp = nsamp,
        seed = 01012025,
        ncat = 10,
        ncell = 10,
        append_collapse = TRUE,
        group_labels = c("Reference","Focal")) {
    
    itm_nms <- names(pop_data_all)[grepl("^x",names(pop_data_all))]
    mod <- paste0("cog = 1-",n_items)
    Theta <- matrix(seq(-6,6,.01))
    i <- 1
    
    while (i <= nsim) {
        tryCatch({
            set.seed(seed)
            itnms <- sample(itm_nms, size = n_items, replace = FALSE)
            info1 <- simulateInfoMG(
                pop_data_all = pop_data_all[,c(itnms,"group")],
                itnms = itnms,
                n_items = n_items,
                nsamp = nsamp,
                seed = seed)
            info1 <- info1 %>% relocate(info,.after = label)
            names(info1)[3] <- paste0(names(info1)[3],"_i",i)
            if (i == 1){
                info <- info1
            } else {
                info <- info %>% left_join(info1, by = c("Theta","label"))
            }
            cat(paste0("Simulated Data Set - ",i,"\n"))
            i <- i + 1
            seed <- seed + 1
        }, error=function(e){i <- i})
    }
    
    info$info_mn <- apply(info[,grepl("info_i",names(info))],1, function(x) mean(x))
    info$info_sd <- apply(info[,grepl("info_i",names(info))],1, function(x) sd(x))
    
    return(info)        
}

#   -------------------- End Test Information Functions ------------------------

#   ******************************** Example Code ******************************

# # with external z variable
# sim <- simNonNormal(N=5000,nitems=30,seed=12725,focal_mn=-0.90,focal_sd=0.75,
#     focal_mn_z = 0, focal_sd_z = 1)
# s <- sim[["data"]] %>% mutate(
#     group = case_when(
#         row_number() <= n()/2 ~ "Reference",
#         TRUE ~ "Focal"
#     )
# )
# sim75 <- simNonNormal(N=50000,nitems=50,seed=12725,focal_mn=-0.90,focal_sd=0.75,
#     focal_mn_z = 0, focal_sd_z = 1)
# s75 <- sim75[["data"]] %>% mutate(
#     group = case_when(
#         row_number() <= n()/2 ~ "Reference",
#         TRUE ~ "Focal"
#     )
# )
# sim100 <- simNonNormal(N=50000,nitems=50,seed=12725,focal_mn=-0.90,focal_sd=01.0,
#     focal_mn_z = 0, focal_sd_z = 1)
# s100 <- sim100[["data"]] %>% mutate(
#     group = case_when(
#         row_number() <= n()/2 ~ "Reference",
#         TRUE ~ "Focal"
#     )
# )
# sim125 <- simNonNormal(N=50000,nitems=50,seed=12725,focal_mn=-0.90,focal_sd=01.25,
#     focal_mn_z = 0, focal_sd_z = 1)
# s125 <- sim125[["data"]] %>% mutate(
#     group = case_when(
#         row_number() <= n()/2 ~ "Reference",
#         TRUE ~ "Focal"
#     )
# )
# # saveRDS(sim,file="Analysis/Results/simulation_population_data.rds")
# # sim <- readRDS("Analysis/Results/simulation_population_data.rds")
# 
# # without external z
# sim <- simNonNormal(N=5000,nitems=30,seed=12725,focal_mn=-0.90,focal_sd=0.75)
# s <- sim[["data"]] %>% mutate(
#     group = case_when(
#         row_number() <= n()/2 ~ "Reference",
#         TRUE ~ "Focal"
#     )
# )
# 
# summary(s[1:5000,])
# summary(s[5001:10000,])
# cor(s[1:5000,!names(s) %in% "group"])
# cor(s[5001:10000,!names(s) %in% "group"])
# sim[["loadings"]]
# sim[["intercepts"]]
# sim[["skewness"]]
# sim[["skewnessf"]]
# sim[["kurtosis"]]
# sim[["kurtosisf"]]
# 
# rm(s,sim)
# 
# simr <- simNonNormal(N=5000,nitems=30,seed=12725,focal_mn=0.45,focal_sd=1.0,
#     focal_mn_z = 0, focal_sd_z = 1)
# sr <- simr[["data"]] %>% mutate(
#     group = case_when(
#         row_number() <= n()/2 ~ "Reference",
#         TRUE ~ "Focal"
#     )
# )
# 
# simf <- simNonNormal(N=5000,nitems=30,seed=12725,focal_mn= -0.45,focal_sd=0.75,
#     focal_mn_z = 0, focal_sd_z = 1)
# sf <- simf[["data"]] %>% mutate(
#     group = case_when(
#         row_number() <= n()/2 ~ "Reference",
#         TRUE ~ "Focal"
#     )
# )
# 
# summary(sr[1:5000,])
# summary(sr[5001:10000,])
# summary(sf[1:5000,])
# summary(sf[5001:10000,])
# 
# summary(sr)
# summary(sf)
# skew(sr[,1:30])
# skew(sf[,1:30])
# 
# cor(sr[,1:31])
# cor(sf[,1:31])
# 
# simr[["loadings"]]
# simr[["intercepts"]]
# simr[["skewness"]]
# simr[["kurtosis"]]
# 
# 
# s <- rbind(sr,sf)
# 
# skew(s[,1:30])
# 
# sim4 <- simNonNormal(N=5000,nitems=30,seed=12725,focal_mn=-0.9,focal_sd=0.75)
# # Reference group cases are the rows 1-N, Focal group N+1-2N
# s4 <- sim4[["data"]] %>% mutate(
#     group = case_when(
#         row_number() <= n()/2 ~ "Reference",
#         TRUE ~ "Focal"
#     )
# )
# 
# cor(data.frame(sim4[["data"]]))
# cor(s4[,!names(s4) %in% "group"])
# cor(s4[s4$group == "Reference",!names(s4) %in% "group"])
# cor(s4[s4$group == "Focal",!names(s4) %in% "group"])
# sim4[["loadings"]]
# sim4[["intercepts"]]
# sim4[["skewness"]]
# sim4[["kurtosis"]]
# 
# skew(s4[s4$group == "Reference",1:30])
# skew(s4[s4$group == "Focal",1:30])
#           
# hist(s4[s4$group == "Reference","x30"],xlim=c(-6,6))
# hist(s4[s4$group == "Focal","x30"],xlim=c(-6,6))
# 
# mean(sim4[["loadings"]])
# mean(sim4[["skewness"]])
# mean(sim4[["intercepts"]])
# hist(sim4[["skewness"]])
# 
# summary(s4[s4$group == "Reference",])
# summary(s4[s4$group == "Focal",])
# skew(s4[,1:30])
# cor(sim4[["skewness"]],skew(s4[,1:30]))
# plot(sim4[["skewness"]],skew(s4[,1:30]))
# 
# n <- 5000
# mean <- sim4[["intercepts"]]
# skew <- sim4[["skewness"]]
# kurt <- sim4[["kurtosis"]]
# sigma <- sim4[["loadings"]] %*% t(sim4[["loadings"]])
# diag(sigma) <- 1
# 
# 
# 
# i4 <- simInfoMGMult(
#     pop_data_all = s,
#     nsim = 25,
#     n_items = 15,
#     nsamp = 1000,
#     seed = 01012025)
# 
# i4 <- i4 %>% mutate(
#     lower = info_mn - (1.96 * info_sd),
#     upper = info_mn + (1.96 * info_sd)
# )
# 
# pl <- ggplot(i4,aes(x = Theta, y = info_mn, color = label)) +
#     geom_line()
# pl
# pl + geom_ribbon(aes(ymin = lower, ymax = upper, fill = label), alpha = 0.2, colour = NA)
# 
# i75 <- simInfoMGMult(
#     pop_data_all = s75,
#     nsim = 50,
#     n_items = 15,
#     nsamp = 1000,
#     seed = 01012025)
# 
# i75 <- i75 %>% mutate(
#     lower = info_mn - (1.96 * info_sd),
#     upper = info_mn + (1.96 * info_sd)
# )
# 
# i100 <- simInfoMGMult(
#     pop_data_all = s100,
#     nsim = 50,
#     n_items = 15,
#     nsamp = 1000,
#     seed = 01012025)
# 
# i100 <- i100 %>% mutate(
#     lower = info_mn - (1.96 * info_sd),
#     upper = info_mn + (1.96 * info_sd)
# )
# 
# i125 <- simInfoMGMult(
#     pop_data_all = s125,
#     nsim = 50,
#     n_items = 15,
#     nsamp = 1000,
#     seed = 01012025)
# 
# i125 <- i125 %>% mutate(
#     lower = info_mn - (1.96 * info_sd),
#     upper = info_mn + (1.96 * info_sd)
# )
# 
# ggplot(i75,aes(x = Theta, y = info_mn, color = label)) +
#     geom_line()
# ggplot(i100,aes(x = Theta, y = info_mn, color = label)) +
#     geom_line()
# ggplot(i125,aes(x = Theta, y = info_mn, color = label)) +
#     geom_line()
# 
# 
# ggplot(s75, aes(x = x3, color = group)) +
#     geom_density()
# ggplot(s100, aes(x = x3, color = group)) +
#     geom_density()
# ggplot(s125, aes(x = x3, color = group)) +
#     geom_density()
# 
# 
# i75a <- simInfoMGMult(
#     pop_data_all = s75,
#     nsim = 20,
#     n_items = 15,
#     nsamp = 2000,
#     seed = 01012025)
# 
# i75a <- i75a %>% mutate(
#     lower = info_mn - (1.96 * info_sd),
#     upper = info_mn + (1.96 * info_sd)
# )
# ggplot(i75a,aes(x = Theta, y = info_mn, color = label)) +
#     geom_line()
# 
# i75b <- simInfoMGMult(
#     pop_data_all = s75,
#     nsim = 20,
#     n_items = 15,
#     nsamp = 500,
#     seed = 01012025)
# 
# i75b <- i75b %>% mutate(
#     lower = info_mn - (1.96 * info_sd),
#     upper = info_mn + (1.96 * info_sd)
# )
# ggplot(i75b,aes(x = Theta, y = info_mn, color = label)) +
#     geom_line()
# 
# i125a <- simInfoMGMult(
#     pop_data_all = s125,
#     nsim = 20,
#     n_items = 15,
#     nsamp = 2000,
#     seed = 01012025)
# 
# i125a <- i125a %>% mutate(
#     lower = info_mn - (1.96 * info_sd),
#     upper = info_mn + (1.96 * info_sd)
# )
# 
# ggplot(i125a,aes(x = Theta, y = info_mn, color = label)) +
#     geom_line()
# 
# 
# i125b <- simInfoMGMult(
#     pop_data_all = s125,
#     nsim = 20,
#     n_items = 15,
#     nsamp = 500,
#     seed = 01012025)
# 
# i125b <- i125b %>% mutate(
#     lower = info_mn - (1.96 * info_sd),
#     upper = info_mn + (1.96 * info_sd)
# )
# 
# ggplot(i125b,aes(x = Theta, y = info_mn, color = label)) +
#     geom_line()
# 



#   --------------------------------- End Example ------------------------------

