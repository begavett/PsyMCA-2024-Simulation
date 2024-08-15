require(sqldf)
require(dplyr)

# Blom formula
# y <- qnorm((r-c)/(N-2c+1))
# where qnorm is the standard normal quantile function
# r is the rank of your variable
# c = 3/8, the Blom constant
# N is the sample size

# ================================ blom transformation ========================================

# recodeBlom recodes original variables using a Blom rank based inverse normal transformation.
# recodeBlom parameters:
#   df - label for the data frame that contains the variables to be recoded.
#       Can be dataframe or dataframe name (in quotes)
#   varlist_orig - List of labels for the original variables to be recoded
#   varlist_tr - List of labels for the recoded variables

recodeBlom <- function(df,varlist_orig,varlist_tr){
    if (is.data.frame(df)) {
        tr <- df
    } else {
        tr <- eval(parse(text=df))
    }
  for (j in 1:length(varlist_orig)){
      varlab <- varlist_tr[j]
      tr[,varlab] <- ifelse(is.na(tr[,varlist_orig[j]]),NA,qnorm((rank(tr[,varlist_orig[j]], ties="average")-0.375)/(sum(!is.na(tr[,varlist_orig[j]])) - 2*0.375 + 1))) 
  }
  return(tr)
}

# ------------------------------------- end blom ----------------------------------------------


# ================================ recode into ordinal variables ==============================

# ------------------------------- Combined Interval and Quantile ------------------------------

# recodeOrdinal recodes original variables into ordinal variables using either equal-quantile 
#   or equal-interval recoding. Equal-quantile divides the original variable into k categories
#   based on k quantiles; the transformed variable will have a uniform distribution. Equal-interval
#   divides the original variable into k categories with the same range the original variable units.
#   Equal-interval will start with the number of specified categories (ncat) but will iteratively
#   decrease this number until the minimum number of observations (nobs) is achieved in each category.
#   Equal-quantile is appropriate for continuous measures where the sample size divided by the
#   number of categories is greater than the minimum number of observations desired for each response
#   category. Equal-quantile does not work well with ordinal scales that have highly skewed
#   distributions; equal-interval is recommended. It is important to check the distributions
#   of the recoded variables to assure that the recoding matches the purpose and goals for recoding.
# recodeOrdinal parameters:
#   df - label for the data frame that contains the variables to be recoded.
#       Can be dataframe or dataframe name (in quotes).
#   varlist_orig - List of labels for the original variables to be recoded
#   varlist_tr - List of labels for the recoded variables
#   type - "interval" for equal-interval recoding, "quantile" for equal quantile
#   ncat - maximum number of categories for recoded variables
#   nobs - miminum number of observed responses in each category of the transformed variable

recodeOrdinal <- function(df,varlist_orig,varlist_tr,type="interval",ncat=10, nobs=10) {
  if (is.data.frame(df)) {
      tro <- df
  } else {
      tro <- eval(parse(text=df))
  }
  for (j in 1:length(varlist_orig)){
    cuts <- {}
    min <- min(tro[,varlist_orig[j]], na.rm=TRUE)
    max <- max(tro[,varlist_orig[j]], na.rm=TRUE)
    ncat1 <- ncat
    if (type == "interval"){
      #  Trims tails by recoding the min and max values to the value
      #  below or above which there are cumulative nobs observations
      for (i in 1:(ncat1+1)){
        cuts[i] <- min + ((i-1) * ((max - min)/ncat1))
      }
      ordt <- cut(tro[,varlist_orig[j]], breaks=cuts, include.lowest=TRUE,
                  labels=c(1:ncat1))
      to <- as.data.frame(table(ordt)) 
      to$ordt <- as.numeric(as.character(to$ordt))
      to$cutpt <- cuts[2:length(cuts)]
      frstrow <- c(0,0,cuts[1])
      to <- rbind(frstrow,to)
      cumn <- 0
      for (k in nrow(to):1){
        cumn <- cumn + to[k,"Freq"]
        if (cumn >= nobs){
          high <- to[k-1,"cutpt"]
          break
        }
      }
      cumn <- 0
      for (k in 1:nrow(to)){
        cumn <- cumn + to[k,"Freq"]
        if (cumn >= nobs){
          low <- to[k-1,"cutpt"]
          break
        }
      }
      tvar <- paste(varlist_orig[j],"t",sep="")
      tro[,tvar] <- ifelse(tro[,varlist_orig[j]] > high,high,
                           tro[,varlist_orig[j]])       
      tro[,tvar] <- ifelse(tro[,tvar] < low,low,
                           tro[,tvar])       
      #  ------------------end trim tails-----------------         
      cuts <- {}
      min <- min(tro[,tvar], na.rm=TRUE)
      max <- max(tro[,tvar], na.rm=TRUE)
      go <- TRUE
      repeat{
        for (i in 1:(ncat1+1)){
          cuts[i] <- min + ((i-1) * ((max - min)/ncat1))
        }
        ordt <- cut(tro[,tvar], breaks=cuts, include.lowest=TRUE,
                    labels=c(1:ncat1))
        #                 ordt <- cut(tro[,varlist_orig[j]], breaks=cuts, include.lowest=TRUE,
        #                       labels=c(1:ncat1))
        if (min(table(ordt)) >= nobs | ncat1 == 2) {go <- FALSE}
        if (go==FALSE) {break}
        ncat1 <- ncat1-1
        cuts <- {}
      }
      varlab <- varlist_tr[j]
      tro[,varlab] <- as.numeric(cut(tro[,tvar], breaks=cuts, 
          include.lowest=TRUE))
      tro <- tro[,!(names(tro) %in% tvar)]
    } else {
      for (i in 1:(ncat1+1)){
        if (i==1){
          cuts[i] <- min
        } else{
          cuts[i] <- quantile(tro[,varlist_orig[j]],(i-1) * (1.0/ncat1), na.rm=TRUE)
        }  
      }
      cuts <- unique(cuts)
      varlab <- varlist_tr[j]
      tro[,varlab] <- as.numeric(cut(tro[,varlist_orig[j]], breaks=cuts, 
          include.lowest=TRUE))
    }
  }
  return(tro)
}

# -------------------------------- end ordinal recode --------------------------------------- 


# =================================== Recode Lookup ========================================= 

# recodeLookup generates a reference table from recoded variables within a dataframe and uses
#   this as a lookup table to recode values that have not yet been recoded. This is applicable
#   where a subset of records are used for the original recode and this subset is then used 
#   as the reference to recode the other records. It is applicable for continuous transformed variables
#   (including but not limited to Blom) and for ordinal transformations.
# recodeLookup parameters:
#   df - label for the data frame that contains the variables to be recoded.
#       Can be dataframe or dataframe name (in quotes).
#   varlist_orig - List of labels for the original variables to be recoded
#   varlist_tr - List of labels for the recoded variables
#   type - "continuous" for continuous (numeric) transformation, "ordinal" for ordinal
#   lookup - NULL if the input table (df) will be used to lookup recode values.
#       Can be dataframe or dataframe name (in quotes). 
#     that contains the variables to be recoded (in quotes) if an external table is used. The external
#     table has to be opened as a dataframe, and must contain all of the columns 
#     in varlist_orig and varlist_tr.
#   lu_type = "lookup" for lookup table. For ordinal level recodes only. Lookup table
#     must have column "recode_score" that represents the ordinal recode values, 
#     "min_<varnm>", and "max_<varnm>" for each original variable (<varnm>). This
#     table is generated by createLookupTable. "data" is default and refers to data
#     frame that contains original and transformed variables.
#     

#   df <- "npsy1"
#   varlist_orig <- blomvarin
#   varlist_tr <- blomvarout
#   type <- "continuous"
#   lookup <- NULL
#   lookup <- npsy1
#   lu_type="data"

recodeLookup <- function(df,varlist_orig,varlist_tr,type="continuous",lookup=NULL,
      lu_type="data") {
  require(sqldf)
  if (is.data.frame(df)){
    rcd <- df
  } else {
    rcd <- eval(parse(text=df))
  }
  if (is.null(lookup)){
    rcdlu <- rcd[,c(varlist_orig,varlist_tr)]
  } else{
    if (is.data.frame(lookup)) {
      if (lu_type == "lookup") {
        rcdlu <- lookup
      } else {
        rcdlu <- lookup[,c(varlist_orig,varlist_tr)]
      }
    } else {
      if (lu_type == "lookup") {
        rcdlu <- eval(parse(text=lookup))
      } else {
        rcdlu <- eval(parse(text=lookup))
        rcdlu <- rcdlu[,c(varlist_orig,varlist_tr)]
      }
    }
  }
  for (j in 1:length(varlist_orig)){
      if (lu_type == "lookup") {
      luv <- rcdlu[,c("recode_score",paste0("min_",varlist_orig[j]),paste0("max_",varlist_orig[j]))]
      names(luv) <- sub("recode_score","tr_min",names(luv))
      luv$tr_max <- luv$tr_min
      names(luv) <- c("tr_min","orig_min","orig_max","tr_max")
      luv <- luv[!is.na(luv$orig_min),]
      for (i in 1:nrow(luv)) {
        if (i < nrow(luv)) {
          luv[i,"orig_max"] <- luv[i+1,"orig_min"]
        } 
      }
      # luv$orig_min <- round(luv$orig_min,4)
      # luv$orig_max <- round(luv$orig_max,4)
    } else {
      t5 <- unique(rcdlu[!is.na(rcdlu[,varlist_tr[j]]),c(varlist_orig[j],varlist_tr[j])])
      t5 <- t5 %>% arrange(.data[[varlist_orig[j]]])
      # t5 <- t5[order(t5[varlist_orig[j]]),]
      luv <- as.data.frame(cbind(t5,rbind(t5[2:nrow(t5),],c(NA,NA))))
      colnames(luv) <- c("orig_min","tr_min","orig_max","tr_max")
    }
    mino <- min(luv[,"orig_min"], na.rm=TRUE)
    mint <- min(luv[,"tr_min"], na.rm=TRUE)
    maxo <- max(luv[,"orig_max"], na.rm=TRUE)
    maxt <- max(luv[,"tr_max"], na.rm=TRUE)
    
    rcd[,varlist_tr[j]] <- ifelse(rcd[,varlist_orig[j]] <= mino,mint,NA)
    rcd[,varlist_tr[j]] <- ifelse(rcd[,varlist_orig[j]] >= maxo,maxt,rcd[,varlist_tr[j]])
    t3 <- rcd[,c(varlist_orig[j],varlist_tr[j])]
    
    if (lu_type == "lookup") {
      sqlcd <- paste("SELECT * FROM t3 AS t3
                       LEFT JOIN
                       (SELECT * FROM luv)
                       AS luv1 ON t3.",varlist_orig[j],
                     " >= luv1.orig_min AND t3.", varlist_orig[j],
                     " < luv1.orig_max",sep="")
    } else {
      sqlcd <- paste("SELECT * FROM t3 AS t3
                       LEFT JOIN
                       (SELECT * FROM luv)
                       AS luv1 ON t3.",varlist_orig[j],
                     " >= luv1.orig_min AND t3.", varlist_orig[j],
                     " < luv1.orig_max",sep="")
    }
    options(sqldf.driver = "SQLite")
    t4 <- sqldf(sqlcd)
    if (type == "continuous") {
      t4[,varlist_tr[j]] <- ifelse(!is.na(t4[,varlist_tr[j]]),t4[,varlist_tr[j]],
          ifelse(t4[,varlist_orig[j]] >= maxo,maxt,
            ( ( (t4[,varlist_orig[j]] - t4$orig_min) / (t4$orig_max - t4$orig_min) ) *
              (t4$tr_max - t4$tr_min) ) + t4$tr_min))
      #       t4[,varlist_tr[j]] <- ifelse(!is.na(t4[,varlist_tr[j]]),t4[,varlist_tr[j]],( ( (t4[,varlist_orig[j]] - t4$orig_min) / (t4$orig_max - t4$orig_min) ) * (t4$tr_max - t4$tr_min) ) + t4$tr_min)
    } else {
      t4[,varlist_tr[j]] <- ifelse(!is.na(t4[,varlist_tr[j]]),t4[,varlist_tr[j]],
        ifelse(t4[,varlist_orig[j]] >= max(t4$orig_max,na.rm=TRUE),max(t4$tr_max,na.rm=TRUE),
        t4$tr_min))
    }
    rcd[,varlist_tr[j]] <- t4[,varlist_tr[j]]
  }
  return(rcd) 
}

# t4[t4[,varlist_orig[j]] %in% 328,]
# 
#   t4[!is.na(t4[,varlist_orig[j]]) & t4[,varlist_orig[j]] > 211,]
# summary(t4$tau)
                                                      
# ------------------------------------ end recode lookup --------------------------------------


# =============================createLookupTable ====================================

#   createLookupTable generates a reference table from recoded variables within a dataframe that
#   can be used as a lookup table to recode vales that have not yet been recoded. 
#   This is applicable for recoded variables that take on integer values.
#   createLookupTable parameters:
#   df - label for the data frame that contains the variables to be recoded.
#       Can be dataframe or dataframe name (in quotes).
#   varlist_orig - List of labels for the original variables that were recoded
#   varlist_tr - List of labels for the recoded variables
#     the input table (df) has to be opened as a dataframe, and must contain  
#     all of the columns in varlist_orig and varlist_tr.
#   lookup - NULL if the input table (df) will be used to lookup recode values.
#       Can be dataframe or dataframe name (in quotes). 
#   This function returns the look-up reference table in the form of a data frame.

createLookupTable <- function(df,varlist_orig,varlist_tr,lookup=NULL) {
    if (is.data.frame(df)) {
        rcd <- df
    } else {
        rcd <- eval(parse(text=df))
    }
    if (is.null(lookup)){
        rcdlu <- rcd[,c(varlist_orig,varlist_tr)]
    } else{
        if (is.data.frame(lookup)) {
            rcdlu <- lookup[,c(varlist_orig,varlist_tr)]
        } else {
            rcdlu <- eval(parse(text=lookup))
            rcdlu <- rcdlu[,c(varlist_orig,varlist_tr)]
        }
    }
  for (j in 1:length(varlist_orig)){
    t5 <- unique(rcdlu[!is.na(rcdlu[,varlist_tr[j]]),c(varlist_orig[j],varlist_tr[j])])
    t5 <- t5 %>% arrange(.data[[varlist_orig[j]]])
    # t5 <- t5[order(t5[varlist_orig[j]]),]
    colnames(t5) <- c("orig","recode_score")
    t6 <- t5 %>% group_by(recode_score) %>%
      summarise(
        min = min(orig),
        max = max(orig)
      )
    colnames(t6) <- gsub("min",paste("min_",varlist_orig[j],sep=""),colnames(t6))
    colnames(t6) <- gsub("max",paste("max_",varlist_orig[j],sep=""),colnames(t6))
    if (j == 1){
      luv <- t6
    } else {
      luv <- merge(luv,t6,by="recode_score",all.x=TRUE,all.y=TRUE)
    }
  }
  return(luv) 
}

# --------------------------end createLookupTable -----------------------------------
