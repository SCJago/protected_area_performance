# Load libraries #
library("foreign")
library("dplyr")
library("xlsx")
library("ggplot2")
library("broom")
library("MatchIt")
library("ggrepel")
library("calibrate")
library("plm")
library("gridExtra")
library("tidyr")
library("foreach")
library("Rmisc")
library("lmtest") #for displaying coef. table
library("sandwich") #vcovCL
library("marginaleffects")


# --------------------- 1) make the function -----------------------------------------#

the_cov_parameters_loop <- function(group, output, distance, replacement, caliper_value) { 
  # Define covariates
  covariates_min <- names(group[,c("Elevation", "Slope", "Temperature", "Precipitation", "Ethnicity1", "Agri_suitability")]) 
  covariates_add <- names(group[,c("Access", "Forest", "Population", "Grassland", "Agriculture")])    
  
  param_combinations <- expand.grid(distance = distance, replacement = replacement, caliper_value = caliper_value)
  
  id <- unlist(
    lapply(1:length(covariates_add),
           function(i) combn(1:length(covariates_add), i, simplify=FALSE)
    ),
    recursive=FALSE
  )
  
  set.seed(12345)
  
  mydbs_full <- lapply(id, function(i) {
    lapply(1:nrow(param_combinations), function(j) {
      c(as.character(param_combinations[j, "distance"]),
        param_combinations[j, "replacement"],
        param_combinations[j, "caliper_value"],
        covariates_min,
        covariates_add[i])
    })
  })
  mydbs_full <- unlist(mydbs_full, recursive = FALSE)
  
  maha <- lapply(mydbs_full, function(cov) {
    cali <- ifelse(cov[1] == "mahalanobis", 0, as.numeric(cov[3]))
    
    data <- group[, c("Elevation", "Slope", "Temperature", "Precipitation", "Ethnicity1", "Agri_suitability", 
                      "Access", "Forest", "Population", "Grassland", "Agriculture", "Type", "Forest_change", 
                      "Grass_change", "Agri_change", "Land", "Ecoregion")]
    
    treated <- group$Type
    treated <- as.factor(treated)
    
    idx <- match(cov, names(data))
    idx <- sort(idx)
    NewDF <- data[, c(idx)]
    
    formulae  <- as.formula(paste("treated ~ ", paste0(names(NewDF), collapse = "+")))
    
    m.out.test <- matchit(formulae, method = "nearest", data = cbind(treated, data), ratio = 1)
    test_data_S <- match.data(m.out.test)
    ps.sd_S <- sd(test_data_S$distance)
    cal <- (as.numeric(cali)) * ps.sd_S 
    
    COVAR.ex <- c("Land", "Ecoregion")
    
    eval_code <- TRUE
    
    m.out <- tryCatch(
      {
        if (as.character(cov[1]) == "glm") {
          matchit(formulae, data = cbind(treated, data), 
                  method = "nearest", 
                  distance = as.character(cov[1]), 
                  replace = as.logical(cov[2]),
                  caliper = cal,
                  exact = COVAR.ex)
        } else {
          matchit(formulae, data = cbind(treated, data), 
                  method = "nearest", 
                  distance = as.character(cov[1]), 
                  replace = as.logical(cov[2]),
                  exact = COVAR.ex)
        }
      }, 
      error = function(e) { eval_code <<- FALSE }
    )
    
    if (eval_code) {
      match_output <- summary(m.out, standardize = TRUE)
      test <- as.data.frame(match_output$sum.matched[, 3])
      colnames(test) <- "Stand_Mean_Diff"
      mean_std_diff <- mean(abs(test$Stand_Mean_Diff))
      max_std_diff <- max(abs(test$Stand_Mean_Diff))
      unmatched <- match_output$nn[5, 2]
      cov_balance <- cbind(unmatched, mean_std_diff, max_std_diff)
      
      m.data <- match.data(m.out)
      m.data$treated <- as.factor(m.data$treated)  # Ensure treated is factor
      
      # Regression formula including covariates and exact match vars
      reg_formula <- as.formula(paste(output, "~ treated +", paste(names(NewDF), collapse = " + "), "+ Land + Ecoregion"))
      
      # Fit weighted regression
      regression_model <- lm(reg_formula, weights = m.data$weights, data = m.data)
      
      # Estimate ATE using marginaleffects
      ATE_result <- avg_comparisons(regression_model,
                                    variables = "treated",
                                    vcov = ~subclass,
                                    weights = m.data$weights)
      
      ATE <- ATE_result$estimate
      se <- ATE_result$std.error
      sig <- ATE_result$p.value
      z = ATE_result$statistic
      
      var <- as.data.frame(t(covariates_add %in% names(NewDF)))
      names(var) <- covariates_add
      
      m1i <- data.frame(
        z = z,
        sig = sig,
        ATE = ATE,
        se = se,
        cov_balance,
        var,
        cal025 = cali == 0.25, cal05 = cali == 0.5, cal1 = cali == 1, cal0 = cali == 0,
        replacement = as.logical(cov[2]),
        maha = cov[1] == "mahalanobis",
        glm = cov[1] == "glm"
      )
    }
  })
  
  df_name <- paste0(output, "_covparam_synthesis_loop")
  assign(df_name, data.frame(do.call("rbind", maha)), envir = .GlobalEnv)
}


##############################################################################################


the_hh_cov_parameters_loop <- function(group, output, distance, replacement, caliper_value) { 
  # A matrix will all possible covariates, and distinguish between essential and additional
  
  
  # We define two set of covariates
  covariates_min <- names(group[,c("Elevation", "Slope", "Temperature", "Precipitation", "Ethnicity1", "Agri_suitability")]) 
  covariates_add <- names(group[,c("Access", "Population", "Agriculture")])    
  
  # Generate all possible combinations of distance, replacement, and caliper_value
  param_combinations <- expand.grid(distance = distance, replacement = replacement, caliper_value = caliper_value)
  
  
  # We determine how many different models can be constructed based on the inclusion of additional covariates
  id <- unlist(
    lapply(1:length(covariates_add),
           function(i)combn(1:length(covariates_add),i,simplify=FALSE)
    )
    ,recursive=FALSE)
  
  # We create all combinations of required + optional covariates
  
  set.seed(12345)
  
  mydbs_full <- lapply(id, function(i) {
    lapply(1:nrow(param_combinations), function(j) {
      c(as.character(param_combinations[j, "distance"]),
        param_combinations[j, "replacement"],
        param_combinations[j, "caliper_value"],
        covariates_min,
        covariates_add[i])
    })
  })
  
  mydbs_full <- unlist(mydbs_full, recursive = FALSE)
  
  # We loop over all possible models based on any combination of 
  
  maha <- lapply(mydbs_full, function(cov) {
    
    
    # We determine the calipers for each variables
    cali <- ifelse(cov[1] == "mahalanobis", 0, as.numeric(cov[3]))
    
    
    # We only keep the ppol of covaiates, treatment index and outcome
    data <- group[, c("Elevation", "Slope", "Temperature", "Precipitation", "Ethnicity1", "Agri_suitability", 
                      "Access", "Population", "Agriculture", "Type", "MAHFP_change", 
                      "HDDS_change", "Asset_change")]
    

    treated <- as.factor(group$Type)
    
    
    idx <- match(cov, names(data))   # We determines the columns of the covariates
    idx <- sort(idx)               # We need to arrange them in increasing order
    NewDF <- data[,c(idx)]           # We create a database that contains only the covariates needed in this iteration of the loop
    
    formulae  <- as.formula(paste("treated ~ ", paste0(names(NewDF),collapse ="+")))  # We create a formula to be inserted in the matching
    
    # Calipers
    m.out.test=matchit(formulae,method="nearest",data= cbind(treated,data),ratio=1)
    test_data_S=match.data(m.out.test)
    ps.sd_S=sd(test_data_S$distance)
    cal <- (as.numeric(cali))*ps.sd_S 
    
    
    # Applying the matching function to 
    eval_code <- TRUE
    
    m.out <- tryCatch(
      {
        if (as.character(cov[1]) == "glm") {
          matchit(formulae, data = cbind(treated, data), 
                  method = "nearest", 
                  distance = as.character(cov[1]), 
                  replace = as.logical(cov[2]),
                  caliper = cal)
        } else {
          matchit(formulae, data = cbind(treated, data), 
                  method = "nearest", 
                  distance = as.character(cov[1]), 
                  replace = as.logical(cov[2]))
        }
      }, error = function(e) { eval_code <<- FALSE})
    
    if (eval_code) {
      match_output <- summary(m.out, standardize = TRUE)
      test <- as.data.frame(match_output$sum.matched[, 3])
      colnames(test) <- "Stand_Mean_Diff"
      mean_std_diff <- mean(abs(test$Stand_Mean_Diff))
      max_std_diff <- max(abs(test$Stand_Mean_Diff))
      unmatched <- match_output$nn[5, 2]
      cov_balance <- cbind(unmatched, mean_std_diff, max_std_diff)
      
      m.data <- match.data(m.out)
      m.data$treated <- as.factor(m.data$treated)  # Ensure treated is factor
      print(names(NewDF))
      
      # Regression formula including covariates and exact match vars
      reg_formula <- as.formula(paste(output, "~ treated +", paste(names(NewDF), collapse = " + ")))
      
      # Fit weighted regression
      regression_model <- lm(reg_formula, weights = m.data$weights, data = m.data)
      
      # Estimate ATE using marginaleffects
      ATE_result <- avg_comparisons(regression_model,
                                    variables = "treated",
                                    vcov = ~subclass,
                                    weights = m.data$weights)
      ATE <- ATE_result$estimate
      se <- ATE_result$std.error
      sig <- ATE_result$p.value
      z = ATE_result$statistic
      
      var <- as.data.frame(t(covariates_add %in% names(NewDF)))
      names(var) <- covariates_add
      
      m1i <- data.frame(
        z = z,
        sig = sig,
        ATE = ATE,
        se = se,
        cov_balance,
        var,
        cal025 = cali == 0.25, cal05 = cali == 0.5, cal1 = cali == 1, cal0 = cali == 0,
        replacement = as.logical(cov[2]),
        maha = cov[1] == "mahalanobis",
        glm = cov[1] == "glm"
      )
    }
  })
    
  df_name <- paste0(output, "_covparam_synthesis_loop")
  
  assign(df_name, data.frame(do.call("rbind", maha)), envir = .GlobalEnv)
}  


##############################################################################################################


# 1 - Load the functions 
the_cov_parameters_loop
the_hh_cov_parameters_loop

# 2 - bring in pre-matched data
STRICT <- read.csv("S_sample.csv")
STRICT <- subset(STRICT, select = c("Forest_change", "Percent_grass_change", "Percent_agriculturechange", 
                                    "Elevation", "Slope", "Temperature", "Precipitation", "Ethnicity1", 
                                    "Agri_suitability","Access", "Forest", "Population", "Grassland", 
                                    "Agriculture", "Land", "Ecoregion", "type"))
                 
names(STRICT) <- c("Forest_change", "Grass_change", "Agri_change", 
                   "Elevation", "Slope", "Temperature", "Precipitation", "Ethnicity1", 
                   "Agri_suitability","Access", "Forest", "Population", "Grassland", 
                   "Agriculture", "Land", "Ecoregion", "Type")

STRICT$ID <- seq(nrow(STRICT))
STRICT$group[STRICT$Type == 1] <- "S"
STRICT$group[STRICT$Type == 0] <- "Cont"
STRICT$Land <- as.factor(STRICT$Land)
STRICT$Ecoregion <- as.factor(STRICT$Ecoregion)
STRICT <- na.omit(STRICT)

LESS_STRICT <- read.csv("LS_sample.csv")
LESS_STRICT <- subset(LESS_STRICT, select = c("Forest_change", "Percent_grass_change", "Percent_agriculturechange", 
                                              "Elevation", "Slope", "Temperature", "Precipitation", "Ethnicity1", 
                                              "Agri_suitability","Access", "Forest", "Population", "Grassland", 
                                              "Agriculture", "Land", "Ecoregion", "type"))
names(LESS_STRICT) <- c("Forest_change", "Grass_change", "Agri_change", 
                        "Elevation", "Slope", "Temperature", "Precipitation", "Ethnicity1", 
                        "Agri_suitability","Access", "Forest", "Population", "Grassland", 
                        "Agriculture", "Land", "Ecoregion", "Type")
LESS_STRICT$ID <- seq(nrow(LESS_STRICT))
LESS_STRICT$group[LESS_STRICT$Type == 1] <- "LS"
LESS_STRICT$group[LESS_STRICT$Type == 0] <- "Cont"
LESS_STRICT$Land <- as.factor(LESS_STRICT$Land)
LESS_STRICT$Ecoregion <- as.factor(LESS_STRICT$Ecoregion)
LESS_STRICT <- na.omit(LESS_STRICT)

HOUSEHOLD <- read.csv("household_sample.csv")
HOUSEHOLD <- subset(HOUSEHOLD, select = c("MINS_change", "HDDS_change", "asset_change",  
                                          "Elevation", "Slope", "Temperature", "Precipitation", "Ethnicity1", 
                                          "Agri_suitability","Access",  "Population", 
                                          "Agriculture", "Land", "Ecoregion", "type"))
names(HOUSEHOLD) <- c("MAHFP_change", "HDDS_change", "Asset_change", 
                      "Elevation", "Slope", "Temperature", "Precipitation", "Ethnicity1", 
                      "Agri_suitability","Access", "Population", 
                      "Agriculture", "Land", "Ecoregion", "Type")

HOUSEHOLD$ID <- seq(nrow(HOUSEHOLD))
HOUSEHOLD$group[HOUSEHOLD$Type == 1] <- "Treat"
HOUSEHOLD$group[HOUSEHOLD$Type == 0] <- "Cont"
HOUSEHOLD <- na.omit(HOUSEHOLD)

# 2 -----------------  Run functions for each match group and each output----------------------------------------------#

#strict forest
the_cov_parameters_loop(STRICT,
                        "Forest_change",
                        caliper_value = c(0.25,0.5,1), 
                        replacement = c(T,F), 
                        distance = c("mahalanobis","glm"))
S_forest_covparam_synthesis_loop <- Forest_change_covparam_synthesis_loop
write.csv(S_forest_covparam_synthesis_loop,"cov_parameter_loop_forest_S.csv")

#strict grass
the_cov_parameters_loop(STRICT,
                        "Grass_change",
                        caliper_value = c(0.25,0.5,1), 
                        replacement = c(T,F), 
                        distance = c("mahalanobis","glm"))
S_grass_covparam_synthesis_loop <- Grass_change_covparam_synthesis_loop
write.csv(S_grass_covparam_synthesis_loop,"cov_parameter_loop_grass_S.csv")

#strict agri
the_cov_parameters_loop(STRICT,
                        "Agri_change",
                        caliper_value = c(0.25,0.5,1), 
                        replacement = c(T,F), 
                        distance = c("mahalanobis","glm"))
S_agri_covparam_synthesis_loop <- Agri_change_covparam_synthesis_loop
write.csv(S_agri_covparam_synthesis_loop,"cov_parameter_loop_agri_S.csv")

#less strict forest
the_cov_parameters_loop(LESS_STRICT,
                        "Forest_change",
                        caliper_value = c(0.25,0.5,1), 
                        replacement = c(T,F), 
                        distance = c("mahalanobis","glm"))
LS_forest_covparam_synthesis_loop <- Forest_change_covparam_synthesis_loop
write.csv(LS_forest_covparam_synthesis_loop,"cov_parameter_loop_forest_LS.csv")

#less strict forest
the_cov_parameters_loop(LESS_STRICT,
                        "Grass_change",
                        caliper_value = c(0.25,0.5,1), 
                        replacement = c(T,F), 
                        distance = c("mahalanobis","glm"))
LS_grass_covparam_synthesis_loop <- Grass_change_covparam_synthesis_loop
write.csv(LS_grass_covparam_synthesis_loop,"cov_parameter_loop_grass_LS.csv")


#less strict Agri
the_cov_parameters_loop(LESS_STRICT,
                        "Agri_change",
                        caliper_value = c(0.25,0.5,1), 
                        replacement = c(T,F), 
                        distance = c("mahalanobis","glm"))
LS_agri_covparam_synthesis_loop <- Agri_change_covparam_synthesis_loop
write.csv(LS_agri_covparam_synthesis_loop,"cov_parameter_loop_agri_LS.csv")

#household MAHFP
the_hh_cov_parameters_loop(HOUSEHOLD,
                           "MAHFP_change",
                           caliper_value = c(0.25,0.5,1), 
                           replacement = c(T,F), 
                           distance = c("mahalanobis","glm"))
MAHFP_covparam_synthesis_loop <- MAHFP_change_covparam_synthesis_loop
write.csv(MAHFP_covparam_synthesis_loop,"cov_parameter_loop_MAHFP.csv")

#household HDDS
the_hh_cov_parameters_loop(HOUSEHOLD,
                           "HDDS_change",
                           caliper_value = c(0.25,0.5,1), 
                           replacement = c(T,F), 
                           distance = c("mahalanobis","glm"))
HDDS_covparam_synthesis_loop <- HDDS_change_covparam_synthesis_loop
write.csv(HDDS_covparam_synthesis_loop,"cov_parameter_loop_HDDS.csv")

#household asset
the_hh_cov_parameters_loop(HOUSEHOLD,
                           "asset_change",
                           caliper_value = c(0.25,0.5,1), 
                           replacement = c(T,F), 
                           distance = c("mahalanobis","glm"))
Asset_covparam_synthesis_loop <- asset_change_covparam_synthesis_loop
write.csv(Asset_covparam_synthesis_loop,"cov_parameter_loop_asset.csv")

#check results for primary match are same as results in paper

# ---------------------------------------3. PLOT OUTPUTS ----------------------------------------------------------#

################################# schart function ##################################################

schart <- function(data, labels=NA, highlight=NA, highlight2=NA, n=1, index.est=1, index.se=2, index.ci=NA,
                   order="asis", ci=.95, ylim=NA, axes=T, heights=c(1,1), leftmargin=11, offset=c(0,0), ylab="ATE", lwd.border=1,
                   lwd.est=4, pch.est=21, lwd.symbol=2, ref=0, lwd.ref=1, lty.ref=2, col.ref="black", band.ref=NA, col.band.ref=NA,length=0,
                   col.est=c("grey60", "red3"), col.est2=c("grey80","lightcoral"), bg.est=c("white", "white"),
                   col.dot=c("grey60","grey95","grey95","black"),
                   bg.dot=c("grey60","grey95","grey95","white"),
                   pch.dot=c(22,22,22,22), fonts=c(2,1), adj=c(1,1),cex=c(1,1)) {
  
  # Authors: Ariel Ortiz-Bobea (ao332@cornell.edu).
  # Version: March 10, 2020
  # If you like this function and use it, please send me a note. It might motivate
  # me to imporove it or write new ones to share.
  
  # Description of arguments
  
  # Data:
  # data: data.frame with data, ideally with columns 1-2 with coef and SE, then logical variables.
  # labels: list of labels by group. Can also be a character vector if no groups. Default is rownames of data.
  # index.est: numeric indicating position of the coefficient column.
  # index.se: numeric indicating position of the SE column.
  # index.ci: numeric vector indicating position of low-high bars for SE. Can take up to 2 CI, so vector can be up to length 4
  
  # Arrangement and basic setup:
  # highlight: numeric indicating position(s) of models (row) to highlight in original dataframe.
  # n: size of model grouping. n=1 removes groupings. A vector yields arbitrary groupings.
  # order: whether models should be sorted or not. Options: "asis", "increasing", "decreasing"
  # ci: numeric indicating level(s) of confidence. 2 values can be indicated.
  # ylim: if one wants to set an arbitrary range for Y-axis
  
  # Figure layout:
  # heights: Ratio of top/bottom panel. Default is c(1,1) for 1 50/50 split
  # leftmargin: amount of space on the left margin
  # offset: vector of numeric with offset for the group and specific labels
  # ylab: Label on the y-axis of top panel. Default is "Coefficient"
  # lwd.border: width of border and other lines
  
  # Line and symbol styles and colors:
  # lwd.est: numeric indicating the width of lines in the top panel
  # ref: numeric vector indicating horizontal reference lines(s) Default is 0.
  # lty.ref: Style of reference lines. Default is dash line (lty=2).
  # lwd.ref. Width of reference lines. Default is 1.
  # col.ref: vector of colors of reference lines. Default is black.
  # band.ref: vector of 2 numerics indicating upper abdn lower height for a band
  # col.band.ref: color of this band
  # col.est: vector of 2 colors indicating for "other" and "highlighted" models
  # col.est2: same for outer confidence interval if more than 1 confidence interval
  # col.dot: vector of 4 colors indicating colors for borders of symbol in bottom panel for "yes", "no", "NA", and "yes for highlighted model"
  # bg.dot : vector of 4 colors indicating colors for background of symbol in bottom panel for "yes", "no", "NA", and "yes for highlighted model"
  # pch.dot: style of symbols in bottom panel for "yes", "no", "NA", and "yes for highlighted model"
  # length: length of the upper notch on th vertical lines. default is 0.
  
  # Letter styles
  # fonts: numeric vector indicating font type for group (first) and other labels (second) (e.g. 1:normal, 2:bold, 3:italic)
  # adj: numeric vector indicating alignment adjustment for text label: 0 is left, .5 is center, 1 is right.
  # cex: numeric vector for size of fonts for top panel (first) and bottom panel (Second)
  
  # 1. Set up
  if (T) {
    # Arrange data
    d <- data
    rownames(d) <- 1:nrow(d)
    
    # Create ordering vector
    if (order=="asis")       o <- 1:length(d[,index.est])
    if (order=="increasing") o <- order(d[,index.est])
    if (order=="decreasing") o <- order(-d[,index.est])
    if (!is.numeric(d[,index.est])) {warning("index.est does not point to a numeric vector.") ; break}
    d <- d[o,]
    est <- d[,index.est] # Estimate
    if (length(index.ci)>1) {
      l1 <- d[,index.ci[1]]
      h1 <- d[,index.ci[2]]
      if (length(index.ci)>2) {
        l2 <- d[,index.ci[3]]
        h2 <- d[,index.ci[4]]
      }
    } else {
      if (!is.numeric(d[,index.se]))  {warning("index.se does not point to a numeric vector.") ; break}
      se  <- d[,index.se] # Std error
      ci <- sort(ci)
      a <- qnorm(1-(1-ci)/2)
      l1 <- est - a[1]*se
      h1 <- est + a[1]*se
      if (length(ci)>1) {
        l2 <- est - a[2]*se
        h2 <- est + a[2]*se
      }
    }
    
    # Table
    if (length(index.ci)>1) remove.index <- c(index.est,index.ci) else remove.index <- c(index.est,index.se)
    remove.index <- remove.index[!is.na(remove.index)]
    tab <- t(d[,-remove.index]) # get only the relevant info for bottom panel
    if (!is.list(labels) & !is.character(labels)) labels <- rownames(tab)
    
    # Double check we have enough labels
    if ( nrow(tab) != length(unlist(labels))) {
      print("Warning: number of labels don't match number of models.")
      labels <- rownames(tab)
    }
    
    # Plotting objects
    xs <- 1:nrow(d) # the Xs for bars and dots
    if (n[1]>1 & length(n)==1) xs <- xs + ceiling(seq_along(xs)/n) - 1 # group models by n
    if (length(n)>1) {
      if (sum(n) != nrow(d) ) {
        warning("Group sizes don't add up.")
      } else {
        idx <- unlist(lapply(1:length(n), function(i) rep(i,n[i])))
        xs <- xs + idx - 1
      }
    }
    h <- nrow(tab) + ifelse(is.list(labels),length(labels),0) # number of rows in table
    # Location of data and labels
    if (is.list(labels)) {
      index <- unlist(lapply(1:length(labels), function(i) rep(i, length(labels[[i]])) ))
      locs <- split(1:length(index),index)
      locs <- lapply(unique(index), function(i) {
        x <- locs[[i]]+i-1
        x <- c(x,max(x)+1)
      })
      yloc  <- unlist(lapply(locs, function(i) i[-1])) # rows where data points are located
      yloc2 <- sapply(locs, function(i) i[1]) # rows where group lables are located
    } else {
      yloc <- 1:length(labels)
    }
    
    # Range
    if (is.na(ylim[1]) | length(ylim)!=2) {
      if (length(index.ci)>2 | length(ci)>1) {
        ylim <- range(c(l2,h2,ref)) # range that includes reference lines
      } else {
        ylim <- range(c(l1,h1,ref))
      }
      ylim <- ylim + diff(ylim)/10*c(-1,1) # and a bit more
    }
    xlim <- range(xs) #+ c(1,-1)
  }
  
  # 2. Plot
  if (T) {
    #par(mfrow=c(2,1), mar=c(0,leftmargin,0,0), oma=oma, xpd=F, family=family)
    layout(t(t(2:1)), height=heights, widths=1)
    par(mar=c(0,leftmargin,0,0), xpd=F)
    
    # Bottom panel (plotted first)
    plot(est, xlab="", ylab="", axes=F, type="n", ylim=c(h,1), xlim=xlim)
    lapply(1:nrow(tab), function(i) {
      # Get colors and point type
      type <- ifelse(is.na(tab[i,]),3,ifelse(tab[i,]==TRUE,1, ifelse(tab[i,]==FALSE,2,NA)))
      type <- ifelse(names(type) %in% paste(highlight2) & type==1,4,type) # replace colors for baseline model
      col <- col.dot[type]
      bg  <- bg.dot[type]
      pch <- as.numeric(pch.dot[type])
      sel <- is.na(pch)
      # Plot points
      points(xs, rep(yloc[i],length(xs)), col=col, bg=bg, pch=pch, lwd=lwd.symbol)
      points(xs[sel], rep(yloc[i],length(xs))[sel], col=col[sel], bg=bg[sel], pch=pch.dot[3]) # symbol for missing value
      
    })
    par(xpd=T)
    if (is.list(labels)) text(-offset[1], yloc2, labels=names(labels), adj=adj[1], font=fonts[1], cex=cex[2])
    # Does not accomodate subscripts
    text(-rev(offset)[1], yloc , labels=unlist(labels), adj=rev(adj)[1], font=fonts[2], cex=cex[2])
    # Accomodates subscripts at the end of each string
    if (F) {
      labels1 <- unlist(labels)
      lapply(1:length(labels1), function(i) {
        a  <- labels1[i]
        a1 <- strsplit(a,"\\[|\\]")[[1]][1]
        a2 <- rev(strsplit(a,"\\[|\\]")[[1]])[1]
        if (identical(a1,a2))  a2 <- NULL
        text(-rev(offset)[1], yloc[i], labels=bquote(.(a1)[.(a2)]), adj=adj[2], font=fonts[2], cex=cex[2])
      })
    }
    par(xpd=F)
    
    # Top panel (plotted second)
    colvec  <- ifelse(colnames(tab) %in% paste(highlight), col.est[2], col.est[1])
    colvec  <- ifelse(colnames(tab) %in% paste(highlight2), col.est[3], colvec)
    bg.colvec  <- ifelse(colnames(tab) %in% paste(highlight), bg.est[2], bg.est[1])
    colvec2 <- ifelse(colnames(tab) %in% paste(highlight),col.est2[2], col.est2[1])
    colvec2 <- ifelse(colnames(tab) %in% paste(highlight2),col.est2[3], colvec2)
    plot(est, xlab="", ylab="", axes=F, type="n", ylim=ylim, xlim=xlim)
    # Band if present
    if (!is.na(band.ref[1])) {
      rect(min(xlim)-diff(xlim)/10, band.ref[1], max(xlim)+diff(xlim)/10, band.ref[2], col=col.band.ref, border=NA)
    }
    # Reference lines
    abline(h=ref, lty=lty.ref, lwd=lwd.ref, col=col.ref)
    # Vertical bars
    if (length(ci)>1 | length(index.ci)>2) arrows(x0=xs, y0=l2, x1=xs, y1=h2, length=length, code=3, lwd=rev(lwd.est)[1], col=colvec2, angle=90)
    arrows(x0=xs, y0=l1, x1=xs, y1=h1, length=length, code=3, lwd=lwd.est[1]     , col=colvec, angle=90)
    points(xs, est, pch=pch.est, lwd=lwd.symbol, col=colvec, bg=bg.colvec)
    # Axes
    if (axes) {
      axis(2, las=2, cex.axis=cex[1], lwd=lwd.border)
      axis(4, labels=NA, lwd=lwd.border)
    }
    mtext(ylab, side=2, line=3.5, cex=cex[1])
    box(lwd=lwd.border)
    
  }
  
}

#########################################################################################################


labels <- list("Balance:" = c(">75% obs. unmatched", "Mean bal. achieved", "Balanced for all cov."),
               "Significance" = c("p value < 0.05"),
               "Additional covariates:" = c("Access", "Forest", "Population", "Agriculture", "Grassland"),
               "Model parameters:" = c("0.25SD caliper","0.5SD caliper", "1SD caliper", "No caliper", "With replacement"),
               "Distance:" = c("Mahalanobis", "glm"))


## PLOT S FOREST ##########################################
S_forest  <- read.csv("cov_parameter_loop_forest_S.csv")
S_forest <- S_forest[,-1]

#convert to the number matched
S_forest$matched <- nrow(STRICT[STRICT$Type == 1, ])-S_forest$unmatched


# Define the thresholds
threshold_replacement_T <- 0.75 * nrow(STRICT[STRICT$Type == 1, ])

# Create a new column more_75pc_matched in S_forest based on the appropriate threshold
S_forest$more_75pc_matched <-
  ifelse(S_forest$matched > threshold_replacement_T,T,F)

#clean dataset
S_forest$bal__mean_achieved  <- ifelse(S_forest$mean_std_diff<0.25,T,F)
S_forest$bal_achieved_all_cov  <- ifelse(S_forest$max_std_diff<0.25,T,F)
S_forest$significance <- ifelse(S_forest$sig<0.05,T,F)
S_forest  <- S_forest[,-c(1:2, 5:7, 20)]
S_forest <- S_forest[,c(1,2,15:18,3:14)]
S_forest <- S_forest[!duplicated(S_forest), ]

#get invalid models
rows_invalid_models = which( 
  S_forest$more_75pc_matched==F | 
    S_forest$bal__mean_achieved==F |
    S_forest$bal_achieved_all_cov==F)  

#remove invalid models
valid_S_forest <- if (length(rows_invalid_models) == 0) {
  S_forest
} else {
  S_forest[-rows_invalid_models, ]
}


#get significant models
rows_sig_models <- which(valid_S_forest$significance == T)

#get main model
row_main_model = which(valid_S_forest$cal05==T & valid_S_forest$replacement==T & valid_S_forest$maha==F & 
                         valid_S_forest$Access==T & valid_S_forest$Forest==T & valid_S_forest$Population==T 
                       &  valid_S_forest$Agriculture == T & valid_S_forest$Grassland == T)  


png(file="fig_forest_S.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_S_forest, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3))


dev.off()

png(file="MAIN_text_fig_forest_S.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_S_forest, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3),
       cex = c(2,1))


dev.off()


## PLOT S GRASS #######################################
S_grass  <- read.csv("cov_parameter_loop_grass_S.csv")
S_grass <- S_grass[,-1]


#convert to the number matched
S_grass$matched <- nrow(STRICT[STRICT$Type == 1, ])-S_grass$unmatched


# Define the thresholds
threshold_replacement_T <- 0.75 * nrow(STRICT[STRICT$Type == 1, ])

# Create a new column more_75pc_matched in S_grass based on the appropriate threshold
S_grass$more_75pc_matched <- 
  ifelse(S_grass$matched > threshold_replacement_T,T,F)

#clean dataset
S_grass$bal__mean_achieved  <- ifelse(S_grass$mean_std_diff<0.25,T,F)
S_grass$bal_achieved_all_cov  <- ifelse(S_grass$max_std_diff<0.25,T,F)
S_grass$significance <- ifelse(S_grass$sig<0.05,T,F)
S_grass  <- S_grass[,-c(1:2, 5:7, 20)]
S_grass <- S_grass[,c(1,2,15:18,3:14)]
S_grass <- S_grass[!duplicated(S_grass), ]

#get invalid models
rows_invalid_models = which( 
  S_grass$more_75pc_matched==F | 
    S_grass$bal__mean_achieved==F |
    S_grass$bal_achieved_all_cov==F)  

#remove invalid models
valid_S_grass <- if (length(rows_invalid_models) == 0) {
  S_grass
} else {
  S_grass[-rows_invalid_models, ]
}

#get significant models
rows_sig_models <- which(valid_S_grass$significance == T)

#get main model
row_main_model = which(valid_S_grass$cal05==T & valid_S_grass$replacement==T & valid_S_grass$maha==F &  
                         valid_S_grass$Access==T & valid_S_grass$Forest==T & valid_S_grass$Population==T & 
                         valid_S_grass$Agriculture == T & valid_S_grass$Grassland == T)  


png(file="fig_grass_S.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_S_grass, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3))

dev.off()


png(file="MAIN_text_fig_grass_S.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_S_grass, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3),
       cex = c(2,1))

dev.off()

## PLOT S AGRI ################################################
S_agri  <- read.csv("cov_parameter_loop_agri_S.csv")
S_agri <- S_agri[,-1]

#convert to the number matched
S_agri$matched <- nrow(STRICT[STRICT$Type == 1, ])-S_agri$unmatched

# Define the thresholds
threshold_replacement_T <- 0.75 * nrow(STRICT[STRICT$Type == 1, ])

# Create a new column more_75pc_matched in S_agri based on the appropriate threshold
S_agri$more_75pc_matched <- 
  ifelse(S_agri$matched > threshold_replacement_T,T,F)

#clean dataset
S_agri$bal__mean_achieved  <- ifelse(S_agri$mean_std_diff<0.25,T,F)
S_agri$bal_achieved_all_cov  <- ifelse(S_agri$max_std_diff<0.25,T,F)
S_agri$significance <- ifelse(S_agri$sig<0.05,T,F)
S_agri  <- S_agri[,-c(1:2, 5:7, 20)]
S_agri <- S_agri[,c(1,2,15:18,3:14)]
S_agri <- S_agri[!duplicated(S_agri), ]

#get invalid models
rows_invalid_models = which( 
  S_agri$more_75pc_matched==F | 
    S_agri$bal__mean_achieved==F |
    S_agri$bal_achieved_all_cov==F)  

#remove invalid models
valid_S_agri <- if (length(rows_invalid_models) == 0) {
  S_agri
} else {
  S_agri[-rows_invalid_models, ]
}

#get significant models
rows_sig_models <- which(valid_S_agri$significance == T)

#get main model
row_main_model = which(valid_S_agri$cal05==T & valid_S_agri$replacement==T & valid_S_agri$maha==F 
                       &  valid_S_agri$Access==T & valid_S_agri$Forest==T & valid_S_agri$Population==T &
                       valid_S_agri$Agriculture == T & valid_S_agri$Grassland == T)  


png(file="fig_agri_S.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_S_agri, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3))

dev.off()

png(file="MAIN_text_fig_agri_S.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_S_agri, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3),
       cex=c(2,1))

dev.off()

## PLOT LS FOREST ###################################################
LS_forest  <- read.csv("cov_parameter_loop_forest_LS.csv")
LS_forest <- LS_forest[,-1]


#convert to the number matched
LS_forest$matched <- nrow(LESS_STRICT[LESS_STRICT$Type == 1,])-LS_forest$unmatched

# Define the thresholds
threshold_replacement_T <- .75 * nrow(LESS_STRICT[LESS_STRICT$Type == 1,])

# Create a new column more_75pc_matched in LS_forest based on the appropriate threshold
LS_forest$more_75pc_matched <- ifelse(LS_forest$matched > threshold_replacement_T,T,F)

#clean dataset
LS_forest$bal__mean_achieved  <- ifelse(LS_forest$mean_std_diff<0.25,T,F)
LS_forest$bal_achieved_all_cov  <- ifelse(LS_forest$max_std_diff<0.25,T,F)
LS_forest$significance <- ifelse(LS_forest$sig<0.05,T,F)
LS_forest  <- LS_forest[,-c(1:2, 5:7, 20)]
LS_forest <- LS_forest[,c(1,2,15:18,3:14)]
LS_forest <- LS_forest[!duplicated(LS_forest), ]

#get invalid models
rows_invalid_models = which( 
  LS_forest$more_75pc_matched==F | 
    LS_forest$bal__mean_achieved==F |
    LS_forest$bal_achieved_all_cov==F)  


#remove invalid models
valid_LS_forest <- if (length(rows_invalid_models) == 0) {
  LS_forest
} else {
  LS_forest[-rows_invalid_models, ]
}

#get significant models
rows_sig_models <- which(valid_LS_forest$significance == T)

#get main model
row_main_model = which(valid_LS_forest$cal0==T & valid_LS_forest$replacement==F & valid_LS_forest$maha==T &  
                         valid_LS_forest$Access==T & valid_LS_forest$Forest==T & valid_LS_forest$Population==T & 
                         valid_LS_forest$Agriculture == T & valid_LS_forest$Grassland == T)  


png(file="fig_forest_LS.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_LS_forest, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3))

dev.off()

png(file="MAIN_text_fig_forest_LS.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_LS_forest, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3),
       cex = c(2,1))

dev.off()

## PLOT LS GRASS ###################################################
LS_grass  <- read.csv("cov_parameter_loop_grass_LS.csv")
LS_grass <- LS_grass[,-1]


#convert to the number matched
LS_grass$matched <- nrow(LESS_STRICT[LESS_STRICT$Type == 1, ])-LS_grass$unmatched


# Define the thresholds
threshold_replacement_T <- .75 * nrow(LESS_STRICT[LESS_STRICT$Type == 1, ])

# Create a new column more_75pc_matched in LS_grass based on the appropriate threshold
LS_grass$more_75pc_matched <- 
  ifelse(LS_grass$matched > threshold_replacement_T,T,F)

#clean dataset
LS_grass$bal__mean_achieved  <- ifelse(LS_grass$mean_std_diff<0.25,T,F)
LS_grass$bal_achieved_all_cov  <- ifelse(LS_grass$max_std_diff<0.25,T,F)
LS_grass$significance <- ifelse(LS_grass$sig<0.05,T,F)
LS_grass  <- LS_grass[,-c(1:2, 5:7, 20)]
LS_grass <- LS_grass[,c(1,2,15:18,3:14)]
LS_grass <- LS_grass[!duplicated(LS_grass), ]

#get invalid models
rows_invalid_models = which( 
  LS_grass$more_75pc_matched==F | 
    LS_grass$bal__mean_achieved==F |
    LS_grass$bal_achieved_all_cov==F)  

#remove invalid models
valid_LS_grass <- if (length(rows_invalid_models) == 0) {
  LS_grass
} else {
  LS_grass[-rows_invalid_models, ]
}

#get significant models
rows_sig_models <- which(valid_LS_grass$significance == T)

#get main model
row_main_model = which(valid_LS_grass$cal0==T & valid_LS_grass$replacement==F & valid_LS_grass$maha==T &  
                         valid_LS_grass$Access==T & valid_LS_grass$Forest==T & valid_LS_grass$Population==T & 
                         valid_LS_grass$Agriculture == T & valid_LS_grass$Grassland == T)  


png(file="fig_grass_LS.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_LS_grass, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3))

dev.off()

png(file="MAIN_text_fig_grass_LS.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_LS_grass, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3),
       cex=c(2,1))

dev.off()


## PLOT LS AGRO ###################################################
LS_agri  <- read.csv("cov_parameter_loop_agri_LS.csv")
LS_agri <- LS_agri[,-1]


#convert to the number matched
LS_agri$matched <- nrow(LESS_STRICT[LESS_STRICT$Type == 1, ])-LS_agri$unmatched


# Define the thresholds
threshold_replacement_T <- .75 * nrow(LESS_STRICT[LESS_STRICT$Type == 1, ])

# Create a new column more_75pc_matched in LS_agri based on the appropriate threshold
LS_agri$more_75pc_matched <-
  ifelse(LS_agri$matched > threshold_replacement_T,T,F)

#clean dataset
LS_agri$bal__mean_achieved  <- ifelse(LS_agri$mean_std_diff<0.25,T,F)
LS_agri$bal_achieved_all_cov  <- ifelse(LS_agri$max_std_diff<0.25,T,F)
LS_agri$significance <- ifelse(LS_agri$sig<0.05,T,F)
LS_agri  <- LS_agri[,-c(1:2, 5:7, 20)]
LS_agri <- LS_agri[,c(1,2,15:18,3:14)]
LS_agri <- LS_agri[!duplicated(LS_agri), ]

#get invalid models
rows_invalid_models = which( 
  LS_agri$more_75pc_matched==F | 
    LS_agri$bal__mean_achieved==F |
    LS_agri$bal_achieved_all_cov==F)  

#remove invalid models
valid_LS_agri <- if (length(rows_invalid_models) == 0) {
  LS_agri
} else {
  LS_agri[-rows_invalid_models, ]
}

#get significant models
rows_sig_models <- which(valid_LS_agri$significance == T)

#get main model
row_main_model = which(valid_LS_agri$cal0==T & valid_LS_agri$replacement==F & valid_LS_agri$maha==T &  
                         valid_LS_agri$Access==T & valid_LS_agri$Forest==T & valid_LS_agri$Population==T & 
                         valid_LS_agri$Agriculture == T & valid_LS_agri$Grassland == T)  


png(file="fig_agri_LS.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_LS_agri, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3))

dev.off()

png(file="MAIN_text_fig_agri_LS.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_LS_agri, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3),
       cex=c(2,1))

dev.off()


##hh labels
labels <- list("Balance:" = c(">75% obs. unmatched", "Mean bal. achieved", "Balanced for all cov."),
               "Significance" = c("p value < 0.05"),
               "Additional covariates:" = c("Population","Access", "Agriculture"),
               "Model parameters:" = c("0.25SD caliper","0.5SD caliper", "1SD caliper", "No caliper", "With replacement"),
               "Distance:" = c("Mahalanobis", "glm"))


## PLOT MAHFP ###################################################
HH_MAHFP  <- read.csv("cov_parameter_loop_MAHFP.csv")
HH_MAHFP <- HH_MAHFP[,-1]

#convert to the number matched
HH_MAHFP$matched <- nrow(HOUSEHOLD[HOUSEHOLD$Type == 1,])-HH_MAHFP$unmatched

# Define the thresholds
threshold_replacement_T <- 0.75 * nrow(HOUSEHOLD[HOUSEHOLD$Type == 1,])

# Create a new column more_75pc_matched in HH_MAHFP based on the appropriate threshold
HH_MAHFP$more_75pc_matched <- 
  ifelse(HH_MAHFP$matched > threshold_replacement_T,T,F)

#clean dataset
HH_MAHFP$bal__mean_achieved  <- ifelse(HH_MAHFP$mean_std_diff<0.25,T,F)
HH_MAHFP$bal_achieved_all_cov  <- ifelse(HH_MAHFP$max_std_diff<0.25,T,F)
HH_MAHFP$significance <- ifelse(HH_MAHFP$sig<0.05,T,F)
HH_MAHFP  <- HH_MAHFP[,-c(1:2, 5:7, 18)]
HH_MAHFP <- HH_MAHFP[,c(1,2,13:16,3:12)]
HH_MAHFP <- HH_MAHFP[!duplicated(HH_MAHFP), ]
HH_MAHFP_unique <- unique(HH_MAHFP)

# Print all duplicated rows (including the first instance)
HH_MAHFP[duplicated(HH_MAHFP) | duplicated(HH_MAHFP, fromLast = TRUE), ]


#get invalid models
rows_invalid_models = which( 
  HH_MAHFP$more_75pc_matched==F | 
    HH_MAHFP$bal__mean_achieved==F |
    HH_MAHFP$bal_achieved_all_cov==F)  

#remove invalid models
valid_HH_MAHFP <- if (length(rows_invalid_models) == 0) {
  HH_MAHFP
} else {
  HH_MAHFP[-rows_invalid_models, ]
}

#get significant models
rows_sig_models <- which(valid_HH_MAHFP$significance == T)

#get main model
row_main_model = which(valid_HH_MAHFP$cal1==T & valid_HH_MAHFP$replacement==F & valid_HH_MAHFP$maha==F & 
                         valid_HH_MAHFP$Access==T & valid_HH_MAHFP$Population==T & valid_HH_MAHFP$Agriculture==T)  


png(file="fig_MAHFP.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_HH_MAHFP, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3))

dev.off()


png(file="MAIN_text_fig_MAHFP.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_HH_MAHFP, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3),
       cex=c(2,1))

dev.off()

## PLOT HDDS ###################################################
HH_HDDS  <- read.csv("cov_parameter_loop_HDDS.csv")
HH_HDDS <- HH_HDDS[,-1]


#convert to the number matched
HH_HDDS$matched <- nrow(HOUSEHOLD[HOUSEHOLD$Type == 1, ])-HH_HDDS$unmatched


# Define the thresholds
threshold_replacement_T <- .75* nrow(HOUSEHOLD[HOUSEHOLD$Type == 1, ])


# Create a new column more_75pc_matched in HH_HDDS based on the appropriate threshold
HH_HDDS$more_75pc_matched <-  ifelse(HH_HDDS$matched > threshold_replacement_T,T,F)


#clean dataset
HH_HDDS$bal__mean_achieved  <- ifelse(HH_HDDS$mean_std_diff<0.25,T,F)
HH_HDDS$bal_achieved_all_cov  <- ifelse(HH_HDDS$max_std_diff<0.25,T,F)
HH_HDDS$significance <- ifelse(HH_HDDS$sig<0.05,T,F)
HH_HDDS  <- HH_HDDS[,-c(1:2, 5:7, 18)]
HH_HDDS <- HH_HDDS[,c(1,2,13:16,3:12)]
HH_HDDS <- HH_HDDS[!duplicated(HH_HDDS), ]

#get invalid models
rows_invalid_models = which( 
  HH_HDDS$more_75pc_matched==F | 
    HH_HDDS$bal__mean_achieved==F |
    HH_HDDS$bal_achieved_all_cov==F)  

#remove invalid models
valid_HH_HDDS <- if (length(rows_invalid_models) == 0) {
  HH_HDDS
} else {
  HH_HDDS[-rows_invalid_models, ]
}

#get significant models
rows_sig_models <- which(valid_HH_HDDS$significance == T)

#get main model
row_main_model = which(valid_HH_HDDS$cal1==T & valid_HH_HDDS$replacement==F & valid_HH_HDDS$maha==F & 
                         valid_HH_HDDS$Access==T & valid_HH_HDDS$Population==T & valid_HH_HDDS$Agriculture==T)  


png(file="fig_HDDS.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_HH_HDDS, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3))

dev.off()

png(file="MAIN_text_fig_HDDS.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_HH_HDDS, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3),
       cex=2,1)

dev.off()
## PLOT asset ###################################################
HH_asset  <- read.csv("cov_parameter_loop_asset.csv")
HH_asset <- HH_asset[,-1]


#convert to the number matched
HH_asset$matched <- nrow(HOUSEHOLD[HOUSEHOLD$Type == 1, ])-HH_asset$unmatched


# Define the thresholds
threshold_replacement_T <- .75 * nrow(HOUSEHOLD[HOUSEHOLD$Type == 1, ])


# Create a new column more_75pc_matched in HH_asset based on the appropriate threshold
HH_asset$more_75pc_matched <- 
  ifelse(HH_asset$matched > threshold_replacement_T,T,F)

#clean dataset
HH_asset$bal__mean_achieved  <- ifelse(HH_asset$mean_std_diff<0.25,T,F)
HH_asset$bal_achieved_all_cov  <- ifelse(HH_asset$max_std_diff<0.25,T,F)
HH_asset$significance <- ifelse(HH_asset$sig<0.05,T,F)
HH_asset  <- HH_asset[,-c(1:2, 5:7, 18)]
HH_asset <- HH_asset[,c(1,2,13:16,3:12)]
HH_asset <- HH_asset[!duplicated(HH_asset), ]

#get invalid models
rows_invalid_models = which( 
  HH_asset$more_75pc_matched==F | 
    HH_asset$bal__mean_achieved==F |
    HH_asset$bal_achieved_all_cov==F)  

#remove invalid models
valid_HH_asset <- if (length(rows_invalid_models) == 0) {
  HH_asset
} else {
  HH_asset[-rows_invalid_models, ]
}

#get significant models
rows_sig_models <- which(valid_HH_asset$significance == T)

#get main model
row_main_model = which(valid_HH_asset$cal1==T & valid_HH_asset$replacement==F & valid_HH_asset$maha==F & 
                         valid_HH_asset$Access==T & valid_HH_asset$Population==T & valid_HH_asset$Agriculture==T)  

png(file="fig_asset.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_HH_asset, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3))

dev.off()

png(file="MAIN_text_fig_asset.png",
    height = 220, width = 380, res = 300, units = "mm")

schart(valid_HH_asset, labels, 
       order = "increasing",
       highlight = rows_sig_models,
       highlight2 = row_main_model,
       col.est = c("grey90","grey60","black"),
       bg.dot=c("grey60","grey95","grey95","grey60"),
       fonts=c(2,3),
       cex=c(2,1))

dev.off()

#summary info
min(valid_S_forest$ATE) #0.017
mean(valid_S_forest$ATE) # 0.11
max(valid_S_forest$ATE) # 0.26
min(valid_S_grass$ATE) #-1.42
mean(valid_S_grass$ATE) # 3.56
max(valid_S_grass$ATE) # 6.54
min(valid_S_agri$ATE) #-2.03
mean(valid_S_agri$ATE) # -0.97
max(valid_S_agri$ATE) # -0.47

min(valid_LS_forest$ATE) #-0.27
mean(valid_LS_forest$ATE) # 0.03
max(valid_LS_forest$ATE) # 0.37
min(valid_LS_grass$ATE) #-1.02
mean(valid_LS_grass$ATE) # 1.004
max(valid_LS_grass$ATE) # 2.92
min(valid_LS_agri$ATE) #-1.78
mean(valid_LS_agri$ATE) # -1.16
max(valid_LS_agri$ATE) # -0.70
count(valid_LS_forest$significance)
count(valid_LS_grass$significance)

min(valid_HH_MAHFP$ATE) #-1.73
mean(valid_HH_MAHFP$ATE) # -1.21
max(valid_HH_MAHFP$ATE) # -0.98
min(valid_HH_HDDS$ATE) #-0.60
mean(valid_HH_HDDS$ATE) # - 0.12
max(valid_HH_HDDS$ATE) # 0.23
min(valid_HH_asset$ATE) #-1.36
mean(valid_HH_asset$ATE) # -0.34
max(valid_HH_asset$ATE) # 1.89
