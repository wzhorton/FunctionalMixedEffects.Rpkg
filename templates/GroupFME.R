#### GroupFME.R ####

# Version 1.4.3

# This script is a template script found within the FunctionalMixedEffects
# package. See github.com/wzhorton/FunctionalMixedEffects.Rpkg for install
# instructions and up to date script versions.


#-------------------------------------------------------------------------------
# Model Overview #
#-------------------------------------------------------------------------------

# Core model: comparison of waveform shapes between groups
# Optional: controlling for covariates (both subject and trial level)
# Optional: repeated measures (sum-to-zero random intercepts) subject effects


#-------------------------------------------------------------------------------
# Data Overview #
#-------------------------------------------------------------------------------

# Functional response curves:
# - Time normalized
# - Same number of trials is NOT required
# - Subject ID and group label included

# Subject level variables (optional):
# - Subject ID
# - Other clinical info
# - Group label (optional*)

# Trial level variables (optional):
# - Subject ID (optional*)
# - Trial number (optional*)
# - Other trial info

# *will be removed

#-------------------------------------------------------------------------------
# Output Overview #
#-------------------------------------------------------------------------------

# Group effects / covariate effects / comparisons
# - CSV file with 3 columns for each group
# - 4 columns: mean estimate, sd estimate, and 95% upper and lower curves
# - Extra columns for Cohen d effect sizes on comparisons

# Significant regions
# - Separate CSV file for each analysis above
# - Contains significant region boundaries and max value


#-------------------------------------------------------------------------------
# Julia Installation #
#-------------------------------------------------------------------------------

# The Julia backend used by this package is fast, easy to maintain, and cross
# platform compatible. This script can automatically install and configure Julia
# if no installation is found.

# If you have a custom Julia setup, set this to FALSE to prevent accidentally
# installing a duplicate version. TRUE is recommended otherwise.
install_julia_if_not_found <- TRUE

#-------------------------------------------------------------------------------
# Data Loading Details #
#-------------------------------------------------------------------------------

# Data should be organized into 3 separate CSV files, one for each of these: trial
# level variables, subject level variables, and response curves. Details for
# each are given below.

#---------------------------#
# Response curves/waveforms #
#---------------------------#
# Structure:
# - Curves correspond with columns
# - First row should contain group labels (no spaces)
# - Second row should contain ID numbers (no spaces)
# - Time normalized; each curve is the same length
# - No missing values anywhere, remove missing columns
# - Aggregated; all groups contained in one file

# Specify name and path of waveform file
waveform_path <- "path/to/waveformfile.csv"

#----------------------------#
# Subject level variables #
#----------------------------#
# Optional; only required for variables beyond group label.

# Structure:
# - Variables correspond with columns
# - 1 required columns: ID values
# - Group label may be supplied, but will be removed
# - ID numbers should uniquely match all values in waveform file
# - First row should contain variable names (no spaces)

# Specify name and path of individual variable file.
# Set to NULL if not including, or comment out.
subj_path <- "path/to/subjvars.csv"

# Provide the file column names for group labels and ID values.
# If path is NULL, these are ignored. Set to NULL if otherwise not present.
group_column_subj <- "Grp"
id_column_subj <- "Sub" #required

#-----------------------#
# Trial level variables #
#-----------------------#
# Structure:
# - Variables correspond with columns
# - 1 required column: (repeating) ID values
# - Order must match the waveform file order
# - First row should contain variable names (no spaces)
# - Trial number, ID, and group may be supplied, but will be removed

# Specify name and path of trial variable file.
# Set to NULL if not including, or comment out.
trial_path <- "path/to/trialvars.csv"

# Provide file column name for ID values and trial numbers.
# If path is NULL, these are ignored. Set to NULL if otherwise not present
id_column_trial <- "Sub"
trial_column_trial <- "Trial"
group_column_trial <- "Grp"


#-------------------------------------------------------------------------------
# Group Comparison Builder #
#-------------------------------------------------------------------------------

# Along with the group analysis, the model will also output inference for
# specific group comparisons, as well as aggregate comparisons. Each comparison
# is composed of two super groups (groups of groups) and ultimately a difference
# curve between groups is computed. Note that the group names must match those
# found in the subject variable file.

# Indicate whether comparisons should be computed and saved. Set TRUE or FALSE.
output_comparisons <- TRUE

# Set the comparisons in the list below, including a name to use for saving
# results. If output_comparisons is FALSE, this will be ignored and does not
# need to be modified.
comparison_list <- list(
  list(c("G1T1"), c("G2T1"), "G1T1-minus-G2T1"),
  list(c("G1T1","G1T2"), c("G2T1","G2T2"), "G1-minus-G2")
)

#-------------------------------------------------------------------------------
# Covariate Baseline Values #
#-------------------------------------------------------------------------------

# When controlling for covariates in a group analysis, results are standardized
# to some baseline covariate value. It is important to select covariate baselines
# with meaningful interpretations, such as known averages, values of interest,
# sample averages, or zero.

# If file paths above are set to NULL, the respective lists below are ignored.

# Provide the subject-level baselines (not including ID or group label) in
# this list with correctly matching names. Note that "avg" can be used to
# automatically compute the sample average.
baselines_subj <- list(
  weight = 80.5,
  height = "avg"
)

# Provide the trial-level baselines (not including ID or trial) in this list
# with correctly matching names.
baselines_trial <- list(
  time = "avg",
  hope = 0
)


#-------------------------------------------------------------------------------
# Include Repeated Measures (Random Intercepts) #
#-------------------------------------------------------------------------------

# The mixed effects structure of the package allows for incorporating
# repeated measures in the form of a random intercepts model. The model output
# will not include random effects estimates, but the fixed effects will be
# adjusted to account for the repeated measures structure.

# Set the repeated measures indicator to either TRUE or FALSE
repeated_measures <- FALSE


#-------------------------------------------------------------------------------
# Output Details #
#-------------------------------------------------------------------------------

# Provide both the folder path and a common name to prefix to each output file.
output_folder_path <- "path/to/folder"
output_prefix <- "VGRF"

# Indicate whether sig-region output should be saved. Note that group averages
# and comparisons are always saved, and covariate analysis is saved if included.
# Set the appropriate flag to TRUE or FALSE.
output_sig_regions <- TRUE


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Main Script (Modify at your own risk) #
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-- Utility functions --#
# Validation
abort <- function(msg){
  stop(msg)
  invokeRestart("abort")
}

validate <- function(expr, msg){
  tryCatch(
    expr = expr,
    error = function(e){
      abort(msg)
    }
  )
}

#Cohen's d calculation
cohend <- function(diff, ns, sds, conf = 0.95){
  poolsd <- sqrt(sum(sds^2*(ns-1))/sum(ns-1))
  effect_size <- diff/poolsd
  effect_se <- sqrt(sum(1/ns) + 0.5*effect_size^2/sum(ns))
  es_vec <- effect_size + c(0,-1,1)*qnorm(conf)*effect_se
  es_vec <- c(es_vec[1], effect_se, es_vec[2:3])
  names(es_vec) <- c("eff_size","std_err","es_lower","es_upper")
  return(es_vec)
}

#Significant area extraction
extract_sig_regions <- function(outmat, lb_col = 3, ub_col = 4,
                                mn_col = 1){
  sig_flags <- as.numeric(outmat[,lb_col]*outmat[,ub_col] > 0)
  change_points <- diff(sig_flags)
  starts <- which(change_points == 1) + 1
  ends <- which(change_points == -1)
  if(head(sig_flags,1) != 0){ starts <- c(1, starts) }
  if(tail(sig_flags,1) != 0){ ends <- c(ends, nrow(outmat)) }

  sig_mat <- cbind(starts, ends)
  if(length(sig_mat) == 0){
    area_mat <- matrix(c(NA,NA,NA,NA,NA,NA), nrow=1)
  } else {
    maxes <- apply(sig_mat, 1, function(r){
      max_loc <- which.max(abs(outmat[r[1]:r[2],mn_col])) + r[1] - 1
      max_val <- outmat[,mn_col][max_loc]
      max_ub <- outmat[,lb_col][max_loc]
      max_lb <- outmat[,ub_col][max_loc]
      return(c(max_loc, max_val, max_lb, max_ub))
    })
    area_mat <- cbind(sig_mat, t(maxes))
  }
  colnames(area_mat) <- c("start","end","max_location","max_value",
                          "max_lower", "max_upper")
  return(area_mat)
}


#-- Load package --#
validate(library(FunctionalMixedEffects),
         "FunctionalMixedEffects not installed.")

#-- Load Curve Data --#
validate(curves <- read.csv(waveform_path, header=FALSE,skip=2, stringsAsFactors=F),
         "Waveform path invalid. File failed to load.")
id_curves <- as.character(read.csv(waveform_path, header=FALSE, nrows=1,
                                   skip=1, stringsAsFactors=F))
grp_curves <- as.character(read.csv(waveform_path, header = FALSE,
                                    nrows = 1, stringsAsFactors=F))
curves <- as.matrix.data.frame(curves)

if(typeof(curves) != "double"){
  abort(paste("Waveform data is not numeric. Type read:", typeof(curves)))
}
if(any(is.na(curves))) {
  abort("Missing values found in waveform file.")
}


#-- Load subject variable data --#
if(exists("subj_path") && !is.null(subj_path)){
  validate(subj_covs <- read.csv(subj_path, stringsAsFactors=F),
           "Individual variable file failed to load.")

  # ID and group label processing
  subj_cov_labels <- colnames(subj_covs)

  if(is.null(id_column_subj) || !(id_column_subj %in% subj_cov_labels)){
    abort("ID column label is not found in individual variables file")
  }
  id_subj <- subj_covs[,id_column_subj]
  if(!setequal(id_subj,id_curves)){
    abort("Subject covariate IDs do not match waveform IDs")
  }
  if(any(duplicated(id_subj)){
    abort("Non-unique subject IDs found in subject covariate file")
  }
  subj_covs[,id_column_subj] <- NULL

  if(!is.null(group_column_subj)){
    if(!(group_column_subj %in% subj_cov_labels)){
      abort("Group column label is not found in subject variables file")
    }
    grp_subj <- subj_covs[, group_column_subj]
    if(!setequal(grp_subj, grp_curves)){
      abort("Subject group labels do not match waveform group labels")
    }
    subj_covs[,group_column_subj] <- NULL
  }

  if(ncol(subj_covs) == 0){
    abort("No other covariates found in covariate file")
  }

  # Baseline processing
  subj_cov_labels <- colnames(subj_covs)
  if(!setequal(names(baselines_subj), subj_cov_labels)){
    abort("Subject covariate names do not match baseline list")
  }
  for(i in seq_along(baselines_subj)){
    cov_base <- baselines_subj[[i]]
    cov_name <- names(baselines_subj)[i]
    if(cov_base == "avg"){
      cov_base <- mean(subj_covs[,cov_name])
    }
    subj_covs[,cov_name] <- subj_covs[,cov_name] - cov_base
  }
}


#-- Load trial variable data --#
if(exists("trial_path") && !is.null(trial_path)){
  validate(trial_covs <- read.csv(trial_path, stringsAsFactors=F),
           "Trial level variable file failed to load.")

  #ID, group, and trial number processing
  trial_cov_labels <- colnames(trial_covs)

  if(!is.null(id_column_trial)){
    if(!(id_column_trial %in% trial_cov_labels)){
      abort("ID column label is not found in trial variables file")
    }
    id_trial <- trial_covs[,id_column_trial]
    if(any(id_trial != id_curves)){
      abort("Trial variable IDs do not match waveform IDs")
    }
    trial_covs[,id_column_trial] <- NULL
  }

  if(!is.null(group_column_trial)){
    if(!(group_column_trial %in% trial_cov_labels)){
      abort("Group column label is not found in trial variables file")
    }
    group_trial <- trial_covs[,group_column_trial]
    if(any(group_trial != grp_curves)){
      abort("Trial variable group labels do not match waveform group labels")
    }
    trial_covs[,group_column_trial] <- NULL
  }

  if(!is.null(trial_column_trial)){
    if(!(trial_column_trial %in% trial_cov_labels)){
      abort("Trial column label is not found in trial variables file")
    }
    trial_trial <- trial_covs[,trial_column_trial]
    trial_covs[,trial_column_trial] <- NULL
  }

  # Baseline processing
  trial_cov_labels <- colnames(trial_covs)
  if(!setequal(names(baselines_trial), trial_cov_labels)){
    abort("Trial covariate names do not match baseline list")
  }
  for(i in seq_along(baselines_trial)){
    cov_base <- baselines_trial[[i]]
    cov_name <- names(baselines_trial)[i]
    if(cov_base == "avg"){
      cov_base <- mean(trial_covs[,cov_name])
    }
    trial_covs[,cov_name] <- trial_covs[,cov_name] - cov_base
  }
}


#-- Extra Input Validation --#
if(output_comparisons){
  for(comp in comparison_list){
    if(length(comp) != 3){
      abort("A comparison list item appears malformed. Expected format: list(c(...),c(...),'...')")
    }
  }
}


#-- Construct Design Matrices --#
grp_design <- setNames(data.frame(matrix(as.numeric(model.matrix(~-1+grp_curves)),
                                  ncol = length(unique(grp_curves)))), levels(as.factor(grp_curves)))
id_design <- setNames(data.frame(matrix(as.numeric(model.matrix(~-1+id_curves)),
                                         ncol = length(unique(id_curves)))), levels(as.factor(id_curves)))

Xfix <- cbind(grp_design)
if(!is.null(trial_path)){
  Xfix <- cbind(grp_design, trial_covs)
}

if(!is.null(subj_path)){
  Xsubj <- NULL
  for(i in seq_along(id_curves)){
    Xsubj <- rbind(Xsubj, subj_covs[which(id_subj == id_curves[i]), , drop = FALSE])
  }
  Xfix <- cbind(Xfix,Xsubj)
}
Xfix <- t(Xfix) #transpose for model compatibility

if(repeated_measures){
  if(all(sapply(split(grp_curves, id_curves), function(v) length(unique(v))) == 1)){ # check for subject-group nesting
    # Implement nested sum constraints
    grouped_ids <- lapply(split(id_curves, grp_curves), function(l) unique(l))
    Xrand <- as.matrix(id_design)
    for(ids in grouped_ids){
      inds <- match(ids, colnames(Xrand))
      Xrand[Xrand[,inds[1]] == 1, inds] <- -1
      Xrand <- Xrand[,-inds[1]]
    }
    Xrand <- t(Xrand) #transpose for model compatibility
    Xcen <- matrix(1, ncol = nrow(Xrand)) #default marginal term
  } else {
    abort("Non-nested subjects not yet implemented. Automatic sum constraints were not determined.")
  }

} else {
  Xrand <- NULL
  Xcen <- NULL
}


#-- MCMC Model Fit --#
setwd(output_folder_path)

num_spline_bases <- 25L
mcmc_iterations <- 25000L
mcmc_burnin <- as.integer(round(mcmc_iterations/10))

if(!exists("julia")){
  julia <- setup_julia(install_julia = install_julia_if_not_found)
}

chains <- fit_model_fme(julia, curves, Xfix, Xrand, Xcen, p = num_spline_bases,
                        niter = mcmc_iterations, nburn = mcmc_burnin)


#-- Compute and Save Group Effects --#
ngrp <- ncol(grp_design)
H <- splines::bs(seq(0,1,len=101),#nrow(curves)),
                 df = num_spline_bases, intercept = TRUE)

for(g in 1:ngrp){ #group effects are first
  grp_bchain <- H%*%chains$Bfix[,g,]
  grp_mean <- rowMeans(grp_bchain)
  grp_sd <- apply(grp_bchain, 1, sd)
  grp_low <- apply(grp_bchain, 1, function(r) quantile(r, 0.025))
  grp_up <- apply(grp_bchain, 1, function(r) quantile(r, 0.975))
  out_grp <- cbind(grp_mean, grp_sd, grp_low, grp_up)
  write.csv(out_grp, paste0(output_prefix,"-", names(grp_design)[g],
                            "-GroupMean.csv"))
  if(output_sig_regions){
    write.csv(extract_sig_regions(out_grp), paste0(output_prefix,"-",
                names(grp_design)[g],"-GroupMean-Areas.csv"),
              row.names = FALSE)
  }
}


#-- Compute and Save Remaining Covariate Effects --#
if((exists("subj_path") && !is.null(subj_path)) || (exists("trial_path") && !is.null(trial_path))){
  for(cc in (ngrp+1):nrow(Xfix)){
    cc_bchain <- H%*%chains$Bfix[,cc,]
    cc_mean <- rowMeans(cc_bchain)
    cc_sd <- apply(cc_bchain, 1, sd)
    cc_low <- apply(cc_bchain, 1, function(r) quantile(r, 0.025))
    cc_up <- apply(cc_bchain, 1, function(r) quantile(r, 0.975))
    out_cc <- cbind(cc_mean, cc_sd, cc_low, cc_up)
    write.csv(out_cc, paste0(output_prefix,"-", rownames(Xfix)[cc],
                              "-Effect.csv"))
    if(output_sig_regions){
      write.csv(extract_sig_regions(out_cc), paste0(output_prefix,"-",
                   rownames(Xfix)[cc],"-Effect-Areas.csv"),
                row.names = FALSE)
    }
  }
}


#-- Compute and Save Difference Effects --#
if(output_comparisons){
  for(i in 1:length(comparison_list)){
    grps1 <- comparison_list[[i]][[1]]
    grp1_inds <- sapply(grps1, function(g) which(rownames(Xfix) == g))
    grp1_bchain <- H%*%apply(chains$Bfix[,grp1_inds,,drop=FALSE], c(1,3), mean)
    es1_subset <- curves[,as.logical(colSums(Xfix[grp1_inds,,drop=FALSE]))]
    es1_mean <- rowMeans(es1_subset)
    es1_sd <- apply(es1_subset, 1, sd)
    n1 <- sum(grp_design[,grps1])

    grps2 <- comparison_list[[i]][[2]]
    grp2_inds <- sapply(grps2, function(g) which(rownames(Xfix) == g))
    grp2_bchain <- H%*%apply(chains$Bfix[,grp2_inds,,drop=FALSE], c(1,3), mean)
    es2_subset <- curves[,as.logical(colSums(Xfix[grp2_inds,,drop=FALSE]))]
    es2_mean <- rowMeans(es2_subset)
    es2_sd <- apply(es2_subset, 1, sd)
    n2 <- sum(grp_design[,grps2])

    ns <- c(n1,n2)

    diff_bchain <- grp1_bchain - grp2_bchain
    diff_mean <- rowMeans(diff_bchain)
    diff_sd <- apply(diff_bchain, 1, sd)
    diff_low <- apply(diff_bchain, 1, function(r) quantile(r, 0.025))
    diff_up <- apply(diff_bchain, 1, function(r) quantile(r, 0.975))
    es <- t(sapply(seq_along(diff_mean), function(j){
        cohend(es1_mean[j] - es2_mean[j], ns, c(es1_sd[j],es2_sd[j]))
    }))

    out_diff <- cbind(diff_mean, diff_sd, diff_low, diff_up, es)
    write.csv(out_diff, paste0(output_prefix,"-", comparison_list[[i]][3],
                             "-DiffMean.csv"))
    if(output_sig_regions){
      sig_areas <- extract_sig_regions(out_diff)
      if(all(is.na(sig_areas))){
        sig_out <- cbind(sig_areas, matrix(NA, nrow = 1, ncol = 4,
               dimnames = list(NULL, c("eff_size","std_err","es_lower","es_upper"))))
      } else {
        sig_out <- cbind(sig_areas,t(sapply(sig_areas[,3], function(max_ind){
          out_diff[max_ind,5:8]
        })))
      }
      write.csv(sig_out, paste0(output_prefix,"-", comparison_list[[i]][3],
                                 "-DiffMean-Areas.csv"), row.names = FALSE)
    }
  }
}


