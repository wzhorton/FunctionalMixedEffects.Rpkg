#### Functions interfacing to the FunctionalMixedEffect.jl moduls ####

#' Set the \code{FME} config variable in the Julia environment.
#'
#' The argument names and values are saved to the \code{julia_obj}
#' environment, then the config object is built.
#'
#' @param julia_obj Julia environment object from \code{julia_setup}.
#' #' @param p number of spline coefficients. Must be strict integer (5L instead of 5).
#' @param save_random_effects Logical indicating whether the random effects
#' chain should be saved. Strict logical is required. Note the centered effects
#' will be saved even if this is false
#' @param save_theta Logical indicating whether the individual coefficients
#' chain should be saved. Strict logical is required.
#' @param niter,nburn,nthin MCMC iteration parameters. Strict integers are required.
#'
#' @return Invisibly returns list of config variables. The config object \code{cfg}
#' is saved in the \code{julia_obj} environment.
#' @export

set_config_fme <- function(julia_obj, p=20L, niter = 5000L, nburn = 1000L, nthin = 1L,
                           save_random_effects = FALSE, save_theta = FALSE){
  if(!all(is.integer(c(p,niter, nburn, nthin)))){
    stop("p, niter, nburn, nthin must all be strict integers
         (i.e. niter=5000L instead of =5000")
  }
  if(!is.logical(save_random_effects)){
    stop("save_random_effects must be strict logical (i.e. TRUE/FALSE)")
  }
  arg_list <- as.list(environment())[-1]
  for(i in seq_along(arg_list)){
    julia_obj$assign(names(arg_list)[i], arg_list[[i]])
  }
  julia_obj$command("cfg = OutputConfigFME(p, niter, nburn, nthin,
                    save_random_effects, save_theta)")
  invisible(arg_list)
}


#' Set the \code{FME} hyperparameters in the Julia environment.
#'
#' The argument names and values are save to the \code{julia_obj}
#' environment, then the hyperparameter object is built. The values must
#' all return \code{TRUE} to \code{is.double}.
#'
#' @param julia_obj Julia environment object from \code{julia_setup}.
#' @param a_sig,b_sig sigma priors. Must have \code{a > 2}
#' @param a_tau,b_tau tau priors. Must have \code{a > 2}
#' @param a_lam,b_lam lambda priors. Must have \code{a > 2}
#' @param v_fix,v_cent fixed and centered effects prior variance.
#'
#' @return Invisibly returns list of hyperparameter variables. The hyperparameter
#' object \code{hyps} is saved in the \code{julia_obj} environment.
#' @export

set_hyperparm_fme <- function(julia_obj, a_sig=3, b_sig=1, a_tau=3, b_tau=1,
                              a_lam=3, b_lam=1, v_fix=1000, v_cent=1000){
  if(!all(sapply(c(a_sig, b_sig, a_tau, b_tau,a_lam, b_lam, v_fix, v_cent), is.double))){
    stop("hyperparameters must all be strict doubles")
  }
  arg_list <- as.list(environment())[-1]
  for(i in seq_along(arg_list)){
    julia_obj$assign(names(arg_list)[i], arg_list[[i]])
  }
  julia_obj$command("hyps = HyperParametersFME(a_sig, b_sig, a_tau, b_tau,
                    a_lam, b_lam, v_fix, v_cent)")
  invisible(arg_list)
}


#' Set the \code{FME} data object in the Julia environment
#'
#' The argument names and values are save to the \code{julia_obj}
#' environment, then the data object is built. The values must
#' all return \code{TRUE} to \code{is.double}.
#'
#' @param julia_obj Julia environment from \code{julia_setup}.
#' @param Yobs matrix containing curves/functions as columns.
#' @param Xfix matrix containing fixed effect covariates, with individuals contained
#' as columns and variables as rows (reverse of usual design matrices). Can be set
#' to \code{NULL} if only random effects are present.
#' @param Xrand matrix containing fixed effect covariates, with individuals contained
#' as columns and variables as rows (reverse of usual design matrices). Can be
#' set to \code{NULL} to perform functional regression without random effects.
#' @param Xcent matrix containing the mapping between random effects and their
#' hierarchical centered effect. Rows contain indicators for which random effects
#' are centered together. Can be set to \code{NULL} when random effects are also
#' missing.
#'
#' @return Invisibly returns list of data variables. The data
#' object \code{data} is saved in the \code{julia_obj} environment.
#' @export

set_data_fme <- function(julia_obj, Yobs, Xfix = NULL, Xrand = NULL, Xcent = NULL){
  if(!is.double(Yobs)){
    stop("Y must be strict double")
  }
  if(is.null(Xfix) && is.null(Xrand)){
    stop("Xfix and Xrand cannot both be NULL")
  }
  if(!is.null(Xfix) && !is.double(Xfix)){
    stop("Xfix must be strict double")
  }
  if(is.null(Xrand) != is.null(Xcent)){
    stop("Xrand and Xcent must either both be provided or both not.")
  }
  if(!is.null(Xrand) && !is.double(Xrand)){
    stop("Xrand must be strict double")
  }
  if(!is.null(Xcent) && !is.double(Xcent)){
    stop("Xcent must be strict double")
  }
  arg_list <- as.list(environment())[-1]
  for(i in seq_along(arg_list)){
    julia_obj$assign(names(arg_list)[i], arg_list[[i]])
  }

  julia_obj$command("data = DataFME(Yobs, Xfix, Xrand, Xcent)")
  invisible(arg_list)
}



#' Functional mixed effects MCMC
#'
#' Fit the functional mixed model to data. \code{HyperParameter}, \code{Data},
#' and \code{OutputConfig} objects must be defined in the Julia environment..
#'
#' @param julia_obj Julia environment from \code{julia_setup}.

#' @return Invisibly returns NULL. The resulting object \code{chains}
#' is saved in the \code{julia_obj} environment.
#' @export

mcmc_fme <- function(julia_obj){
  if(!julia_obj$exists("cfg")){
    stop("Config object not found")
  }
  if(!julia_obj$exists("hyps")){
    stop("Hyperparameter object not found")
  }
  if(!julia_obj$exists("data")){
    stop("Data object not found")
  }

  julia_obj$command("chains = mcmc_fme(data, hyps, cfg);")
  invisible(NULL)
}


#' Extract \code{FME} model chain objects into an R list
#'
#' Custom structs like the returned chains object don't return properly in R.
#' This function extracts the elements into a list.
#'
#' @param julia_obj Julia environment object from \code{julia_setup} which has
#' also run the \code{FME} MCMC code.
#'
#' @return List of chain objects. Note the special unicode characters in the
#' Julia object are renamed using standard characters. Also note that parameters
#' associated with NULL data objects will be returned as empty arrays.
#' @export

extract_chains_fme <- function(julia_obj){
  if(!julia$exists("chains")){
    stop("No chains object found")
  }
  chains <- list()
  chains$sigma <- julia_obj$eval("chains.σ")
  chains$tau <- julia_obj$eval("chains.τ")
  chains$lambda <- julia_obj$eval("chains.λ")
  chains$Theta <- julia_obj$eval("chains.θ")
  chains$Bfix <- julia_obj$eval("chains.Bfix")
  chains$Brand <- julia_obj$eval("chains.Brand")
  chains$Bcent <- julia_obj$eval("chains.Bcent")
  return(chains)
}


#' Functional mixed effects model wrapper function
#'
#' Fits the functional mixed effects model. This allows the user to provide all
#' parameters at once rather than sequentially calling the interface functions.
#'
#' @param julia_obj A Julia environment from \code{julia_setup}.
#' @param Yobs,Xfix,Xrand,Xcent see documentation for \code{set_data_fme}
#' @param p,niter,nburn,nthin,save_random_effects see documentation for \code{set_config_fme}
#' @param a_sig,b_sig,a_tau,b_tau,a_lam,b_lam,v_fix,v_cent see documentation for \code{set_hyperparm_fme}
#'
#' @return MCMC chains object, see documentation for \code{mcmc_fme}.
#' @export

fit_model_fme <- function(julia_obj, Yobs, Xfix = NULL, Xrand = NULL, Xcent = NULL,
                          p = 20L, niter = 15000L, nburn = 5000L, nthin = 1L,
                          save_random_effects = FALSE, save_theta = FALSE,
                          a_sig = 3, b_sig = 1, a_tau = 3, b_tau = 1,
                          a_lam = 3, b_lam = 1, v_fix = 1000, v_cent=1000){
  set_config_fme(julia_obj, p, niter, nburn, nthin, save_random_effects, save_theta)
  set_hyperparm_fme(julia_obj, a_sig, b_sig, a_tau, b_tau, a_lam, b_lam, v_fix, v_cent)
  set_data_fme(julia_obj, Yobs, Xfix, Xrand, Xcent)
  mcmc_fme(julia_obj)
  extract_chains_fme(julia_obj)
}

