#### Julia Setup Functions ####

#' Build Julia environment and load \code{FunctionalMixedEffects.jl} module.
#'
#' @param install_julia Logical indicating if Julia should be installed in the
#' event that an install either is not installed or cannot be found on the path.
#'
#' @return A \code{JuliaCall} object which interfaces to a Julia session.
#' @export

setup_julia <- function(install_julia = FALSE){
  julia <- JuliaCall::julia_setup(installJulia = install_julia, rebuild = TRUE)
  julia$command("using Pkg")
  julia$call("Pkg.activate", system.file("./julia/FunctionalMixedEffects.jl/",
                                         package = "FunctionalMixedEffects"))
  julia$command("Pkg.instantiate()")
  julia$source(system.file("./julia/FunctionalMixedEffects.jl/src/FunctionalMixedEffects.jl",
                           package = "FunctionalMixedEffects"))
  julia$command("using .FunctionalMixedEffects")
  return(julia)
}
