# FunctionalMixedEffects.Rpkg

`FunctionalMixedEffects` is an R package designed to perform linear mixed effect regression on functional response variables. It 
functions mainly as a wrapper for the [FunctionalMixedEffects.jl](github.com/wzhorton/FunctionalMixedEffects.jl) julia package. The README 
file there contains more in-depth information about the model construction and sampling scheme. Pre-configured scripts can be found under the `templates` folder.

---
## Install Information

Package installation can be done using the `devtools` package:`devtools::install_github("wzhorton/FunctionalMixedEffects.Rpkg")`. On Windows you may be required to install git.

Note that although the GitHub repository name contains `.Rpkg`, the package itself is named `FunctionalMixedEffects` (i.e. no `.Rpkg` and 
is properly loaded by `library(FunctionalMixedEffects)`).

---
## Update Information

For the sake of simplicity (and not requiring users to know anything about julia), `FunctionalMixedEffect.Rpkg` comes bundled with 
julia code, which is done by linking to the main julia repo as a submodule. Proper updating is handled most easily through
the above `devtools` command. However, it is possible to maintain a local copy and install from source.

- To properly clone the repo, use the `--recurse-submodules` flag
- To update, first `pull` to grab the new submodule commit pointer, then run `git submodules update --remote` to pull submodule content.
