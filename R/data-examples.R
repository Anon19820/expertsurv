##' Breast cancer survival data
##' 
##' Survival times of 686 patients with primary node positive breast cancer. See \pkg{flexsurv} for details.
##' @keywords datasets
"bc"


##' Bronchiolitis obliterans syndrome after lung transplants
##' 
##' See \pkg{flexsurv} for details.
##' @name bos
##' @aliases bosms3 bosms4
##' @docType data
##' @format A data frame containing a sequence of observed or censored
##' transitions to the next stage of severity or death.  It is grouped
##' by patient and includes histories of 204 patients.  All patients
##' start in state 1 (no BOS) at six months after transplant, and may
##' subsequently develop BOS or die.
##' @references Heng. D. et al. (1998).  Bronchiolitis Obliterans Syndrome:
##' Incidence, Natural History, Prognosis, and Risk Factors.  Journal of Heart
##' and Lung Transplantation 17(12)1255--1263.
##' @source Papworth Hospital, U.K.
##' @keywords datasets
NULL

#' Object with Compiled Stan Code
#'
#' Pre-compiled stan models that can be accessed by \code{expertsurv::compiled_models_saved}. These compiled stan models can be supplied to the `fit.models.expert` function with the argument `compile_mods`. This greatly facilitates the estimation of the parametric survival models run by stan.
#'
#' @source Generated from \code{expertsurv::compile_stan()}
#' @name compiled_models_saved
#' @docType data
#' @keywords datasets
"compiled_models_saved"