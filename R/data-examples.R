##' Breast cancer survival data
##' 
##' Survival times of 686 patients with primary node positive breast cancer.
##' 
##' 
##' @format A data frame with 686 rows.  \tabular{rll}{ \code{censrec} \tab
##' (numeric) \tab 1=dead, 0=censored \cr \code{rectime} \tab (numeric) \tab
##' Time of death or censoring in days\cr \code{group} \tab (numeric) \tab
##' Prognostic group: \code{"Good"},\code{"Medium"} or \code{"Poor"}, \cr \tab
##' \tab from a regression model developed by Sauerbrei and Royston (1999).\cr
##' }
##' @seealso \code{\link{flexsurvspline}}
##' @references Royston, P. and Parmar, M. (2002).  Flexible parametric
##' proportional-hazards and proportional-odds models for censored survival
##' data, with application to prognostic modelling and estimation of treatment
##' effects. Statistics in Medicine 21(1):2175-2197.
##' 
##' Sauerbrei, W. and Royston, P. (1999). Building multivariable prognostic and
##' diagnostic models: transformation of the predictors using fractional
##' polynomials.  Journal of the Royal Statistical Society, Series A 162:71-94.
##' @source German Breast Cancer Study Group, 1984-1989.  Used as a reference
##' dataset for the spline-based survival model of Royston and Parmar (2002),
##' implemented here in \code{\link{flexsurvspline}}.  Originally provided with
##' the \code{stpm} (Royston 2001, 2004) and \code{stpm2} (Lambert 2009, 2010)
##' Stata modules.
##' @keywords datasets
"bc"


##' Bronchiolitis obliterans syndrome after lung transplants
##' 
##' A dataset containing histories of bronchiolitis obliterans syndrome (BOS)
##' from lung transplant recipients. BOS is a chronic decline in lung function,
##' often observed after lung transplantation.
##' 
##' The entry time of each patient into each stage of BOS was estimated by
##' clinicians, based on their history of lung function measurements and acute
##' rejection and infection episodes.  BOS is only assumed to occur beyond six
##' months after transplant.  In the first six months the function of each
##' patient's new lung stabilises.  Subsequently BOS is diagnosed by comparing
##' the lung function against the "baseline" value.
##' 
##' The same data are provided in the \pkg{msm} package, but in the
##' native format of \pkg{msm} to allow Markov models to be fitted.
##' In \pkg{flexsurv}, much more flexible models can be fitted.
##' @name bos
##' @aliases bosms3 bosms4
##' @docType data
##' @format A data frame containing a sequence of observed or censored
##' transitions to the next stage of severity or death.  It is grouped
##' by patient and includes histories of 204 patients.  All patients
##' start in state 1 (no BOS) at six months after transplant, and may
##' subsequently develop BOS or die.
##' 
##' \code{bosms3} contains the data for a three-state model: no BOS, BOS or
##' death. \code{bosms4} uses a four-state representation: no BOS, mild BOS,
##' moderate/severe BOS or death.  \tabular{rll}{ \code{id} \tab (numeric) \tab
##' Patient identification number \cr \code{from} \tab (numeric) \tab Observed
##' starting state of the transition \cr \code{to} \tab (numeric) \tab Observed
##' or potential ending state of the transition \cr \code{Tstart} \tab
##' (numeric) \tab Time at the start of the interval \cr \code{Tstop} \tab
##' (numeric) \tab Time at the end of the interval \cr \code{time} \tab
##' (numeric) \tab Time difference \code{Tstart}-\code{Tstop} \cr \code{status}
##' \tab (numeric) \tab 1 if the transition to state \code{to} was observed, or
##' 0 if the transition to state \code{to} was censored (for example, if the
##' patient was observed to move to a competing state) \cr \code{trans} \tab
##' (factor) \tab Number of the transition \code{from}-\code{to} in the set of
##' all \code{ntrans} allowed transitions, numbered from 1 to \code{ntrans}.  }
##' @references Heng. D. et al. (1998).  Bronchiolitis Obliterans Syndrome:
##' Incidence, Natural History, Prognosis, and Risk Factors.  Journal of Heart
##' and Lung Transplantation 17(12)1255--1263.
##' @source Papworth Hospital, U.K.
##' @keywords datasets
NULL

#' Object with Compiled Stan Code
#'
#' Pre-compiled stan models that can be made accessed expertsurv::compiled_models_saved
#'
#' @source Generated from \code{expertsurv::compile_stan()}
"compiled_stan"