#' Plot survival curves for the models fitted using \code{fit.models}
#' 
#' @param ... 
#' Must include at least one result object saved as the call to the \code{fit.models} function. 
#' May include other optional parameters. These include:
#' 
#' * \code{add.km}: Whether the KM curve should be added.
#' * \code{newdata}: Specifies a profile of covariates (in the list \code{newdata}). Other possibilities are additional (mainly graphical) options:
#'     * \code{xlab}: A string with the label for the x-axis (default = "time").
#'     * \code{ylab}: A string with the label for the y-axis (default = "Survival").
#'     * \code{lab.profile}: A (vector of) string(s) indicating the labels associated with the strata defining the different survival curves to plot. Defaults to the value used by the Kaplan Meier estimate given in \code{fit.models}.
#'     * \code{xlim}: A vector determining the limits for the x-axis.
#'     * \code{colors}: A vector of characters defining the colours in which to plot the different survival curves.
#'     * \code{lab.profile}: A vector of characters defining the names of the models fitted.
#'     * \code{add.km}: TRUE (whether to also add the Kaplan Meier estimates of the data).
#'     * \code{annotate}: FALSE (whether to also add text to highlight the observed vs extrapolated data).
#'     * \code{legend.position}: A vector of proportions to place the legend. Default to 'c(.75,.9)', which means 75% across the x-axis and 90% across the y-axis.
#'     * \code{legend.title}: Suitable instructions to format the title of the legend; defaults to 'element_text(size=15,face="bold")' but other arguments can be added (using 'ggplot' facilities).
#'     * \code{legend.text}: Suitable instructions to format the text of the legend; defaults to 'element_text(colour="black", size=14, face="plain")' but other arguments can be added (using 'ggplot' facilities).
#' * \code{plot_opinion}: TRUE will provide an illustration of the expert opinion at each time-point.
#' * \code{plot_ci}: Statistical uncertainty can be plotted using the \code{plot_ci = TRUE} argument and by specifying \code{nsim} equal to the number of desired simulations (for Bayesian models, this must be less than the total number of simulations from the posterior). By default, the confidence/credible intervals are plotted as dashed lines. If an area/ribbon plot is preferred, set \code{ci_plot_ribbon = TRUE}.
#' * \code{nsim}: Even if statistical uncertainty is not required in the plots, it is recommended that \code{nsim} is set to a reasonable number. If \code{nsim = 1} by default, the maximum likelihood estimates or the posterior mean of the parameters will be used to plot the results. In most cases, this should suffice (particularly for maximum likelihood). However, the expected survival estimated by the full sampling distribution may be different from the estimate at its expectation/maximum likelihood estimate.
#' @return A ggplot2 object of the survival curves.
#' @author Gianluca Baio
#' @seealso \code{fit.models}, \code{write.surv}
#' @keywords Parametric survival models
#' @examples
#' require("dplyr")
#' param_expert_example1 <- list()
#' param_expert_example1[[1]] <- data.frame(dist = c("norm","t"),
#'                                          wi = c(0.5,0.5), # Ensure Weights sum to 1
#'                                          param1 = c(0.1,0.12),
#'                                          param2 = c(0.15,0.5),
#'                                          param3 = c(NA,3))
#' 
#' timepoint_expert <- 14
#' data2 <- data %>% rename(status = censored) %>% mutate(time2 = ifelse(time > 10, 10, time),
#'                                                        status2 = ifelse(time> 10, 0, status))
#' example1_mle <- fit.models.expert(formula=Surv(time2,status2)~1,data=data2,
#'                                   distr=c("wph", "exp"),
#'                                   method="mle",
#'                                   pool_type = "log pool",
#'                                   opinion_type = "survival",
#'                                   times_expert = timepoint_expert,
#'                                   param_expert = param_expert_example1)
#' 
#' 
#' #plot(example1_mle, add.km = TRUE, t = 0:30,plot_opinion = TRUE)
#'
#' @references 
#' \insertRef{Baio.2020}{expertsurv}
#' 
#' @export
plot.expertsurv <- function(...) {
  
  # Collects all the extra arguments
  exArgs=list(...)
  
  # Finds out whether there are objects with no name (if so, they will be 'survHE' objects!)
  # If there are any, then needs to rename them to make the rest of the function work
  if(length(names(exArgs))==0) {
    # This is the case where the only argument(s) is/are unnamed 'survHE' object(s)
    names(exArgs)=paste0("Object",1:length(exArgs))
  }
  if(length(which(names(exArgs)==""))>0){
    names(exArgs)[which(names(exArgs)=="")] = paste0("Object",1:length(which(names(exArgs)=="")))
  }
  
  # The default is to go with the 'ggplot' version of the graph. 
  if (exists("graph",exArgs)) {graph=exArgs$graph} else {graph="ggplot"}

  # If so, then call the function 'plot_ggplot_survHE
  if(graph=="ggplot") {
    return(plot_ggplot_expertsurv(exArgs))
  }

  # If the user selects 'base' (only for back-compatibility), then runs the old code
  ### NB: Do I want this? (probably not...)
  #if(graph=="base") {
  #  do.call(plot_base_expertsurv,exArgs)
  #}
}
