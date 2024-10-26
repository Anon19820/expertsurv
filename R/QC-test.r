#' Title of the Function
#'
#' Detailed description of what the function does. 
#'
#' @examples
#' \dontrun{
#' library("expertsurv")
#' require("dplyr")
#' 
#' param_expert_example1 <- list()
#' param_expert_example1[[1]] <- data.frame(dist = c("norm"),
#'                                          wi = c(1), # Ensure Weights sum to 1
#'                                          param1 = c(0.1),
#'                                          param2 = c(0.01),
#'                                          param3 = c(NA))
#' timepoint_expert <- 14
#' data2 <- data %>% rename(status = censored) %>% mutate(time2 = ifelse(time > 10, 10, time),
#'                                                        status2 = ifelse(time > 10, 0, status))
#' example1 <- fit.models.expert(formula = Surv(time2, status2) ~ 1, data = data2,
#'                               distr = c("wph", "exp", "gomp"),
#'                               method = "mle",
#'                               pool_type = "log pool",
#'                               opinion_type = "survival",
#'                               times_expert = timepoint_expert,
#'                               param_expert = param_expert_example1)
#' 
#' #plot(example1, add.km = TRUE, t = 0:30, plot_opinion = TRUE)
#' 
#' pars <- example1$models$`Weibull (PH)`$res[, 1]
#' LL_data1 <- sum(hweibullPH(data2$time2[data2$status2 == 1], shape = pars[1], scale = pars[2], log = TRUE)) +
#'             sum(pweibullPH(data2$time2, pars[1], pars[2], log = TRUE, lower.tail = FALSE))
#' St_expert <- pweibullPH(timepoint_expert, pars[1], pars[2], log = FALSE, lower.tail = FALSE)
#' LL_expert1 <- dnorm(St_expert, mean = 0.1, sd = 0.01, log = TRUE)
#' LL_data1 + LL_expert1
#' 
#' pars <- example1$models$Gompertz$res[, 1]
#' LL_data2 <- sum(hgompertz(data2$time2[data2$status2 == 1], pars[1], pars[2], log = TRUE)) +
#'             sum(pgompertz(data2$time2, pars[1], pars[2], log = TRUE, lower.tail = FALSE))
#' St_expert <- pgompertz(timepoint_expert, pars[1], pars[2], log = FALSE, lower.tail = FALSE)
#' LL_expert2 <- dnorm(St_expert, mean = 0.1, sd = 0.01, log = TRUE)
#' LL_data2 + LL_expert2
#' }
#' @noRd
my_function <- function() {
  # Function code
}
