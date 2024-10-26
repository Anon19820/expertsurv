utils::globalVariables(c("(Intercept)", "expert",'times_expert','weights','sum_weights','v_array','p_mat','distr3',".",
                         'S','model_name','strata','object_name',
                         'low','upp', 'lower', 'uppper','value',
                         'object_name','model_name',
                         'xmin','xmax','ymin','ymax','plot_base_expertsurv','xlim','strata',
                         'event','upper', 'lab','event','model',
                         'Survival', 'Time', 'X2.5.', 'X50.', 'X97.5.', 'actionButton', 'arm', 'checkboxInput', 'column',
                         'downloadButton', 'downloadHandler', 'eventReactive', 'fileInput', 'fluidPage',
                         'fluidRow', 'ftype', 'fx', 'h3', 'head', 'helpText', 'hr', 'mainPanel', 'need', 'numericInput',
                         'observeEvent', 'p', 'plotOutput', 'reactive', 'reactiveValues', 'renderPlot',
                         'selectInput', 'sidebarPanel', 'tabPanel', 'tabsetPanel', 'textInput', 'titlePanel',
                         'updateSelectInput', 'validate', 'varSelectInput', 'wellPanel', 'x', 'y','stats', '.checkMFClasses', '.getXlevels', 'BIC',
						 'as.formula', 'coef', 'confint', 'contrasts<-', 'dbeta',
						'delete.response', 'dgamma', 'dlnorm', 'dnorm', 'dt',
						'dweibull', 'formula', 'integrate', 'median',
						'model.extract', 'model.frame', 'model.matrix', 'na.pass',
						'optim', 'pbeta', 'pexp', 'pgamma', 'plnorm', 'plogis',
						'pnorm', 'pt', 'pweibull', 'qbeta', 'qexp', 'qf', 'qgamma',
						'qlnorm', 'qlogis', 'qnorm', 'qt', 'quantile', 'qweibull',
						'reformulate', 'rf', 'rgamma', 'rlnorm', 'runif',
						'rweibull', 'sd', 'setNames', 'terms', 'time', 'uniroot',
						'update', 'var', 'vcov', 'weighted.mean','scaling_fac','fx_final'))
						 



#' Title
#'
#' @param x 
#'
#' @return x
#' @importFrom Rdpack reprompt
#' @noRd
#' @import broom
cran_req <- function(x){
  return(x)
}
						 
						   

