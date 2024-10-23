#' Fitting Parametric Survival models with Expert Opinion
#'
#' Implementation of survival models with expert opinion on the survival probabilities or expected difference in survival.
#' Function is equivalent to the `fit.models` in `survHE` except for the inclusion of the "expert_type" and "param_expert" arguments. 
#' Worked examples can be found in the [README](https://github.com/Philip-Cooney/expertsurv/blob/master/README.md) file.
#' Note that the default method is "bayes", however, the user may use "mle"  (method "inla" is not included).
#'
#' @param formula As per `fit.models` in `survHE`
#' @param data As per `fit.models` in `survHE`
#' @param distr As per `fit.models` in `survHE`. Note Generalized F model is not available for method = "bayes" nor Royston-Parmar available with opinion on the mean survival.
#' @param method As per `fit.models` in `survHE`. (except for the inla method). It should be noted that serval models are fit using stan and others by JAGS, therefore the method argument for a Bayesian analysis is "bayes" and not "hmc" as in survHE.
#' @param expert_type Either "survival", which indicates expert opinion on the survival function or "mean" (actually anything that does not contain "survival") which represents a belief on difference in survival.
#' @param param_expert A list containing a dataframe for each timepoint (if applicable). Each dataframe should have columns with the following names and each row representing an expert:
#'  \itemize{
#'   \item \strong{dist}: Names of the distribution assigned to each expert which may be "norm", "t", "lnorm", "gamma", "beta".
#'   \item \strong{wi}: Weight of the expert, must sum to 1.
#'   \item \strong{param1}: First parameter of the distribution (e.g. mean for norm distribution). Parameters as per `SHELF` package. 
#'   \item \strong{param2}: Second parameter of the distribution.
#'   \item \strong{param3}: Third parameter of the distribution (NA expect for degrees of freedom for t distribution)
#' }
#' @param id_St This is only necessary if the model includes covariates (usually treatments) and expert opinion on survival probabilities is being integrated into the analysis. It specifies the row number in the data frame, indicating that the covariate pattern for this row represents the group for which the expert opinion is provided.
#' @param id_trt This is only necessary if including expert opinion about differences in expected survival (i.e. area under the curve) between a treatment and a comparator and specifies the row number in the data frame corresponding to a treated patient. Similarly the location of the reference id_comp (i.e. comparator) needs to be specified.  
#' @param ... Other arguments may be required depending on the example. See \href{../README.md}{README} for details and further examples.
#' @return An object of class ``expertsurv`` which contains the parameters of the models estimated with expert opinion.
#' @importFrom magrittr %>%
#' @keywords models
#' @export
#' @md
#' 
#' @examples 
#' require("dplyr")
#' #Expert Opinion as a normal distribution centered on 0.1 with sd 0.005
#' param_expert_example1 <- list()
#' param_expert_example1[[1]] <- data.frame(dist = c("norm"),
#'                                          wi = c(1), # Ensure Weights sum to 1
#'                                          param1 = c(0.1),
#'                                          param2 = c(0.05),
#'                                         param3 = c(NA))
#' timepoint_expert <- 14 # Expert opinion at t = 14
#' 
#' data2 <- data %>% rename(status = censored) %>% mutate(time2 = ifelse(time > 10, 10, time),
#' status2 = ifelse(time> 10, 0, status))
#' 
#' example1  <- fit.models.expert(formula=Surv(time2,status2)~1,data=data2,
#'                               distr=c("wei", "gomp"),
#'                               method="mle",
#'                               opinion_type = "survival",
#'                               times_expert = timepoint_expert, 
#'                               param_expert = param_expert_example1)
#'                               
#' plot(example1, add.km = TRUE, t = seq(0:20)) #Plot Survival
#' model.fit.plot(example1, type = "aic")  #Plot AIC 
#'
#' example1_bayes  <- fit.models.expert(formula=Surv(time,status)~as.factor(arm),data=data2,
#'                                     distr=c("wph", "lnorm","exp"),
#'                                     method="bayes",
#'                                     pool_type = "log pool", 
#'                                     opinion_type = "survival",
#'                                     id_St = min(which(data2$arm ==0)), #We want our opinion to refer to the treatment called "0"
#'                                     times_expert = timepoint_expert, 
#'                                     param_expert = param_expert_example1)
#'
#'
                               
fit.models.expert <- function(formula = NULL, data, distr = NULL, method = "bayes", 
                              expert_type = "survival", param_expert = NULL, ...){
  exArgs <- list(...)
  

  #will need to be modified
  exArgs$formula <- formula
  exArgs$data = data
  exArgs$param_expert <- param_expert
  if(!is.null(expert_type) && method == "inla"){
    warning("Expert Opinion is not implemented with the inla method")
    stop()
  }
  
  if(!is.null(expert_type) && is.null(param_expert)){
    warning("You have not specified any expert opinions using the param_expert argument - Evaluating survival model without expert opinion")
    exArgs$times_expert <- 1
    param_expert_vague <- list()
    param_expert_vague[[1]] <- data.frame(dist = "beta", wi = 1, param1 = 1, param2 = 1, param2 = NA)
    param_expert <- param_expert_vague
    exArgs$param_expert <- param_expert
    exArgs$opinion_type <- "survival"
  }
  
  if(!is.null(expert_type) && expert_type != "survival" && any(distr == "rps")){
    warning("Mean Difference is not implemented for RPS models")
    stop()
  }
  
  if(method == "bayes" && any(distr == "genf")){
    warning("Generalized F models are implemented")
    stop()
  }
  
  
  if(is.null(exArgs$pool_type)){
    nrow_vec <- rep(NA,length(param_expert))
    for(i in 1:length(param_expert)){
      nrow_vec[i] <- nrow(param_expert[[i]])
    }
    
    if(any(nrow_vec>1)){
      warning("Assuming Linear pooling for the multiple expert opinions")
    }
    exArgs$pool_type <-"linear pool"
  }
  
  
  
  
  fit.models(formula = formula, data = data, distr = distr, method = method, exArgs = exArgs)
}


fit.models <- function (formula = NULL, data, distr = NULL, method = "mle", exArgs, 
                        ...){
  
  if (is.null(formula)) {
    stop("You need to specify a model 'formula', e.g. 'formula=Surv(time,event)~treat'")
  }
  method <- tolower(method)
  if (!method %in% c("bayes", "inla", "mle")) {
    stop("Methods available for use are 'mle', 'bayes' or 'inla'")
  }
  check_distributions(method, distr)
  if (method == "mle") {
    
    res <- format_output_fit.models(lapply(distr, function(x)runMLE(x, 
                                                                    exArgs)), method, distr, formula, data)
  }
  if (method == "inla") {
    
    stop("INLA is not implemented in expertsurv")
    
  }
  if (method == "bayes") {
    if(any(distr %in% c("gam", "gomp", "gga"))){
      
      if (!isTRUE(requireNamespace("rjags", quietly = TRUE))|!isTRUE(requireNamespace("R2jags", quietly = TRUE))) {
        stop("You need to install the R packages 'rjags' and 'R2jags' along with JAGS")
      }
    }else{
      
      if (!isTRUE(requireNamespace("rstan", quietly = TRUE))|!isTRUE(requireNamespace("rstantools", quietly = TRUE))) {
        stop("You need to install the R packages 'rstan' and 'rstantools' to run the models evaluated by Stan")
      }
      
    }
    res <- format_output_fit.models(lapply(distr, function(x) runBAYES(x, 
                                                                     exArgs)), method, distr, formula, data)
  }
  
  return(res)
}

#' Helper function to run the survival models using Bayesian inference (rstan or JAGS)
#' for a given formula and dataset
#' 
#' @param x a (vector of) string(s) containing the name(s) of the model(s)
#' to be fitted
#' @param exArgs a list of extra arguments passed from the main 'fit.models' 
#' function
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Hamiltonian Monte Carlo
#' @noRd 
runBAYES <- function (x, exArgs){
  if (!isTRUE(requireNamespace("rstan", quietly = TRUE))) {
    stop("You need to install the R package 'rstan'. Please run in your R terminal:\n install.packages('rstan')")
  }
  formula <- exArgs$formula
  data = exArgs$data
  availables <- load_availables()
  d3 <- manipulate_distributions(x)$distr3
  method <- "bayes"
  if (exists("chains", where = exArgs)) {
    chains <- exArgs$chains
  }
  else {
    chains <- 2
  }
  if (exists("iter", where = exArgs)) {
    iter <- exArgs$iter
  }
  else {
    iter <- 2000
  }
  if (exists("warmup", where = exArgs)) {
    warmup <- exArgs$warmup
  }
  else {
    warmup <- floor(iter/2)
  }
  if (exists("thin", where = exArgs)) {
    thin <- exArgs$thin
  }
  else {
    thin <- 1
  }
  if (exists("control", where = exArgs)) {
    control <- exArgs$control
  }
  else {
    control <- list(NULL)
  }
  if (exists("seed", where = exArgs)) {
    seed <- exArgs$seed
  }
  else {
    seed <- sample.int(.Machine$integer.max, 1)
  }
  if (exists("pars", where = exArgs)) {
    pars <- exArgs$pars
  }
  else {
    pars <- c("lambda_cens", "lambda_obs", "cens", 
              "d", "lp__", "loglambda_cens", 
              "loglambda_obs", "mu", "logP", 
              "linpred")
  }
  if (exists("include", where = exArgs)) {
    include <- exArgs$include
  }
  else {
    include <- FALSE
  }
  if (exists("k", where = exArgs)) {
    k <- exArgs$k
  }
  else {
    k <- 0
  }
  if (exists("cores", where = exArgs)) {
    cores <- exArgs$cores
  }
  else {
    cores <- 1
  }
  
  if (exists("iter_jags", where = exArgs)) {
    iter_jags <- exArgs$iter_jags
  }
  else {
    iter_jags <- iter*5
  }
  if (exists("init", where = exArgs)) {
    if(d3%in% names(exArgs$init) ){
      init <- exArgs$init[[d3]] 
    }else{
      init = "random"
    }
  }
  else {
    init = "random"
  }
  if (exists("save.stan", where = exArgs)) {
    save.stan <- exArgs$save.stan
    if(is.null(exArgs$save.stan_path)){
      stop("You must specify the path to save the model files using the save.stan_path argument")
    }else{
      save.stan.path <- exArgs$save.stan_path
    }
    
  }
  else {
    save.stan = FALSE
  }
  if (exists("refresh", where = exArgs)) {
    refresh = exArgs$refresh
  }else {
    refresh = max(iter/10, 1)
  }
  
  if(exists("compile_mods", where = exArgs)){
    compile_mods = exArgs$compile_mods
  }else{
    compile_mods = NULL
  }
  
  
  d <- names(availables[[method]][match(d3, availables[[method]])])
  data.stan <- make_data_stan(formula, data, d3, exArgs)
  
  tic <- proc.time()
  
  if (d3 %in% c("gam", "gga", "gom")){
    data.jags <- data.stan
    if(d3 %in% c( "gom")){
      parameters.to.save_jags = c("alpha","beta", "rate")
      
      #Inits as per flexsurvreg (reparameterized)
      modelinits <- function(){
        beta = c(log(1/mean(data.jags$t)*stats::runif(1,0.8,1.5)),rep(0,data.jags$H -1))
        list(alpha1 = stats::runif(1,0.001,0.003),alpha2 = stats::runif(1,0.001,0.003), beta = beta) 
      }
      
    }else if(d3 == "gga"){ #(d3 == "gga")
      parameters.to.save_jags = c("Q","sigma", "beta", "r", "b","mu")
      tinits1 <-data.jags$t + max(data.jags$t)
      is.na(tinits1)<-data.jags$d ==1
      data.jags$is.censored <- ifelse(data.jags$d==0, 1, 0)
      data.jags$t_jags <- ifelse(data.jags$is.censored ==1, NA, data.jags$t) 
      data.jags$t_cen <- data.jags$t+data.jags$d
      modelinits <- function(){list(t_jags = tinits1)}
      #Stop JAGS Warning messages
      data.jags <- data.jags[names(data.jags) %!in% c("t", "d", "a0")]
      
      
    }else{ #"gam",
      parameters.to.save_jags = c("alpha","beta", "rate")
      modelinits <- NULL
    }
    data.jags <- data.jags[names(data.jags) %!in% "max_param"]
    
    message(paste0(" \n SAMPLING FOR MODEL '",d,"_expert' NOW.  \n"))
    suppressWarnings({
      model <-R2jags::jags(model.file = textConnection(get(paste0(d,".jags"))),
                           data=data.jags,
                           n.chains=chains,
                           inits=modelinits,
                           parameters.to.save = c(parameters.to.save_jags,"St_expert"),
                           n.iter = iter_jags,
                           n.thin = thin,
                           n.burnin = iter,
                           jags.module = c("glm","dic"))
    })
    
    
  }else{
    # dso <- textConnection(get(paste0(d, "_expert")))
    # model <- rstan::sampling(dso, data.stan, chains = chains, 
    #                          iter = iter, warmup = warmup, thin = thin, seed = seed, 
    #                          control = control, pars = pars, include = include, cores = cores, 
    #                          init = init, refresh = refresh)
    browser()
    if(is.na(match(paste0(d, "_expert"), names(compile_mods)))){
      stan_code <- get(paste0(d, "_expert"))
      stan_model <- rstan::stan_model(model_code = stan_code, model_name  = paste0(d, "_expert"))
    }else{
      stan_model <- compile_mods[[paste0(d, "_expert")]]
    }

    model <- rstan::sampling(stan_model, data.stan, chains = chains, 
                             iter = iter, warmup = warmup, thin = thin, seed = seed, 
                             control = control, pars = pars, include = include, cores = cores, 
                             init = init, refresh = refresh)
    
    time_stan <- sum(rstan::get_elapsed_time(model))
    
  }
  
  toc <- proc.time() - tic
  time_survHE <- toc[3]
  ics <- compute_ICs_stan(model, d3, data.stan)
  
  if (save.stan) {
    if(d3 %in% c("gam", "gga", "gom")){
      
      model_code <- get(paste0(d,".jags"))
      con <- paste0(save.stan.path, d, ".txt")
    }else{
      model_code <- attr(model@stanmodel, "model_code")
      con <- paste0(save.stan.path,d, ".stan")
      
    }
    
    writeLines(model_code, con = con)
    message(paste0("Model code saved to the file: ", con, 
                   "\n"))
    
    ## Add in for Jags
  }
  model_name <- d3
  
  list(model = model, aic = ics$aic, bic = ics$bic, dic = ics$dic, 
       dic2 = ics$dic2,waic = ics$waic, pml = ics$pml,  time2run = time_survHE, 
       data.stan = data.stan, save.stan = save.stan, model_name = model_name)
}




compile_stan <- function(dist_stan = c("exp","wei","wph","rps","llo","lno")){
  availables_all <- expertsurv:::load_availables()[["bayes"]]
  availables_stan <- availables_all[match(dist_stan,availables_all)]
  names_stan <- names(availables_stan)
  list_stan_model <- list()
  for(i in 1:length(names_stan)){
    mode_name_temp <- paste0(names_stan[i],"_expert")
    stan_code_temp <- get(mode_name_temp)
    list_stan_model[[mode_name_temp]] <- rstan::stan_model(model_code = stan_code_temp, model_name = mode_name_temp)
  }
  
  return(list_stan_model)
  
}


#' Helper function to create data in the correct format for rstan
#' 
#' @param formula a formula specifying the model to be used, in the form
#' \code{Surv(time,event)~treatment[+covariates]} in flexsurv terms, or
#' \code{inla.surv(time,event)~treatment[+covariates]} in INLA terms.
#' @param data A data frame containing the data to be used for the analysis.
#' This must contain data for the 'event' variable. In case there is no
#' censoring, then \code{event} is a column of 1s.
#' @return \item{data.stan}{A list containing the variables needed to pass
#' to 'stan' when calling \code{fit.models} with \code{method="bayes"}}.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Bayesian inference Frequentist inference Expert Opinion
#' @import tibble
#' @importFrom Rdpack reprompt
#' @importFrom utils getFromNamespace
#' @import dplyr
#' @import stats 
#' @import survival
#' @import graphics
#' @importFrom stringr str_replace_all
#' @noRd 
make_data_stan <- function (formula, data, distr3, exArgs = globalenv()){
  availables <- load_availables()
  method <- "bayes"
  formula_temp <- stats::update(formula, paste(all.vars(formula, data)[1], 
                                               "~", all.vars(formula, data)[2], "+."))
  mf <- tibble::as_tibble(stats::model.frame(formula_temp, data)) %>% 
    dplyr::rename(time = 1,event = 2) %>% dplyr::rename_if(is.factor, .funs = ~gsub("as.factor[( )]","", .x)) %>% 
    dplyr::rename_if(is.factor, .funs = ~gsub("[( )]","", .x)) %>% 
    dplyr::bind_cols(tibble::as_tibble(stats::model.matrix(formula_temp,data)) %>% dplyr::select(contains("Intercept"))) %>%
    dplyr::select(time,event, contains("Intercept"), everything()) %>% tibble::rownames_to_column("ID")
  
  ####Code Change Here
  ######
  if (distr3 %!in% c("rps")) {
    data.stan <- list(t = (mf$time), d = mf$event, n = nrow(mf), 
                      X = matrix(stats::model.matrix(formula, data), nrow = nrow(mf)), 
                      H = ncol(stats::model.matrix(formula, data)))
    if (data.stan$H == 1) {
      data.stan$X <- cbind(data.stan$X, rep(0, data.stan$n))
      data.stan$H <- ncol(data.stan$X)
    }
  }
  if (distr3 == "rps") {
    if (exists("k", where = exArgs)) {
      
      k <- exArgs$k
    }
    else {
      k <- 0
    }
    knots <- quantile(log((mf %>% filter(event == 1))$time), 
                      seq(0, 1, length = k + 2))
    B <- flexsurv::basis(knots, log(mf$time))
    B_expert <- flexsurv::basis(knots, log(exArgs$times_expert))
    DB <- flexsurv::dbasis(knots, log(mf$time))
    mm <- stats::model.matrix(formula, data)[, -1]
    if (length(mm) < 1) {
      mm <- matrix(rep(0, nrow(mf)), nrow = nrow(mf), ncol = 2)
    }
    if (is.null(dim(mm))) {
      mm <- cbind(mm, rep(0, length(mm)))
    }
    data.stan <- list(t = mf$time, d = mf$event, n = nrow(mf), 
                      M = k, X = mm, H = ncol(mm), B = B, DB = DB, mu_gamma = rep(0,k + 2),
                      sigma_gamma = rep(5, k + 2), knots = knots, B_expert = B_expert)
  }
  data.stan$mu_beta = rep(0, data.stan$H)
  if (distr3 %in% c("lno")) {
    
    data.stan$sigma_beta <- rep(100, data.stan$H)
  }
  data.stan$sigma_beta <- rep(5, data.stan$H)
  
  if (distr3 %in% c("gam","gom", "gga", "llo", "wei", 
                    "wph")) {
    data.stan$a_alpha = data.stan$b_alpha = 0.1
  }else if(distr3 %in% c("lno")){
    data.stan$a_alpha = 0
    data.stan$b_alpha = 5
  }
  d <- names(availables[[method]][match(distr3, availables[[method]])])
  priors <- list()
  if (exists("priors", where = exArgs)) {
    abbrs = manipulate_distributions(names(exArgs$priors))$distr3
    pos = grep(distr3, abbrs)
    if (length(pos) > 0) {
      priors = exArgs$priors[[pos]]
    }
  }
  if (!is.null(priors$mu_beta)) {
    data.stan$mu_beta = priors$mu_beta
  }
  if (!is.null(priors$sigma_beta)) {
    data.stan$sigma_beta <- priors$sigma_beta
  }
  if (!is.null(priors$mu_gamma) & distr3 == "rps") {
    data.stan$mu_gamma <- priors$mu_gamma
  }
  if (!is.null(priors$sigma_gamma) & distr3 == "rps") {
    data.stan$sigma_gamma <- priors$sigma_gamma
  }
  if (!is.null(priors$a_sigma)) {
    data.stan$a_sigma = priors$a_sigma
  }
  if (!is.null(priors$b_sigma)) {
    data.stan$b_sigma = priors$b_sigma
  }
  if (!is.null(priors$mu_P)) {
    data.stan$mu_P = priors$mu_P
  }
  if (!is.null(priors$sigma_P)) {
    data.stan$sigma_P = priors$sigma_P
  }
  if (!is.null(priors$mu_Q)) {
    data.stan$mu_Q = priors$mu_Q
  }
  if (!is.null(priors$sigma_Q)) {
    data.stan$sigma_Q = priors$sigma_Q
  }
  if (!is.null(priors$a_alpha)) {
    data.stan$a_alpha = priors$a_alpha
  }
  if (!is.null(priors$b_alpha)) {
    data.stan$b_alpha = priors$b_alpha
  }
  
  
  if(exArgs$opinion_type == "survival"){
    data.stan$St_indic <- 1
    #even if survival need to define these (just put as 1)
    data.stan$id_comp <- 1
    data.stan$id_trt <- 1
    
    if(is.null(exArgs$id_St)){
      data.stan$id_St <- 1
    }else{
      data.stan$id_St <- exArgs$id_St
    }
    
  }else{
    data.stan$St_indic <- 0
    #even if survival need to define these (just put as 1)
    data.stan$id_St <- 1
    
    data.stan$id_trt <- exArgs$id_trt
    data.stan$id_comp <- exArgs$id_comp
    
    if(is.null(exArgs$id_trt|exArgs$id_comp)){
      message("You need to supply the location within the dataframe row number of a treatment and a comparator arm to arguments id_trt and id_comp")
      stop()
    }
    
    
  }
  
  # if(ncol(mf) == 4){
  #   #No covariates
  #   # Has to be opinion_type survival 
  #   data.stan$id_St <- 1
  #   
  # }else if(ncol(mf) == 5){
  #   
  #   if(exArgs$opinion_type == "survival"){
  #     data.stan$id_St <- min(which(mf[,5] == exArgs$id_St))
  #   }else{# Survival Difference
  #     data.stan$id_trt <- min(which(mf[,5] == exArgs$id_trt)) 
  #     if(length(unique(mf[,5] %>% pull()))==2){
  #       data.stan$id_comp <- min(which(mf[,5] != exArgs$id_trt)) 
  #     }else{
  #       data.stan$id_comp <- min(which(mf[,5] == exArgs$id_comp))  
  #     }
  #     
  #   }
  #   #put the number in  could put in a combination of numbers
  # }else{
  #   message("We do not allow more than one covariate (i.e. treatment) in the analysis - although it is technically possible")
  #   stop()
  # }
  
  
  
  
  param_expert <- exArgs$param_expert
  n.experts <- c()
  
  for(i in 1:length(param_expert)){
    n.experts <- c(n.experts, nrow(param_expert[[i]])) 
  }
  
  data_dist_ind <- num_param <- matrix(-999.2,nrow = max(n.experts), ncol =  length(param_expert))
  expert.array <- array(-999.2,dim = c(max(n.experts),5,length(param_expert))) 
  
  for(i in 1:length(param_expert)){
    lk_up_dist <- c("norm", "t", "gamma", "lnorm","beta")
    dist_fit <- param_expert[[i]][,1]
    if(length(dist_fit) - length(expert.array[,1,i])){
      dist_fit <- c(dist_fit, rep(-999.2,length(dist_fit) - length(expert.array[,1,i])))
    }
    expert.array[,1,i] <- as.numeric(sapply(dist_fit, function(x){which(x==lk_up_dist)}))
    weight_vec <- param_expert[[i]][,2]
    expert.array[1:length(weight_vec),2,i] <- weight_vec
    expert.array[1:nrow(param_expert[[i]][,3:5]),3:5,i] <- as.matrix(param_expert[[i]][,3:5])
  }
  
  
  #Stan does not allow NA
  expert.array[is.na(expert.array)] <- -999.2
  
  if(!is.null(exArgs$times_expert)){
    data.stan$n_time_expert <- length(exArgs$times_expert)
    data.stan$time_expert <- as.array(exArgs$times_expert)
  }else{
    data.stan$n_time_expert <- 1
    data.stan$time_expert <- numeric(0) #This produces an array of size 0
    #https://dev.to/martinmodrak/optional-parametersdata-in-stan-4o33
    if (distr3 %in% c("gam", "gga", "gom")){
      data.stan$time_expert <- 1 # Has to be defined for JAGS
      
    }
    
    
  }
  
  data.stan$param_expert <-expert.array
  data.stan$n_experts <- as.array(n.experts)  
  
  if(is.null(exArgs$pool_type)){
    
    data.stan$pool_type <- 1
    
  }else{
    data.stan$pool_type <- as.numeric(grepl("line", exArgs$pool_type)) 
    
  }
  
  if(data.stan$pool_type == 0){
    k_norm <- rep(NA,data.stan$n_time_expert)
	
   for(i in 1:data.stan$n_time_expert){
	
	if(data.stan$n_experts[i] ==1){
	k_norm[i] <- 1
	}else{
      
      param_expert[[i]]$dist <- stringr::str_replace_all(param_expert[[i]]$dist, "normal", "norm") 
      param_expert[[i]]$dist <- stringr::str_replace_all(param_expert[[i]]$dist, "lognorm", "lnorm") 
      
      if(data.stan$St_indic==1){
        min_quant <- 0
        max_quant <- 1
      }else{
        quant.vec <- t(apply(param_expert[[i]], 1, function(x){get_quant_val(
          dist = x["dist"],
          param1 = x["param1"],
          param2 = x["param2"],
          param3 = x["param3"],
          probs = c(0.001,0.025,0.5,0.975,0.999))}))
        
        central.cauchy <- mean(quant.vec[,3])#mean
        sd.cauchy <- max(apply(quant.vec,1, function(x){(x[4]-x[2])/4})) #sd
        min_quant <- min(quant.vec)
        max_quant <- max(quant.vec)
      }
      
      
      
      x.eval <- seq(min_quant, max_quant, length.out = 100)
      dens.eval <- eval_dens_pool(x.eval,param_expert[[i]],pool_type = "log pool",St_indic =data.stan$St_indic)
      k_norm[i] <- integrate.xy(x = x.eval,fx = dens.eval)
    }
	
	}
    data.stan$k_norm <- k_norm
    
  }
  
  #Power prior
  
  if(!is.null(exArgs$a0)){
    data.stan$a0 <- exArgs$a0
  }else{
    data.stan$a0 <- rep(1, nrow(data))
  }
  
  data.stan	
}


#' Helper function to compute the information criteria statistics
#' when using bayes as the inferential engine. 'rstan' does not do
#' DIC automatically and AIC/BIC are also not standard for Bayesian
#' models, so can compute them post-hoc by manipulating the 
#' likelihood functions.
#' 
#' @param model The 'rstan' object with the model fit
#' @param distr3 The 'rstan' object with the model fit
#' @importFrom loo loo
#' @return \item{list}{A list containing the modified name of the 
#' distribution, the acronym (3-letters abbreviation), or the
#' labels (humane-readable name)}.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Bayesian inference via Hamiltonian
#' Monte Carlo Bayesian inference via Integrated Nested Laplace Approximation
#' @noRd 
compute_ICs_stan <- function (model, distr3, data.stan){
  if (distr3 %!in% c("gam", "gga", "gom")) {
    beta <- rstan::extract(model)$beta
  }
  else {
    beta <- model$BUGSoutput$sims.matrix[, grep("beta", 
                                                colnames(model$BUGSoutput$sims.matrix))]
  }
  beta.hat <- apply(beta, 2, stats::median)
  linpred <- beta %*% t(data.stan$X)
  linpred.hat <- beta.hat %*% t(data.stan$X)
  model.eval <- paste0("lik_", distr3)
  out = do.call(what = eval(parse(text = model.eval)), args = list(distr3, 
                                                                   linpred, linpred.hat, model, data.stan))
  logf = out$logf
  logf.hat = out$logf.hat
  npars = out$npars
  logf_comb <- matrix(nrow = nrow(logf), ncol = ncol(logf))
  for (i in 1:nrow(logf)) {
    logf_comb[i, ] <- logf[i, ] + out$logf.expert[i]/ncol(logf)
  }
  tryCatch(suppressWarnings(WAIC <- loo::loo(logf_comb)[["estimates"]][grep("looic", 
                                                                            rownames(loo::loo(logf_comb)[["estimates"]])), "Estimate"]), 
           error = function(e) message("Cannot Evaluate WAIC"))
  if (!exists("WAIC")) {
    WAIC <- Inf
  }
  PML <- -2 * sum(log(nrow(logf_comb)/colSums(1/exp(logf_comb))))
  loglik <- apply(logf, 1, sum) + out$logf.expert
  loglik.bar <- apply(logf.hat, 1, sum) + out$logf.hat.expert
  D.theta <- -2 * loglik
  D.bar <- -2 * loglik.bar
  pD <- mean(D.theta) - D.bar
  if (is.nan(pD)) {
    warning(paste0("pD is not defined for ", distr3, "; DIC estimates invalid, use WAIC or PML."))
    pD <- 0
  }
  if (pD < 0) {
    warning(paste0("pD is ", round(pD), " for ", distr3, 
                   "; DIC estimates unreliable, use WAIC or PML."))
  }
  pV <- 0.5 * stats::var(D.theta)
  dic <- mean(D.theta) + pD
  dic2 <- mean(D.theta) + pV
  aic <- D.bar + 2 * npars
  bic <- D.bar + npars * log(data.stan$n)
  log_lik_data_avg <- mean(apply(logf, 1, sum))
  log_lik_expert_avg <- mean(out$logf.expert)
  list(aic = aic, bic = bic, dic = dic, dic2 = dic2, waic = WAIC, 
       pml = PML, log_lik_data_avg = log_lik_data_avg,log_lik_expert_avg= log_lik_expert_avg)
}


#' Helper function to format the output of the modelling (produced either
#' by running 'runMLE', or 'runINLA', 'runBAYES'), in a way that is consistent
#' with the architecture of 'survHE'
#' 
#' @param output The output of one of the helper functions used to run the
#' models.
#' @param method The method used to do the estimation
#' @param distr The abbreviated name for the distribution to be used
#' @param formula The model formula
#' @param data The dataset used
#' @return \item{res}{A 'survHE' object containing all the relevant output
#' conveniently formatted}.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Bayesian inference via Hamiltonian
#' Monte Carlo Bayesian inference via Integrated Nested Laplace Approximation
#' @noRd 
format_output_fit.models <- function (output, method, distr, formula, data){
  labs <- manipulate_distributions(distr)$labs
  models <- lapply(output, function(x) x$model)
  model.fitting <- list(aic = unlist(lapply(output, function(x) x$aic)), 
                        bic = unlist(lapply(output, function(x) x$bic)), dic = unlist(lapply(output, 
                                                                                             function(x) x$dic)))
  misc <- list(time2run = unlist(lapply(output, function(x) x$time2run)), 
               formula = formula, data = data, model_name = unlist(lapply(output, 
                                                                          function(x) x$model_name)))
  if (any(distr == "polyweibull")) {
    misc$km = lapply(formula, function(f) make_KM(f, data))
  }
  else {
    misc$km = make_KM(formula, data)
  }
  if (method == "bayes") {
    misc$data.stan <- lapply(output, function(x) x$data.stan)
    model.fitting$dic2 <- unlist(lapply(output, function(x) x$dic2))
    model.fitting$waic <- unlist(lapply(output, function(x) x$waic))
    model.fitting$pml <- unlist(lapply(output, function(x) x$pml))
    
  }
  names(models) <- labs
  res <- list(models = models, model.fitting = model.fitting, 
              method = method, misc = misc)
  class(res) <- "expertsurv"
  return(res)
}

make_sim_bayes <- function (m, t, X, nsim, newdata, dist, summary_stat, ...){
  iter_stan <- m@stan_args[[1]]$iter
  beta = rstan::extract(m)$beta
  if (ncol(X) == 1) {
    beta = beta[, 1,drop = F]
  }
  
  linpred <- beta %*% t(X)
  sim <- lapply(1:nrow(X), function(x) {
    do.call(paste0("rescale_bayes_", dist), args = list(m, 
                                                      X, linpred[, x]))
  })
  if (nsim > iter_stan) {
    stop("Please select a value for 'nsim' that is less than or equal to the value set in the call to 'fit.models'")
  }
  if (nsim == 1) {
    sim <- lapply(sim, function(x) as.matrix(tibble::as_tibble(x) %>% 
                                               dplyr::summarise_all(summary_stat), nrow = 1, ncol = ncol(x)))
  }
  if (nsim > 1 & nsim < iter_stan) {
    sim <- lapply(sim, function(x) as.matrix(tibble::as_tibble(x) %>% 
                                               dplyr::sample_n(nsim, replace = FALSE), nrow = nsim, ncol = ncol(x)))
  }
  return(sim)
}






get_Surv <- function(dist, time, param1 = NULL, param2 = NULL, param3 = NULL, log = F, data.stan = NULL){
  
  if(dist == "wei"){
    return(stats::pweibull(time, 
                           shape = param1, scale = param2, lower.tail = F, log = log))
  }
  
  if(dist == "wph"){
    return(flexsurv::pweibullPH(time, 
                                shape = param1, scale = param2, lower.tail = F, log = log))
  }
  
  if(dist == "exp"){
    return(stats::pexp(time, rate = param1, lower.tail = F, log = log))
    
  }
  
  if(dist == "gam"){
    return(stats::pgamma(time, shape = param1, rate = param2, lower.tail =  F, log = log))
  }
  
  if(dist == "gga"){
    return(flexsurv::pgengamma(time, mu = param1, sigma = param2, Q = param3, lower.tail =  F, log = log))
  }
  
  if(dist == "gom"){
    return(flexsurv::pgompertz(time, shape = param1, rate = param2, lower.tail =  F, log = log))
  }
  
  if(dist == "lno"){
    return(stats::plnorm(time, 
                         meanlog  = param1, sdlog  = param2, lower.tail = F, log = log))
  }
  
  if(dist == "llo"){
    return(flexsurv::pllogis(time, 
                             shape  = param1, scale  = param2, lower.tail = F, log = log))
  }
  
  if(dist == "rps"){
    #2 ways to do it -- need to check if it is valid
    #eta <-  param1*data.stan$B_expert[which(time == data.stan$time_expert),] + param2
    #return(exp(-exp(eta)))
    return(flexsurv::psurvspline(time, gamma = param1, knots= data.stan$knots,lower.tail = F, log = log, offset = param2 ))
    
  }
  
  
}


get_mean_diff <- function(dist, time, param1 = NULL, param2 = NULL, param3 = NULL, log = F, data.stan = NULL){
  
  if(dist == "wei"){
    return(flexsurv::mean_weibull(shape = param1, scale = param2[1])-flexsurv::mean_weibull(shape = param1, scale = param2[2]))
  }
  
  if(dist == "wph"){
    return(flexsurv::mean_weibullPH(shape = param1, scale = param2[1])-flexsurv::mean_weibullPH(shape = param1, scale = param2[2]))
  }
  
  if(dist == "exp"){
    return(flexsurv::mean_exp(rate = param1[1])-flexsurv::mean_exp(rate = param1[2]))
    
  }
  
  if(dist == "gam"){
    return(flexsurv::mean_gamma(shape = param1, rate = param2[1])-flexsurv::mean_gamma(shape = param1, rate = param2[2]))
  }
  
  if(dist == "gga"){
    return(flexsurv::mean_gengamma(mu = param1[1], sigma = param2, Q = param3)-flexsurv::mean_gengamma(mu = param1[2], sigma = param2, Q = param3))
  }
  
  if(dist == "gom"){
    return(flexsurv::mean_gompertz(shape = param1, rate = param2[1])-flexsurv::mean_gompertz(shape = param1, rate = param2[2]))
  }
  
  if(dist == "lno"){
    return(flexsurv::mean_lnorm(meanlog  = param1[1], sdlog  = param2)-flexsurv::mean_lnorm(meanlog  = param1[2], sdlog  = param2))
  }
  
  if(dist == "llo"){
    return(flexsurv::mean_llogis(shape  = param1[1], scale  = param2)-flexsurv::mean_llogis(shape  = param1[2], scale  = param2))
  }
  
  # if(dist == "rps"){
  # 
  # }
  
  
}



expert_like <- function(data.stan, dist_surv, param1, param2 =NULL, param3= NULL){
  
  log_lik <- rep(NA, length(data.stan$n_time_expert))
  
  for(i in 1:length(data.stan$n_time_expert)){
    
    n_experts <-  dim(data.stan$param_expert[,,i, drop = F])[1] 
    if(data.stan$St_indic ==1){ #Survival
      
      output <- get_Surv(dist_surv, data.stan$time_expert[i], 
                         param1 =param1, param2 = param2, param3 = param3, data.stan = data.stan)
    }else{# Add code for mean
      
      output <- get_mean_diff(dist_surv,param1=param1, param2=param2, param3=param3, data.stan = data.stan)
      
      
    }
    
    if(n_experts == 1){
      param_expert_curr <- matrix(nrow = 1, data.stan$param_expert[,,i])
    }else{
      param_expert_curr <-  data.stan$param_expert[,,i]
    }
    if(data.stan$pool_type == 0){
      log_lik[i]  <-   expert_log_dens(x = output, df = param_expert_curr, pool_type = data.stan$pool_type, k_norm = data.stan$k_norm[i], St_indic = data.stan$St_indic)
    }else{
      log_lik[i]  <- expert_log_dens(x = output, df = param_expert_curr, pool_type = data.stan$pool_type,St_indic = data.stan$St_indic)
      
    }
    
  }
  
  return(sum(log_lik))
  
}


lik_rps <- function (x, linpred, linpred.hat, model, data.stan){
  
  dist <- "rps"
  gamma <- rstan::extract(model)$gamma
  gamma.hat <- apply(gamma, 2, stats::median)
  linpred.hat <- as.numeric(linpred.hat)
  
  # LL<- apply(gamma_iters, 1, function(x){data.stan$d*log(hsurvspline(data.stan$t, gamma = x, knots = m.all$misc$data.stan[[1]]$knots))+
  #     psurvspline(q = data.stan$t, gamma = x, knots =  m.all$misc$data.stan[[1]]$knots, lower.tail = F, log.p =T)})
  # 
  # 
  # 
  # logf <- data.stan$d * (-log(data.stan$t) + log(gamma %*% 
  #                                                  t(data.stan$DB)) + gamma %*% t(data.stan$B) + linpred) - 
  #   exp(gamma %*% t(data.stan$B) + linpred)
  # logf.hat <- t(data.stan$d * (-log(data.stan$t) + log(data.stan$DB %*% 
  #                                                        gamma.hat) + data.stan$B %*% gamma.hat + linpred.hat) - 
  #                 exp(data.stan$B %*% gamma.hat + linpred.hat))
  
  
  logf.hat <- array(dim =  c(1,dim(linpred)[2]))
  
  
  if(all(data.stan$X ==0)){
    
    logf<- apply(gamma, 1, function(x){data.stan$d*log(flexsurv::hsurvspline(data.stan$t, gamma = x, knots = data.stan$knots))+
        flexsurv::psurvspline(q = data.stan$t, gamma = x, knots =  data.stan$knots, lower.tail = F, log.p =T)})
    logf <- t(logf)
  }else{
    logf <- array(dim = dim(linpred))
    #probably can be optimized
    for(i in 1:nrow(logf)){
      
      for(j in 1:ncol(logf)){
        logf[i,j] <- data.stan$d[j]*log(flexsurv::hsurvspline(data.stan$t[j], gamma = gamma[i,], knots = data.stan$knots, offset = linpred[i,j]))+
          flexsurv::psurvspline(q = data.stan$t[j], gamma = gamma[i,], knots =  data.stan$knots, lower.tail = F, log.p =T, offset = linpred[i,j])
        
      }
      
    }
  }
  
  
  for(i in 1:ncol(logf.hat)){
    logf.hat[i] <- data.stan$d[i]*log(flexsurv::hsurvspline(data.stan$t[i], gamma = gamma.hat, knots = data.stan$knots, offset = linpred.hat[i]))+
      flexsurv::psurvspline(q = data.stan$t[i], gamma = gamma.hat, knots =  data.stan$knots, lower.tail = F, log.p =T, offset = linpred.hat[i])
    
  }
  
  logf.expert <- rep(NA, nrow(linpred))
  
  if(data.stan$St_indic == 1){
    index_vec  <- data.stan$id_St
    for(i in 1:nrow(linpred)){
      logf.expert[i] <-  expert_like(data.stan, dist_surv = dist,param1 = gamma[i,], param2 = linpred[index_vec])
    }
    logf.hat.expert <- expert_like(data.stan, dist_surv = dist, param1 = gamma.hat, param2 = linpred.hat[index_vec])
  }else{
    index_vec <- c(data.stan$id_trt,data.stan$id_comp)
    #Enter code for Difference in survival
  }
  
  
  npars <- length(gamma.hat) + sum(apply(data.stan$X, 2, function(x) 1 - 
                                           all(x == 0)))
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, logf.hat.expert = logf.hat.expert)
}

lik_exp <- function (x, linpred, linpred.hat, model, data.stan){
  dist = "exp"
  
  logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
    data.stan$d * log(flexsurv::hexp(data.stan$t, exp(linpred[i, ]))) + 
      log(1 - stats::pexp(data.stan$t, exp(linpred[i, ])))
  })), nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(data.stan$d * log(flexsurv::hexp(data.stan$t, exp(linpred.hat))) + 
                       log(1 - stats::pexp(data.stan$t, exp(linpred.hat))), nrow = 1)
  
  logf.expert <- rep(NA, nrow(linpred))
  
  
  if(data.stan$St_indic == 1){
    index_vec  <- data.stan$id_St
  }else{
    index_vec <- c(data.stan$id_trt,data.stan$id_comp)
    
  }
  
  for(i in 1:nrow(linpred)){
    logf.expert[i] <-  expert_like(data.stan, dist_surv = dist,param1 = exp(linpred[i,index_vec]))
  }
  logf.hat.expert <- expert_like(data.stan, dist_surv = dist, param1 = exp(linpred.hat[1,index_vec]))
  
  npars <- 1 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, logf.hat.expert = logf.hat.expert)
}

lik_wei <- function (x, linpred, linpred.hat, model, data.stan ){
  
  dist = "wei"
  shape <- alpha <- as.numeric(rstan::extract(model)$alpha)
  shape.hat <- stats::median(shape)
  logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
    data.stan$d * log(flexsurv::hweibull(data.stan$t, shape[i], exp(linpred[i, 
    ]))) + log(1 - stats::pweibull(data.stan$t, shape[i], exp(linpred[i, 
    ])))
  })), nrow = nrow(linpred), byrow = T)
  
  
  logf.hat <- matrix(data.stan$d * log(flexsurv::hweibull(data.stan$t, 
                                                          shape.hat, exp(linpred.hat))) + log(1 - stats::pweibull(data.stan$t, 
                                                                                                                  shape.hat, exp(linpred.hat))), nrow = 1)
  logf.expert <- rep(NA, nrow(linpred))
  
  if(data.stan$St_indic == 1){
    index_vec  <- data.stan$id_St
  }else{
    index_vec <- c(data.stan$id_trt,data.stan$id_comp)
    
  }
  for(i in 1:nrow(linpred)){
    logf.expert[i] <-  expert_like(data.stan, dist_surv = dist, param1 = shape[i], param2 = exp(linpred[i,index_vec]))
  }
  
  logf.hat.expert <- expert_like(data.stan, dist_surv = dist, param1 = shape.hat, param2 = exp(linpred.hat[1,index_vec]))
  
  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  
  
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, logf.hat.expert = logf.hat.expert)
}

lik_lno <- function (x, linpred, linpred.hat, model, data.stan){
  dist = "lno"
  sigma = as.numeric(rstan::extract(model)$alpha)
  sigma.hat = stats::median(sigma)
  logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
    data.stan$d * log(flexsurv::hlnorm(data.stan$t, (linpred[i, 
    ]), sigma[i])) + log(1 - stats::plnorm(data.stan$t, 
                                           (linpred[i, ]), sigma[i]))
  })), nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(data.stan$d * log(flexsurv::hlnorm(data.stan$t, 
                                                        (linpred.hat), sigma.hat)) + log(1 - stats::plnorm(data.stan$t, 
                                                                                                           (linpred.hat), sigma.hat)), nrow = 1)
  logf.expert <- rep(NA, nrow(linpred))
  if (data.stan$St_indic == 1) {
    index_vec <- data.stan$id_St
  }
  else {
    index_vec <- c(data.stan$id_trt, data.stan$id_comp)
  }
  for (i in 1:nrow(linpred)) {
    logf.expert[i] <- expert_like(data.stan, dist_surv = dist, 
                                  param1 = linpred[i, index_vec],
                                  param2 = sigma[i])
  }
  logf.hat.expert <- expert_like(data.stan, dist_surv = dist,
                                 param1 = linpred.hat[1, index_vec],
                                 param2 = sigma.hat)
  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 
                                                               0)))
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, 
       logf.hat.expert = logf.hat.expert)
}


lik_llo <- function (x, linpred, linpred.hat, model, data.stan){
  dist = "llo"
  sigma = as.numeric(rstan::extract(model)$alpha)
  sigma.hat = stats::median(sigma)
  logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
    data.stan$d * log(flexsurv::hllogis(data.stan$t, sigma[i], exp(linpred[i,]))) +
      log(1 - flexsurv::pllogis(data.stan$t, sigma[i], exp(linpred[i, 
      ])))
  })), nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(data.stan$d * log(flexsurv::hllogis(data.stan$t, 
                                                         sigma.hat, exp(linpred.hat))) + log(1 - flexsurv::pllogis(data.stan$t, 
                                                                                                                   sigma.hat, exp(linpred.hat))), nrow = 1)
  
  logf.expert <- rep(NA, nrow(linpred))
  
  if(data.stan$St_indic == 1){
    index_vec  <- data.stan$id_St
  }else{
    index_vec <- c(data.stan$id_trt,data.stan$id_comp)
    
  }
  
  for(i in 1:nrow(linpred)){
    logf.expert[i] <-  expert_like(data.stan, dist_surv = dist, param1 = sigma[i], param2 = exp(linpred[i,index_vec]))
  }
  logf.hat.expert <- expert_like(data.stan, dist_surv = dist, param1 = sigma.hat, param2 = exp(linpred.hat[1,index_vec]))
  
  
  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 
                                                               0)))
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, logf.hat.expert = logf.hat.expert)
}




lik_wph <- function (x, linpred, linpred.hat, model, data.stan){
  dist = "wph"
  shape <- alpha <- as.numeric(rstan::extract(model)$alpha)
  shape.hat = stats::median(shape)
  logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
    data.stan$d * log(flexsurv::hweibullPH(data.stan$t, shape[i], exp(linpred[i, 
    ]))) + log(1 - flexsurv::pweibullPH(data.stan$t, shape[i], 
                                        exp(linpred[i, ])))
  })), nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(data.stan$d * log(flexsurv::hweibullPH(data.stan$t, 
                                                            shape.hat, exp(linpred.hat))) + log(1 - flexsurv::pweibullPH(data.stan$t, 
                                                                                                                         shape.hat, exp(linpred.hat))), nrow = 1)
  
  logf.expert <- rep(NA, nrow(linpred))
  
  if(data.stan$St_indic == 1){
    index_vec  <- data.stan$id_St
  }else{
    index_vec <- c(data.stan$id_trt,data.stan$id_comp)
    
  }
  
  for(i in 1:nrow(linpred)){
    logf.expert[i] <-  expert_like(data.stan, dist_surv = dist, param1 = shape[i], param2 = exp(linpred[i,index_vec]))
  }
  logf.hat.expert <- expert_like(data.stan, dist_surv = dist, param1 = shape.hat, param2 = exp(linpred.hat[1,index_vec]))
  
  
  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 
                                                               0)))
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, logf.hat.expert = logf.hat.expert)
}



lik_gam <- function (x, linpred, linpred.hat, model, data.stan){
  dist = "gam"
  
  shape <- alpha <-  as.numeric(model$BUGSoutput$sims.matrix[ , grep("alpha",colnames(model$BUGSoutput$sims.matrix))])
  shape.hat <- stats::median(shape)
  
  logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
    data.stan$d * log(flexsurv::hgamma(data.stan$t, shape = shape[i], 
                                       rate = exp(linpred[i, ]))) + 
      stats::pgamma(q = data.stan$t,shape[i], rate = exp(linpred[i, ]), lower.tail = F, log = T)
  })), nrow = nrow(linpred), byrow = T)
  
  logf.hat <- matrix(data.stan$d * log(flexsurv::hgamma(data.stan$t, 
                                                        shape.hat, exp(linpred.hat))) + 
                       stats::pgamma(data.stan$t,shape.hat, exp(linpred.hat),lower.tail = F,
                                     log = T), nrow = 1)  
  
  logf.expert <- rep(NA, nrow(linpred))
  
  if(data.stan$St_indic == 1){
    index_vec  <- data.stan$id_St
  }else{
    index_vec <- c(data.stan$id_trt,data.stan$id_comp)
    
  }
  
  for(i in 1:nrow(linpred)){
    logf.expert[i] <-  expert_like(data.stan, dist_surv = dist, param1 = shape[i], param2 = exp(linpred[i,index_vec]))
  }
  logf.hat.expert <- expert_like(data.stan, dist_surv = dist, param1 = shape.hat, param2 = exp(linpred.hat[1,index_vec]))
  
  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, logf.hat.expert = logf.hat.expert)
}


lik_gom <- function (x, linpred, linpred.hat, model, data.stan){
  dist = "gom"
  shape <- alpha <- as.numeric(model$BUGSoutput$sims.matrix[ , grep("alpha",colnames(model$BUGSoutput$sims.matrix))])
  shape.hat = stats::median(shape)
  logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
    data.stan$d * log(flexsurv::hgompertz(data.stan$t, shape = shape[i], 
                                          rate = exp(linpred[i, ]))) + 
      flexsurv::pgompertz(data.stan$t,shape[i], rate = exp(linpred[i, ]), lower.tail = F, log = T)
  })), nrow = nrow(linpred), byrow = T)
  logf.hat <- matrix(data.stan$d * log(flexsurv::hgompertz(data.stan$t,shape.hat, exp(linpred.hat))) + 
                       flexsurv::pgompertz(data.stan$t,shape.hat, exp(linpred.hat), lower.tail = F, log = T), nrow = 1)
  
  
  if(data.stan$St_indic == 1){
    index_vec  <- data.stan$id_St
  }else{
    index_vec <- c(data.stan$id_trt,data.stan$id_comp)
    
  }
  logf.expert <- rep(NA, nrow(linpred))
  
  for(i in 1:nrow(linpred)){
    logf.expert[i] <-  expert_like(data.stan, dist_surv = dist, param1 = shape[i], param2 = exp(linpred[i,index_vec]))
  }
  logf.hat.expert <- expert_like(data.stan, dist_surv = dist, param1 = shape.hat, param2 = exp(linpred.hat[1,index_vec]))
  
  npars <- 2 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, logf.hat.expert = logf.hat.expert)
}



lik_gga <- function (x, linpred, linpred.hat, model, data.stan){
  dist = "gga"
  q = as.numeric(model$BUGSoutput$sims.matrix[ , grep("Q",colnames(model$BUGSoutput$sims.matrix))])
  q.bar = stats::median(q)
  scale = as.numeric(model$BUGSoutput$sims.matrix[ , grep("sigma",colnames(model$BUGSoutput$sims.matrix))])
  scale.bar = stats::median(scale)
  
  d2 <- sapply(data.stan$d,function(x){ifelse(x == 1, 0,1)})
  
  logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
    data.stan$d*flexsurv::dgengamma(data.stan$t, 
                                    linpred[i, ], scale[i], q[i], log = T) +
      d2*flexsurv::pgengamma(data.stan$t,linpred[i, ], scale[i], q[i], log = T, lower.tail = F)})),
    nrow = nrow(linpred), byrow = T)
  
  
  logf.hat <- matrix(data.stan$d*flexsurv::dgengamma(data.stan$t,linpred.hat,scale.bar, q.bar, log = T) +
                       d2*flexsurv::pgengamma(data.stan$t, linpred.hat, scale.bar, q.bar, log = T, lower.tail = F),
                     nrow = 1)
  
  
  if(data.stan$St_indic == 1){
    index_vec  <- data.stan$id_St
  }else{
    index_vec <- c(data.stan$id_trt,data.stan$id_comp)
  }
  
  logf.expert <- rep(NA, nrow(linpred))
  
  
  for(i in 1:nrow(linpred)){
    logf.expert[i] <-  expert_like(data.stan, dist_surv = dist, param1 = linpred[i,index_vec], param2 =scale[i], param3 = q[i])
  }
  logf.hat.expert <- expert_like(data.stan, dist_surv = dist, param1 = linpred.hat[1,index_vec], param2 = scale.bar, 
                                 param3 = q.bar[index_vec])
  
  
  npars <- 3 + sum(1 - apply(data.stan$X, 2, function(x) all(x == 0)))
  list(logf = logf, logf.hat = logf.hat, npars = npars, f = NULL, 
       f.bar = NULL, s = NULL, s.bar = NULL, logf.expert = logf.expert, logf.hat.expert = logf.hat.expert)
  
}


get_stats_bayes <- function (x, mod){
  
  if(names(x[["models"]])[mod] %in% c("Gamma", "Gen. Gamma","Gompertz")){
    table = x$models[[mod]]$BUGSoutput$summary[, c("mean", 
                                                   "sd", "2.5%", "97.5%")]
  }else{
    table = rstan::summary(x$models[[mod]])$summary[, c("mean", 
                                                        "sd", "2.5%", "97.5%")]
    table = table[-grep("lp__", rownames(table)), ]
  }
  
  
  
  if (any(apply(x$misc$data.stan[[1]]$X, 2, function(x) all(x == 
                                                            0)))) {
    table = table[-grep("beta\\[2\\]", rownames(table)), ]
  }
  
  
  res = do.call(eval(parse(text=paste0("rescale_stats_bayes_", x$misc$model_name[mod]))), 
                args = list(table = table, x = x))
  
  
  return(res)
}


load_availables <- function(){
  availables = list(mle = c(genf = "gef", genf.orig = "gof", 
                            gengamma = "gga", gengamma.orig = "ggo", 
                            exp = "exp", weibull = "wei", weibullPH = "wph", 
                            lnorm = "lno", gamma = "gam", gompertz = "gom", 
                            llogis = "llo", lognormal = "lno", rps = "rps"), 
                    inla = c(exponential = "exp", weibull = "wei", 
                             weibullPH = "wph", lognormal = "lno", 
                             loglogistic = "llo", rps = "rps"), 
                    bayes = c(Exponential = "exp",Gamma = "gam", GenGamma = "gga", 
                            Gompertz = "gom", 
                            RP = "rps", WeibullAF = "wei", WeibullPH = "wph", 
                            logLogistic = "llo", logNormal = "lno"))
  return(availables)
}

#' Print a summary of the survival model(s) fitted by \code{fit.models}
#' 
#' Prints the summary table for the model(s) fitted, with the estimate of the
#' parameters - ported from ``survHE``.
#' 
#' 
#' @param x the \code{expertsurv} object (the output of the call to
#' \code{fit.models})
#' @param mod is the index of the model. Default value is 1, but the user can
#' choose which model fit to visualise, if the call to fit.models has a vector
#' argument for distr (so many models are fitted & stored in the same object)
#' @param \dots additional options, including: \code{digits} = number of
#' significant digits to be shown in the summary table (default = 6)
#' \code{original} = a flag to say whether the *original* table
#' from either \code{flexsurv} or \code{rstan/JAGS} should be printed
#' @return Printed message (no object returned) providing estimates of the survival models.
#' @author Gianluca Baio
#' @keywords Parametric survival models
#' @examples
#' require("dplyr")
#' param_expert_example1 <- list()
#' param_expert_example1[[1]] <- data.frame(dist = c("norm","t"),
#'                                          wi = c(0.5,0.5), # Ensure Weights sum to 1
#'                                        param1 = c(0.1,0.12),
#'                                       param2 = c(0.15,0.5),
#'                                        param3 = c(NA,3))
#' timepoint_expert <- 14
#' data2 <- data %>% rename(status = censored) %>% mutate(time2 = ifelse(time > 10, 10, time),
#'                                                        status2 = ifelse(time> 10, 0, status))
#' mle = example1 <- fit.models.expert(formula=Surv(time2,status2)~1,data=data2,
#'                    distr=c("wph", "gomp"),
#'                    method="mle",
#'                    pool_type = "log pool",
#'                    opinion_type = "survival",
#'                    times_expert = timepoint_expert,
#'                    param_expert = param_expert_example1)
#' print(mle)
#' @references 
#' \insertRef{Baio.2020}{expertsurv}
#' 
#' @export
print.expertsurv <-function (x, mod = 1, ...) 
{
  
  
  exArgs <- list(...)
  availables <- load_availables()
  if (!exists("digits", where = exArgs)) {
    digits = 6
  }
  else {
    digits = exArgs$digits
  }
  if (!exists("original", where = exArgs)) {
    original = FALSE
  }
  else {
    original = exArgs$original
  }
  if (exists("orig", exArgs)) {
    original = exArgs$orig
  }
  if (original == TRUE) {
    do.call(paste0("original_table_", x$method), args = list(x, 
                                                             mod, digits))
  }
  else {
    method_eval <- paste0("get_stats_", x$method)
    res = do.call(method_eval, args = list(x,mod))
    format_table(x, mod, res, digits)
  }
}


#### Adjusts SurvHE functions:

#error in a rps function 



make_sim_bayes <- function (m, t, X, nsim, newdata, dist, summary_stat, ...){
  
  if(inherits(m,"rjags")){
    iter_stan <- m[["n.iter"]]
    beta <- m$BUGSoutput$sims.matrix[, grep("beta",colnames(m$BUGSoutput$sims.matrix))]
  }else{
    
    iter_stan <- m@stan_args[[1]]$iter
    beta = rstan::extract(m)$beta
  }
  if (ncol(X) == 1) {
    beta = beta[, 1, drop = F] #PC: add drop to stop it coreceing to a vector
  }
  if (dist == "rps" & any(grepl("Intercept", colnames(X)))) {
    X <- as.matrix(tibble::as_tibble(X) %>% dplyr::select(-`(Intercept)`))
    beta = beta[, -ncol(beta)]
  }
  linpred <- beta %*% t(X)
  sim <- lapply(1:nrow(X), function(x) {
    do.call(paste0("rescale_bayes_", dist), args = list(m, 
                                                      X, linpred[, x]))
  })
  if (nsim > iter_stan) {
    stop("Please select a value for 'nsim' that is less than or equal to the value set in the call to 'fit.models'")
  }
  if (nsim == 1) {
    sim <- lapply(sim, function(x) as.matrix(tibble::as_tibble(x) %>% 
                                               dplyr::summarise_all(summary_stat), nrow = 1, ncol = ncol(x)))
  }
  if (nsim > 1 & nsim < iter_stan) {
    sim <- lapply(sim, function(x) as.matrix(tibble::as_tibble(x) %>% 
                                               dplyr::sample_n(nsim, replace = FALSE), nrow = nsim, ncol = ncol(x)))
  }
  return(sim)
}


rescale_bayes_gam <- function (m, X, linpred){
  if(inherits(m,"rjags")){
    shape <- as.numeric(m$BUGSoutput$sims.matrix[,"alpha"])
  }else{
    shape <- as.numeric(rstan::extract(m)$alpha)
  }
  
  rate <- exp(linpred)
  sim <- cbind(shape, rate)
  colnames(sim) <- c("shape", "rate")
  return(sim)
}



rescale_bayes_gom <- function (m, X, linpred){
  if(inherits(m,"rjags")){
    shape <- as.numeric(m$BUGSoutput$sims.matrix[,"alpha"])
  }else{
    shape <- as.numeric(rstan::extract(m)$alpha)
  }
  rate <- exp(linpred)
  sim <- cbind(shape, rate)
  colnames(sim) <- c("shape", "rate")
  return(sim)
}


rescale_bayes_gga<- function (m, X, linpred){
  
  if(inherits(m,"rjags")){
    Q <- as.numeric(m$BUGSoutput$sims.matrix[,"Q"])
    sigma <- as.numeric(m$BUGSoutput$sims.matrix[,"sigma"])
  }else{
    Q <- as.numeric(rstan::extract(m)$Q)
    sigma <- as.numeric(rstan::extract(m)$sigma)
  }
  mu <- linpred
  sim <- cbind(mu, sigma, Q)
  colnames(sim) <- c("mu", "sigma", "Q")
  return(sim)
}

get_stats_bayes <- function(x, mod){
  
  if(inherits(x$models[[mod]],"rjags")) {
    table =  x$models[[mod]]$BUGSoutput$summary[,c("mean", 
                                                   "sd", "2.5%", "97.5%")]
  }else{
    table = rstan::summary(x$models[[mod]])$summary[, c("mean", 
                                                        "sd", "2.5%", "97.5%")]
    table = table[-grep("lp__", rownames(table)), ]
    
  }
  
  if ("X_obs" %in% names(x$misc$data.stan[[1]])) {
    if (any(apply(x$misc$data.stan[[1]]$X_obs, 2, function(x) all(x == 
                                                                  0)))) {
      table = table[-grep("beta\\[2\\]", rownames(table)), 
      ]
    }
  }
  else {
    if (any(apply(x$misc$data.stan[[1]]$X, 2, function(x) all(x == 
                                                              0)))) {
      table = table[-grep("beta\\[2\\]", rownames(table)), 
      ]
    }
  }
  res = do.call(paste0("rescale_stats_bayes_", x$misc$model_name[mod]), 
                args = list(table = table, x = x))
  return(res)
}



rescale_stats_bayes_gam <- function (table, x){
  rate <- matrix(table[grep("rate", rownames(table)),], ncol = 4)
  rownames(rate) <- "rate"
  shape <- matrix(table[grep("alpha", rownames(table)), 
  ], ncol = 4)
  rownames(shape) <- "shape"
  effects = add_effects_bayes(table, x) #typo in this function
  res <- rbind(shape, rate, effects)
  if (is.null(dim(res))) {
    names(res) <- c("mean", "se", "L95%", 
                    "U95%")
  }
  else {
    colnames(res) <- c("mean", "se", "L95%", 
                       "U95%")
  }
  return(res)
}

get_stats_bayes <- function (x, mod){
  if (inherits(x$models[[mod]],"rjags")) {
    table = x$models[[mod]]$BUGSoutput$summary[, c("mean", 
                                                   "sd", "2.5%", "97.5%")]
  }
  else {
    table = rstan::summary(x$models[[mod]])$summary[, c("mean", 
                                                        "sd", "2.5%", "97.5%")]
    table = table[-grep("lp__", rownames(table)), ]
  }
  if ("X_obs" %in% names(x$misc$data.stan[[1]])) {
    if (any(apply(x$misc$data.stan[[1]]$X_obs, 2, function(x) all(x == 0)))) {
      #Error (mod instead of 1)
      #for RPS X matrix can be 0 for both columns
      beta_drop <- which(apply(x$misc$data.stan[[mod]]$X, 2, function(x) all(x == 0)) == TRUE)
      beta_drop <- paste0("beta\\[", beta_drop,"\\]")
      if(length(beta_drop)>1){
        beta_drop <- paste(beta_drop, collapse = "|")
      }
      table <-  table[-grep(beta_drop, rownames(table)), ]
    }
  }
  else {
    if (any(apply(x$misc$data.stan[[1]]$X, 2, function(x) all(x == 0)))) {
      beta_drop <- which(apply(x$misc$data.stan[[mod]]$X, 2, function(x) all(x == 0)) == TRUE)
      beta_drop <- paste0("beta\\[", beta_drop,"\\]")
      if(length(beta_drop)>1){
        beta_drop <- paste(beta_drop, collapse = "|")
      }
      table <-table[-grep(beta_drop, rownames(table)), ]
    }
  }
  res = do.call(paste0("rescale_stats_bayes_", x$misc$model_name[mod]), 
                args = list(table = table, x = x))
  return(res)
}



#' Graphical representation of the measures of model fitting based on Information Criteria
#'
#' Plots a summary of the model fit for all the models fitted.
#' 
#' @param ... Optional inputs. Must include an \code{expertsurv} object.
#' @param type should the DIC, WAIC, PML be plotted (AIC, BIC also allowed but only valid for frequentist approach). 
#'
#' @import dplyr
#' @import ggplot2
#' 
#'
#' @return A plot with the relevant model fitting statistics plotted in order of fit.
#' @export
#'
#'
#' @examples
#' require("dplyr")
#' param_expert_example1 <- list()
#' param_expert_example1[[1]] <- data.frame(dist = c("norm"),
#'                                          wi = c(1), # Ensure Weights sum to 1
#'                                          param1 = c(0.1),
#'                                          param2 = c(0.05),
#'                                          param3 = c(NA))
#' timepoint_expert <- 14 # Expert opinion at t = 14
#' 
#' 
#' data2 <- expertsurv::data %>% rename(status = censored) %>% 
#' mutate(time2 = ifelse(time > 10, 10, time),
#' status2 = ifelse(time> 10, 0, status))
#'example1  <- fit.models.expert(formula=Surv(time2,status2)~1,data=data2,
#'                               distr=c("wei", "gomp"),
#'                               method="mle",
#'                               pool_type = "linear pool", 
#'                               opinion_type = "survival",
#'                               times_expert = timepoint_expert, 
#'                               param_expert = param_expert_example1)
#'
#'
#' model.fit.plot(example1, type = "aic")
#' 
#' 
model.fit.plot<- function (..., type = "dic"){
  exArgs = list(...)
  if (length(names(exArgs)) == 0) {
    names(exArgs) = paste0("Object", 1:length(exArgs))
  }
  if (length(which(names(exArgs) == "")) > 0) {
    names(exArgs)[which(names(exArgs) == "")] = paste0("Object", 
                                                       1:length(which(names(exArgs) == "")))
  }
  w <- which(unlist(lapply(1:length(exArgs), function(i) class(exArgs[[i]]))) == 
               "expertsurv")
  if (length(w) == 0) {
    stop("Please give at least one 'expertsurv' object, generated by a call to 'fit.models(...)")
  }
  else {
    survHE_objs = lapply(1:length(w), function(i) exArgs[[w[i]]])
  }
  names(survHE_objs) = names(exArgs)[w]
  if (!exists("mods", exArgs)) {
    mods = 1:sum(unlist(lapply(survHE_objs, function(x) length(x$models))))
  }
  else {
    mods = exArgs$mods
  }
  if (type %in% c("aic", "AIC", "a", "A")) {
    type = "AIC"
  }
  if (type %in% c("bic", "BIC", "b", "B")) {
    type = "BIC"
  }
  if (type %in% c("dic", "DIC", "d", "D")) {
    type = "DIC"
  }
  if (type %in% c("dic2", "DIC2")) {
    type = "DIC2"
  }
  if (type %in% c("waic", "WAIC", "w", "W")) {
    type = "WAIC"
  }
  if (type %in% c("pml", "PML", "p", "P")) {
    type = "PML"
  }
  type2 <- tolower(type) 
  
  
  
  toplot = lapply(1:length(survHE_objs), function(x) survHE_objs[[x]]$model.fitting %>% 
                    bind_rows %>% mutate(object_name = as.factor(names(survHE_objs)[x]), 
                                         model_name = names(survHE_objs[[x]]$models))) %>% bind_rows %>% 
    mutate(lab = paste0(model_name, ":", object_name)) %>% 
    dplyr::select(object_name, model_name, lab, everything()) %>% 
    slice(mods) %>% arrange(desc(!!as.symbol(type2)))
  
  if (exists("xlim", exArgs)) {
    yl = exArgs$xlim
  }else{
    type_vals  = toplot %>% pull(type2)
    yl = c(min(type_vals)*.9, max(type_vals)*1.1)
    #range(pretty(range(toplot %>% pull(type2))))
  }
  
  toplot$model_name <- factor(toplot$model_name, levels =  toplot$model_name) 
  
  mfp = ggplot(data = toplot, aes(x = model_name, y = get(type2), fill = object_name)) + 
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(x = model_name, y = get(type2), label = get(type2) %>% round(digits = 1.5)), hjust = 1.05, 
              color = "white", size = 5.5, position = position_dodge(0.9)) + 
    coord_flip(ylim = yl)
  
  mfp = mfp + theme_bw() + theme(axis.text.x = element_text(color = "black", 
                                                            size = 12, angle = 0, hjust = 0.5, vjust = 0.5), 
                                 axis.text.y = element_text(color = "black", size = 12, 
                                                            angle = 0, hjust = 0.5, vjust = 0.5), axis.title.x = element_text(color = "black", 
                                                                                                                              size = 14, angle = 0, hjust = 0.5, vjust = 0.5), 
                                 axis.title.y = element_text(color = "black", size = 14, 
                                                             angle = 90, hjust = 0.5, vjust = 0.5)) + theme(axis.line = element_line(colour = "black"), 
                                                                                                            panel.background = element_blank(), panel.border = element_blank(), 
                                                                                                            plot.title = element_text(size = 18, face = "bold")) + 
    labs(y = toupper(type), x = "", title = paste0("Model comparison based on ", 
                                                   toupper(type)), fill = "survHE object") + scale_fill_brewer(palette = "Paired") + 
    theme(legend.position = "bottom")
  
  if (exists("col", exArgs)) {
    mfp = mfp + scale_fill_manual(values = exArgs$col)
  }
  if (exists("colour", exArgs)) {
    mfp = mfp + scale_fill_manual(values = exArgs$colour)
  }
  if (exists("color", exArgs)) {
    mfp = mfp + scale_fill_manual(values = exArgs$color)
  }
  if (exists("name_legend", exArgs)) {
    mfp = mfp + labs(fill = exArgs$name_legend)
  }
  
  mfp+ theme(legend.position = "none")
} 
