
make_data_expert <- function(param_expert, times_expert = NULL){
  n.experts <- c()
  
  for(i in 1:length(param_expert)){
    n.experts <- c(n.experts, nrow(param_expert[[i]]))
  }
  
  data_dist_ind <- num_param <- matrix(-999.2,nrow = max(n.experts), ncol =  length(times_expert))
  expert.array <- array(-999.2,dim = c(max(n.experts),5,length(times_expert)))
  
  for(i in 1:length(times_expert)){
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
  
  if(!is.null(times_expert)){
    n_time_expert <- length(times_expert)
    time_expert <- as.array(times_expert)
  }else{
    n_time_expert <- 1
    time_expert <- numeric(0) #This produces an array of size 0
    #https://dev.to/martinmodrak/optional-parametersdata-in-stan-4o33
    
    if (distr3 %in% c("gam", "gga", "gom")){
      time_expert <- 1 # Has to be defined for JAGS
      
    }
    
    
  }
  
  expert.array
  
}

get_k_norm <- function(opinion_list, St_indic = 1){ # Only required if log-pooling
  k_norm <- rep(NA, length(opinion_list))
  
  if(St_indic ==1){
    a <- 0
    b <- 1
  }else{
    a <- -Inf
    b <- +Inf
  }
  
  for(i in 1:length(opinion_list)){
    
    opinion_list[[i]]$dist <- stringr::str_replace_all( opinion_list[[i]]$dist, "normal", "norm") 
    opinion_list[[i]]$dist <- stringr::str_replace_all( opinion_list[[i]]$dist, "lognorm", "lnorm") 
    
    pool.df <- opinion_list[[i]]
    
    if(St_indic ==1){
      min_quant <- 0
      max_quant <- 1
    }else{
      quant.vec <- t(apply(pool.df, 1, function(x){get_quant_val(
        dist = x["dist"],
        param1 = x["param1"],
        param2 = x["param2"],
        param3 = x["param3"],
        probs = c(0.001,0.025,0.5,0.975,0.999))}))
      # central.cauchy <- mean(quant.vec[,3])#mean
      # sd.cauchy <- max(apply(quant.vec,1, function(x){(x[4]-x[2])/4})) #sd
      min_quant <- min(quant.vec)
      max_quant <- max(quant.vec)
    }
    x.eval <- seq(min_quant, max_quant, length.out = 100)
    #Find the range of x.eval to integrate over
    # x.eval <- seq(0,1, by = 0.001)
    dens_eval<- eval_dens_pool(x.eval,opinion_list[[i]],  pool_type = "log pool", St_indic =St_indic)
    k_norm[i] <- integrate.xy(x = x.eval,fx = dens_eval)
   }
  k_norm
}

#' Helper function to run the survival models using MLE and flexsurv
#' 
#' @param x a string containing the name of the model
#' to be fitted
#' #' @param exArgs a list of extra arguments passed from the main 'fit.models' 
#' function
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Maximum likelihood estimation
#' @import tibble
#' @import dplyr
#' @import stats 
#' @import flexsurv
#' @noRd 

runMLE <- function (x, exArgs){
  formula <- exArgs$formula
  data = exArgs$data
  availables <- load_availables()
  d3 <- manipulate_distributions(x)$distr3
  x <- manipulate_distributions(x)$distr
  expert_opinion_flex <- list()
  
 
  if (exArgs$opinion_type != "survival") {
    times <- 99999
    expert_opinion_flex$St_indic <- 0
  }
  
  if (exArgs$opinion_type == "survival") {
    expert_opinion_flex$St_indic <- 1
    times <- exArgs$times_expert
    if(is.null(exArgs$id_St)){
      expert_opinion_flex$id_St <- 1
    }else{
      expert_opinion_flex$id_St <- exArgs$id_St
    }
  }
  if(exArgs$opinion_type == "mean"){
  if(is.null(exArgs$id_trt|exArgs$id_comp)){
    message("You need to supply the location within the dataframe row number of a treatment and a comparator arm to arguments id_trt and id_comp")
    stop()
  }
  expert_opinion_flex$id_trt <- exArgs$id_trt
  expert_opinion_flex$id_comp <- exArgs$id_comp
  }
  
  expert_opinion_flex$param_expert <- make_data_expert(exArgs$param_expert, 
                                                       times)
  expert_opinion_flex$times <- times
  if (exArgs$pool_type == "linear pool") {
    expert_opinion_flex$pool <- 1
  }
  else {
    expert_opinion_flex$pool <- 0
  }
  if (expert_opinion_flex$pool == 0) {
    expert_opinion_flex$k_norm <- get_k_norm(exArgs$param_expert)
  }
  else {
    expert_opinion_flex$k_norm <- NULL
  }
  if (exists("method_mle", where = exArgs)) {
    method_mle <- exArgs$method_mle
  }
  else {
    method_mle <- "BFGS"
  }
  tic <- proc.time()
  if (x == "survspline") {
    if (exists("bhazard", where = exArgs)) {
      bhazard <- exArgs$bhazard
    }
    else {
      bhazard <- NULL
    }
    if (exists("weights", where = exArgs)) {
      weights <- exArgs$weights
    }
    else {
      weights <- NULL
    }
    if (exists("subset", where = exArgs)) {
      subset <- exArgs$subset
    }
    else {
      subset <- NULL
    }
    if (exists("knots", where = exArgs)) {
      knots <- exArgs$knots
    }
    else {
      knots <- NULL
    }
    if (exists("k", where = exArgs)) {
      k <- exArgs$k
    }
    else {
      k <- 0
    }
    if (exists("bknots", where = exArgs)) {
      bknots <- exArgs$bknots
    }
    else {
      bknots <- NULL
    }
    if (exists("scale", where = exArgs)) {
      scale <- exArgs$scale
    }
    else {
      scale <- "hazard"
    }
    if (exists("timescale", where = exArgs)) {
      timescale <- exArgs$scale
    }
    else {
      timescale <- "log"
    }
    suppressWarnings({
      #browser()
      model_mle <- flexsurvspline(formula = formula, data = data, 
                                  k = k, knots = knots, bknots = bknots, scale = scale, 
                                  timescale = timescale, expert_opinion = NULL, 
                                  method = method_mle)
      model <- flexsurvspline(formula = formula, data = data, 
                              k = k, knots = knots, bknots = bknots, scale = scale, 
                              timescale = timescale, expert_opinion = expert_opinion_flex, 
                              method = method_mle, inits = model_mle$res[,1])
    })
  }
  else {
    suppressWarnings({
      model_mle <- flexsurvreg(formula = formula, data = data, 
                               dist = x, expert_opinion = NULL, method = method_mle)
      model <- flexsurvreg(formula = formula, data = data, 
                           dist = x, expert_opinion = expert_opinion_flex, 
                           method = method_mle, inits = model_mle$res[,1])
    })
  }
  toc <- proc.time() - tic
  model_name <- d3
  list(model = model, aic = model$AIC, bic = -2 * model$loglik + 
         model$npars * log(model$N), dic = NULL, time2run = toc[3], 
       model_name = model_name)
}
