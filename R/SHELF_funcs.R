# FUNCTION ----
`%!in%` = Negate(`%in%`)

#NAMESPACE HACK FOR CRAN; won't let me use SHELF 3 times :
logt.error <-utils::getFromNamespace("logt.error", "SHELF")
gamma.error<-utils::getFromNamespace("gamma.error", "SHELF")
lognormal.error<-utils::getFromNamespace("lognormal.error", "SHELF")
logt.error<-utils::getFromNamespace("logt.error", "SHELF")
makeGroupPlot<-utils::getFromNamespace("makeGroupPlot", "SHELF")
makeLinearPoolPlot<-utils::getFromNamespace("makeLinearPoolPlot", "SHELF")
makeSingleExpertPlot<-utils::getFromNamespace("makeSingleExpertPlot", "SHELF")
expertdensity<-utils::getFromNamespace("expertdensity", "SHELF")

normal.error_mod <- function (parameters, values, probabilities, weights, mode,trunc =FALSE){
  
  if(trunc){ #Survival Trunc
    Fx <- pnorm(values, parameters[1], exp(parameters[2]))
    Fa <- pnorm(0, parameters[1], exp(parameters[2]))
    Fb <- pnorm(1, parameters[1], exp(parameters[2]))
    F_final <- (Fx - Fa)/(Fb-Fa)
  }else{
    F_final <- pnorm(values, parameters[1], exp(parameters[2]))
  }
  
  
  res1 <- sum(weights * (F_final -  probabilities)^2)
  if (!is.null(mode)) {
    res1 <- res1 + (parameters[1] - mode)^2
  }
  return(res1)
}

t_error_mod <- function (parameters, values, probabilities, weights, degreesfreedom, 
                         mode,trunc=FALSE){ 
  
  if(trunc){ #Survival Trunc
    Fx <- stats::pt((values - parameters[1])/exp(parameters[2]), 
                    degreesfreedom)
    Fa <- stats::pt((0 - parameters[1])/exp(parameters[2]), 
                    degreesfreedom)
    Fb <- stats::pt((1 - parameters[1])/exp(parameters[2]), 
                    degreesfreedom)
    F_final <- (Fx - Fa)/(Fb-Fa)
  }else{
    F_final <- stats::pt((values - parameters[1])/exp(parameters[2]), 
                         degreesfreedom)
  }
  
  
  res1 <- sum(weights * (F_final - probabilities)^2)
  if (!is.null(mode)) {
    res1 <- res1 + (parameters[1] - mode)^2
  }
  return(res1)
}

gamma.error_mod <- function (parameters, values, probabilities, weights, mode,trunc=FALSE){
  
  if(trunc){ #Survival Trunc
    Fx <- stats::pgamma(values, exp(parameters[1]),exp(parameters[2]))
    Fa <- stats::pgamma(0, exp(parameters[1]),exp(parameters[2]))
    Fb <- stats::pgamma(1, exp(parameters[1]),exp(parameters[2]))
    F_final <- (Fx - Fa)/(Fb-Fa)
  }else{
    F_final <- stats::pgamma(values, exp(parameters[1]),exp(parameters[2]))
  }
  
  res1 <- sum(weights * (F_final - probabilities)^2)
  if (!is.null(mode)) {
    res1 <- res1 + ((exp(parameters[1]) - 1)/exp(parameters[2]) - 
                      mode)^2
  }
  return(res1)
}

lognormal.error_mod <-function (parameters, values, probabilities, weights, mode,trunc =FALSE){
  
  if(trunc){ #Survival Trunc
    Fx <- stats::plnorm(values, parameters[1],exp(parameters[2]))
    Fa <- stats::plnorm(0, parameters[1],exp(parameters[2]))
    Fb <- stats::plnorm(1, parameters[1],exp(parameters[2]))
    F_final <- (Fx - Fa)/(Fb-Fa)
  }else{
    F_final <- stats::plnorm(values, parameters[1],exp(parameters[2]))
  }
  res1 <- sum(weights * ( F_final- probabilities)^2)
  
  if (!is.null(mode)) {
    res1 <- res1 + (exp(parameters[1] - exp(parameters[2])^2) - 
                      mode)^2
  }
  return(res1)
}
beta.error_mod <- function (parameters, values, probabilities, weights, mode ){
  
  
  res1 <- sum(weights * (stats::pbeta(values, exp(parameters[1]), 
                                      exp(parameters[2])) - probabilities)^2)
  if (!is.null(mode)) {
    res1 <- res1 + ((exp(parameters[1]) - 1)/(exp(parameters[1]) + 
                                                exp(parameters[2]) - 2) - mode)^2
  }
  return(res1)
}

dt.scaled <- function (x, df, mean = 0, sd = 1, ncp, log = FALSE){
  if (!log) 
    stats::dt((x - mean)/sd, df, ncp = ncp, log = FALSE)/sd
  else stats::dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - 
    log(sd)
}


expert_log_dens <- function(x, df, pool_type, k_norm = NULL, St_indic){
  if(is.null(dim(df))){ #corced to vector
    df <-matrix(df, nrow = 1, ncol = length(df))
  }
  if(St_indic ==1){
    a <- 0
    b <- 1
  }else{
    a <- -Inf
    b <- +Inf
  }
    
  like <- rep(NA,nrow(df)) 
  for(i in 1:nrow(df)){
    
    if(df[i,1] == 1){ # 1 equal normal
      like[i] <- stats::dnorm(x, df[i,3], df[i,4], log = F) 
      k_trunc <- stats::pnorm(b,df[i,3], df[i,4])-stats::pnorm(a,df[i,3], df[i,4])
    }
    
    
    if(df[i,1] == 2){ # 2 equal t
      like[i] <- dt.scaled(x, df[i,5], df[i,3], df[i,4], log = F)
      k_trunc <- pt.scaled(b, df[i,5], df[i,3], df[i,4])-pt.scaled(a, df[i,5], df[i,3], df[i,4])
      
    }
    
    if(df[i,1] == 3){ # 3 equal gamma
      like[i] <- stats::dgamma(x, df[i,3], df[i,4],  log = F)   
      k_trunc <- stats::pgamma(b,  df[i,3], df[i,4])-stats::pgamma(a,  df[i,3], df[i,4])
      
    }
    
    if(df[i,1] == 4){ # 4 equal lnorm
        like[i] <- stats::dlnorm(x,  df[i,3], df[i,4],log = F)  
       k_trunc <- stats::plnorm(b,  df[i,3], df[i,4])-stats::plnorm(a,  df[i,3], df[i,4])
      
    }
    
    
    if(df[i,1] == 5){# 5 = beta
        like[i] <- stats::dbeta(x,  df[i,3], df[i,4],   log = F)
        k_trunc <- stats::pbeta(b,  df[i,3], df[i,4])-stats::pbeta(a,  df[i,3], df[i,4])
      
    }
    
    like[i] <- like[i]/k_trunc
    if(pool_type ==1){
      like[i] <- like[i]*df[i,2] 
    }else{
      like[i] <- like[i]^df[i,2] 
    }
    
    
    
  }  
  if(pool_type == 1){
    return(log(sum(like)))
  }else{
    return(log(prod(like)/k_norm))
  }
  
  
  
  
}





fitdist_mod <- function (vals, probs, lower = -Inf, upper = Inf, weights = 1, 
                         tdf = 3, expertnames = NULL, excludelog.mirror = TRUE, mode = NULL, trunc = FALSE){
  logt.error <- utils::getFromNamespace("logt.error", "SHELF")
  gamma.error <- utils::getFromNamespace("gamma.error", "SHELF")
  lognormal.error <- utils::getFromNamespace("lognormal.error", 
                                             "SHELF")
  logt.error <- utils::getFromNamespace("logt.error", "SHELF")
  makeGroupPlot <- utils::getFromNamespace("makeGroupPlot", 
                                           "SHELF")
  makeLinearPoolPlot <- utils::getFromNamespace("makeLinearPoolPlot", 
                                                "SHELF")
  makeSingleExpertPlot <- utils::getFromNamespace("makeSingleExpertPlot", 
                                                  "SHELF")
  expertdensity <- utils::getFromNamespace("expertdensity", 
                                           "SHELF")
  if (is.matrix(vals) == F) {
    vals <- matrix(vals, nrow = length(vals), ncol = 1)
  }
  if (is.matrix(probs) == F) {
    probs <- matrix(probs, nrow = nrow(vals), ncol = ncol(vals))
  }
  if (is.matrix(weights) == F) {
    weights <- matrix(weights, nrow = nrow(vals), ncol = ncol(vals))
  }
  if (length(lower) == 1) {
    lower <- rep(lower, ncol(vals))
  }
  if (length(upper) == 1) {
    upper <- rep(upper, ncol(vals))
  }
  if (length(tdf) == 1) {
    tdf <- rep(tdf, ncol(vals))
  }
  n.experts <- ncol(vals)
  normal.parameters <- matrix(NA, n.experts, 2)
  t.parameters <- matrix(NA, n.experts, 3)
  mirrorgamma.parameters <- gamma.parameters <- matrix(NA, 
                                                       n.experts, 2)
  mirrorlognormal.parameters <- lognormal.parameters <- matrix(NA, 
                                                               n.experts, 2)
  mirrorlogt.parameters <- logt.parameters <- matrix(NA, n.experts, 
                                                     3)
  beta.parameters <- matrix(NA, n.experts, 2)
  ssq <- matrix(NA, n.experts, 9)
  colnames(ssq) <- c("normal", "t", "gamma", "lognormal", "logt", 
                     "beta", "mirrorgamma", "mirrorlognormal", "mirrorlogt")
  if (n.experts > 1 & n.experts < 27 & is.null(expertnames)) {
    expertnames <- paste("expert.", LETTERS[1:n.experts], 
                         sep = "")
  }
  if (n.experts > 27 & is.null(expertnames)) {
    expertnames <- paste("expert.", 1:n.experts, sep = "")
  }
 for (i in 1:n.experts) {
    # if (min(probs[, i]) > 0.4) {
    #  stop("smallest elicited probability must be less than 0.4")
    # }
    if (min(probs[, i]) < 0 | max(probs[, i]) > 1) {
      stop("probabilities must be between 0 and 1")
    }
    #  if (max(probs[, i]) < 0.6) {
    #    stop("largest elicited probability must be greater than 0.6")
    #  }
    if (min(vals[, i]) < lower[i]) {
      stop("elicited parameter values cannot be smaller than lower parameter limit")
    }
    if (max(vals[, i]) > upper[i]) {
      stop("elicited parameter values cannot be greater than upper parameter limit")
    }
    if (tdf[i] <= 0) {
      stop("Student-t degrees of freedom must be greater than 0")
    }
    if (min(probs[-1, i] - probs[-nrow(probs), i]) < 0) {
      stop("probabilities must be specified in ascending order")
    }
    if (min(vals[-1, i] - vals[-nrow(vals), i]) <= 0) {
      stop("parameter values must be specified in ascending order")
    }
    inc <- (probs[, i] > 0) & (probs[, i] < 1)
    minprob <- min(probs[inc, i])
    maxprob <- max(probs[inc, i])
    minvals <- min(vals[inc, i])
    maxvals <- max(vals[inc, i])
    
    q.fit <- stats::approx(x = probs[inc, i], y = vals[inc, 
                                                       i], xout = c(0.4, 0.5, 0.6))$y
    l <- q.fit[1]
    u <- q.fit[3]
    minq <- stats::qnorm(minprob)
    maxq <- stats::qnorm(maxprob)
    
    m <- (minvals * maxq - maxvals * minq)/(maxq - minq)
    v <- ((maxvals - minvals)/(maxq - minq))^2
    #browser()
    normal.fit <- stats::optim(c(m, 0.5 * log(v)), normal.error_mod, 
                               values = vals[inc, i], probabilities = probs[inc, 
                                                                            i], weights = weights[inc, i], mode = mode[i],trunc = trunc)
    normal.parameters[i, ] <- c(normal.fit$par[1], exp(normal.fit$par[2]))
    ssq[i, "normal"] <- normal.fit$value
    
    lprob <- 0.000001
    if(is.infinite(lower[i])){
      lower[i] <- stats::qnorm(lprob, normal.parameters[i,1],normal.parameters[i,2])
      upper[i] <- stats::qnorm(1-lprob, normal.parameters[i,1],normal.parameters[i,2])
    }
    
    
    t.fit <- stats::optim(c(m, 0.5 * log(v)), t_error_mod, 
                          values = vals[inc, i], probabilities = probs[inc, 
                                                                       i], weights = weights[inc, i], degreesfreedom = tdf[i], 
                          mode = mode[i],trunc = trunc)
    t.parameters[i, 1:2] <- c(t.fit$par[1], exp(t.fit$par[2]))
    t.parameters[i, 3] <- tdf[i]
    ssq[i, "t"] <- t.fit$value
    if (lower[i] == 0) { #Can't use the distribtuions as they are shifted distributions if lower not equal to 0
      vals.scaled1 <- vals[inc, i] - lower[i]
      m.scaled1 <- m - lower[i]
     # browser()
      gamma.fit <- stats::optim(c(log(m.scaled1^2/v), log(m.scaled1/v)), 
                                gamma.error_mod, values = vals.scaled1, probabilities = probs[inc, 
                                                                                              i], weights = weights[inc, i], mode = mode[i],trunc = trunc)
      gamma.parameters[i, ] <- exp(gamma.fit$par)
      ssq[i, "gamma"] <- gamma.fit$value
      std <- ((log(u - lower[i]) - log(l - lower[i]))/1.35)
      mlog <- (log(minvals - lower[i]) * maxq - log(maxvals - 
                                                      lower[i]) * minq)/(maxq - minq)
      lognormal.fit <- stats::optim(c(mlog, log(std)), 
                                    lognormal.error_mod, values = vals.scaled1, probabilities = probs[inc, 
                                                                                                      i], weights = weights[inc, i], mode = mode[i],trunc = trunc)
      lognormal.parameters[i, 1:2] <- c(lognormal.fit$par[1], 
                                        exp(lognormal.fit$par[2]))
      ssq[i, "lognormal"] <- lognormal.fit$value
      logt.fit <- stats::optim(c(log(m.scaled1), log(std)), 
                               logt.error, values = vals.scaled1, probabilities = probs[inc, 
                                                                                        i], weights = weights[inc, i], degreesfreedom = tdf[i])
      logt.parameters[i, 1:2] <- c(logt.fit$par[1], exp(logt.fit$par[2]))
      logt.parameters[i, 3] <- tdf[i]
      ssq[i, "logt"] <- Inf#logt.fit$value
    }
    if ((lower[i] ==0) & (upper[i] < Inf)) {#Can't use the distribtuions as they are shifted distributions if lower not equal to 0
      vals.scaled2 <- (vals[inc, i] - lower[i])/(upper[i] - 
                                                   lower[i])
      m.scaled2 <- (m - lower[i])/(upper[i] - lower[i])
      v.scaled2 <- v/(upper[i] - lower[i])^2
      alp <- abs(m.scaled2^3/v.scaled2 * (1/m.scaled2 - 
                                            1) - m.scaled2)
      bet <- abs(alp/m.scaled2 - alp)
      if (identical(probs[inc, i], (vals[inc, i] - lower[i])/(upper[i] - 
                                                              lower[i]))) {
        alp <- bet <- 1
      }
      beta.fit <- stats::optim(c(log(alp), log(bet)), beta.error_mod, 
                               values = vals.scaled2, probabilities = probs[inc, 
                                                                            i], weights = weights[inc, i], mode = mode[i], lower = lower[i], upper = upper[i])
      beta.parameters[i, ] <- exp(beta.fit$par)
      

      ssq[i, "beta"] <- beta.fit$value

    }
    if (upper[i] < Inf) {
      valsMirrored <- upper[i] - vals[inc, i]
      probsMirrored <- 1 - probs[inc, i]
      mMirrored <- upper[i] - m
      mirrorgamma.fit <- stats::optim(c(log(mMirrored^2/v), 
                                        log(mMirrored/v)), gamma.error, values = valsMirrored, 
                                      probabilities = probsMirrored, weights = weights[inc, 
                                                                                       i])
      mirrorgamma.parameters[i, ] <- exp(mirrorgamma.fit$par)
      ssq[i, "mirrorgamma"] <- Inf #mirrorgamma.fit$value
      mlogMirror <- (log(upper[i] - maxvals) * (1 - minq) - 
                       log(upper[i] - minvals) * (1 - maxq))/(maxq - 
                                                                minq)
      stdMirror <- ((log(upper[i] - l) - log(upper[i] - 
                                               u))/1.35)
      mirrorlognormal.fit <- optim(c(mlogMirror, log(stdMirror)), 
                                   lognormal.error, values = valsMirrored, probabilities = probsMirrored, 
                                   weights = weights[inc, i])
      mirrorlognormal.parameters[i, 1:2] <- c(mirrorlognormal.fit$par[1], 
                                              exp(mirrorlognormal.fit$par[2]))
      ssq[i, "mirrorlognormal"] <- mirrorlognormal.fit$value
      mirrorlogt.fit <- stats::optim(c(log(mMirrored), 
                                       log(stdMirror)), logt.error, values = valsMirrored, 
                                     probabilities = probsMirrored, weights = weights[inc, 
                                                                                      i], degreesfreedom = tdf[i])
      mirrorlogt.parameters[i, 1:2] <- c(mirrorlogt.fit$par[1], 
                                         exp(mirrorlogt.fit$par[2]))
      mirrorlogt.parameters[i, 3] <- tdf[i]
      ssq[i, "mirrorlogt"] <- Inf#mirrorlogt.fit$value
    }
 }
  
  limits <- data.frame(lower = lower, upper = upper)
  row.names(limits) <- expertnames
  
  dfn <- data.frame(normal.parameters)
  names(dfn) <- c("mean", "sd")
  row.names(dfn) <- expertnames
  dft <- data.frame(t.parameters)
  names(dft) <- c("location", "scale", "df")
  row.names(dft) <- expertnames
  dfg <- data.frame(gamma.parameters)
  names(dfg) <- c("shape", "rate")
  row.names(dfg) <- expertnames
  dfmirrorg <- data.frame(mirrorgamma.parameters)
  names(dfmirrorg) <- c("shape", "rate")
  row.names(dfmirrorg) <- expertnames
  dfln <- data.frame(lognormal.parameters)
  names(dfln) <- c("mean.log.X", "sd.log.X")
  row.names(dfln) <- expertnames
  dfmirrorln <- data.frame(mirrorlognormal.parameters)
  names(dfmirrorln) <- c("mean.log.X", "sd.log.X")
  row.names(dfmirrorln) <- expertnames
  dflt <- data.frame(logt.parameters)
  names(dflt) <- c("location.log.X", "scale.log.X", "df.log.X")
  row.names(dflt) <- expertnames
  dfmirrorlt <- data.frame(mirrorlogt.parameters)
  names(dfmirrorlt) <- c("location.log.X", "scale.log.X", "df.log.X")
  row.names(dfmirrorlt) <- expertnames
  dfb <- data.frame(beta.parameters)
  names(dfb) <- c("shape1", "shape2")
  row.names(dfb) <- expertnames
  ssq <- data.frame(ssq)
  row.names(ssq) <- expertnames
  if (excludelog.mirror) {
    reducedssq <- ssq[, c("normal", "t", "gamma", "lognormal", 
                          "beta")]
    index <- apply(reducedssq, 1, which.min)
    best.fitting <- data.frame(best.fit = names(reducedssq)[index])
  }
  else {
    index <- apply(ssq, 1, which.min)
    best.fitting <- data.frame(best.fit = names(ssq)[index])
  }
  row.names(best.fitting) <- expertnames
  vals <- data.frame(vals)
  names(vals) <- expertnames
  probs <- data.frame(probs)
  names(probs) <- expertnames
  fit <- list(Normal = dfn, Student.t = dft, Gamma = dfg, Log.normal = dfln, 
              Log.Student.t = dflt, Beta = dfb, mirrorgamma = dfmirrorg, 
              mirrorlognormal = dfmirrorln, mirrorlogt = dfmirrorlt, 
              ssq = ssq, best.fitting = best.fitting, vals = t(vals), 
              probs = t(probs), limits = limits)
  class(fit) <- "elicitation"
  fit
}



plotfit <- function (fit, d = "best", xl = -Inf, xu = Inf, ql = NA, qu = NA, 
          lp = FALSE, ex = NA, sf = 3, ind = TRUE, lpw = 1, fs = 12, 
          lwd = 1, xlab = "x", ylab = expression(f[X](x)), legend_full = TRUE, 
          percentages = FALSE, returnPlot = FALSE){
#NAMESPACE HACK FOR CRAN; three times :
logt.error <-utils::getFromNamespace("logt.error", "SHELF")
gamma.error<-utils::getFromNamespace("gamma.error", "SHELF")
lognormal.error<-utils::getFromNamespace("lognormal.error", "SHELF")
logt.error<-utils::getFromNamespace("logt.error", "SHELF")
makeGroupPlot<-utils::getFromNamespace("makeGroupPlot", "SHELF")
makeLinearPoolPlot<-utils::getFromNamespace("makeLinearPoolPlot", "SHELF")
makeSingleExpertPlot<-utils::getFromNamespace("makeSingleExpertPlot", "SHELF")
expertdensity<-utils::getFromNamespace("expertdensity", "SHELF")
	  
		  
  if (d == "beta" & (min(fit$limits) == -Inf | max(fit$limits) == 
                     Inf)) {
    stop("Parameter limits must be finite to fit a beta distribution")
  }
  if (d == "gamma" & min(fit$limits) == -Inf) {
    stop("Lower parameter limit must be finite to fit a (shifted) gamma distribution")
  }
  if (d == "lognormal" & min(fit$limits) == -Inf) {
    stop("Lower parameter limit must be finite to fit a (shifted) log normal distribution")
  }
  if (d == "logt" & min(fit$limits) == -Inf) {
    stop("Lower parameter limit must be finite to fit a (shifted) log t distribution")
  }
  if (is.na(ql) == F & (ql < 0 | ql > 1)) {
    stop("Lower feedback quantile must be between 0 and 1")
  }
  if (is.na(qu) == F & (qu < 0 | qu > 1)) {
    stop("Upper feedback quantile must be between 0 and 1")
  }
  ggplot2::theme_set(ggplot2::theme_grey(base_size = fs))
  ggplot2::theme_update(plot.title = ggplot2::element_text(hjust = 0.5))
  if (nrow(fit$vals) > 1 & is.na(ex) == T & lp == F) {
    if (xl == -Inf & min(fit$limits[, 1]) > -Inf) {
      xl <- min(fit$limits[, 1])
    }
    if (xu == Inf & max(fit$limits[, 2]) < Inf) {
      xu <- max(fit$limits[, 2])
    }
    p1 <- suppressWarnings(makeGroupPlot(fit, xl, xu, d, 
                                         lwd, xlab, ylab, expertnames = rownames(fit$Normal)))

    if (returnPlot) {
      return(p1)
    }
  }
  if (nrow(fit$vals) > 1 & lp == T) {
    if (xl == -Inf & min(fit$limits[, 1]) > -Inf) {
      xl <- min(fit$limits[, 1])
    }
    if (xl == -Inf & min(fit$limits[, 1]) == -Inf) {
      f1 <- SHELF::feedback(fit, quantiles = 0.01, dist = d)
      xl <- min(f1$expert.quantiles)
    }
    if (xu == Inf & max(fit$limits[, 2]) < Inf) {
      xu <- max(fit$limits[, 2])
    }
    if (xu == Inf & max(fit$limits[, 2]) == Inf) {
      f2 <- SHELF::feedback(fit, quantiles = 0.99, dist = d)
      xu <- max(f2$expert.quantiles)
    }
    p1 <- makeLinearPoolPlot(fit, xl, xu, d, lpw, lwd, xlab, 
                             ylab, legend_full, expertnames = rownames(fit$Normal))
   
    if (returnPlot) {
      return(p1)
    }
  }
  if (nrow(fit$vals) > 1 & is.na(ex) == F) {
    if (xl == -Inf & fit$limits[ex, 1] > -Inf) {
      xl <- fit$limits[ex, 1]
    }
    if (xu == Inf & fit$limits[ex, 2] < Inf) {
      xu <- fit$limits[ex, 2]
    }
    p1 <- suppressWarnings(makeSingleExpertPlot(fit, d, 
                                                xl, xu, ql, qu, sf, ex = ex, lwd, xlab, ylab, percentages))
   
    if (returnPlot) {
      return(p1)
    }
  }
  if (nrow(fit$vals) == 1) {
    p1 <- suppressWarnings(makeSingleExpertPlot(fit, d, 
                                                xl, xu, ql, qu, sf, ex = 1, lwd, xlab, ylab, percentages))
   
    if (returnPlot) {
      return(p1)
    }
  }
}



dt.scaled <- function (x, df, mean = 0, sd = 1, ncp, log = FALSE){
  if (!log) 
    stats::dt((x - mean)/sd, df, ncp = ncp, log = FALSE)/sd
  else stats::dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - 
    log(sd)
}

qt.scaled <- function (p, df, mean = 0, sd = 1, ncp, lower.tail = TRUE, log.p = FALSE){
  mean + sd * stats::qt(p, df, ncp = ncp, log.p = log.p)
}


rt.scaled <- function (n, df, mean = 0, sd = 1, ncp) {
  mean + sd * stats::rt(n, df, ncp = ncp)
}

eval_dens_pool <- function(x.eval, pool.df, pool_type, St_indic){

  #Ensure weights sum to 1
    pool.df$wi <- pool.df$wi/sum(pool.df$wi)

  dens.vec <- apply(pool.df, 1, function(x){get_density(
    dist = x["dist"],
    param1 = x["param1"],
    param2 = x["param2"],
    param3 = x["param3"],
    x = x.eval, St_indic =St_indic)})

  if(pool_type == "log pool"){
    
    if(is.matrix(dens.vec)){
      return(apply(dens.vec, 1, function(x){prod(x^pool.df$wi)}))
    }else{
      #scaled by an arbitary constant which we don't need to know for candidate density evaluation
      return(prod(dens.vec^pool.df$wi))
    }
    
  }else{
    if(is.matrix(dens.vec)){
      return(apply(dens.vec, 1, function(x){sum(x*pool.df$wi)}))
      
    }else{
      
      return(sum(dens.vec*pool.df$wi))
    }
    
  }
  
}


ssq_mix <- function(object, values, probs){
  df_ssq <- data.frame(pi = object$pi, mu = object$mu, sd = object$sd)
  
  #Evaluate the pnorm individually
  pnorm_eval  <- apply(df_ssq,1, FUN = function(x){stats::pnorm(values,x["mu"],
                                                         x["sd"])})
  pnorm_eval_weighted <- t(pnorm_eval)*df_ssq$pi
  
  #Sum the pnorm then subtract
  return(sum((colSums(pnorm_eval_weighted) - probs)^2))
  
}


expert_dens <- function(expert_df, probs =  seq(0.01, 0.98, by = 0.002)){
  
if(length(unique(expert_df$expert)) !=1){ #Only one expert, Don't need to anything
  
  
  if(is.null(expert_df$weights) && is.null(expert_df$wi)){
    warning("No weights given.. assuming equally weighted expert opinion")
    expert_df$weights <- 1
  }
  
  if(!is.null(expert_df$wi)){
    expert_df$weights <- expert_df$wi
  }
  
  expert_df_sum <- expert_df %>% dplyr::group_by(times_expert) %>% dplyr::arrange(times_expert) %>%
    dplyr::summarize(sum_weights = sum(weights))
  
  if(any(expert_df_sum$sum_weights != 1)&& any(expert_df$weights != 1)){
    warning("Some weights don't sum to 1.. reweighting")
    
  }
  expert_df <-expert_df %>% dplyr::left_join(expert_df_sum,"times_expert") %>%
    dplyr::mutate(weights = weights/sum_weights)
  
  }
  
  expert_density <- apply(expert_df, 1, function(x){get_quant_val(
    dist = x["dist"],
    param1 = x["param1"],
    param2 = x["param2"],
    param3 = x["param3"],
    probs = probs)})
  
  rownames(expert_density) <- probs
  
  list(expert_df = expert_df  %>% dplyr::select(-sum_weights),
       expert_density = expert_density)
  
}


#Will modify for the SHINY APP
# expert_pooling <- function(expert_quant_list = NULL,
#                            lower_bound = -Inf, upper_bound = Inf, St_indic = 0){
# 
# dfs_expert <- list() 
# plts_pool <- list()
# dfs_pool <- list()
# 
# 
# #if(!is.null(expert_quant_list)){ # If a list of quantiles and probabilities
# 
#  if(is.null(expert_quant_list$weights_mat)){
#    weights_mat <- NULL
#  }
#   suppressMessages(attach(expert_quant_list))
#   
#     
# max.timepoints  <- length(times)
# 
# for(i in 1:max.timepoints){
#   
#   timepoint <- paste0("Time ",times[i])
#   
#   fit.eval <- SHELF::fitdist(vals = na.omit(v_array[,,i]),
#                       probs = na.omit(p_mat[,i]), lower = lower_bound, upper = upper_bound)
#   
#   weights <- na.omit(weights_mat[,i])
#   
#   if(is.null(weights_mat) && ncol(stats::na.omit(v_array[,,i])) == 1){
#     weights <- 1 #Only one expert so weights should be 1
#   }else if(is.null(weights_mat)){
#     warning("No weights assigned assuming equal weights")
#     weights <- 1
#   }
#   
#   best_fit_index  <- apply(fit.eval$ssq[,c("normal","t","gamma", "lognormal", "beta")], 1, which.min)
#   best_fit <- names(fit.eval$ssq[,c("normal","t","gamma", "lognormal", "beta")])[best_fit_index]
#   
#   best_fit_loc <- sapply(best_fit, function(x){which(x  == names(fit.eval$ssq))})
#   fit.eval.dist  <- fit.eval[best_fit_loc]
#   
#   pool.df_output <- matrix(nrow = length(best_fit_loc),ncol = 3)
#   colnames(pool.df_output) <- c("param1", "param2", "param3")
#   
#   for(i in 1:length(best_fit_loc)){
#     pool.df_output[i,1:length(fit.eval.dist[[i]][i,])] <-  as.numeric(as.vector(fit.eval.dist[[i]][i,]))
#   }
#   dfs_expert[[timepoint]] <- data.frame(dist = best_fit, wi = weights, pool.df_output)
# 
#   plts_pool[[timepoint]] <- makePoolPlot(fit  = fit.eval,
#                                          xl =lower_bound,
#                                          xu =upper_bound,
#                                          d = "best",
#                                          w = weights,
#                                          lwd =1,
#                                          xlab = "x",
#                                          ylab =expression(f[X](x)),
#                                          legend_full = TRUE,
#                                          ql = NULL,
#                                          qu = NULL,
#                                          nx = 200,
#                                          addquantile = FALSE,
#                                          fs = 12,
#                                          expertnames = NULL,
#                                          St_indic =St_indic
#                                          )
# 
#   dfs_pool[[timepoint]] <-  plts_pool[[timepoint]][["data"]]
# 
#   }
# 
# # }else{ #This isn't really needed will remove
# # 
# #   times <- unique(expert_density$expert_df[,"times_expert"])
# #   probs <- as.numeric(rownames(expert_density$expert_density))
# #   
# # for(i in 1:length(times)){
# #   
# #   timepoint <- paste0("Time ",times[i])
# #   
# #   index.loc <- which(expert_density$expert_df$times_expert == times[i])
# #   temp_df <- expert_density$expert_df[index.loc, ]
# #   temp_dens <- expert_density$expert_density[,index.loc]
# #   
# #   v <- temp_dens
# #   p <- matrix(rep(probs, ncol(temp_dens)), nrow = length(probs), ncol = ncol(temp_dens))
# #   
# #   # Need to consider upper and lower bounds
# #   fit.eval <- SHELF::fitdist(v, p, lower= lower_bound, upper = upper_bound)
# #   
# #   if(!is.null(temp_df$weights)){
# #     weights <- temp_df$weights
# #   }else{
# #     weights <- 1
# #   } 
# #   
# #   
# #   best_fit_index  <- apply(fit.eval$ssq[,c("normal","t","gamma", "lognormal", "beta")], 1, which.min)
# #   best_fit <- names(fit.eval$ssq[,c("normal","t","gamma", "lognormal", "beta")])[best_fit_index]
# #   
# #   best_fit_loc <- sapply(best_fit, function(x){which(x  == names(fit.eval$ssq))})
# #   fit.eval.dist  <- fit.eval[best_fit_loc]
# #   
# #   pool.df_output <- matrix(nrow = length(best_fit_loc),ncol = 3)
# #   colnames(pool.df_output) <- c("param1", "param2", "param3")
# #   
# #   for(i in 1:length(best_fit_loc)){
# #     pool.df_output[i,1:length(fit.eval.dist[[i]][i,])] <-  as.numeric(as.vector(fit.eval.dist[[i]][i,]))
# #   }
# #   dfs_expert[[timepoint]] <- data.frame(dist = best_fit, wi = weights, pool.df_output)
# #   
# #   
# #   plts_pool[[timepoint]] <- makePoolPlot(fit  = fit.eval,
# #                             xl =lower_bound,
# #                             xu =upper_bound,
# #                             d = "best",
# #                             w = weights,
# #                             lwd =1,
# #                             xlab = "x",
# #                             ylab =expression(f[X](x)),
# #                             legend_full = TRUE, 
# #                             ql = NULL,
# #                             qu = NULL,
# #                             nx = 200,
# #                             addquantile = FALSE,
# #                             fs = 12, 
# #                             expertnames = NULL,St_indic = St_indic)
# #   
# #   dfs_pool[[timepoint]] <-  plts_pool[[timepoint]][["data"]]
# #   
# #   
# #   }
# #   
# # }
#  list(dfs_expert =dfs_expert,
#       plts_pool =plts_pool,
#       dfs_pool = dfs_pool)
# 
# }

#myfit <- fitdist(expert_density, probs)
#plotfit(myfit,lp = T)
# Reweights the weights if they don't sum to 1 anyway



get_quant_val <- function(dist,param1, param2, param3 = NULL, probs = seq(0.01, 0.98, by = 0.01)){
  if(dist == "t"){
    probs_eval <- as.numeric(param1) + as.numeric(param2)*stats::qt(as.numeric(probs),as.numeric(param3))
    return(probs_eval)
     
  }else{
    probs <- paste0(probs, collapse = ",")
    
    probs_eval <-  paste0("q",dist,
                          "(c(",probs,"),", param1,
                          ",",param2,")")
    probs_eval <- eval(parse(text = probs_eval))
    return(probs_eval)
  }
 
}

pt.scaled <-function (q, df, mean = 0, sd = 1, ncp, lower.tail = TRUE, log.p = FALSE){
  stats::pt((q - mean)/sd, df, ncp = ncp, log.p = log.p)
}

get_cdf_val <- function(dist,param1, param2, param3 = NULL, vals = seq(0.01, 0.98, by = 0.01)){
  if(dist == "t"){
    probs_eval <- stats::pt((vals -  as.numeric(param1))/as.numeric(param2), as.numeric(param3), log.p = F)
    return(probs_eval)
    
  }else{
    vals <- paste0(vals, collapse = ",")
    
    probs_eval <-  paste0("p",dist,
                          "(c(",vals,"),", param1,
                          ",",param2,")")
    probs_eval <- eval(parse(text = probs_eval))
    return(probs_eval)
  }
  
}


get_density <- function(dist, param1, param2, param3 = NULL, x = seq(0.01, 0.98, by = 0.01), St_indic){
  x <- paste0(x, collapse = ",")
 
    if(St_indic ==1){
      a <- 0
      b <- 1
    }else{
      a <- -Inf
      b <- +Inf
    }

   if(dist == "t"){
    #From SHELF reference Student.t Parameters of the fitted t distributions. 
    #Note that (X - location) / scale has a standard t distribution 
    dens_x <-  paste0("d",dist,
                      "((c(",x,")-",param1,")/",param2,",", param3,")/",param2)
    cdf_a_b <-paste0("p",dist,
                     "((c(",a,",",b,")-",param1,")/",param2,",", param3,")")


     }else{ #
    dens_x <-  paste0("d",dist,
                      "(c(",x,"),", param1,
                      ",",param2,")")
    cdf_a_b <-  paste0("p",dist,
                        "(c(",a,",",b,"),", param1,
                          ",",param2,")")
    
     }
  k_trunc <- diff(eval(parse(text = cdf_a_b)))

  dens_eval <- eval(parse(text = dens_x))/k_trunc
  
  return(dens_eval)
}


#' Credible interval for pooled distribution
#'
#' Returns the interval based on defined quantiles. 
#' The approach used only provides an approximate (although quite accurate) integral.  
#' @param plt_obj A plot object from `plot_expert_opinion`.
#' @param val The name of the opinion for which the interval will be generated.
#' @param interval A vector of the upper and lower probabilities. Default is the standard 95% interval 
#' @keywords Expert
#' @return Credible interval based on the pooled distribution
#' @export
#' @examples
#' param_expert_example1 <- list()
#' param_expert_example1[[1]] <- data.frame(dist = c("norm","t"),
#'     wi = c(0.5,0.5), # Ensure Weights sum to 1
#'     param1 = c(0.1,0.12),
#'     param2 = c(0.005,0.005),
#'     param3 = c(NA,3))
#' plot_opinion1<- plot_expert_opinion(param_expert_example1[[1]], 
#'               weights = param_expert_example1[[1]]$wi)
#' cred_int(plot_opinion1,val = "linear pool", interval = c(0.025, 0.975))
#' 
#' 
cred_int <- function(plt_obj, val = "linear pool",interval = c(0.025, 0.975)){
  
  plt_df <- plt_obj$data %>% dplyr::filter(expert == val) %>% data.frame()
  
  total_integral <- integrate.xy(plt_df$x, plt_df$fx)
  partial_integral <- rep(NA, nrow(plt_df))
  partial_integral[1] <- 0
  for(i in 2:nrow(plt_df)){
    partial_integral[i] <- integrate.xy(plt_df$x[1:i], plt_df$fx[1:i])/total_integral
  }
  
  plt_df$cdf <- partial_integral
  
  limits <- c(plt_df$x[which.min(abs(plt_df$cdf - interval[1]))],plt_df$x[which.min(abs(plt_df$cdf - interval[2]))])
  names(limits) <- c("lower", "upper")
  return(limits)
  
}


#' makePoolPlot
#'
#' @param fit 
#' @param xl 
#' @param xu 
#' @param d 
#' @param w 
#' @param lwd 
#' @param xlab 
#' @param ylab 
#' @param legend_full 
#' @param ql 
#' @param qu 
#' @param nx 
#' @param addquantile 
#' @param fs 
#' @param expertnames 
#' @param St_indic 
#'
#' @import  ggplot2
#' @importFrom  scales hue_pal
#' @noRd
#' 
makePoolPlot <- function (fit, xl, xu, d = "best", w = 1, lwd = 1, xlab = "x", 
                          ylab = expression(f[X](x)), legend_full = TRUE, ql = NULL, 
                          qu = NULL, nx = 500, addquantile = FALSE, fs = 12, expertnames = NULL, 
                          St_indic){
  logt.error <- utils::getFromNamespace("logt.error", "SHELF")
  gamma.error <- utils::getFromNamespace("gamma.error", "SHELF")
  lognormal.error <- utils::getFromNamespace("lognormal.error", 
                                             "SHELF")
  logt.error <- utils::getFromNamespace("logt.error", "SHELF")
  makeGroupPlot <- utils::getFromNamespace("makeGroupPlot", 
                                           "SHELF")
  makeLinearPoolPlot <- utils::getFromNamespace("makeLinearPoolPlot", 
                                                "SHELF")
  makeSingleExpertPlot <- utils::getFromNamespace("makeSingleExpertPlot", 
                                                  "SHELF")
  expertdensity <- utils::getFromNamespace("expertdensity", 
                                           "SHELF")
  lpname <- c("linear pool", "log pool")
  expert <- ftype <- NULL
  n.experts <- nrow(fit$vals)
  if (length(d) == 1) {
    d <- rep(d, n.experts)
  }
  if (is.null(expertnames)) {
    if (n.experts < 27) {
      expertnames <- LETTERS[1:n.experts]
    }
    if (n.experts > 26) {
      expertnames <- 1:n.experts
    }
  }
  nxTotal <- nx + length(c(ql, qu))
  x <- matrix(0, nxTotal, n.experts)
  fx <- x
  if (min(w) < 0 | max(w) <= 0) {
    stop("expert weights must be non-negative, and at least one weight must be greater than 0.")
  }
  if (length(w) == 1) {
    w <- rep(w, n.experts)
  }
  weight <- matrix(w/sum(w), nxTotal, n.experts, byrow = T)
  sd.norm <- rep(NA, n.experts)
  for (i in 1:n.experts) {
  }
  if (is.infinite(xl) || is.infinite(xu)) {
    if (St_indic == 1) {
      xl <- 0
      xu <- 1
    }
    else {
      min.mean.index <- which.min(fit$Normal$mean)
      min.sd.index <- which.min(fit$Normal$sd)
      
      max.mean.index <- which.max(fit$Normal$mean)
      max.sd.index <- which.max(fit$Normal$sd)
      xl <- qnorm(0.001, fit$Normal[min.mean.index, 1], 
                  fit$Normal[min.sd.index, 2])
      xu <- qnorm(0.999, fit$Normal[max.mean.index, 1], 
                  fit$Normal[max.sd.index, 2])
    }
  }
  for (i in 1:n.experts) {
    densitydata <- expertdensity(fit, d[i], ex = i, xl, 
                                 xu, ql, qu, nx)
    x[, i] <- densitydata$x
    if (St_indic == 1) {
      k_trunc <- integrate.xy(x = x[, 1], fx = densitydata$fx)
    }
    else {
      k_trunc <- 1
    }
    fx[, i] <- densitydata$fx/k_trunc
  }
  fx.lp <- apply(fx * weight, 1, sum)
  if (any(is.infinite(fx^weight))) {
    warning("Print Non finite density for log pooling - Results invalid")
  }
  fx.logp <- apply(fx^weight, 1, prod)
  k_norm <- integrate.xy(x = x[, 1], fx = fx.logp)
  fx.logp <- fx.logp/k_norm
  df1 <- data.frame(x = rep(x[, 1], n.experts + 2), fx = c(as.numeric(fx), 
                                                           fx.lp, fx.logp), expert = factor(c(rep(expertnames, 
                                                                                                  each = nxTotal), rep("linear pool", nxTotal), rep("log pool", 
                                                                                                                                                    nxTotal)), levels = c(expertnames, "linear pool", "log pool")), 
                    ftype = factor(c(rep("individual", nxTotal * n.experts), 
                                     rep("linear pool", nxTotal), rep("log pool", nxTotal)), 
                                   levels = c("individual", "linear pool", "log pool")))
  df1$expert <- factor(df1$expert, levels = c(expertnames, 
                                              "linear pool", "log pool"))
  if (legend_full) {
    cols <- (scales::hue_pal())(n.experts + 2)
    linetypes <- c(rep("dashed", n.experts), "solid", "solid")
    sizes <- lwd * c(rep(0.5, n.experts), 1.5, 1.5)
    names(cols) <- names(linetypes) <- names(sizes) <- c(expertnames, 
                                                         lpname)
    p1 <- ggplot(df1, aes(x = x, y = fx, colour = expert, 
                          linetype = expert, size = expert)) + scale_colour_manual(values = cols, 
                                                                                   breaks = c(expertnames, lpname)) + scale_linetype_manual(values = linetypes, 
                                                                                                                                            breaks = c(expertnames, lpname)) + scale_size_manual(values = sizes, 
                                                                                                                                                                                                 breaks = c(expertnames, lpname))
  }
  else {
    p1 <- ggplot(df1, aes(x = x, y = fx, colour = ftype, 
                          linetype = ftype, size = ftype)) + scale_linetype_manual(name = "distribution", 
                                                                                   values = c("dashed", "solid", "solid")) + scale_size_manual(name = "distribution", 
                                                                                                                                               values = lwd * c(0.5, 1.5, 1.5)) + scale_color_manual(name = "distribution", 
                                                                                                                                                                                                     values = c("black", "red", "blue"))
  }
  if (legend_full) {
    for (i in 1:n.experts) {
      if (d[i] == "hist") {
        p1 <- p1 + geom_step(data = subset(df1, expert == 
                                             expertnames[i]), aes(colour = expert))
      }
      else {
        p1 <- p1 + geom_line(data = subset(df1, expert == 
                                             expertnames[i]), aes(colour = expert))
      }
    }
  }
  else {
    for (i in 1:n.experts) {
      if (d[i] == "hist") {
        p1 <- p1 + geom_step(data = subset(df1, expert == 
                                             expertnames[i]), aes(colour = ftype))
      }
      else {
        p1 <- p1 + geom_line(data = subset(df1, expert == 
                                             expertnames[i]), aes(colour = ftype))
      }
    }
  }
  if (length(unique(d)) == 1 & d[1] == "hist") {
    p1 <- p1 + geom_step(data = subset(df1, expert == lpname), 
                         aes(colour = expert))
  }
  else {
    p1 <- p1 + geom_line(data = subset(df1, expert == lpname[1]), 
                         aes(colour = expert))
    p1 <- p1 + geom_line(data = subset(df1, expert == lpname[2]), 
                         aes(colour = expert))
  }
  p1 <- p1 + labs(x = xlab, y = ylab)
  if ((!is.null(ql)) & (!is.null(qu)) & addquantile) {
    if (legend_full) {
      ribbon_col <- (scales::hue_pal())(n.experts + 2)[n.experts + 
                                                         2]
    }
    else {
      ribbon_col <- "red"
    }
    p1 <- p1 + geom_ribbon(data = with(df1, subset(df1, 
                                                   x <= ql & expert == lpname[1])), aes(ymax = fx, 
                                                                                        ymin = 0), alpha = 0.2, show.legend = FALSE, colour = NA, 
                           fill = ribbon_col) + geom_ribbon(data = with(df1, 
                                                                        subset(df1, x >= qu & expert == lpname[2])), aes(ymax = fx, 
                                                                                                                         ymin = 0), alpha = 0.2, show.legend = FALSE, colour = NA, 
                                                            fill = ribbon_col)
    p1 <- p1 + geom_ribbon(data = with(df1, subset(df1, 
                                                   x <= ql & expert == lpname[2])), aes(ymax = fx, 
                                                                                        ymin = 0), alpha = 0.2, show.legend = FALSE, colour = NA, 
                           fill = ribbon_col) + geom_ribbon(data = with(df1, 
                                                                        subset(df1, x >= qu & expert == lpname[2])), aes(ymax = fx, 
                                                                                                                         ymin = 0), alpha = 0.2, show.legend = FALSE, colour = NA, 
                                                            fill = ribbon_col)
  }
  if (lpname[1] == "marginal") {
    p1 <- p1 + theme(legend.title = element_blank())
  }
  p1 + theme(text = element_text(size = fs))
}




#' Plotting Pooled Expert Opinion
#'
#' Returns a ggplot with the individual expert opinions along with the pooled distributions (both linear and logarithmic).
#'
#' @param object Either a object of class elicitation (from `SHELF`) or a dataframe with parameters of the distribution (see Example below).
#' @param xl_plt Optionally set the lower bound for the plot
#' @param xu_plt Optionally set the upper bound for the plot
#' @param weights A vector with the weight of each expert. If omitted, set to equal weights.
#' @param St_indic Set to 1 if you want to truncate the distributions to be between 0 and 1.
#' @return A ggplot with pooled distributions.
#' @export
#' @examples 
#'  expert_df <- data.frame(dist = c("norm","t"), #Distribution Name
#'                          wi = c(1/3,2/3), #Expert weights
#'                          param1 = c(0.3,0.40), #Parameter 1
#'                          param2 = c(0.05,0.05),# Parameter 2
#'                          param3 = c(NA,3)) #Parameter 3: Only t-distribution
#'  plot_expert_opinion(expert_df , weights = expert_df$wi)
#'
#'                                                         
plot_expert_opinion <- function(object, xl_plt = NULL, xu_plt = NULL, weights = NULL, St_indic =0){
  
  
  if(is.null(weights)){
    weights <- 1
  }
  
  
  
  if(inherits(object,"elicitation")){

      if(is.null(xl_plt)){
        xl_plt <- min(object$limits["lower"])
      }
  


    if(is.null(xu_plt)){
      xu_plt <- max(object$limits["upper"])
      
    }
    
    if(St_indic ==1){
      xl_plt <- max(0, xl_plt)
      xu_plt  <- min(1,xu_plt)
    }
    
    plt <- makePoolPlot(fit= object,
                                     xl =xl_plt,
                                     xu =xu_plt,
                                     d = "best",
                                     w = weights,
                                     lwd =1,
                                     xlab = "x",
                                     ylab =expression(f[X](x)),
                                     legend_full = TRUE,
                                     ql = NULL,
                                     qu = NULL,
                                     nx = 200,
                                     addquantile = FALSE,
                                     fs = 12,
                                     expertnames = NULL,
                                     St_indic =  St_indic)
    
    
  }else{
    
    object$times_expert <- 2 #Just for compatibility
    
    expert_dens_list <- expert_dens(object, probs =  seq(0.001, 0.99, by = 0.005))
    
    lower <- as.numeric(utils::head(expert_dens_list$expert_density, n = 1)-0.1)
    upper <- as.numeric(utils::tail(expert_dens_list$expert_density, n = 1)+0.1)
    
    # if(is.null(lower) || is.null(upper)){
    #   stop("Upper and lower bounds required for distributions")
    # }
    
    if(is.null(xl_plt)){
      xl_plt <- min(lower)
    }
    if(is.null(xu_plt)){
      xu_plt <- max(upper)
      
    }
    if(St_indic ==1){
      xl_plt <- max(0, xl_plt)
      xu_plt  <- min(1,xu_plt)
    }
    
    
    probs_mat <- matrix(as.numeric(rep(rownames(expert_dens_list$expert_density), 
                                       dim(expert_dens_list$expert_density)[2])),
                        ncol = dim(expert_dens_list$expert_density)[2])
    
    fit_shelf  <- SHELF::fitdist(vals = expert_dens_list$expert_density,
                          probs_mat, lower = lower, upper = upper)
    
    plt <- makePoolPlot(fit= fit_shelf,
                                     xl = xl_plt,
                                     xu = xu_plt,
                                     d = "best",
                                     w = weights,
                                     lwd =1,
                                     xlab = "x",
                                     ylab =expression(f[X](x)),
                                     legend_full = TRUE,
                                     ql = NULL,
                                     qu = NULL,
                                     nx = 200,
                                     addquantile = FALSE,
                                     fs = 12,
                                     expertnames = NULL,
                                     St_indic =  St_indic)
    
    
  }
  
  return(plt+theme_bw())
}

