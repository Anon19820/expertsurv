library(expertsurv)
require("dplyr")
param_expert_example1 <- list()
#1 timepoint and 2 experts with equal weight,
#first a normal distribution, second a non-standard t-distribution with
#3 degrees of freedom

param_expert_example1[[1]] <- data.frame(dist = c("norm","t"),
                                         wi = c(0.5,0.5), # Ensure Weights sum to 1
                                         param1 = c(0.1,0.12),
                                         param2 = c(0.005,0.005),
                                         param3 = c(NA,3))


timepoint_expert <- 14

data2 <- data %>% rename(status = censored) 

data2$time2 <- runif(length(data2$time),max = 1)
data2$status2 <- 0

modelinits_gom <- function(){
  #beta[1] is eqivalent to log(rate) in flexsurv
  #beta[2:n_param] is eqivalent beta coeffiecient in a PH model
  beta = c(log(1.60557e-03), 0)
  
  #alpha[1],alpha[2] is eqivalent to shape in flexsurv
  #alpha[2] is constrained to be positive in the case that opinion is on 
  #mean survival as the expected survival would not be finite
  alpha1 = 1.40515e+02
  alpha2 = 1.40515e+02
  
  list(alpha1=alpha1,alpha2=alpha2, beta = beta) 
}

init_eval <- list(gom = modelinits_gom)
example1  <- fit.models.expert(formula=Surv(time2,status2)~1,data=data2,
                               distr=c("gomp"),
                               method="bayes",
                               pool_type = "linear pool", 
                               opinion_type = "survival",
                               times_expert = timepoint_expert, 
                               param_expert = param_expert_example1,
                               inits = init_eval)


example1  <- fit.models.expert(formula=Surv(time2,status2)~1,data=data2,
                               distr=c("gga"),
                               method="bayes",
                               pool_type = "linear pool", 
                               opinion_type = "survival",
                               times_expert = timepoint_expert, 
                               param_expert = param_expert_example1,
                               inits = init_eval)


example1  <- fit.models.expert(formula=Surv(time2,status2)~1,data=data2,
                               distr=c("wph"),
                               method="bayes",
                               pool_type = "linear pool", 
                               opinion_type = "survival",
                               times_expert = timepoint_expert, 
                               param_expert = param_expert_example1,
                               inits = init_eval,
                               iter = 2000)



# \dontrun{


# \dontrun{
model_code <-  "functions{
  real log_density_dist(array[ , ] real params,
                        real x,int num_expert, int pool_type){

    // Evaluates the log density for a range of distributions

    array[num_expert] real dens;

    for(i in 1:num_expert){
    if(params[i,1] == 1){
      if(pool_type == 1){
        dens[i] = exp(normal_lpdf(x|params[i,3], params[i,4]))*params[i,2]; /// Only require the log density is correct to a constant of proportionality
      }else{
        dens[i] = exp(normal_lpdf(x|params[i,3], params[i,4]))^params[i,2]; /// Only require the log density is correct to a constant of proportionality
      }

    }else if(params[i,1] == 2){
       if(pool_type == 1){
          dens[i] = exp(student_t_lpdf(x|params[i,5],params[i,3], params[i,4]))*params[i,2];
        }else{
          dens[i] = exp(student_t_lpdf(x|params[i,5],params[i,3], params[i,4]))^params[i,2];
       }

    }else if(params[i,1] == 3){
        if(pool_type == 1){
            dens[i] = exp(gamma_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
            dens[i] = exp(gamma_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }

    }else if(params[i,1] == 4){

        if(pool_type == 1){
            dens[i] = exp(lognormal_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
             dens[i] = exp(lognormal_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }

    }else if(params[i,1] == 5){
     if(pool_type == 1){
            dens[i] = exp(beta_lpdf(x|params[i,3], params[i,4]))*params[i,2];
        }else{
            dens[i] = exp(beta_lpdf(x|params[i,3], params[i,4]))^params[i,2];
        }


      }

    }


      if(pool_type == 1){
      return(log(sum(dens)));
      }else{
      return(log(prod(dens)));
      }

  } }"




expose_stan_functions(stanc(model_code = model_code))


list_of_objs <-   lapply(1:nrow(example1$misc$data.stan[[1]]$param_expert), FUN = function(i){example1$misc$data.stan[[1]]$param_expert[i,,]})

x <- seq(0, 1, by =0.001)
fx <- rep(NA, length(x))
for(i in 1:length(fx)){
  fx[i] <- log_density_dist(list_of_objs, x[i], 2, 1) 
}


plot(example1, add.km = TRUE, t = 0:30, plot_opinion = TRUE,plot_ci = TRUE, nsim  = 1000)
gg_plot_expert <- plot_expert_opinion(param_expert_example1[[1]],weights = param_expert_example1[[1]]$wi, St_indic = 1)

gg_data <- gg_plot_expert$data %>% filter(expert == "linear pool")

#plot(density(example1$models$`Gen. Gamma`$BUGSoutput$sims.matrix[, "St_expert"]))
lines(gg_data$x , gg_data$fx, col = "red")


# Extracting specific parameter 'theta'
St_samples <- rstan::extract(example1$models$`Weibull (PH)`, pars = "St_expert")
plot(density(St_samples$St_expert), xlim = c(0,0.3))
lines(gg_data$x , gg_data$fx, col = "red")
lines(y = exp(fx), x, col = "blue")


#Maybe if I did the sampling on the log scale. Check when all has been implemented.
library(ggmcmc)
library(coda)
ggmcmc(ggs(coda::as.mcmc(example1$models$`Weibull (PH)`)), file = "Gengamma.pdf")

ggmcmc(ggs(rstan::As.mcmc.list(example1$models$`Weibull (PH)`)), file = "weibull.pdf")





# Define the expressions
weibull_ph <- expression(exp(-scale * t^shape))
weibull_aft <- expression(exp(-(t/scale)^shape))
loglogistic <- expression(1 / (1 + (t / b)^a))

# Find the derivatives with respect to 'a'
derivative_weibull_ph_shape <- D(weibull_ph, "shape")
derivative_weibull_ph_scale <- D(weibull_ph, "scale")
derivative_weibull_aft <- D(weibull_aft, "a")
derivative_loglogistic <- D(loglogistic, "a")




Gompertz.jags <- "
data{
for(i in 1:n){
zeros[i] <- 0
}
for(i in 1:n_time_expert){
zero2[i] <- 0
}
zero_prior <- 0 
C <- 1000 
}

model{


for(i in 1:n){
zeros[i] ~ dpois(zero.mean[i])

log_h[i] = log(mu[i]) + (alpha*t[i]);
log_S[i] = -mu[i]/alpha * (exp(alpha*t[i]) - 1);

LL[i] <- log_h[i]*d[i] + log_S[i] + log(a0[i])
zero.mean[i] <- -LL[i] + C

}

alpha1 ~ dnorm(mu_beta[1],1/(sigma_beta[1])^2)
alpha2 ~ dgamma(a_alpha, b_alpha)
alpha <-   alpha1*ifelse(St_indic == 1, 1,0) +alpha2*ifelse(St_indic == 0, 1,0)


for(i in 1:H){
 prec_beta[i] <- 1/(sigma_beta[i])^2
 beta[i] ~ dnorm(mu_beta[i],prec_beta[i])
}

linpred <- X %*%beta
for( i in 1:n){
mu[i] <- exp(linpred[i])
}


for (i in 1:n_time_expert){
  zero2[i] ~ dpois(phi_con[i])
 
 
    St_expert_temp[i,1] = exp(-mu[id_St]/alpha * (exp(alpha*time_expert[i]) - 1))*St_indic
    mean_surv_trt[i] <- (1/alpha)*exp(mu[id_trt]/alpha)*(exp(loggam(s))*(1-pgamma(mu[id_trt]/max(alpha,0.000000001),s,1)))
    mean_surv_comp[i] <- (1/alpha)*exp(mu[id_comp]/alpha)*(exp(loggam(s))*(1-pgamma(mu[id_comp]/max(alpha,0.0000001),s,1)))
    St_expert_temp[i,2] <- (mean_surv_trt[i] - mean_surv_comp[i])*(ifelse(St_indic == 1, 0,1))
    
    St_expert[i] <- sum(St_expert_temp[i,])

    for(j in 1:n_experts[i]){
    expert_dens[j,1,i] <-  dnorm(St_expert[i], param_expert[j,3,i],pow(param_expert[j,4,i],-2))
    expert_dens[j,2,i] <-  dt(St_expert[i],    param_expert[j,3,i],pow(param_expert[j,4,i],-2),max(param_expert[j,5,i],1)) 
    expert_dens[j,3,i] <-  dgamma(St_expert[i], max(param_expert[j,3,i],0.001),param_expert[j,4,i])
    expert_dens[j,4,i] <-  dlnorm(St_expert[i], param_expert[j,3,i],param_expert[j,4,i])
    expert_dens[j,5,i] <-  dbeta(St_expert[i], max(param_expert[j,3,i], 0.01),param_expert[j,4,i])
    phi_temp[j,i] <- equals(pool_type,1)*(expert_dens[j,param_expert[j,1,i],i]*param_expert[j,2,i])+equals(pool_type,0)*(expert_dens[j,param_expert[j,1,i],i]^param_expert[j,2,i])
    }
 
  phi_con[i] <- -log(sum(phi_temp[,i]))*equals(pool_type,1) +  -log(prod(phi_temp[,i]))*equals(pool_type,0) + C #+test_St_expert[i]*C

 } 
 
  deriv1 <- abs(exp(-(mu[id_St]/alpha)*(exp(alpha*time_expert[1]) - 1)) * (mu[id_St]/(alpha^2)*(exp(alpha*time_expert[1])-1) - (mu[id_St]/alpha)*(exp(alpha*time_expert[1]) * time_expert[1])))
  deriv2 <- abs(-(exp(-(mu[id_St]/alpha)*(exp(alpha*time_expert[1]) - 1))*((1/alpha)*(exp(alpha*time_expert[1]) - 1))))
  
  zero.mean_prior <- -log(deriv1+deriv2) + C
  zero_prior ~ dpois(zero.mean_prior)
 
 
rate = exp(beta[1])

s <- 0.0001

}"






gompertz.jags <- "data{
  for(i in 1:N){
    zero[i] <- 0
  }
  zero_prior <- 0
}
model{
  C <- 10000
  for(i in 1:N){
    logHaz[i] <- (log(b) + a * time[i]) * status[i]
    logSurv[i] <- (-b / a) * (exp(a * time[i]) - 1)
    LL[i] <- logHaz[i] + logSurv[i]
    Like[i] <- exp(LL[i])
    #invLik[i] <- 1/Like[i] Unstable for some datasets (Will calculate outside JAGS)
    zero[i] ~ dpois(zero.mean[i])
    zero.mean[i] <- -logHaz[i] - logSurv[i] + C
  }
  for(i in 1:length(t_pred)){
    St_pred[i] <- exp((-b / a) * (exp(a * t_pred[i]) - 1))
  }
  
  deriv1 <- abs(exp(-(b/a)*(exp(a*time_expert[1]) - 1)) * (b/(a^2)*(exp(a*time_expert[1])-1) - (b/a)*(exp(a*time_expert[1]) * time_expert[1])))
  deriv2 <- abs(-(exp(-(b/a)*(exp(a*time_expert[1]) - 1))*((1/a)*(exp(a*time_expert[1]) - 1))))
  # 
  zero.mean_prior <- 0#-log(deriv1+deriv2) + C
  zero_prior ~ dpois(zero.mean_prior)
  
  
  a ~ dunif(0, 1000)
  b ~ dunif(0, 1000)
  total_LLik <- sum(log(Like))
}"



inits_list <- function(mod, n.chains = 2) {
  list_return <- list()
  for (i in 1:n.chains) {
    list_inits <- list()
    list_inits$t <- tinits1 + runif(1)
    if (mod == "exp") {
      #   list_inits$lambda = 1/mean(df$time)
      list_inits$lambda =0.5
    }
    if (mod == "weibull") {
      # lt <- log(df$time[df$time > 0])
      # shape <- 1.64/var(lt)
      # scale <- exp(mean(lt) + 0.572)
      # list_inits$v <- shape
      # list_inits$lambda <- scale^{
      #   -shape
      # }
      
      list_inits$v <- 1
      list_inits$lambda <- -log(0.5)/10
      
    }
    if (mod == "gompertz") {
      list_inits$a = 0.001
      list_inits$b = 1/mean(df$time)
      list_inits <- list_inits[names(list_inits) %!in%
                                 c("t")]
    }
    if (mod == "lnorm") {
      lt <- log(df$time[df$time > 0])
      list_inits$mu <- mean(lt)
      list_inits$sd <- sd(lt)
    }
    if (mod == "llogis") {
      lt <- log(df$time[df$time > 0])
      list_inits$mu <- mean(lt)
      list_inits$scale <- 3 * var(lt)/(pi^2)
      list_inits$t.log <- log(tinits1 + runif(1))
      list_inits <- list_inits[names(list_inits) %!in%
                                 c("t")]
    }
    if (mod == "gengamma") {
      list_inits$r <- 1
      list_inits$lambda <- 1/mean(df$time)
      list_inits$b <- 1
    }
    if (mod == "gamma") {
      list_inits$lambda = sum(df$time)
      list_inits$shape = sum(df$status)
    }
    list_return[[i]] <- list_inits
  }
  return(list_return)
}


t_pred <- 10
df <- data.frame(time = c(0, 0,0.01), status = c(0,0,0))


n.burnin.jags <- 100
n.iter.jags <- 100000
n.thin.jags <- 1
data_new <- list()
df_jags <- df[, c("time", "status")]
df_jags$t <- df$time
tinits1 <- df_jags$t + max(df$time)
is.na(tinits1) <- df_jags$status == 1
tinits2 <- tinits1 + 5
is.na(df_jags$t) <- df_jags$status == 0
df_jags$is.censored <- 1 - df_jags$status
df_jags$t.cen <- df_jags$time + df_jags$status
data_jags <- list(N = nrow(df_jags), t.cen = df_jags$t.cen,
                  is.censored = df_jags$is.censored, t = df_jags$t)
data_jags$t_pred <- t_pred
data_jags_llogis <- data_jags
data_jags_llogis$t.log <- log(data_jags$t)
data_jags_llogis$t.cen.log <- log(data_jags$t.cen)
`%!in%` <- Negate(`%in%`)
data_jags_llogis <- data_jags_llogis[names(data_jags_llogis) %!in%
                                       c("t", "t.cen")]
data_gomp <- list()
data_gomp$time <- df$time
data_gomp$status <- df$status
data_gomp$N <- nrow(df)
data_gomp$t_pred <- data_jags$t_pred
n.chains = 2
data_gomp$time_expert <- data_jags$t_pred
cat("Exponential Model \n")


# Its the derivative with respect to the parameter in the model
# the contribution is d'S(t) with respect to lambda. So it is t*exp(-lambda*t) 

expo.mod <- R2jags::jags(model.file = textConnection(gompertz.jags),
                         data = data_gomp, 
                         n.chains = n.chains, parameters.to.save = c("Like",
                                                                     "lambda", "St_pred", "total_LLik","zero.mean_prior"), n.iter = n.iter.jags,
                         n.thin = n.thin.jags, n.burnin = n.burnin.jags)



