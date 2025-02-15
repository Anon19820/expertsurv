if(FALSE){
library("truncnorm")
expo <- "
data{
    zero_prior <- 0
    C <- 10000
    #t_expert <- 10
    #mu_expert <- 0.1
    #sd_expert <- 0.05
    tau_expert <- pow(sd_expert,-2)
}

model{

  for(i in 1:N){
    is.censored[i] ~ dinterval(t[i], t.cen[i])
    t[i] ~ dexp(lambda)
    Like[i] <- ifelse(is.censored[i], 1 - pexp(t.cen[i], lambda), dexp(t[i], lambda))
    #invLik[i] <- 1/Like[i] Unstable for some datasets (Will calculate outside JAGS)
  }
  for(i in 1:length(t_pred)){
    St_pred[i] <- 1 - pexp(t_pred[i], lambda)
  }
  lambda ~ dunif(0, 10000)
  St_expert <- 1 - pexp(t_expert, lambda)
  zero.mean_prior <- -log(t_expert*exp(-lambda*t_expert)) - logdensity.norm(St_expert, mu_expert, tau_expert) +  C #derivative of St for exponential is just the pdf
  zero_prior ~ dpois(zero.mean_prior)
  total_LLik <- 0#sum(log(Like))
}"

weibull <- "
data{
    zero_prior <- 0
    C <- 1000000
    tau_expert <- pow(sd_expert,-2)
}

model{

  for(i in 1:N){
    is.censored[i] ~ dinterval(t[i], t.cen[i])
    t[i] ~ dweib(v, lambda)
    Like[i] <- ifelse(is.censored[i], 1 - pweib(t.cen[i], v, lambda), dweib(t[i], v, lambda))
    #invLik[i] <- 1/Like[i] Unstable for some datasets (Will calculate outside JAGS)
    #x[i] <- 1
  }
  for(i in 1:length(t_pred)){
    St_pred[i] <- 1 - pweib(t_pred[i], v, lambda)
  }
  
  lambda ~ dunif(0, 10)
  v ~  dunif(0, 10)
  #beta ~ dnorm(0, pow(5,-2))
  #lambda <- exp(beta)
  #v ~  dgamma(0.1, 0.10)
  St_expert <- 1 - pweib(t_expert, v, lambda)
  
  
  deriv1 <- abs(-pow(t_expert,v)*log(t_expert)*lambda*exp(-pow(t_expert,v)*lambda))
  deriv2 <- abs(-t_expert^v*exp(-pow(t_expert,v)*lambda))
  
  zero.mean_prior <- -log(deriv1+deriv2) - logdensity.norm(St_expert, mu_expert, tau_expert) + C
  expert_log_dens <- logdensity.norm(St_expert, mu_expert, tau_expert)
  zero_prior ~ dpois(zero.mean_prior)

  total_LLik <- 0#sum(log(Like))
}"

gamma.jags <- "data{
    zero_prior <- 0
    C <- 1000000
    tau_expert <- pow(sd_expert,-2)
    epsilon <- 0.001
}

model{
  for(i in 1:N){
    is.censored[i] ~ dinterval(t[i], t.cen[i])
    t[i] ~ dgamma(shape, lambda)
    Like[i] <- ifelse(is.censored[i], 1 - pgamma(t.cen[i], shape, lambda), dgamma(t[i], shape, lambda))
    #invLik[i] <- 1/Like[i] Unstable for some datasets (Will calculate outside JAGS)
  }
  for(i in 1:length(t_pred)){
    St_pred[i] <- 1 - pgamma(t_pred[i], shape, lambda)
  }
  
  St_expert <- 1 - pgamma(t_expert, shape, lambda)
  
  lambda ~ dunif(epsilon*2, 100)
  shape ~ dunif(epsilon*2, 100)
  
  deriv1 <- abs((1- pgamma(t_expert[1],shape+epsilon,lambda)) - (1 - pgamma(t_expert[1],shape-epsilon,lambda)))/(2*epsilon)
  deriv2 <- abs((1- pgamma(t_expert[1],shape,lambda+epsilon)) - (1 - pgamma(t_expert[1],shape,lambda-epsilon)))/(2*epsilon)

  zero.mean_prior <- -log(deriv1+deriv2) - logdensity.norm(St_expert, mu_expert, tau_expert) + C
  zero_prior ~ dpois(zero.mean_prior)
 
  total_LLik <- 0#sum(log(Like))
}"

lnorm.jags <- "
data{
    zero_prior <- 0
    C <- 1000000
    tau_expert <- pow(sd_expert,-2)
    epsilon <- 0.001
}

model{
  for(i in 1:N){
    is.censored[i] ~ dinterval(t[i], t.cen[i])
    t[i] ~ dlnorm(mu, tau)
    Like[i] <- ifelse(is.censored[i], 1 - plnorm(t.cen[i], mu, tau), dlnorm(t[i], mu, tau))
    #invLik[i] <- 1/Like[i] Unstable for some datasets (Will calculate outside JAGS)
  }
  for(i in 1:length(t_pred)){
    St_pred[i] <- 1 - plnorm(t_pred[i], mu, tau)
  }
  St_expert <- 1 - plnorm(t_expert, mu, tau)
  
  deriv1 <- abs((1- plnorm(t_expert[1],mu+epsilon,tau)) - (1 - plnorm(t_expert[1],mu-epsilon,tau)))/(2*epsilon)
  deriv2 <- abs((1- plnorm(t_expert[1],mu,tau+epsilon)) - (1 - plnorm(t_expert[1],mu,tau-epsilon)))/(2*epsilon)

  zero.mean_prior <- -log(deriv1+deriv2) - logdensity.norm(St_expert, mu_expert, tau_expert) + C
  zero_prior ~ dpois(zero.mean_prior)
  
  mu ~ dnorm(0, 0.001)
  sd ~ dunif(0.01, 10)
  tau <- pow(sd, -2)

  total_LLik <- 0#sum(log(Like))
}"

llogis.jags <- "
data{
    zero_prior <- 0
    C <- 1000000
    tau_expert <- pow(sd_expert,-2)
      epsilon <- 0.001
}

model{
  for(i in 1:N){
    #is.censored[i] ~ dinterval(t.log[i], t.cen.log[i])
    #t.log[i] ~ dlogis(mu, tau)
    #Like[i] <- ifelse(is.censored[i], 1 / (1 + pow(exp(t.cen.log[i]) / beta, alpha)), 
    #                  (alpha / beta) * pow(exp(t.log[i]) / beta, alpha - 1) / pow(1 + pow(exp(t.log[i]) / beta, alpha), 2))
    #invLik[i] <- 1/Like[i] Unstable for some datasets (Will calculate outside JAGS)
    x[i] <- 1
  }
  for(i in 1:length(t_pred)){
    St_pred[i] <- 1 / (1 + pow(t_pred[i] / beta, alpha))
  }
   #St_expert <- 1- plogis(log(t_expert),mu, tau)
   St_expert <-  1 / (1 + pow(t_expert / beta, alpha))
 
  
	#deriv1 = (alpha*pow(t_expert[1]/beta,alpha))/(beta*pow(pow(t_expert[1]/beta,alpha)+1,2));
	#deriv2 = -(pow(t_expert[1]/beta,alpha)*log(t_expert[1]/beta))/pow(pow(t_expert[1]/beta,alpha)+1,2);
  # The derivative needed to be with respect to the mu and tau not the log-logistic parameters.
  deriv1 <- ((1-plogis(log(t_expert[1]),mu+epsilon,tau)) - (1-plogis(log(t_expert[1]),mu-epsilon,tau)))/(2*epsilon)
  deriv2 <- ((1-plogis(log(t_expert[1]),mu,tau+epsilon)) - (1-plogis(log(t_expert[1]),mu,tau-epsilon)))/(2*epsilon)
  #deriv1 <- 1
  #deriv2 <- 1 
  zero.mean_prior <- -logdensity.norm(St_expert, mu_expert, tau_expert)  + C -log(abs(deriv1)+abs(deriv2))
  zero_prior ~ dpois(zero.mean_prior)
  
  mu ~ dnorm(0, 0.001)
  scale ~ dgamma(0.001, 0.001) # Dont use this anymore dgamma gave weird results
  #tau <- pow(scale, -1) # Inverse of scale which is beta on the log-logistic dist
  tau ~ dunif(0.001*3, 10) # Try the prior on the tau
  beta <- exp(mu)
  alpha <- tau
  total_LLik <- 0#sum(log(Like))
}"

gompertz.jags <- "data{
  for(i in 1:N){
    zero[i] <- 0
  }
    zero_prior <- 0
    C <- 1000000
    tau_expert <- pow(sd_expert,-2)
    epsilon <- 0.001
}
model{

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
    St_expert <- exp((-b / a) * (exp(a * t_expert) - 1))
 
  deriv1 <- abs(exp(-(b/a)*(exp(a*t_expert[1]) - 1)) * (b/(a^2)*(exp(a*t_expert[1])-1) - (b/a)*(exp(a*t_expert[1]) * t_expert[1])))
  deriv2 <- abs(-(exp(-(b/a)*(exp(a*t_expert[1]) - 1))*((1/a)*(exp(a*t_expert[1]) - 1))))

  zero.mean_prior <- -log(deriv1+deriv2) - logdensity.norm(St_expert, mu_expert, tau_expert) + C
  zero_prior ~ dpois(zero.mean_prior)

  a ~ dnorm(0, 0.001)
  b ~ dunif(0, 10)
  total_LLik <- 0#sum(log(Like))
}"

gen.gamma.jags <- "
data{
    zero_prior <- 0
    C <- 1000000
    tau_expert <- pow(sd_expert,-2)
    epsilon <- 0.001
}


model{
  for(i in 1:N){
    is.censored[i] ~ dinterval(t[i], t.cen[i])
    t[i] ~ dgen.gamma(r, lambda, b)
    Like[i] <- ifelse(is.censored[i], 1 - pgen.gamma(t.cen[i], r, lambda, b), dgen.gamma(t[i], r, lambda, b))
    #invLik[i] <- 1/Like[i] Unstable for some datasets (Will calculate outside JAGS)
  }
  for(i in 1:length(t_pred)){
    St_pred[i] <- 1 - pgen.gamma(t_pred[i], r, lambda, b)
  }
  #r ~ dgamma(0.001, 0.001)T(,10)
  #lambda ~ dgamma(0.001, 0.001)T(,10)
  #b ~ dgamma(0.001, 0.001)T(,10)
  
  r ~ dunif(0.01, 10)
  lambda ~ dunif(0.01, 10)
  b ~ dunif(0.01, 10)


  total_LLik <- 0#sum(log(Like))
  
  St_expert <-  1 - pgen.gamma(t_expert, r, lambda, b)
 
  deriv1 <- ((1-pgen.gamma(t_expert[1],r+epsilon,lambda,b)) - (1-pgen.gamma(t_expert[1],r-epsilon,lambda,b)))/(2*epsilon)
  deriv2 <- ((1-pgen.gamma(t_expert[1],r,lambda+epsilon,b)) - (1-pgen.gamma(t_expert[1],r,lambda-epsilon,b)))/(2*epsilon)
  deriv3 <- ((1-pgen.gamma(t_expert[1],r,lambda,b+epsilon)) - (1-pgen.gamma(t_expert[1],r,lambda,b-epsilon)))/(2*epsilon)

  zero.mean_prior <- -log(abs(deriv1)+abs(deriv2)+abs(deriv3)) - logdensity.norm(St_expert, mu_expert, tau_expert) +  C
  zero_prior ~ dpois(zero.mean_prior)
}"



t_expert <- t_pred <- 10
mu_expert <- 0.1
sd_expert <- 0.05

inits_list <- function(mod, n.chains = 2) {
  list_return <- list()
  for (i in 1:n.chains) {
    list_inits <- list()
    list_inits$t <- tinits1 + runif(1)
    if (mod == "exp") {
      #   list_inits$lambda = 1/mean(df$time)
      list_inits$lambda =-log(mu_expert)/t_expert
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
      list_inits$lambda <- -log(mu_expert)/t_expert

    }
    if (mod == "gompertz") {
      list_inits$a = -log(mu_expert)/t_expert
      list_inits$b = 1
      list_inits <- list_inits[names(list_inits) %!in%
                                 c("t")]
    }
    if (mod == "lnorm") {
      lt <- log(df$time[df$time > 0])
      list_inits$mu <- mean(lt)
      list_inits$sd <- sd(lt)
    }
    if (mod == "llogis") {
      # lt <- log(df$time[df$time > 0])
      list_inits$mu <- 1#mean(lt)
      list_inits$scale <- 1#3 * var(lt)/(pi^2)
      list_inits$t.log <- 1#log(tinits1 + runif(1))
      list_inits <- list_inits[names(list_inits) %!in%
                                 c("t")]
    }
    if (mod == "gengamma") {
      list_inits$r <- 1
      list_inits$lambda <-  -log(mu_expert)/t_expert #1/mean(df$time)
      list_inits$b <- 1
    }
    if (mod == "gamma") {
      list_inits$lambda = -log(mu_expert)/t_expert
      list_inits$shape = 1#sum(df$status)
    }
    list_return[[i]] <- list_inits
  }
  return(list_return)
}



df <- data.frame(time = c(0.001, 0.0001), status = c(0,0))


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
data_jags$mu_expert <- mu_expert
data_jags$sd_expert <- sd_expert
data_jags$t_expert <- t_expert

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
data_gomp$N <- 0#nrow(df)
data_gomp$t_pred <- data_jags$t_pred
data_gomp$mu_expert <- mu_expert
data_gomp$sd_expert <- sd_expert
data_gomp$t_expert <- t_expert

n.chains = 2
cat("Exponential Model \n")


# Its the derivative with respect to the parameter in the model
# the contribution is d'S(t) with respect to lambda. So it is t*exp(-lambda*t) 

expo.mod <- R2jags::jags(model.file = textConnection(expo),
                         data = data_jags, inits = inits_list("exp", n.chains),
                         n.chains = n.chains, parameters.to.save = c("Like",
                                                                     "lambda", "St_pred", "total_LLik","zero.mean_prior"), n.iter = n.iter.jags,
                         n.thin = n.thin.jags, n.burnin = n.burnin.jags)

weib.mod <- R2jags::jags(model.file = textConnection(weibull),
                         data = data_jags, inits = inits_list("weibull", n.chains),
                         n.chains = n.chains, parameters.to.save = c("lambda",
                                                                     "v", "Like", "St_pred", "total_LLik", "deriv1","deriv2","zero.mean_prior","expert_log_dens"), n.iter = n.iter.jags,
                         n.thin = n.thin.jags, n.burnin = n.burnin.jags)

if(FALSE){
plot(weib.mod$BUGSoutput$sims.matrix[,"lambda"], y = weib.mod$BUGSoutput$sims.matrix[,"v"])
df <- weib.mod$BUGSoutput$sims.matrix[,c("lambda","v")]
library(ggplot2)
 
 # Create a sample data frame
 set.seed(123)
 df <- data.frame(
   x = rnorm(100),
   y = rnorm(100)
 )  
 
# Create the 2D density contour plot
ggplot(df, aes(x = lambda, y = v)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_viridis_c() +
  labs(title = "2D Density Contour Plot", x = "X-Axis", y = "Y-Axis") +
  theme_minimal()
}

cat("Gamma Model \n")
gamma.mod <- R2jags::jags(model.file = textConnection(gamma.jags),
                          data = data_jags, inits = inits_list("gamma", n.chains),
                          n.chains = n.chains, parameters.to.save = c("lambda",
                                                                      "shape", "Like", "St_pred", "total_LLik"), n.iter = n.iter.jags,
                          n.thin = n.thin.jags, n.burnin = n.burnin.jags)

cat("LogNormal Model \n")
lnorm.mod <- R2jags::jags(model.file = textConnection(lnorm.jags),
                          data = data_jags, inits = inits_list("lnorm", n.chains),
                          n.chains = n.chains, parameters.to.save = c("mu", "sd",
                                                                      "Like", "St_pred", "total_LLik"), n.iter = n.iter.jags,
                          n.thin = n.thin.jags, n.burnin = n.burnin.jags)

cat("LogLogistic Model \n")
llogis.mod <- R2jags::jags(model.file = textConnection(llogis.jags),
                           data = data_jags_llogis, #inits = inits_list("llogis",n.chains), 
                           n.chains = n.chains, parameters.to.save = c("alpha","beta", "Like", "St_pred", "total_LLik","deriv1","deriv2"), n.iter = n.iter.jags,
                           n.thin = n.thin.jags, n.burnin = n.burnin.jags)


cat("Gompertz Model \n")
gomp.mod <- R2jags::jags(model.file = textConnection(gompertz.jags),
                         data = data_gomp, inits = inits_list("gompertz", n.chains),
                         n.chains = n.chains, parameters.to.save = c("a", "b",
                                                                     "Like", "St_pred", "total_LLik"), n.iter = n.iter.jags,
                         n.thin = n.thin.jags, n.burnin = n.burnin.jags)

cat("Generalized Gamma Model \n")
gen.gamma.mod <- R2jags::jags(model.file = textConnection(gen.gamma.jags),
                              data = data_jags, inits = inits_list("gengamma", n.chains),
                              n.chains = n.chains, parameters.to.save = c("r", "lambda",
                                                                          "b", "Like", "St_pred", "total_LLik"), n.iter = n.iter.jags,
                              n.thin = n.thin.jags, n.burnin = n.burnin.jags)


# Define parameters

a <- 0
b <- 1

# Generate random samples
dtrunc <-dtruncnorm(seq(0,1, by = 0.001), a = a, b = b, mean = mu_expert, sd = sd_expert)
# Estimate and plot the density

St_expert_expo <- expo.mod$BUGSoutput$sims.matrix[,c("St_pred")]
St_expert_weib <- weib.mod$BUGSoutput$sims.matrix[,c("St_pred")]
St_expert_gamma <- gamma.mod$BUGSoutput$sims.matrix[,c("St_pred")]
St_expert_lnorm <- lnorm.mod$BUGSoutput$sims.matrix[,c("St_pred")]
St_expert_llogis <- llogis.mod$BUGSoutput$sims.matrix[,c("St_pred")]
St_expert_gomp <- gomp.mod$BUGSoutput$sims.matrix[,c("St_pred")]
St_expert_gengam<- gen.gamma.mod$BUGSoutput$sims.matrix[,c("St_pred")]

plot(density(St_expert_llogis))
lines(y = dtrunc, x=seq(0,1, by = 0.001),  main = "Density of Truncated Normal Distribution",
      xlab = "Value", ylab = "Density", col = "red")


# Install and load the required package
library(numDeriv)

# Define a wrapper function for pgamma
pgamma_wrapper <- function(x, shape, rate, lower.tail = FALSE) {
  pgamma(x, shape = shape, rate = rate,lower.tail= lower.tail)
}

# Example parameters
shape <- 3
rate <- 2
x <- 1

# Load the rstan package
#install.packages("rstan")
library(rstan)

# Define the Stan model code
stan_model_code <- "
functions {
  real numerical_derivative(real x, real epsilon, real alpha, real beta, int param) {
    real f_plus;
    real f_minus;
    
    if(param==1){
    f_plus = gamma_cdf(x, alpha+ epsilon, beta);
    f_minus = gamma_cdf(x , alpha- epsilon, beta);
    }else{
        f_plus = gamma_cdf(x, alpha, beta + epsilon);
        f_minus = gamma_cdf(x , alpha, beta- epsilon);
    }
    
    return (f_plus - f_minus) / (2 * epsilon);
  }
}

data {
  real<lower=0> x; // Point at which to evaluate the derivative
  real<lower=0> alpha; // Shape parameter for gamma
  real<lower=0> beta; // Inverse scale parameter for gamma
  real<lower=0> epsilon;
  
}


generated quantities {
  real derivative1;
  real derivative2;
  derivative1 = numerical_derivative(x, epsilon, alpha, beta, 1);
  derivative2 = numerical_derivative(x, epsilon, alpha, beta, 2);
}
"

# Define the data
stan_data <- list(
  x = 1.0,
  alpha = 2.0,
  beta = 3.0,
  epsilon = 0.001
)

# Fit the model
fit <- stan(
  model_code = stan_model_code, 
  data = stan_data,
  chains = 1,
  iter = 1,
  warmup = 0,
  algorithm = "Fixed_param"
)

# Print the results
print(fit, pars = "derivative1")

as.matrix(fit)
# Install and load the required package
library(numDeriv)
# Example parameters
shape <- 2
rate <-3
x <- 1

# Compute the numerical derivative of pgamma with respect to x
grad(func = function(shape) pgamma_wrapper(x, shape, rate), x = shape)
derivative_pgamma <- grad(func = function(rate) pgamma_wrapper(x, shape, rate), x = rate)
grad(func = function(rate) pgamma_wrapper(x, shape, rate), x = rate)

gomp_surv <- expression(exp(-(b/a)*(exp(a*t) - 1)))

D(gomp_surv, "a")
D(gomp_surv, "b")

# Define the expressions
weibull_ph <- expression(exp(-lambda * t^a))
weibull_aft <- expression(exp(-(t/scale)^shape))
loglogistic <- expression(1 / (1 + (t / b)^a))

# Find the derivatives with respect to 'a'
derivative_weibull_ph <- D(weibull_ph, "a")
derivative_weibull_aft <- D(weibull_aft, "a")
derivative_loglogistic <- D(loglogistic, "a")

# Print the derivatives
print(derivative_weibull_ph)
print(derivative_weibull_aft)
print(derivative_loglogistic)


}


