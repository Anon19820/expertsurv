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

#alpha1 ~ dunif(-100,100)
#alpha2 ~ dunif(0,100)
alpha1 ~ dnorm(mu_beta[1],1/(sigma_beta[1])^2)
alpha2 ~ dgamma(a_alpha, b_alpha)
alpha <-   alpha1*ifelse(St_indic == 1, 1,0) +alpha2*ifelse(St_indic == 0, 1,0)


for(i in 1:H){
 prec_beta[i] <- 1/(sigma_beta[i])^2
 beta[i] ~ dnorm(mu_beta[i],prec_beta[i])
 #beta[i] ~ dunif(-100,100)
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
 
  phi_con[i] <- -log(sum(phi_temp[,i]))*equals(pool_type,1) +  -log(prod(phi_temp[,i]))*equals(pool_type,0) + C 

 } 
 
  #deriv1 <- abs(exp(-(mu[id_St]/alpha)*(exp(alpha*time_expert[1]) - 1)) * (mu[id_St]/(alpha^2)*(exp(alpha*time_expert[1])-1) - (mu[id_St]/alpha)*(exp(alpha*time_expert[1]) * time_expert[1])))
  #deriv2 <- abs(-(exp(-(mu[id_St]/alpha)*(exp(alpha*time_expert[1]) - 1))*((1/alpha)*(exp(alpha*time_expert[1]) - 1))))
  #zero.mean_prior <- -log(deriv1+deriv2) + C
  #zero_prior ~ dpois(zero.mean_prior)
 
 
rate = exp(beta[1])

s <- 0.0001

}"

Gamma.jags <- "
data{
for(i in 1:n){
    zeros[i] <- 0
    d2[i] <- ifelse(d[i] == 1, 0,1)
  }
for(i in 1:n_time_expert){
    zero2[i] <- 0
}
 
  zero_prior <- 0 
  C <- 1000 
  epsilon <- 0.0001
 
}

model{

for(i in 1:H){
 prec_beta[i] <- 1/(sigma_beta[i])^2
 beta[i] ~ dnorm(mu_beta[i],prec_beta[i])
 #beta[i] ~ dunif(-100,100)
}

linpred <- X %*%beta
for( i in 1:n){
lambda[i] <- exp(linpred[i])
}

for(i in 1:n){

log_d[i]  <- log(dgamma(t[i],alpha,lambda[i]))
log_S[i]  <- log(1-pgamma(t[i],alpha, lambda[i]))
LL[i] <- log_d[i]*d[i] + log_S[i]*d2[i] + log(a0[i])
zeros[i] ~ dpois(zero.mean[i])
zero.mean[i] <- -LL[i] + C

}

for (i in 1:n_time_expert){
  zero2[i] ~ dpois(phi_con[i])

    St_expert_temp[i,1] =  (1-pgamma(time_expert[i],alpha,lambda[id_St]))*St_indic
    mean_surv_trt[i] <- alpha/lambda[id_trt]
    mean_surv_comp[i] <- alpha/lambda[id_comp]
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
 
 phi_con[i] <- -log(sum(phi_temp[,i]))*equals(pool_type,1) +  -log(prod(phi_temp[,i]))*equals(pool_type,0) + C 

 } 

  #deriv1 <- ((1- pgamma(time_expert[1],alpha+epsilon,lambda[id_St])) - (1 - pgamma(time_expert[1],alpha-epsilon,lambda[id_St])))/(2*epsilon)
  #deriv2 <- ((1- pgamma(time_expert[1],alpha,lambda[id_St]+epsilon)) - (1 - pgamma(time_expert[1],alpha,lambda[id_St]-epsilon)))/(2*epsilon)
  #zero.mean_prior <- -log(abs(deriv1)+abs(deriv2)) + C
  #zero_prior ~ dpois(zero.mean_prior)

  alpha ~ dgamma(a_alpha,b_alpha);
  #alpha ~ dunif(0,100);
  rate <- exp(beta[1]);


}"
GenGamma.jags <- "

data{

    for(i in 1:n){
    zero[i] <- 0
    }
    
  zero_prior <- 0 
  C <- 1000 
  epsilon <- 0.0001
}


model{


for(i in 1:H){
 prec_beta[i] <- 1/(sigma_beta[i])^2
 beta_jags[i] ~ dnorm(mu_beta[i],prec_beta[i])
  #beta_jags[i] ~ dunif(-100,100)

 
}

linpred <- X %*%beta_jags

for( i in 1:n){
lambda[i] <- exp(linpred[i])
}


for(i in 1:n){
is.censored[i]~dinterval(t_jags[i],t_cen[i])
t_jags[i] ~ dgen.gamma(r,lambda[i],b)

}

for (i in 1:n_time_expert){
  zero[i] ~ dpois(phi_con[i])

    St_expert_temp[i,1] =  (1-pgen.gamma(time_expert[i],r,lambda[id_St],b))*St_indic
    mean_surv_trt[i] <- exp(loggam((b*r +1)/b) - loggam(r))/lambda[id_comp]
    mean_surv_comp[i] <- exp(loggam((b*r +1)/b) - loggam(r))/lambda[id_trt]
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
  
  phi_con[i] <- -log(sum(phi_temp[,i]))*equals(pool_type,1) +  -log(prod(phi_temp[,i]))*equals(pool_type,0) + C 
 } 


r ~ dgamma(a_alpha,b_alpha);
b ~ dgamma(a_alpha,b_alpha);
#r ~ dunif(0,100);
#b ~ dunif(0,100);

sigma <- 1/(b*pow(r,0.5))
Q <- pow(r,-0.5)
mu <- -beta_jags[1] + (log(r)/b)

beta[1] <- mu

for(i in 2:H){
beta[i] <- beta_jags[i]

}
  #deriv1 <- ((1-pgen.gamma(time_expert[1],r+epsilon,lambda[id_St],b)) - (1-pgen.gamma(time_expert[1],r-epsilon,lambda[id_St],b)))/(2*epsilon)
  #deriv2 <- ((1-pgen.gamma(time_expert[1],r,lambda[id_St]+epsilon,b)) - (1-pgen.gamma(time_expert[1],r,lambda[id_St]-epsilon,b)))/(2*epsilon)
  #deriv3 <- ((1-pgen.gamma(time_expert[1],r,lambda[id_St],b+epsilon)) - (1-pgen.gamma(time_expert[1],r,lambda[id_St],b-epsilon)))/(2*epsilon)
  #zero.mean_prior <- -log(abs(deriv1)+abs(deriv2)+abs(deriv3)) + C
  #zero_prior ~ dpois(zero.mean_prior)


}"
