library("flexsurv")
library("rjags")
library("dplyr")
n <- 2000
times <- rgenf(n,  1.3, 2, 0.5, 0.1)
df <- data.frame(time_event = times)

censoring_time <-5
df <- mutate(df, time = ifelse(time_event > censoring_time, censoring_time,time_event))
df <- mutate(df, status = ifelse(time_event > censoring_time, 0, 1))

obs_events <- df[which(df$status==1),"time"]
obs_censors <- df[which(df$status==0),"time"]
km.df <- survfit(Surv(time, status)~1,data = df)

plot(km.df)


data_new <- list()
df_jags <- df[,c("time","status")]
df_jags$t <- df$time

tinits1<-df_jags$t + max(df$time)
is.na(tinits1)<-df_jags$status==1
tinits2<-tinits1 + 5

is.na(df_jags$t)<-df_jags$status==0
df_jags$is.censored<-1-df_jags$status
df_jags$t.cen<-df_jags$time+df_jags$status


modelinits <- list(list(t = tinits1),
                   list(t = tinits2))
data_jags <- list(N = nrow(df_jags),
                  t.cen = df_jags$t.cen,
                  is.censored = df_jags$is.censored,
                  t = df_jags$t)


mu <- 1
sigma <- 2
Q <- 3
P <- 0.5
x <- time <- 2

pgenf(time, mu, sigma, Q,P, lower.tail = F)

tmp = Q * Q + 2 * P
delta = sqrt(tmp)
s1 = 2 / (tmp + Q*delta)
s2 = 2 / (tmp - Q*delta)
#expw = std::pow(x, delta/sigma) * std::exp(-mu*delta/sigma)
expw <- x^(delta/sigma)*exp(-mu*delta/sigma)


ldensity_calc <- log(delta) + s1/sigma*delta*(log(x) - mu) +
  s1*(log(s1) - log(s2)) - log(sigma*x)-
  (s1+s2)*log(1 + s1*expw/s2) - lgamma(s1)- lgamma(s2)+ lgamma(s1+s2)
dgenf(time, mu, sigma, Q,P, log = T)
qb <- c()
qb[1] = s2/(s2 + s1*expw)
qb[2] = s1*expw/(s2 + s1*expw)
qb_indicator <- ifelse(qb[1] > 0.99, 2, 1)

Surv_calc <- ifelse(qb_indicator == 1, pbeta(qb[1],s2,s1,lower.tail = T,log.p = F),
                    pbeta(qb[2],s1,s2,lower.tail = F,log.p = F))

pgenf(time, mu, sigma, Q,P, lower.tail = F)

genf.jags <- "
data{
for(i in 1:N){
zero[i] <- 0}
}

model{

C <- 10000
for(i in 1:N){

logdens[i] <- (log(delta) + (s1/sigma)*delta*(log(time[i]) - mu) + 
              s1*(log(s1) - log(s2)) - log(sigma*time[i])-
              (s1+s2)*log(1 + s1*expw[i]/s2) - lbeta)

logSurv[i] <- log(Surv[i])


zero[i] ~ dpois(zero.mean[i])
zero.mean[i] <- -logdens[i]*(status[i])-logSurv[i]*(1-status[i]) + C


expw[i] <- pow(time[i],(delta/sigma))*exp(-mu*delta/sigma)
Surv[i] <- pbeta(qb[i],s2,s1)

qb[i] = s2/(s2 + s1*expw[i])

#For some reason the overflow protection in flexsurvreg gives the wrong vals 

}

#Calculate the parameters
tmp <- Q * Q + 2 * P
delta <- sqrt(tmp)
s1 <- 2 / (tmp + Q*delta)
s2 <- 2 / (tmp - Q*delta)

#Need Log Beta lbeta(s1, s2)
lbeta <- loggam(s1) + loggam(s2) - loggam(s1+s2)

#mu ~ dunif(-5,5)
#sigma ~ dunif(0.2,5)
#Q ~ dunif(0.2,5)
#P ~ dunif(0.2,5)

mu ~ dnt(prior_mean[1], prior_tau[1], 1)
sigma ~ dnt(prior_mean[2], prior_tau[2], 1)T(0.2,)
Q ~ dnt(prior_mean[3], prior_tau[3], 1)
P ~ dnt(prior_mean[4], prior_tau[4], 1)T(0.2,)
}"

genf.mle <- flexsurvreg(Surv(time,status)~1,data = df, dist = "genf")
genf.mle.res <- data.frame(genf.mle$res[,c(1,4)])
genf.mle.res$tau <- as.vector(sapply(genf.mle.res[,c(2)],function(x){(2*x)^(-2)}))

data_surv <- list()
data_surv$time <- df$time
data_surv$status <- df$status
data_surv$N <- nrow(df)

data_surv$prior_mean <- genf.mle.res$est
data_surv$prior_tau <- genf.mle.res$tau



genf.mod <-R2jags::jags(model.file = textConnection(genf.jags),
                        data=data_surv,
                        #inits=list(list(a = 1, b = 1), list(a = 2,b = 2)),
                        n.chains=4,
                        parameters.to.save = c("mu","sigma","Q","P", "Surv", "logdens"),
                        #parameters.to.save = c("mu","sigma","Q","P","logdens", "Surv"),
                        n.iter = 2000,
                        n.thin = 5,
                        n.burnin = 10)


genf.mod$BUGSoutput$summary[rownames(genf.mod$BUGSoutput$summary) %in% c("mu", "sigma", "Q", "P"),]

list.modules()

time_seq <- seq(0,5, by = 0.01)
surv.array <- apply(genf.mod$BUGSoutput$sims.matrix[,c("mu", "sigma", "Q", "P")], 1, function(x){pgenf(time_seq, x[1],x[2], x[3], x[4], lower.tail = F)})

plot(genf.mle)
lines(time_seq,rowMeans(surv.array))
lines(time_seq,apply(surv.array, 1,quantile, prob = c(0.025)))
lines(time_seq,apply(surv.array, 1,quantile, prob = c(0.975)))