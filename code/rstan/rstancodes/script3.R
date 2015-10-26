library(rstan)
setwd("~/Google Drive/courses/Bayesian/Bayesian software/rstan")

schools_data <- list(J = 8, 
                     y = c(28, 8, -3, 7, -1, 1, 18, 12), 
                     sigma = c(15, 10, 16, 11, 9, 11, 10, 18))

fit <- stan(file="schools.stan", data=schools_data,iter=100, chains=4)
fit



fit2 <- stan(file="schools.stan", data=schools_data,iter=100, chains=4, 
             warmup = 20, thin=2)
fit2

print(fit2, pars=c("theta", "mu", "tau", "lp__"),probs=c(.1,.5,.9))

summary(fit2)

plot(fit2)

traceplot(fit2)

extract(fit2, permuted = FALSE) 
