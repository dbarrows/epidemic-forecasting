library(rstan)
setwd("~/Google Drive/courses/Bayesian/Bayesian software/rstan")

rt <- stanc(file = "schools.stan", model_name = '8schools')

sm <- stan_model(stanc_ret = rt, verbose = TRUE)

schools_data <- list(J = 8, 
                     y = c(28, 8, -3, 7, -1, 1, 18, 12), 
                     sigma = c(15, 10, 16, 11, 9, 11, 10, 18))

fit <- sampling(sm, data=schools_data, chains=4)
fit


