library(rstan)
library(shinystan)

rmod = stanc(file = "schools.stan", model_name = "school_model")

smod = stan_model(stanc_ret = rmod, verbose = TRUE)

schools_data <- list(J = 8,
                     y = c(28, 8, -3, 7, -1, 1, 18, 12),
                     sigma = c(15, 10, 16, 11, 9, 11, 10, 18))

fit <- sampling (smod, data = schools_data, chains = 100)

sso <- as.shinystan(fit)
sso <- launch_shinystan(sso)