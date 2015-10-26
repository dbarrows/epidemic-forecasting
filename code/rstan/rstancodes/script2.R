library(rstan)

modelString = "data {
  int<lower=0> J; // number of schools
  real y[J]; // estimated treatment effects: 
             //  y is a vector of integers with length J
  real<lower=0> sigma[J]; // s.e. of effect estimates: 
                          // sigma is a vector of integers with length J
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[J] eta;
}
transformed parameters {
  vector[J] theta;
  theta <- mu + tau * eta;
}
model {
  eta ~ normal(0, 1); // vectorization allowed
  y ~ normal(theta, sigma); // vectorization allowed
}
"

sm <- stan_model(model_code = modelString)

schools_data <- list(J = 8, 
                     y = c(28, 8, -3, 7, -1, 1, 18, 12), 
                     sigma = c(15, 10, 16, 11, 9, 11, 10, 18))
fit <- sampling(sm, data=schools_data, chains=4)
fit