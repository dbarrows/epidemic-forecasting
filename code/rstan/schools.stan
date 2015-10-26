data {
    int<lower=0> J;             // number of schools
    real y[J];                  // estimated treatment effects: y is a vector of integers with length J
    real<lower=0> sigma[J];     // s.e. of effect estimates: sigma is a vector of integers with length J
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
    eta ~ normal(0, 1);         // vectorization allowed
    y ~ normal(theta, sigma);   // vectorization allowed
}