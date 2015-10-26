data {
    int     <lower=1>   T;      // number of time steps
    real    <lower=0>   y[T];   // observed number of cases at each time step, must be positive
    int     <lower=1>   N;      // population size
}


parameters {

    real    <lower=0>               r;
    real    <lower=0>               R0;
    real    <lower=0>               sigma;

}


transformed parameters {
}


model {

    vector[T] S;
    vector[T] I;
    vector[T] R;

    S[1] <- N - y[1];
    I[1] <- y[1];
    R[1] <- 0;

    y[1] ~ normal(I[1],sigma);    

    for (t in 2:T) {

        S[t] <- S[t-1] - R0*r/N*S[t-1]*I[t-1];
        I[t] <- I[t-1] + R0*r/N*S[t-1]*I[t-1] - r*I[t-1];
        R[t] <- R[t-1] + r*I[t-1];

        y[t] ~ normal( I[t], sigma );

    }

}