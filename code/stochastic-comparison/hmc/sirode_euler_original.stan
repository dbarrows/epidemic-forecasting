data {

    int     <lower=1>   T;      // total integration steps
    real                y[T];   // observed number of cases
    int     <lower=1>   N;      // population size
    real                h;      // step size

}

parameters {

    real <lower=0, upper=10>    R0;     // R0
    real <lower=0, upper=10>    r;      // recovery rate
    real <lower=0, upper=20>    sigma;  // observation error
    real <lower=0, upper=500>   y0[3];  // initial conditions

}

model {

    real S[T];
    real I[T];
    real R[T];

    S[1] <- y0[1];
    I[1] <- y0[2];
    R[1] <- y0[3];
    
    y[1] ~ normal(y0[2], sigma);

    for (t in 2:T) {

        S[t] <- S[t-1] + h*( - S[t-1]*I[t-1]*R0*r/N );
        I[t] <- I[t-1] + h*( S[t-1]*I[t-1]*R0*r/N  - I[t-1]*r );
        R[t] <- R[t-1] + h*( I[t-1]*r );
        
        if (y[t] > 0) {
            y[t] ~ normal( I[t], sigma );
        }

    }
    
    y0[1] ~ normal(N - y[1], sigma);
    y0[2] ~ normal(y[1], sigma);
    
    R0      ~ lognormal(1,1);
    r       ~ lognormal(1,1);
    sigma   ~ lognormal(1,1);
        
}