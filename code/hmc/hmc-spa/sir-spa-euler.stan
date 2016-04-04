## Dexter Barrows
## dbarrows.github.io
## McMaster University
## 2016

data {

    int     <lower=1>   T;      // total integration steps
    int     <lower=1>   nloc;   // number of locations
    real                y[nloc, T];   // observed number of cases
    int     <lower=1>   N;      // population size
    real                h;      // step size
    int     <lower=0>   neinum[nloc];        // number of neighbors each location has
    int                 neibmat[nloc, nloc]; // neighbor list for each location

}

parameters {

    real <lower=0, upper=10>        R0;     // R0
    real <lower=0, upper=10>        r;      // recovery rate
    real <lower=0, upper=20>        sigma;  // observation error
    real <lower=0, upper=30>        Iinit[nloc];    // initial infected for each location
    real <lower=0, upper=1>         eta;    // geometric walk attraction strength
    real <lower=0, upper=1>         berr;   // beta walk noise
    real <lower=-1.5, upper=1.5>    Bnoise[nloc,T];   // Beta vector
    real <lower=0, upper=1>         phi;    // interconnectivity strength

}

model {

    real S[nloc, T];
    real I[nloc, T];
    real R[nloc, T];
    real B[nloc, T];
    real B0;

    real BSI[nloc, T];
    real rI[nloc, T];
    int  n;
    real sphi;
    real ophi;
    real nBIsum;

    B0 <- R0 * r / N;

    for (loc in 1:nloc) {
        S[loc, 1] <- N - Iinit[loc];
        I[loc, 1] <- Iinit[loc];
        R[loc, 1] <- 0.0;
        B[loc, 1] <- B0;
    }

    for (t in 2:T) {
        for (loc in 1:nloc) {

            Bnoise[loc, t] ~ normal(0,berr);
            B[loc, t] <- exp( log(B[loc, t-1]) + eta * ( log(B0) - log(B[loc, t-1]) ) + Bnoise[loc, t] );

            n <- neinum[loc];
            sphi <- 1.0 - phi*( n/(n+1.0) );
            ophi <- phi/(n+1.0);

            nBIsum <- 0.0;
            for (j in 1:n)
                nBIsum <- nBIsum + B[neibmat[loc, j], t-1] * I[neibmat[loc, j], t-1];

            BSI[loc, t] <- S[loc, t-1]*( sphi*B[loc, t-1]*I[loc, t-1] + ophi*nBIsum );
            rI[loc, t]  <- r*I[loc, t-1];

            S[loc, t] <- S[loc, t-1] + h*( - BSI[loc, t] );
            I[loc, t] <- I[loc, t-1] + h*( BSI[loc, t] - rI[loc, t] );
            R[loc, t] <- R[loc, t-1] + h*( rI[loc, t] );
            
            if (y[loc, t] > 0) {
                y[loc, t] ~ normal( I[loc, t], sigma );
            }

        }
    }
    
    R0      ~ lognormal(1,1);
    r       ~ lognormal(1,1);
    sigma   ~ lognormal(1,1);
    for (loc in 1:nloc) {
        Iinit[loc] ~ normal(y[loc, 1], sigma);
    }

}