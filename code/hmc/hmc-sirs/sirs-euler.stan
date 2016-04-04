## Dexter Barrows
## dbarrows.github.io
## McMaster University
## 2016

data {

    int     <lower=1>   T;      // total integration steps
    real                y[T];   // observed number of cases
    int     <lower=1>   N;      // population size
    real                h;      // step size

}

parameters {

    real <lower=0, upper=10>        R0;     // R0
    real <lower=0, upper=10>        r;      // recovery rate
    real <lower=0, upper=10>        re;     // resusceptibility rate
    real <lower=0, upper=20>        sigma;  // observation error
    real <lower=0, upper=30>        Iinit;    // initial infected
    real <lower=0, upper=1>         eta;    // geometric walk attraction strength
    real <lower=0, upper=1>         berr;   // beta walk noise
    real <lower=-1.5, upper=1.5>    Bnoise[T];   // Beta vector

}

model {

    real S[T];
    real I[T];
    real R[T];
    real B[T];
    real B0;

    real pi;
    real Bfac;

    pi <- 3.1415926535;

    B0 <- R0 * r / N;

    B[1] <- B0;

    S[1] <- N - Iinit;
    I[1] <- Iinit;
    R[1] <- 0.0;

    for (t in 2:T) {

        Bnoise[t] ~ normal(0,berr);
        Bfac <- exp(2*cos((2*pi/365)*t) - 2);
        B[t] <- exp( log(B0) + eta * ( log(B[t-1]) - log(B0) ) + Bnoise[t] );

        S[t] <- S[t-1] + h*( - Bfac*B[t]*S[t-1]*I[t-1] + re*R[t-1] );
        I[t] <- I[t-1] + h*( Bfac*B[t]*S[t-1]*I[t-1] - I[t-1]*r );
        R[t] <- R[t-1] + h*( I[t-1]*r - re*R[t-1] );
        
        if (y[t] > 0) {
            y[t] ~ normal( I[t], sigma );
        }

    }
    
    R0      ~ lognormal(1,1);
    r       ~ lognormal(1,1);
    sigma   ~ lognormal(1,1);
    re      ~ lognormal(1,1);
    Iinit   ~ normal(y[1], sigma);
        
}