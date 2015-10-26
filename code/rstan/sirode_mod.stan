functions {

    real[] sir( real t,
                real[] y,
                real[] theta,
                real[] x_r,
                int[] x_i) {

        real dydt[3];

        dydt[1] <- - theta[1]*theta[2]*y[1]*y[2]/x_i[1];
        dydt[2] <- theta[1]*theta[2]*y[1]*y[2]/x_i[1] - theta[2]*y[2];
        dydt[3] <- theta[2]*y[2];

        return dydt;

    }
}

data {

    int     <lower=1>   T;        // total integration time
    real    <lower=0>   y[T+1];   // observed number of cases at each time step, must be positive
    int     <lower=1>   N;        // population size
    real    <lower=0>   ts[2];    // integration sample times

}

transformed data {

    real    x_r [0];
    int     x_i [1];
    x_i[1] <- N;

}

parameters {

    real <lower=0,upper=20> sigma;
    real <lower=0,upper=50> theta [2];
    real <lower=0,upper=500> y0[3];

}


model {

    real yhat[2,3];
    real yhat_old[3];

    real S[T];
    real I[T];
    real R[T];

    theta[1]    ~ lognormal(1, 1);
    theta[2]    ~ lognormal(1, 1);
    sigma       ~ lognormal(1, 1);

    yhat <- integrate_ode(sir, y0, 0, ts, theta, x_r, x_i);

    S[1] <- yhat[2,1];
    I[1] <- yhat[2,2];
    R[1] <- yhat[2,3];

    y[1] ~ normal( y0[2], sigma );

    for (t in 2:T) {

        yhat_old[1] <- S[t-1];
        yhat_old[2] <- I[t-1];
        yhat_old[3] <- R[t-1];

        yhat <- integrate_ode(sir, yhat_old, 0, ts, theta, x_r, x_i);

        S[t] <- yhat[2,1];
        I[t] <- yhat[2,2];
        R[t] <- yhat[2,3];

        y[t] ~ normal( I[t], sigma );

    }

}