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
    real    <lower=0>   ts[T];    // integration sample times

}

transformed data {

    real    x_r [0];
    int     x_i [1];
    x_i[1] <- N;

}

parameters {

    real <lower=0,upper=10> sigma;
    real <lower=0,upper=10> theta [2];
    real <lower=0,upper=500> y0[3];

}


model {

    real yhat[T,3];

    y0[1] ~ normal(N - y[1], sigma);
    y0[2] ~ normal(y[1], sigma);

    theta[1]    ~ lognormal(1,1);
    theta[2]    ~ lognormal(1,1);
    sigma       ~ lognormal(1,1);

    yhat <- integrate_ode(sir, y0, 0, ts, theta, x_r, x_i);

    y[1] ~ normal( y0[2], sigma );

    for (t in 1:T)
        y[t+1] ~ normal( yhat[t,2], sigma );

}