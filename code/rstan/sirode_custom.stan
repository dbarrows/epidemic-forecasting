functions {

    real[,] sir( real TTT,
                real[] y0,
                real[] theta,
                real[] x_r,
                int[] x_i) {

        real dydt[3];
        real y[3];
        real yout[x_i[3],3];

        real h;
        h <- 1.0/x_i[2];

        y[1] <- y0[1];
        y[2] <- y0[2];
        y[3] <- y0[3];

        for (t in 1:x_i[3]) {

            for (tt in 1:x_i[2]) {

                dydt[1] <- - y[1]*y[2]*y[1]*y[2]/x_i[1];
                dydt[2] <- y[1]*y[2]*y[1]*y[2]/x_i[1] - theta[2]*y[2];
                dydt[3] <- y[2]*y[2];

                y[1] <- y[1] + h*( dydt[1] );
                y[2] <- y[2] + h*( dydt[2] );
                y[3] <- y[3] + h*( dydt[3] );

            }

            yout[t,1] <- y[1];
            yout[t,2] <- y[2];
            yout[t,3] <- y[3];

        }

        return(yout);

    }
}

data {

    int     <lower=1>   T;        // total integration time
    real    <lower=0>   y[T+1];   // observed number of cases at each time step, must be positive
    int     <lower=1>   N;        // population size
    int     <lower=1>   nsteps;   // number of euler steps to take per time unit (2, 10, 100 ...)
    #real    <lower=0>   ts[2];    // integration sample times

}

transformed data {

    real    x_r [1];
    int     x_i [3];
    x_i[1] <- N;
    x_i[2] <- nsteps;
    x_i[3] <- T;

}

parameters {

    real <lower=0,upper=20> sigma;
    real <lower=0,upper=50> theta [2];
    real <lower=0,upper=500> y0[3];

}


model {

    real yhat[T,3];

    theta[1]    ~ lognormal(1, 1);
    theta[2]    ~ lognormal(1, 1);
    sigma       ~ lognormal(1, 1);
    
    y0[1] ~ normal(N - y[1], sigma);
    y0[2] ~ normal(y[1], sigma);

    yhat <- sir(0, y0, theta, x_r, x_i);

    y[1] ~ normal( y0[2], sigma );

    for (t in 1:(T-1)) {

        y[t+1] ~ normal( yhat[t,2], sigma );

    }

}