functions {

    row_vector sir( row_vector y0,
                real[] theta ) {

        real dydt[3];
        row_vector[3] y;

        real h;
        h <- 1.0/1.0;

        y[1] <- y0[1];
        y[2] <- y0[2];
        y[3] <- y0[3];

        for (t in 1:10) {

            dydt[1] <- - y[1]*y[2]*y[1]*y[2]/500;
            dydt[2] <- y[1]*y[2]*y[1]*y[2]/500 - theta[2]*y[2];
            dydt[3] <- y[2]*y[2];

            y[1] <- y[1] + h*( dydt[1] );
            y[2] <- y[2] + h*( dydt[2] );
            y[3] <- y[3] + h*( dydt[3] );


        }

        return(y);

    }
}

data {

    int     <lower=1>   T;        // total integration time
    real    <lower=0>   y[T+1];   // observed number of cases at each time step
    int     <lower=1>   N;        // population size

}

parameters {

    real <lower=0,upper=7> sigma;
    real <lower=0,upper=7> theta [2];
    real <lower=0,upper=500> y0[3];

}


model {

    row_vector[3] Y[T];
    row_vector[3] init;

    y[1] ~ normal( y0[2], sigma );

    theta[1]    ~ lognormal(1, 1);
    theta[2]    ~ lognormal(1, 1);
    sigma       ~ lognormal(1, 1);
    
    y0[1] ~ normal(N - y[1], sigma);
    y0[2] ~ normal(y[1], sigma);

    init[1] <- y0[1];
    init[2] <- y0[2];
    init[3] <- y0[3];

    Y[1] <- sir(init, theta);

    for (t in 2:T) {

        Y[t] <- sir(Y[t-1], theta);

        y[t+1] ~ normal( Y[t], sigma );

    }

}