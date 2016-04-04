## Dexter Barrows
## dbarrows.github.com
## McMaster University
## 2016

StocSIR <- function(y, pars, T, steps) {

    out <- matrix(NA, nrow = (T+1), ncol = 4)

    R0 <- pars[['R0']]
    r <- pars[['r']]
    N <- pars[['N']]
    eta <- pars[['eta']]
    berr <- pars[['berr']]

    S <- y[['S']]
    I <- y[['I']]
    R <- y[['R']]

    B0 <- R0 * r / N
    B <- B0

    out[1,] <- c(S,I,R,B)

    h <- 1 / steps

    for ( i in 1:(T*steps) ) {

        B <- exp( log(B) + eta*(log(B0) - log(B)) + rnorm(1, 0, berr) )

        BSI <- B*S*I
        rI <- r*I

        dS <- -BSI
        dI <- BSI - rI
        dR <- rI

        S <- S + h*dS
        I <- I + h*dI
        R <- R + h*dR

        if (i %% steps == 0)
            out[i/steps+1,] <- c(S,I,R,B)

    }

    return(out)

}

### Suggested parameters
#
# T       <- 60
# i_infec <- 5
# steps   <- 7
# N       <- 500
# sigma   <- 10
#
# pars <- c(R0 = 3.0,     # new infected people per infected person
#           r = 0.1,      # recovery rate
#           N = 500,      # population size
#           eta = 0.5,    # geometric random walk
#           berr = 0.5)   # Beta geometric walk noise