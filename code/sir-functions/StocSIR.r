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

        S <- S + h*dS  #newInf
        I <- I + h*dI  #newInf - h*dR
        R <- R + h*dR  #h*dR

        if (i %% steps == 0)
            out[i/steps+1,] <- c(S,I,R,B)

    }

    return(out)

}