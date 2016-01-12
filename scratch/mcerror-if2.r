library(foreach)
library(parallel)
library(doParallel)
library(Rcpp)

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


set.seed(1004)

T       <- 60
i_infec <- 5
steps   <- 7
N       <- 500
sigma   <- 10
Tlim    <- T

## Generate true trajecory and synthetic data
##

pars_true <- c(R0 = 3.0,    # new infected people per infected person
          r = 0.1,      # recovery rate
          N = 500,      # population size
          eta = 0.5,    # geometric random walk
          berr = 0.5)   # Beta geometric walk noise

true_init_cond <- c(S = N - i_infec,
                    I = i_infec,
                    R = 0)

# setup cluster
numCores <- detectCores()
numCores
cl <- makeCluster(numCores)
registerDoParallel(cl)


## Re-run PF stuff with target number of particles
##

NP          <- 3000
nPasses     <- 50
coolrate    <- 0.975

if2file <- paste(getwd(),"if2-d.cpp",sep="/")

#if2vals <- numeric(10)
#if2res <- numeric(10)

if2mat <- foreach (i = 1:10, .combine = rbind, .packages = "Rcpp") %dopar% {

    sdeout_true <- StocSIR(true_init_cond, pars_true, T, steps)
    colnames(sdeout_true) <- c('S','I','R','B')

    infec_counts_raw <- sdeout_true[,'I'] + rnorm(T+1, 0, sigma)
    infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

    sourceCpp(if2file)
    invisible( if2time <- system.time( if2data <- if2(infec_counts[1:(Tlim+1)],
                                                        Tlim+1,
                                                        N,
                                                        NP,
                                                        nPasses,
                                                        coolrate)
                                    )
                )

    paramdata <- data.frame( if2data$paramdata )
    names(paramdata) <- c("R0", "r", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")

    parnames <- c("R0","r","sigma","eta","berr","Iinit")
    parvars <- var(paramdata[parnames])
    parmeans <- colMeans(paramdata[parnames])

    res <- t(1 / as.matrix(parmeans)) %*% parvars %*% (1 / as.matrix(parmeans))
    val <- sqrt( res / ( (NP) * length(parmeans)) )

    #if2vals[i] <- val
    #if2res[i] <- res
    return( c(val,res) )

}

colnames(if2mat) <- c("vals","res")
if2vals <- if2mat[,'vals']
if2res  <- if2mat[,'res']

if2_err <- mean(if2vals)
print(if2_err)

save.image(file = "mcerror-if2.RData")