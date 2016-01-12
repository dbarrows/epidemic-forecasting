# Dexter Barrows
# run with:
#
# Rscript sc1-multi.r > sc1-multi.log 2>&1
#
# or
#
# Rscript sc1-multi.r > sc1-multi-time.log 2>&1


library(reshape2)
library(foreach)
library(parallel)
library(doParallel)
library(rstan)

## Stochastic SIR function
##############################################################################################

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

## Setup
##############################################################################################

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

nTrials <- 10

numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)

## external files

if2file <- paste(getwd(),"../../../code/stochastic-comparison/if2/if2-d.cpp",sep="/")
hmcfile <- paste(getwd(), "../../../code/stochastic-comparison/hmc", "sirode_euler.stan", sep="/")

## IF2 settings

NP          <- 3000
nPasses     <- 50
coolrate    <- 0.975


## fit once to compile
##

# get raw data
sdeout_true <- StocSIR(true_init_cond, pars_true, T, steps)
colnames(sdeout_true) <- c('S','I','R','B')

infec_counts_raw <- sdeout_true[,'I'] + rnorm(T+1, 0, sigma)
infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

# conform data to stan's required format
datlen <- T*7 + 1

data <- matrix(data = -1, nrow = T+1, ncol = steps)
data[,1] <- infec_counts
standata <- as.vector(t(data))[1:datlen]

#options
rstan_options(auto_write = TRUE)

stan_options <- list(   chains = 1,         # number of chains
                        iter   = 3000,      # iterations per chain
                        warmup = 1000,      # warmup interations
                        thin   = 5)         # thinning number

# data
sir_data <- list( T = datlen,       # simulation time
                  y = standata,     # infection count data
                  N = 500,          # population size
                  h = 1/steps )     # step size per day

# fit!
initialfit <- with(stan_options,
                stan(file   = hmcfile,
                    data    = sir_data,
                    chains  = chains,
                    iter    = iter,
                    warmup  = warmup,
                    thin    = thin)
                )


## Run multiple trajectories
#########################################################################################################

estmat <- foreach ( i = 1:nTrials, .combine = rbind, .packages = c("Rcpp","rstan","reshape2") ) %dopar% {

    sdeout_true <- StocSIR(true_init_cond, pars_true, T, steps)
    colnames(sdeout_true) <- c('S','I','R','B')

    infec_counts_raw <- sdeout_true[,'I'] + rnorm(T+1, 0, sigma)
    infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

    # IF2

    sourceCpp(if2file)
    if2time <- system.time( if2data <- if2(infec_counts[1:(Tlim+1)], Tlim+1, N, NP, nPasses, coolrate) )

    paramdata <- data.frame( if2data$paramdata )
    names(paramdata) <- c("R0", "r", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")

    parnames <- c("R0","r","sigma","eta","berr","Iinit")
    parvars <- var(paramdata[parnames])
    if2parmeans <- colMeans(paramdata[parnames])
    names(if2parmeans) <- c("if2.R0","if2.r","if2.sigma","if2.eta","if2.berr","if2.Iinit")


    ## HMC fitting with stan

    datlen <- T*7 + 1

    data <- matrix(data = -1, nrow = T+1, ncol = steps)
    data[,1] <- infec_counts
    standata <- as.vector(t(data))[1:datlen]

    sir_data <- list( T = datlen,       # simulation time
                      y = standata,     # infection count data
                      N = 500,          # population size
                      h = 1/steps )     # step size per day

    hmctime <- system.time(fit <- with(stan_options,
                            stan(#file      = hmcfile,
                                 fit = initialfit,
                                data    = sir_data,
                                chains  = chains,
                                iter    = iter,
                                warmup  = warmup,
                                thin    = thin)
                            )
                )

    exfit <- extract(fit, permuted = FALSE, inc_warmup = FALSE)

    paramdata <- data.frame(R0 = melt(exfit[,,'R0'])$value,
                            r = melt(exfit[,,'r'])$value,
                            sigma = melt(exfit[,,'sigma'])$value,
                            eta = melt(exfit[,,'eta'])$value,
                            berr = melt(exfit[,,'berr'])$value,
                            Sinit = melt(exfit[,,'y0[1]'])$value,
                            Iinit = melt(exfit[,,'y0[2]'])$value,
                            Rinit = melt(exfit[,,'y0[3]'])$value )

    parnames <- c("R0","r","sigma","eta","berr","Iinit")
    parvars <- var(paramdata[parnames])
    hmcparmeans <- colMeans(paramdata[parnames])
    names(hmcparmeans) <- c("hmc.R0","hmc.r","hmc.sigma","hmc.eta","hmc.berr","hmc.Iinit")

    return( c(if2parmeans, if2times = if2time[['user.self']], hmcparmeans, hmctimes = hmctime[['user.self']]) )

}


stopCluster(cl)

save.image(file = "sc1-multi-time.RData")