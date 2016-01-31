# Dexter Barrows
#
# run: Rscript if2-paraboot.r > if2-paraboot.log 2>&1

library(foreach)
library(parallel)
library(doParallel)
library(ggplot2)
library(Rcpp)

source(paste(getwd(),"../../sir-functions", "StocSIR.r", sep = "/"))

set.seed(112358)


## data params
T       <- 60
i_infec <- 5
steps   <- 7
N       <- 500
sigma   <- 10
Tlim    <- round(T / 2)

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

# get true trajectory
sdeout_true <- StocSIR(true_init_cond, pars_true, T, steps)
colnames(sdeout_true) <- c('S','I','R','B')

# get initial data
infec_counts_raw <- sdeout_true[,'I'] + rnorm(T+1, 0, sigma)
infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

## parallel setup
##

nTrials <- 100

numCores <- detectCores()
cl <- makeCluster(min(numCores, nTrials))
registerDoParallel(cl)

## if2 setup
##

NP          <- 3000
nPasses     <- 50
coolrate    <- 0.975

## initial fit
##
if2file <- paste(getwd(),"../../if2", "if2.cpp", sep="/")

#if (FALSE) {
sourceCpp(if2file)
if2time <- system.time( if2data <- if2(infec_counts[1:(Tlim+1)], Tlim+1, N, NP, nPasses, coolrate) )

# get parameter estimates

paramdata <- data.frame( if2data$paramdata )
names(paramdata) <- c("R0", "r", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")
parmeans <- colMeans(paramdata)
names(parmeans) <- c("R0", "r", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")

# get fit state means
fitdatafull <- data.frame(if2data$statemeans)
names(fitdatafull) <- c("S","I","R")
fitdatainfec <- fitdatafull$I
#}


## use parametric bootstrapping to generate forcasts
##
trajectories <- foreach( i = 1:nTrials, .combine = rbind, .packages = "Rcpp") %dopar% {

	# draw new data from fitted means + normal sigma noise (obs noise)
	est_counts_raw <- fitdatainfec + rnorm(length(fitdatainfec), 0, parmeans[['sigma']])
	#est_counts_raw     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)
    est_counts     <- ifelse(est_counts_raw < 0, 0, est_counts_raw)

    # refit using new data
    rm(if2) # because stupid things get done in packages
    sourceCpp(if2file)
    if2time <- system.time( if2data <- if2(est_counts, Tlim+1, N, NP, nPasses, coolrate) )

    paramdata <- data.frame( if2data$paramdata )
	names(paramdata) <- c("R0", "r", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")
	parmeans <- colMeans(paramdata)
	names(parmeans) <- c("R0", "r", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")

	## generate the rest of the trajectory
	##

	# pack new parameter estimates

	pars <- with( 	as.list(parmeans),
	              c(R0 = R0,
          			r = r,
          			N = N,
          			eta = eta,
          			berr = berr) )

	init_cond <- c(S = if2data$statemeans[dim(if2data$statemeans)[1],1],
                   I = if2data$statemeans[dim(if2data$statemeans)[1],2],
                   R = if2data$statemeans[dim(if2data$statemeans)[1],3])

	# generate trajectory

	sdeout <- StocSIR(init_cond, pars, T-Tlim-1, steps)
	colnames(sdeout) <- c('S','I','R','B')

	return( c( counts = sdeout[,'I'],
	       		parmeans,
	       		time = if2time[['user.self']]) )


}


stopCluster(cl)

save.image(file = "if2-paraboot.RData")