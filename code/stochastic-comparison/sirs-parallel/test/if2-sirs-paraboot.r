# Dexter Barrows
#
# IF2 parametric bootstrapping function

library(foreach)
library(parallel)
library(doParallel)
library(Rcpp)

if2_sirs_paraboot <- function(if2data, T, Tlim, steps, N, nTrials, if2file, stoc_sir_file, NP, nPasses, coolrate) {

	source(stoc_sir_file)
	#print(stoc_sir_file)

	if(nTrials < 2) {
		nTrials <- 2
	}

	# unpack if2 first fit data
	# ...parameters
	paramdata <- data.frame( if2data$paramdata )
	names(paramdata) <- c("R0", "r", "re", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")
	parmeans <- colMeans(paramdata)
	names(parmeans) <- c("R0", "r", "re", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")
	# ...states
	statedata <- data.frame( if2data$statedata )
	names(statedata) <- c("S","I","R","B")
	statemeans <- colMeans(statedata)
	names(statemeans) <- c("S","I","R","B")

	#print(1:nTrials)

	## use parametric bootstrapping to generate forcasts
	##
	trajectories <- foreach( i = 1:nTrials, .combine = rbind, .packages = "Rcpp") %dopar% {

		source(stoc_sir_file)

		## draw new data
		##

		pars <- with( as.list(parmeans),
		              c(R0 = R0,
	          			r = r,
	          			re = re,
	          			N = N,
	          			eta = eta,
	          			berr = berr) )

		init_cond <- with( as.list(parmeans),
		                   c(S = Sinit,
	                   		 I = Iinit,
	                   		 R = Rinit) )

		# generate trajectory
		sdeout <- StocSIRS(init_cond, pars, Tlim + 1, steps)
		colnames(sdeout) <- c('S','I','R','B')

		# add noise
		counts_raw <- sdeout[,'I'] + rnorm(dim(sdeout)[1], 0, parmeans[['sigma']])
	    counts     <- ifelse(counts_raw < 0, 0, counts_raw)

	    ## refit using new data
	    ##

	    rm(if2) # because stupid things get done in packages
	    sourceCpp(if2file)
	    if2time <- system.time( if2data <- if2_sirs(counts, Tlim+1, N, NP, nPasses, coolrate) )

	    paramdata <- data.frame( if2data$paramdata )
		names(paramdata) <- c("R0", "r", "re", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")
		parmeans <- colMeans(paramdata)
		names(parmeans) <- c("R0", "r", "re", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")

		## generate the rest of the trajectory
		##

		# pack new parameter estimates
		pars <- with( 	as.list(parmeans),
		              c(R0 = R0,
	          			r = r,
	          			re = re,
	          			N = N,
	          			eta = eta,
	          			berr = berr) )
		init_cond <- c(S = statemeans[['S']],
	                   I = statemeans[['I']],
	                   R = statemeans[['R']])

		# generate remaining trajectory part
		sdeout_future <- StocSIRS(init_cond, pars, T-Tlim, steps)
		colnames(sdeout_future) <- c('S','I','R','B')

		return ( c( counts = unname(sdeout_future[,'I']),
		       	    parmeans,
		       	    time = if2time[['user.self']]) )


	}

	return(trajectories)

}
