# Dexter Barrows
#
# IF2 parametric bootstrapping function

library(foreach)
library(parallel)
library(doParallel)
library(Rcpp)

if2_spa_paraboot <- function(if2data_parent, T, Tlim, steps, N, nTrials, if2file, stoc_sir_file, NP, nPasses, coolrate, neinum, neibmat) {

	nloc <- length(neinum)

	source(stoc_sir_file)
	#print(stoc_sir_file)

	if(nTrials < 2) {
		nTrials <- 2
	}

	# unpack if2 first fit data

	# ...parameters
	paramdata_parent <- data.frame( if2data_parent$paramdata )
	names(paramdata_parent) <- c("R0", "r", "sigma", "eta", "berr", "phi")
	parmeans_parent <- colMeans(paramdata_parent)
	names(parmeans_parent) <- names(paramdata_parent)

	# ...intital infected
	Iinitdata_parent <- data.frame( if2data$initInfec )
	Iinitmeans_parent <- rowMeans(Iinitdata_parent)

	# ...states
	statedata_parent <- data.frame( if2data_parent$finalstate )
	names(statedata_parent) <- c("S","I","R","B")
	#statemeans_parent <- colMeans(statedata_parent)
	#names(statemeans_parent) <- names(statedata_parent)


	## use parametric bootstrapping to generate forcasts
	##
	trajectories <- foreach( i = 1:nTrials, .combine = list, .packages = "Rcpp") %dopar% {

		source(stoc_sir_file)

		## draw new data
		##

		pars <- with( as.list(parmeans_parent),
		              c(R0 = R0,
	          			r = r,
	          			N = N,
	          			eta = eta,
	          			berr = berr,
	          			phi = phi) )

		init_cond <- with( as.list(parmeans_parent),
		                   matrix(c(500 - Iinitmeans_parent, Iinitmeans_parent, rep(0, nloc)), nloc, 3) )
		colnames(init_cond) <- c("S", "I", "R")

		# generate trajectory
		ssdeout <- StocSSIR(init_cond, pars, Tlim, steps, neinum, neibmat)
		colnames(ssdeout) <- c('S','I','R','B')

		# add noise
	    counts_raw 	<- ssdeout[,'I',] + matrix( rnorm(nloc*(Tlim+1), 0, parmeans_parent[['sigma']]), nloc, Tlim+1)
		counts 		<- ifelse(counts_raw < 0, 0, counts_raw)

	    ## refit using new data
	    ##

	    rm(if2_spa) # because stupid things get done in packages
	    sourceCpp(if2file)
	    if2time <- system.time( if2data <- if2_spa(counts, Tlim+1, N, NP, nPasses, coolrate, neinum, neibmat, nloc) )

	    paramdata <- data.frame( if2data$paramdata )
		names(paramdata) <- c("R0", "r", "sigma", "eta", "berr", "phi")
		parmeans <- colMeans(paramdata)
		names(parmeans) <- names(paramdata)

		## generate the rest of the trajectory
		##

		# pack new parameter estimates
		pars <- with( 	as.list(parmeans),
		              c(R0 = R0,
	          			r = r,
	          			N = N,
	          			eta = eta,
	          			berr = berr,
	          			phi = phi) )

		init_cond <- as.matrix(statedata_parent[,c("S", "I", "R")])

		# generate remaining trajectory part
		ssdeout_future <- StocSSIR(init_cond, pars, T-Tlim, steps, neinum, neibmat)
		colnames(ssdeout_future) <- c('S','I','R','B')

		return ( list(  counts = unname(ssdeout_future[,'I',]),
		       	    	parmeans = parmeans,
		       	    	time = if2time[['user.self']]) )

	}

	return(trajectories)

}
