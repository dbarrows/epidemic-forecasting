library(reshape2)
library(rstan)
#library(shinystan)
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

StocSIRstan <- function(y, pars, T, steps, berrvec) {

	out <- matrix(NA, nrow = (T+1), ncol = 4)

	R0 <- pars[['R0']]
	r <- pars[['r']]
	N <- pars[['N']]
	eta <- pars[['eta']]
	#berr <- pars[['berr']]

	S <- y[['S']]
	I <- y[['I']]
	R <- y[['R']]

	B0 <- R0 * r / N
	B <- B0

	out[1,] <- c(S,I,R,B)

	h <- 1 / steps

	for ( i in 1:(T*steps) ) {

		B <- exp( log(B) + eta*(log(B0) - log(B)) + berrvec[i])

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

## Generate initial trajectory and data
##

set.seed(1004)

T 		<- 60
i_infec <- 5
steps 	<- 7
N 		<- 500
sigma 	<- 10
Tlim 	<- T

## Generate true trajecory and synthetic data
##

pars_true <- c(R0 = 3.0, 	# new infected people per infected person
          r = 0.1, 		# recovery rate
		  N = 500, 		# population size
		  eta = 0.5, 	# geometric random walk
		  berr = 0.5) 	# Beta geometric walk noise

true_init_cond <- c(S = N - i_infec,
					I = i_infec,
					R = 0)

# setup cluster
numCores <- detectCores()
numCores
cl <- makeCluster(numCores)
registerDoParallel(cl)

## STAN 
##################################################################################################


# options
stan_options <- list(   chains = 1,  # number of chains
                        iter   = 2000,      # iterations per chain, including warmup
                        warmup = 1000,      # warmup interations
                        thin   = 1)         # thinning number


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
sir_data <- list( T = datlen,       # simulation time
                  y = standata,     # infection count data
                  N = 500,          # population size
                  h = 1/steps )     # step size per day

hmcfile <- "sirode_euler.stan"

# fit!
initialfit <- with(stan_options,
                stan(file   = hmcfile,
                    data    = sir_data,
                    chains  = chains,
                    iter    = iter,
                    warmup  = warmup,
                    thin    = thin)
                )

#hmcvals <- numeric(10)

## get new data and refit in parallel
##
hmcvals <- foreach (i = 1:10, .combine = cbind, .packages = c("rstan","reshape2")) %dopar% {

	sdeout_true <- StocSIR(true_init_cond, pars_true, T, steps)
	colnames(sdeout_true) <- c('S','I','R','B')

	infec_counts_raw <- sdeout_true[,'I'] + rnorm(T+1, 0, sigma)
	infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)


	## HMC fitting with stan
	##

	datlen <- T*7 + 1

	data <- matrix(data = -1, nrow = T+1, ncol = steps)
	data[,1] <- infec_counts
	standata <- as.vector(t(data))[1:datlen]

	sir_data <- list( T = datlen,   	# simulation time
	                  y = standata, 	# infection count data
	                  N = 500,      	# population size
	                  h = 1/steps )   	# step size per day 
	                    
	#rstan_options(auto_write = TRUE)
	#options(mc.cores = numCores)

	hmctime <- system.time(fit <- with(stan_options,
			            	stan(#file  	= hmcfile,
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

	#for (j in 1:datlen) {
	#	varname <- paste('Bnoise[', j, ']', sep = '')
	#	paramdata[[varname]] <- melt( exfit[,i,varname] )$value
	#}

	parnames <- c("R0","r","sigma","eta","berr","Iinit")
	parvars <- var(paramdata[parnames])
	parmeans <- colMeans(paramdata[parnames])

	res <- t(1 / as.matrix(parmeans)) %*% parvars %*% (1 / as.matrix(parmeans))
	val <- sqrt( res / ( (stan_options$iter - stan_options$warmup) / stan_options$thin * length(parmeans)) )

	#hmcvals[i] <- val

    return(val)

}

print("hmcvals")
hmcvals

target_err <- mean(hmcvals)

print("target_error")
target_err

##################################################################################################



## IF2 
##################################################################################################

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
	if2time <- system.time( if2data <- if2(infec_counts[1:(Tlim+1)], Tlim+1, N, NP, nPasses, coolrate) )

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

target_particles <- mean(if2res) / (target_err^2 * 6)
target_particles


## Re-run PF stuff with target number of particles
##

NP          <- target_particles
nPasses     <- 15
coolrate    <- 0.90

if2file <- paste(getwd(),"if2-d.cpp",sep="/")

#if2vals <- numeric(10)
#if2res <- numeric(10)

if2mat <- foreach (i = 1:10, .combine = rbind, .packages = "Rcpp") %dopar% {

	sdeout_true <- StocSIR(true_init_cond, pars_true, T, steps)
	colnames(sdeout_true) <- c('S','I','R','B')

	infec_counts_raw <- sdeout_true[,'I'] + rnorm(T+1, 0, sigma)
	infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

    sourceCpp(if2file)
	if2time <- system.time( if2data <- if2(infec_counts[1:(Tlim+1)], Tlim+1, N, NP, nPasses, coolrate) )

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
if2_err


##################################################################################################

# save everything
save.image(file = "mchammererror.RData")