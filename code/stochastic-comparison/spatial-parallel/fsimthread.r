library(reshape2)
library(rstan)
library(Rcpp)

## expects 'outdir' as argument
trial <- TRIAL

## get command line arg
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (output dir).\n")
} else if (length(args)==1) {
  outdir <- args[1]
}

print(outdir)
imfile <- paste(outdir, "/", "thread-", trial, ".RData", sep = "")
print(imfile)

##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------

## number of forcast trajecotries to draw for each method
## (if2: parametric bootstrap, hmcmc: bootstrap)
nTraj <- 200

##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------

T 		<- 60
N 		<- 500
i_infec <- 5
nloc  	<- 5
steps   <- 7
sigma   <- 10
Tlim 	<- 50

pars <- c(	R0 	 = 3.0, # new infected people per infected person
          	r 	 = 0.1, # recovery rate
          	N 	 = 500, # population size
          	eta  = 0.5, # geometric random walk
          	berr = 0.5, # Beta geometric walk noise
          	phi  = 0.1) # degree of interconnectivity

initcond <- matrix(NA, nloc, 3)
colnames(initcond) <- c("S","I","R")

## default (everyone susceptible)
initcond[,"S"] <- N
initcond[,"I"] <- 0
initcond[,"R"] <- 0

## one region infected
initcond[2,] <- c(N - i_infec, i_infec, 0);

## set up number of neighbors vector and neighbor index matrix

## ring topology:
neinum <- rep(2,nloc)
neibmat <- matrix(-1, nloc, nloc)
neibmat[,1] <- c( ((0:(nloc-1)-1)%%nloc)+1 )
neibmat[,2] <- c( ((0:(nloc-1)+1)%%nloc)+1 )

## external files

# R source files
stoc_ssir_file 			<- paste(getwd(),"../../sir-functions", "StocSSIR.r", sep = "/")
stoc_ssir_stan_file 	<- paste(getwd(),"../../sir-functions", "StocSSIRstan.r", sep = "/")
if2_ssir_paraboot_file 	<- paste(getwd(), "if2-spa-paraboot.r", sep = "/")
dewdrop_file 				<- paste(getwd(), "../../smap/dewdrop", "dewdrop.r", sep = "/")
source(stoc_ssir_file)
source(stoc_ssir_stan_file)
source(if2_ssir_paraboot_file)
source(dewdrop_file)

# C++ source files
if2file 	<- paste(getwd(), "../../if2/if2-spa", "if2-spa.cpp", sep="/")
hmcfile 	<- paste(getwd(), "../../hmc/hmc-spa", "sir-spa-euler.stan", sep="/")

## IF2 settings

NP          <- 2500*2
nPasses     <- 50
coolrate    <- 0.975

## Smap settings

E 			<- 14
theta 		<- 3

# get raw data
ssdeout_true <- StocSSIR(initcond, pars, T, steps, neinum, neibmat)
colnames(ssdeout_true) <- c('S','I','R','B')

infec_counts_raw <- ssdeout_true[,'I',] + matrix(rnorm(nloc*(T+1), 0, sigma), nloc, T+1)
infec_counts <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

## fit once to compile
##

# conform data to stan's required format, short data (not full)
samsize <- 10
datlen <- (samsize-1)*steps + 1

standata <- matrix(NA, nloc, datlen)

for (loc in 1:nloc) {
	data <- matrix(data = -1, nrow = T+1, ncol = steps)
	data[,1] <- infec_counts[loc,]
	standata[loc,] <- as.vector(t(data))[1:datlen]
}

# data
sir_data <- list( T = datlen,   		# simulation time
                  nloc = nloc, 			# number of locations
                  y = standata, 		# infection count data
                  N = 500,      		# population size
                  h = 1/steps,    		# step size per day 
                  neinum = neinum, 		# number of neibours for each location
                  neibmat = neibmat)	# neighbour list for each location
                    
rstan_options(auto_write = TRUE)

# options
stan_options <- list(   chains = 1,    		# number of chains
                        iter   = 4000, 		# iterations per chain
                        warmup = 2000, 		# warmup interations
                        thin   = 1,   		# thinning number
                        verbose = FALSE,
                        refresh = 10)

# fit!
initialfit <- with(stan_options,
            stan(file  	= hmcfile,
	            data    = sir_data,
	            chains  = chains,
	            iter    = iter,
	            warmup  = warmup,
	            thin    = thin,
	            verbose = verbose,
	            refresh = refresh )
        	)


## Run trajectory
#########################################################################################################

datapart <- infec_counts[,1:(Tlim+1)]

## IF2
################################################################################################
################################################################################################

# initial fit
if2time1 <- system.time( if2data <- if2_spa(datapart, Tlim+1, N, NP, nPasses, coolrate, neinum, neibmat, nloc) )

# IF2 parametric bootstrap
source(if2_ssir_paraboot_file) ## temp
if2time2 <- system.time(if2_paraboot_data <- if2_spa_paraboot(if2data,
                          T, Tlim, steps, N, nTraj,
                          if2file, stoc_ssir_file,
                          NP, nPasses, coolrate,
                          neinum, neibmat)
)

if2time <- if2time1 + if2time2

countmeans <- matrix(NA, nloc, T-Tlim+1)
parabootdata <- array(NA, c(nTraj, nloc, T-Tlim+1))

# get mean of trajectories
for (i in 1:nTraj) {
	datalist <- if2_paraboot_data[[i]]
	trajectories <- datalist$counts
	parabootdata[i,,] <- trajectories
}
for (i in 1:nloc) {
	locslice <- parabootdata[,i,]
	countmeans[i,] <- colMeans(locslice)
}

# get SSE, save
truefuture  <- ssdeout_true[,'I',(Tlim+2):(T+1)]
estfuture   <- countmeans[,-1]
err <- estfuture - truefuture
if2sse <- sum(err^2)
if2normsse <- if2sse / (nloc*T-Tlim)

if2pred <- estfuture


## HMC fitting with stan
################################################################################################
################################################################################################

datlen <- Tlim*steps + 1

standata <- matrix(NA, nloc, datlen)

for (loc in 1:nloc) {
	data <- matrix(data = -1, nrow = T+1, ncol = steps)
	data[,1] <- infec_counts[loc,]
	standata[loc,] <- as.vector(t(data))[1:datlen]
}

sir_data <- list( T = datlen,   		# simulation time
                  nloc = nloc, 			# number of locations
                  y = standata, 		# infection count data
                  N = 500,      		# population size
                  h = 1/steps,    		# step size per day 
                  neinum = neinum, 		# number of neibours for each location
                  neibmat = neibmat)	# neighbour list for each location

hmctime1 <- system.time(fit <- with(stan_options,
			                        stan( fit = initialfit,
			                              data    = sir_data,
			                              chains  = chains,
			                              iter    = iter,
			                              warmup  = warmup,
			                              thin    = thin
			                            )
			                        )
            			)


exfit <- extract(fit, permuted = FALSE, inc_warmup = FALSE)

paramdata <- data.frame(R0 		= melt(exfit[,,'R0'])$value,
               			r 		= melt(exfit[,,'r'])$value,
               			sigma 	= melt(exfit[,,'sigma'])$value,
               			eta 	= melt(exfit[,,'eta'])$value,
               			berr 	= melt(exfit[,,'berr'])$value,
               			phi 	= melt(exfit[,,'phi'])$value )

# Iinit extraction
for (j in 1:nloc) {
	varname <- paste('Iinit[', j, ']', sep = '')
	paramdata[[varname]] <- melt( exfit[,,varname] )$value
}

# Bnoise extraction
for (loc in 1:nloc) {
	for (t in 1:datlen) {
		varname <- paste('Bnoise[', loc,',',t,']', sep = '')
		paramdata[[varname]] <- melt( exfit[,,varname] )$value
	}
}

bootstrapdata <- array(NA, c(nTraj, nloc, T+1))

pardatlen 	<- dim(paramdata)[1]
inds 	    <- sample.int(pardatlen,nTraj,replace = TRUE)
params 	    <- paramdata[inds,]

hmctime2 <- system.time(
for (i in 1:nTraj) {

	paramset <- params[i,]

	init_cond <- matrix(NA, nloc, 3)
	for (loc in 1:nloc) {
		varname <- paste("Iinit[", loc, "]", sep = "")
		I0 <- paramset[[varname]]
		init_cond[loc,] <- c(N - I0, I0, 0)
	}
	colnames(init_cond) <- c("S", "I", "R")

	pars <- c(R0 = paramset$R0,
	          r = paramset$r,
	          N = 500.0,
	          eta = paramset$eta,
	          berr = paramset$berr,
	          phi = paramset$phi)

	berrmat <- matrix(NA, nloc, datlen)
	for (loc in 1:nloc) {
		for (t in 1:datlen) {
			varname <- paste("Bnoise[", loc, ",", t, "]", sep = "")
			berrmat[loc,t] <- paramset[[varname]]
		}
	}

	ssdeout <- StocSSIRstan(init_cond, pars, T, steps, neinum, neibmat, berrmat, datlen)
	colnames(ssdeout) <- c('S','I','R','B')

	bootstrapdata[,i,] <- ssdeout[,'I',]

}
)

hmctime <- hmctime1 + hmctime2

meanTraj <- matrix(NA, nloc, T+1)

# in case of explosion
for (loc in 1:nloc) {
	locdata <- bootstrapdata[loc,,]
	locdata <- locdata[complete.cases(locdata),]
	meanTraj[loc,] <- colMeans(locdata)
}

meanTraj 	<- colMeans(bootstrapdata)

truefuture  <- ssdeout_true[,'I',(Tlim+2):(T+1)]
estfuture   <- meanTraj[,(Tlim+2):(T+1)]
err <- estfuture - truefuture
hmcsse <- sum(err^2)
hmcnormsse <- hmcsse / (nloc*(T-Tlim))

hmcpred <- estfuture


## Smapping
################################################################################################
################################################################################################

stepsAhead <- T - Tlim
smaptime <- system.time( predictions <- dewdrop(datapart, E, theta, stepsAhead) )

truefuture  <- ssdeout_true[,'I',(Tlim+2):(T+1)]
estfuture   <- predictions
err <- estfuture - truefuture
smapsse <- sum(err^2)
smapnormsse <- smapsse / (nloc*(T-Tlim))

smappred <- predictions

## Save results
#########################################################################################################

save.image(imfile)