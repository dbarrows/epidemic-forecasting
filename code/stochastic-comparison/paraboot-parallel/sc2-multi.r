# Dexter Barrows
#
# run: Rscript sc2-multi.r > sc2-multi.log 2>&1


library(reshape2)
library(foreach)
library(parallel)
library(doParallel)
library(rstan)

## Stochastic SIR functions
##############################################################################################

source(paste(getwd(),"../../sir-functions", "StocSIR.r", sep = "/"))
source(paste(getwd(),"../../sir-functions", "StocSIRstan.r", sep = "/"))



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
maxTrunc <- 40
nTraj <- 200

numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)

## external files

# R source files
stoc_sir_file 		<- paste(getwd(),"../../sir-functions", "StocSIR.r", sep = "/")
if2_paraboot_file 	<- paste(getwd(), "if2-paraboot.r", sep = "/")
source(stoc_sir_file)
source(if2_paraboot_file)

# C++ source files
if2file 	<- paste(getwd(), "../../if2", "if2.cpp", sep="/")
if2_s_file 	<- paste(getwd(), "../../if2", "if2-s.cpp", sep="/")
hmcfile 	<- paste(getwd(), "../../hmc", "sirode_euler.stan", sep="/")

## IF2 settings

NP          <- 2500
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

stan_options <- list(   chains = 1,    		# number of chains
                        iter   = 2000, 		# iterations per chain
                        warmup = 1000, 		# warmup interations
                        thin   = 1)   		# thinning number

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

SSEmat <- foreach ( trunc = seq(4,maxTrunc,4), .combine = rbind, .packages = c("Rcpp","rstan","reshape2","foreach") ) %dopar% {

	Tlim <- T - trunc

	if2SSEs <- numeric(nTrials)
	hmcSSEs <- numeric(nTrials)

	for(trial in 1:nTrials) {

		# get true trajectory
	    sdeout_true <- StocSIR(true_init_cond, pars_true, T, steps)
	    colnames(sdeout_true) <- c('S','I','R','B')

	    # get full data
	    infec_counts_raw <- sdeout_true[,'I'] + rnorm(T+1, 0, sigma)
	    infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

	    # truncate data
	    datapart <- c(infec_counts[1:(Tlim+1)],rep(NA,T-Tlim))

	    ## IF2
	    ##

	    # initial fit
	    sourceCpp(if2file)
	    if2time <- system.time( if2data <- if2(datapart, Tlim+1, N, NP, nPasses, coolrate) )

	    # IF2 parametric bootstrap
	    if2_paraboot_data <- if2_paraboot(if2data,
                                  T, Tlim, steps, N, nTraj,
                                  if2file, if2_s_file, stoc_sir_file,
                                  NP, nPasses, coolrate)

	    # get mean of trajectories
		parabootdata <- data.frame(if2_paraboot_data)
		countnames <- paste("counts",1:(dim(parabootdata)[2]-9), sep = "")
		countdata <- parabootdata[,countnames]
		countmeans <- colMeans(countdata)

		# get SSE, save
		truefuture  <- sdeout_true[(Tlim+2):(T+1),'I']
    	estfuture   <- countmeans[-1]
    	err <- estfuture - truefuture
    	sse <- sum(err^2)
    	if2SSEs[trial] <- sse


	    ## HMC fitting with stan
	    ##

	    datlen <- Tlim*7 + 1

	    data <- matrix(data = -1, nrow = T+1, ncol = steps)
	    data[,1] <- datapart
	    standata <- as.vector(t(data))[1:datlen]

	    sir_data <- list( T = datlen,       # simulation time
	                      y = standata,     # infection count data
	                      N = 500,          # population size
	                      h = 1/steps )     # step size per day

	    hmctime <- system.time(fit <- with(stan_options,
	                            stan( fit = initialfit,
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

	    for (j in 1:datlen) {
			varname <- paste('Bnoise[', j, ']', sep = '')
			paramdata[[varname]] <- melt( exfit[,,varname] )$value
		}

	    parnames <- c("R0","r","sigma","eta","berr","Iinit")
	    parvars <- var(paramdata[parnames])
	    hmcparmeans <- colMeans(paramdata[parnames])
	    names(hmcparmeans) <- c("hmc.R0","hmc.r","hmc.sigma","hmc.eta","hmc.berr","hmc.Iinit")

	    bootstrapdata <- matrix(NA, nrow = nTraj, ncol = (T+1) )
    
    	pardatlen 	<- dim(paramdata)[1]
    	inds 	    <- sample.int(pardatlen,nTraj,replace = TRUE)
    	params 	    <- paramdata[inds,]
    
    	for (i in 1:nTraj) {
    
    		paramset <- params[i,]

        	init_cond <- c(S = paramset$Sinit,
        	               I = paramset$Iinit,
        	               R = paramset$Rinit)
        	pars <- c(R0 = paramset$R0,
        	          r = paramset$r,
        	          N = 500.0,
        	          eta = paramset$eta,
        	          berr = paramset$berr)
        
        	berrvec <- numeric(datlen)
        	for (j in 1:datlen) {
        		varname <- paste("Bnoise[", j, "]", sep = "")
        		berrvec[j] <- paramset[[varname]]
        	}
        
        	sdeout <- StocSIRstan(init_cond, pars, T, steps, berrvec, datlen)
        	colnames(sdeout) <- c('S','I','R','B')
        
        	bootstrapdata[i,] <- sdeout[,'I']
    
    	}
    
    	# in case of explosion
    	bootstrapdata <- bootstrapdata[complete.cases(bootstrapdata),]
    	goodbootsize <- dim(bootstrapdata)[1]
    
    	meanTraj 	<- colMeans(bootstrapdata)
    
    	truefuture  <- sdeout_true[(Tlim+2):(T+1),'I']
		estfuture   <- meanTraj[(Tlim+2):(T+1)]
    	err <- estfuture - truefuture
    	sse <- sum(err^2)
    	hmcSSEs[trial] <- sse

	}

	return( c(if2SSEs, hmcSSEs) )

}


stopCluster(cl)

save.image(file = "sc2-multi.RData")
