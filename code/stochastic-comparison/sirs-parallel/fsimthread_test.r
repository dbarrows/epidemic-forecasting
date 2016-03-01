## expects 'outdir' as argument

outdir <- getwd()
trunc <- 1
trial <- 1

## get command line arg
#args = commandArgs(trailingOnly=TRUE)
#if (length(args)==0) {
#  stop("At least one argument must be supplied (output dir).\n")
#} else if (length(args)==1) {
#  outdir <- args[1]
#}

print(outdir)
rdsfile <- paste(outdir, "/", "thread-", trunc, "-", trial, ".rds", sep = "")
print(rdsfile)

##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------

## number of forcast trajecotries to draw for each method
## (if2: parametric bootstrap, hmcmc: bootstrap)
nTraj <- 2

##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------

library(reshape2)
library(rstan)
library(Rcpp)

T       <- 100
i_infec <- 10
steps   <- 7
N       <- 500
sigma   <- 3

## Generate true trajecory and synthetic data
##

pars_true <- c(R0 = 3.0,    # new infected people per infected person
          r = 0.1,      # recovery rate
          N = 500,      # population size
          eta = 0.5,    # geometric random walk
          berr = 0.5,	# Beta geometric walk noise
          re = 1 )    	# resusceptibility

true_init_cond <- c(S = N - i_infec,
                    I = i_infec,
                    R = 0)

## external files

# R source files
stoc_sirs_file 			<- paste(getwd(),"../../sir-functions", "StocSIRS.r", sep = "/")
stoc_sirs_stan_file 	<- paste(getwd(),"../../sir-functions", "StocSIRSstan.r", sep = "/")
if2_sirs_paraboot_file 	<- paste(getwd(), "if2-sirs-paraboot.r", sep = "/")
smap_file 				<- paste(getwd(), "../../smap", "smap.r", sep = "/")
source(stoc_sirs_file)
source(stoc_sirs_stan_file)
source(if2_sirs_paraboot_file)
source(smap_file)

# C++ source files
if2file 	<- paste(getwd(), "../../if2/if2-sirs", "if2-sirs.cpp", sep="/")
hmcfile 	<- paste(getwd(), "../../hmc/hmc-sirs", "sirs-euler.stan", sep="/")

## IF2 settings

NP          <- 2500
nPasses     <- 50
coolrate    <- 0.975

## Smap settings

E 			<- 10
theta 		<- 10

# get raw data
sdeout_true <- StocSIRS(true_init_cond, pars_true, T, steps)
colnames(sdeout_true) <- c('S','I','R','B')

infec_counts_raw <- sdeout_true[,'I'] + rnorm(T+1, 0, sigma)
infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

## fit once to compile
##

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
sirs_data <- list(T = datlen,       # simulation time
                  y = standata,     # infection count data
                  N = 500,          # population size
                  h = 1/steps )     # step size per day

# fit!
initialfit <- with(stan_options,
                stan(file   = hmcfile,
                    data    = sirs_data,
                    chains  = chains,
                    iter    = iter,
                    warmup  = warmup,
                    thin    = thin)
                )


## Run multiple trajectories
#########################################################################################################

Tlim <- T - trunc

#if2SSEs  <- numeric(nTrials)
#hmcSSEs  <- numeric(nTrials)
#smapSSEs <- numeric(nTrials)

# get true trajectory
sdeout_true <- StocSIRS(true_init_cond, pars_true, T, steps)
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
if2time1 <- system.time( if2data <- if2_sirs(datapart, Tlim+1, N, NP, nPasses, coolrate) )
print(if2time1)

# IF2 parametric bootstrap
if2time2 <- system.time(if2_paraboot_data <- if2_sirs_paraboot(if2data,
                          T, Tlim, steps, N, nTraj,
                          if2file, stoc_sirs_file,
                          NP, nPasses, coolrate)
)

if2time <- if2time1 + if2time2
#print(if2time1)

# get mean of trajectories
parabootdata <- data.frame(if2_paraboot_data)
countnames <- paste("counts",1:(dim(parabootdata)[2]-10), sep = "")
countdata <- parabootdata[,countnames]
countmeans <- colMeans(countdata)

# get SSE, save
truefuture  <- sdeout_true[(Tlim+2):(T+1),'I']
estfuture   <- countmeans[-1]
err <- estfuture - truefuture
if2sse <- sum(err^2)


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

hmctime1 <- system.time( fit <- with(stan_options,
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
                        re = melt(exfit[,,'re'])$value,
                        sigma = melt(exfit[,,'sigma'])$value,
                        eta = melt(exfit[,,'eta'])$value,
                        berr = melt(exfit[,,'berr'])$value,
                        Iinit = melt(exfit[,,'Iinit'])$value )

for (j in 1:datlen) {
	varname <- paste('Bnoise[', j, ']', sep = '')
	paramdata[[varname]] <- melt( exfit[,,varname] )$value
}

parnames <- c("R0","r","re","sigma","eta","berr","Iinit")
parvars <- var(paramdata[parnames])
hmcparmeans <- colMeans(paramdata[parnames])
names(hmcparmeans) <- paste("hmc.", parnames,sep="")

bootstrapdata <- matrix(NA, nrow = nTraj, ncol = (T+1) )

pardatlen 	<- dim(paramdata)[1]
inds 	    <- sample.int(pardatlen,nTraj,replace = TRUE)
params 	    <- paramdata[inds,]

hmctime2 <- system.time(
for (i in 1:nTraj) {

	paramset <- params[i,]

	init_cond <- c(S = N - paramset$Iinit,
	               I = paramset$Iinit,
	               R = 0.0 )
	pars <- c(R0 = paramset$R0,
	          r = paramset$r,
	          re = paramset$re,
	          N = 500.0,
	          eta = paramset$eta,
	          berr = paramset$berr)

	berrvec <- numeric(datlen)
	for (j in 1:datlen) {
		varname <- paste("Bnoise[", j, "]", sep = "")
		berrvec[j] <- paramset[[varname]]
	}

	sdeout <- StocSIRSstan(init_cond, pars, T, steps, berrvec, datlen)
	colnames(sdeout) <- c('S','I','R','B')

	bootstrapdata[i,] <- sdeout[,'I']

}
)

hmctime <- hmctime1 + hmctime2

# in case of explosion
bootstrapdata <- bootstrapdata[complete.cases(bootstrapdata),]
goodbootsize <- dim(bootstrapdata)[1]

meanTraj 	<- colMeans(bootstrapdata)

truefuture  <- sdeout_true[(Tlim+2):(T+1),'I']
estfuture   <- meanTraj[(Tlim+2):(T+1)]
err <- estfuture - truefuture
hmcsse <- sum(err^2)

## Smapping
##

stepsAhead 	<- trunc
smaptime <- system.time( predictions <- smap(datapart[1:(Tlim+1)], E, theta, stepsAhead) )

truefuture  <- sdeout_true[(Tlim+2):(T+1),'I']
estfuture   <- predictions
err <- estfuture - truefuture
smapsse <- sum(err^2)


##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------

vecout <- c(trunc = trunc, trial = trial, if2sse = if2sse, hmcsse = hmcsse, smapsse = smapsse,
            if2time = if2time[['user.self']], hmctime = hmctime[['user.self']], smaptime = smaptime[['user.self']])

print(vecout)

saveRDS(vecout, rdsfile)