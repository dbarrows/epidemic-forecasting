library(Rcpp)
library(rstan)
library(reshape2)

source(paste(getwd(), "../../sir-functions", "StocSSIR.r", sep = "/"))

T 		<- 60
N 		<- 500
i_infec <- 5
nloc  	<- 5
steps   <- 7
sigma   <- 10

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

ssdeout <- StocSSIR(initcond, pars, T, steps, neinum, neibmat)

infec_counts <- ssdeout[,'I',] + matrix(rnorm(nloc*(T+1), 0, sigma), nloc, T+1)
infec_counts <- ifelse(infec_counts < 0, 0, infec_counts)

datlen <- T*7 + 1

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
                    
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stan_options <- list(   chains = 1,    		# number of chains
                        iter   = 4000, 		# iterations per chain
                        warmup = 2000, 		# warmup interations
                        thin   = 1,   		# thinning number
                        verbose = FALSE,
                        refresh = 50)

fit <- with(stan_options,
            stan(file  	= "sir-spa-euler.stan",
	            data    = sir_data,
	            chains  = chains,
	            iter    = iter,
	            warmup  = warmup,
	            thin    = thin,
	            verbose = verbose,
	            refresh = refresh )
        	)

exfit <- extract(fit, permuted = FALSE, inc_warmup = FALSE)
paramdata <- data.frame(R0 = melt(exfit[,,'R0'])$value,
               			r = melt(exfit[,,'r'])$value,
               			sigma = melt(exfit[,,'sigma'])$value,
               			eta = melt(exfit[,,'eta'])$value,
               			berr = melt(exfit[,,'berr'])$value,
               			phi = melt(exfit[,,'phi'])$value )

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