library(deSolve)
library(ggplot2)
library(RColorBrewer)
library(pracma)

set.seed(1010)

## external files
##
stoc_sirs_file 	<- paste(getwd(), "../sir-functions", "StocSIRS.r", sep = "/")
smap_file 		<- paste(getwd(), "smap.r", sep = "/")
source(stoc_sirs_file)
source(smap_file)



## parameters
##
T 		<- 6*52
Tlim 	<- T - 52
i_infec <- 10
steps 	<- 7
N 		<- 500
sigma 	<- 5

true_pars <- c(	R0 = 3.0, 	# new infected people per infected person
          		r = 0.1, 	# recovery rate
		    	N = 500,    # population size
		    	eta = 0.5, 	# geometric random walk
		    	berr = 0.5, # Beta geometric walk noise
          		re = 1)  	# resuceptibility rate

true_init_cond <- c(S = N - i_infec,
                    I = i_infec,
                    R = 0)

## trial parameter values to check.options
##
Elist <- 1:20
thetalist <- 10*exp(-(seq(0,9.5,0.5)))
nTrials <- 100

ssemat <- matrix(NA, 20, 20)

for (i in 1:length(Elist)) {
	for (j in 1:length(thetalist)) {

		ssemean <- 0

		for (k in 1:nTrials) {

			E <- Elist[i]
			theta <- thetalist[j]

			## get true trajectory
			##
			sdeout <- StocSIRS(true_init_cond, true_pars, T, steps)

			## perturb to get data
			##
			infec_counts_raw <- sdeout[1:(Tlim+1),'I'] + rnorm(Tlim+1,0,sigma)
			infec_counts 	 <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

			predictions <- smap(infec_counts, E, theta, 52)

			err <- sdeout[(Tlim+2):dim(sdeout)[1],'I'] - predictions
			sse <- sum(err^2)

			ssemean <- ssemean + (sse / nTrials)

		}

		ssemat[i,j] <- ssemean

		
	}
}

quartz()
image(-ssemat)
quartz()
filled.contour(-ssemat)

#print(ssemat)
#cms <- colMeans(ssemat)
#rms <- rowMeans(ssemat)

#Emin <- Elist[which.min(rms)]
#thetamin <- thetalist[which.min(cms)]
#print(Emin)
#print(thetamin)

mininds <- which(ssemat==min(ssemat),arr.ind=TRUE)

Emin <- Elist[mininds[,'row']]
thetamin <- thetalist[mininds[,'col']]

print(Emin)
print(thetamin)