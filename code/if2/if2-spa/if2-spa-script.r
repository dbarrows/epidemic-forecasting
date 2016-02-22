library(Rcpp)

source(paste(getwd(), "../../sir-functions", "StocSSIR.r", sep = "/"))

T 		<- 60
N 		<- 500
i_infec <- 5
nloc  	<- 10
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
neibmat <- matrix(NA, nloc, 3)
neibmat[,1] <- c( ((0:(nloc-1)-1)%%nloc)+1 )
neibmat[,2] <- c( ((0:(nloc-1)+1)%%nloc)+1 )

ssdeout <- StocSSIR(initcond, pars, T, steps, neinum, neibmat)

infec_counts <- ssdeout[,'I',] + matrix(rnorm(nloc*(T+1), 0, sigma), nloc, T+1)
infec_counts <- ifelse(infec_counts < 0, 0, infec_counts)

NP <- 2500*nloc
nPasses <- 50
coolrate <- 0.975

if2_spa_file <- paste(getwd(), "if2-spa.cpp", sep = "/")
sourceCpp(if2_spa_file)
if2data <- if2_spa(infec_counts, T, N, NP, nPasses, coolrate, neinum, neibmat, nloc)