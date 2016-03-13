## ymat: 	Contains the initial conditions where:
#  			- rows are locations
#  			- columns are S, I, R
## pars: 	Contains the parameters: global values for R0, r, N, eta, berr
## T: 		The stop time. Since 0 in included, there should be T+1 time steps in the simulation
## neinum:	Number of neighbors for each location, in order
## neibmat: Contains lists of neighbors for each location
#			- rows are parent locations (nodes)
# 			- columns are locations each parent is attached to (edges)
StocSSIRstan <- function(ymat, pars, T, steps, neinum, neibmat, berrmat, bmatlim) {

	## number of locations
    nloc <- dim(ymat)[1]

    ## storage
    ## dims are locations, (S,I,R,B), times
    # output array
    out <- array(NA, c(nloc, 4, T+1), dimnames = list(NULL, c("S","I","R","B"), NULL))
    # temp storage
    BSI <- numeric(nloc)
    rI <- numeric(nloc)

    ## extract parameters
    R0 <- pars[['R0']]
    r <- pars[['r']]
    N <- pars[['N']]
    eta <- pars[['eta']]
    berr <- pars[['berr']]
    phi <- pars[['phi']]

    B0 <- rep(R0*r/N, nloc)

    ## state vectors
    S <- ymat[,'S']
    I <- ymat[,'I']
    R <- ymat[,'R']
    B <- B0

    ## assign starting to output matrix
    out[,,1] <- cbind(ymat, B0)

    h <- 1 / steps

    for ( i in 1:(T*steps) ) {

    	if (i <= bmatlim) {
		    B <- exp( log(B) + eta*(log(B0) - log(B)) + berrmat[,i])
	    } else {
	        B <- exp( log(B) + eta*(log(B0) - log(B)) + rnorm(nloc, 0, berr) )
	    }
        

        for (loc in 1:nloc) {
        	n <- neinum[loc]
        	sphi <- 1 - phi*(n/(n+1))
        	ophi <- phi/(n+1)
        	nBIsum <- B[neibmat[loc,1:n]] %*% I[neibmat[loc,1:n]]
        	BSI[loc] <- S[loc]*( sphi*B[loc]*I[loc] + ophi*nBIsum )
        }

        #if(i == 1)
        #	print(BSI)

        rI <- r*I

        dS <- -BSI
        dI <- BSI - rI
        dR <- rI

        S <- S + h*dS
        I <- I + h*dI
        R <- R + h*dR

        if (i %% steps == 0)
            out[,,i/steps+1] <- cbind(S,I,R,B)

    }

    #out[,,2] <- cbind(S,I,R,B)

	return(out)

}

### Suggested parameters
#
# T       <- 60
# i_infec <- 5
# steps   <- 7
# N       <- 500
# sigma   <- 10
#
# pars <- c(R0 = 3.0,     # new infected people per infected person
#           r = 0.1,      # recovery rate
#           N = 500,      # population size
#           eta = 0.5,    # geometric random walk
#           berr = 0.5,   # Beta geometric walk noise
# 			phi = 0.5 )	  # interconnectivity degree