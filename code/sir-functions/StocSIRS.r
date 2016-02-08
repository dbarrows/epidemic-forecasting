StocSIRS <- function(y, pars, T, steps) {

	out <- matrix(NA, nrow = (T+1), ncol = 4)

	R0 <- pars[['R0']]
	r <- pars[['r']]
	N <- pars[['N']]
	eta <- pars[['eta']]
	berr <- pars[['berr']]
    re <- pars[['re']]

	S <- y[['S']]
	I <- y[['I']]
	R <- y[['R']]

	B0 <- R0 * r / N
	B <- B0

	out[1,] <- c(S,I,R,B)

	h <- 1 / steps

	for ( i in 1:(T*steps) ) {

        #Bfac <- 1/2 - cos((2*pi/365)*i)/2
        Bfac <- exp(2*cos((2*pi/365)*i) - 2)

		B <- exp( log(B) + eta*(log(B0) - log(B)) + rnorm(1, 0, berr) )

		BSI <- Bfac*B*S*I
		rI <- r*I
        reR <- re*R 

		dS <- -BSI + reR
		dI <- BSI - rI
		dR <- rI - reR

		S <- S + h*dS  #newInf
		I <- I + h*dI  #newInf - h*dR
		R <- R + h*dR  #h*dR

		if (i %% steps == 0)
			out[i/steps+1,] <- c(S,I,R,B)

	}

	return(out)

}

### suggested parameters
#
# T 		<- 200
# i_infec 	<- 5
# steps 	<- 7
# N 		<- 500
# sigma 	<- 5
#
# pars <- c(R0 = 3.0, 	# new infected people per infected person
#           r = 0.1, 	# recovery rate
#           gam = 2, 	# new infected shock intensity
#		    N = 500,    # population size
#		    eta = 0.5, 	# geometric random walk
#		    berr = 0.5, # Beta geometric walk noise
#           re = 0.05)  # resuceptibility rate