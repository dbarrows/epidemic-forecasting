StocSIR <- function(y, pars, T, steps) {

	out <- matrix(NA, nrow = T, ncol = 3)

	R0 <- pars[['R0']]
	r <- pars[['r']]
	gam <- pars[['gam']]
	N <- pars[['N']]

	S <- y[['S']]
	I <- y[['I']]
	R <- y[['R']]

	B <- R0 * r / N

	h <- 1 / steps

	for ( i in 1:(T*steps) ) {

		BSI <- B*S*I
		newInf <- rpois(1,h*BSI)
		rI <- r*I

		dS <- -BSI
		dI <- BSI - rI
		dR <- rI

		S <- S - newInf
		I <- I + newInf - h*dR
		R <- R + h*dR

		if (i %% steps == 0)
			out[i/steps,] <- c(S,I,R)

	}

	return(out)

}

T 		<- 100
i_infec <- 5
steps 	<- 1

pars <- c(R0 = 3.0, # new infected people per infected person
          r = 0.1, # recovery rate
          gam = 2, # new infected shock intensity
		  N = 500) # population size

true_init_cond <- c(S = N - i_infec,
					I = i_infec,
					R = 0)

sdeout <- StocSIR(true_init_cond, pars, T, steps)