library(ggplot2)

StocSIR <- function(y, pars, T, steps) {

	out <- matrix(NA, nrow = (T+1), ncol = 4)

	R0 <- pars[['R0']]
	r <- pars[['r']]
	gam <- pars[['gam']]
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

		B <- exp( log(B0) + eta*(log(B) - log(B0)) + rnorm(1, 0, berr) )
		if (B < 0)
			B <- 0.0

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

set.seed(1001)

T 		<- 60
i_infec <- 5
steps 	<- 7
N 		<- 500
sigma 	<- 10

pars <- c(R0 = 3.0, 	# new infected people per infected person
          r = 0.1, 		# recovery rate
          gam = 2, 		# new infected shock intensity
		  N = 500, 		# population size
		  eta = 0.5, 	# geometric random walk
		  berr = 0.5) 	# Beta geometric walk noise

true_init_cond <- c(S = N - i_infec,
					I = i_infec,
					R = 0)

sdeout <- StocSIR(true_init_cond, pars, T, steps)
colnames(sdeout) <- c('S','I','R','B')

quartz()
plot(1:(T+1), sdeout[,'B'])

infec_counts_raw <- sdeout[,'I'] + rnorm(T+1, 0, sigma)
infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

plotdata <- data.frame(times = 0:T, true = sdeout[,'I'], data = infec_counts)

g <- ggplot(plotdata, aes(times)) +
        geom_line(aes(y = true, colour = "True")) + 
        geom_point(aes(y = data, colour = "Data")) +
        labs(x = "Time", y = "Infection count", color = "") +
        scale_color_brewer(palette="Paired") +
        theme(panel.background = element_rect(fill = "#F0F0F0"))

quartz()
print(g)