library(ggplot2)
library(reshape2)

StocSIRS <- function(y, pars, T, steps) {

	out <- matrix(NA, nrow = (T+1), ncol = 4)

	R0 <- pars[['R0']]
	r <- pars[['r']]
	gam <- pars[['gam']]
	N <- pars[['N']]
	eta <- pars[['eta']]
	berr <- pars[['berr']]
    fr <- pars[['fr']]

	S <- y[['S']]
	I <- y[['I']]
	R <- y[['R']]

	B0 <- R0 * r / N
	B <- B0

	out[1,] <- c(S,I,R,B)

	h <- 1 / steps

	for ( i in 1:(T*steps) ) {

        Bfac <- 1/2 - cos((2*pi/365)*i)/2

		B <- exp( log(B) + eta*(log(B0) - log(B)) + rnorm(1, 0, berr) )

		BSI <- Bfac * B*S*I
		rI <- r*I
        frR <- fr*R 

		dS <- -BSI + frR
		dI <- BSI - rI
		dR <- rI - frR

		S <- S + h*dS  #newInf
		I <- I + h*dI  #newInf - h*dR
		R <- R + h*dR  #h*dR

		if (i %% steps == 0)
			out[i/steps+1,] <- c(S,I,R,B)

	}

	return(out)

}

set.seed(1001)

T 		<- 200
i_infec <- 5
steps 	<- 7
N 		<- 500
sigma 	<- 10

pars <- c(R0 = 3.0, 	  # new infected people per infected person
          r = 0.1, 		  # recovery rate
          gam = 2, 		  # new infected shock intensity
		  N = 500,        # population size
		  eta = 0.5, 	  # geometric random walk
		  berr = 0.5,  # Beta geometric walk noise
          fr = 0.1)       # resuceptibility rate

true_init_cond <- c(S = N - i_infec,
					I = i_infec,
					R = 0)

sdeout <- StocSIRS(true_init_cond, pars, T, steps)
colnames(sdeout) <- c('S','I','R','B')

infec_counts_raw <- sdeout[,'I'] + rnorm(T+1, 0, sigma)
infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

smapfile <- paste(getwd(),"../code/smap/smap.r", sep = "/")
smapfile
source(smapfile)

E <- 10
theta <- 10
stepsAhead <- 100

predictions <- smap(infec_counts, E, theta, stepsAhead)

dlen <- length(infec_counts)
predictee <- infec_counts[(dlen-E+1):dlen]

plotdata <- data.frame(times = 1:(length(infec_counts)+stepsAhead),                                         			  # master index of times
                       tseries = c(infec_counts, rep(NA,stepsAhead) ),                                      			  # time series data, padded
                       predictions = c( rep(NA,length(infec_counts)-1), infec_counts[length(infec_counts)], predictions), # predictions
                       predictee = c(rep(NA,length(infec_counts)-E),predictee,rep(NA,stepsAhead)) )    					  # last predictee vector


p <- ggplot(plotdata, aes(times)) +
        geom_line(aes(y = tseries), linetype = "solid") +
        geom_line(aes(y = predictions), linetype = "dashed") +
        geom_point(aes(y = predictee)) +
        labs(x = "Time", y = "Infection count", color = "") +
        theme_bw()

print(p)