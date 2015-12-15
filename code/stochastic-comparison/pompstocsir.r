library(pomp)
library(ggplot2)

## Fake data

StocSIR <- function(y, pars, T, steps) {

	out <- matrix(NA, nrow = (T+1), ncol = 4)

	R0 <- pars[['R0']]
	r <- pars[['r']]
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

		B <- exp( log(B) + eta*(log(B0) - log(B)) + rnorm(1, 0, berr) )

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

set.seed(1234)

T 		<- 60
i_infec <- 5
steps 	<- 7
N 		<- 500
sigma 	<- 10
Tlim 	<- T

## Generate true trajecory and synthetic data
##

pars_true <- c(R0 = 3.0, 	# new infected people per infected person
          r = 0.1, 		# recovery rate
		  N = 500, 		# population size
		  eta = 0.5, 	# geometric random walk
		  berr = 0.5) 	# Beta geometric walk noise

true_init_cond <- c(S = N - i_infec,
					I = i_infec,
					R = 0)

sdeout <- StocSIR(true_init_cond, pars_true, T, steps)
colnames(sdeout) <- c('S','I','R','B')

infec_counts_raw <- sdeout[,'I'] + rnorm(T+1, 0, sigma)
infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

data <- data.frame(time = 0:T, y = infec_counts)

sir_step <- Csnippet("
    B = exp( log(B) + eta*(log(R0*r/500.0) - log(B)) + rnorm(0, berr) );
	double dS = - B*S*I;
	double dI = B*S*I - r*I;
	double dR = r*I;
	S += dt*dS;
	I += dt*dI;
	R += dt*dR;
")

sir_init <- Csnippet("
	S = 500.0 - I0;
	I = I0;
	R = 0;
	B = R0*r/500.0;
")

dmeas <- Csnippet('
  	lik = dnorm(y, I, sigma, give_log);
  	//Rprintf("%lf %lf, %lf, %lf, %d , %lf, %lf, %lf", lik, y, I, sigma, give_log, S, I, R);
')

rmeas <- Csnippet('
  	y = rnorm(I, sigma);
')

logtrans <- Csnippet("
	TR0 = log(R0);
	Tr = log(r);
	Tsigma = log(sigma);
	Teta = log(eta);
	Tberr = log(berr);
	TI0 = log(I0);
")

exptrans <- Csnippet("
	TR0 = exp(R0);
	Tr = exp(r);
	Tsigma = exp(sigma);
	Teta = exp(eta);
	Tberr = exp(berr);
	TI0 = exp(I0);
")

options(verbose=FALSE) 
sir <- pomp(data = data,
            time = "time",
            t0 = 0,
            initializer = sir_init,
            rprocess = euler.sim(step.fun = sir_step,
                                 delta.t = 1.0/steps),
            dmeasure = dmeas,
            rmeasure = rmeas,
            statenames = c("S","I","R","B"),
            paramnames = c("R0","r","I0","sigma","eta","berr"),
            toEstimationScale = logtrans,
            fromEstimationScale = exptrans)

#simStates <- simulate(sir,nsim=1,params=c(R0 = 3.0, r = 0.1, I0 = 5.0, sigma = 10.0, eta = 0.5, berr = 0.3),
#                      states = TRUE, obvs = TRUE, include = FALSE, as = TRUE)

transfun <- function(params,...) {
    params[c("R0","r","I0","sigma","eta", "berr")] <- exp(params[c("R0","r","I0","sigma","eta", "berr")])
    return(params)
}

transfun_inv <- function(params,...) {
    params[c("R0","r","I0","sigma","eta", "berr")] <- log(params[c("R0","r","I0","sigma","eta", "berr")])
    return(params)
}


m1 <- mif2(sir,
	      Nmif = 300,
	      start = c(R0 = 3.0, r = 0.1, I0 = 5.0, sigma = 10.0, eta = 0.5, berr = 0.3),
	      rw.sd = rw.sd(R0 = 0.3, r = 0.01, I0 = 0.5, sigma = 1.0, eta = 0.05, berr = 0.03),
	      cooling.fraction.50 = 0.90,
	      Np = 3000,
	      toEstimationScale = transfun,
          fromEstimationScale = transfun_inv,
	      transform = TRUE
	      )

#plot(m1)
coef(m1)

if2pars <- coef(m1)[-c(grep('sigma', names(coef(m1))),grep('I0', names(coef(m1))))]
if2pars[['N']] <- 500

numtraj <- 100

trajData <- matrix(NA, nrow = numtraj, ncol = (T+1))
for(i in 1:numtraj) {
	sdeout <- StocSIR(true_init_cond, if2pars, T, steps)
	colnames(sdeout) <- c('S','I','R','B')
	trajData[i,] <- sdeout[,'I']
}

meanTraj 	<- colMeans(trajData)
quantTraj 	<- t(apply(trajData, 2, quantile, probs = c(0.025,0.975)))
colnames(quantTraj) <- c("025","975")

qplot(1:(T+1), meanTraj, geom = "line", xlab = "Time", ylab = "Infection count") +
	geom_ribbon(aes(ymin = quantTraj[,'025'], ymax=quantTraj[,'975']), alpha=0.1) +
	geom_line(aes(y = infec_counts)) +
	theme_bw()