library(pomp)

## Fake data

N <- 500
T <- 60
steps <- 7
h <- 1 / steps
sigma <- 10

pars_true <- c(R0 = 3.0, 	# new infected people per infected person
          		r = 0.1, 	# recovery rate
          		sigma = 5,	# observation error intensity
          		I0 = 5)		# initial condition, number of infected people

B <- pars_true['R0'] * pars_true['r'] / N
     
odeout <- matrix(NA, nrow = (T+1), ncol = 3)
odeout[1,] <- c(N - pars_true['I0'], pars_true['I0'],0)
colnames(odeout) <- c("S", "I", "R")

S <- odeout[1,1]
I <- odeout[1,2]
R <- odeout[1,3]
R0 <- pars_true['R0']
r <- pars_true['r']

for (i in 1:(T*steps)) {

    dS <- - B*S*I
    dI <- B*S*I - r*I
    dR <- r*I

    S <- S + h*dS;
    I <- I + h*dI;
    R <- R + h*dR; 

    if (i %% 7 == 0)
    	odeout[(i/7) + 1,] = c(S, I, R)

}

infec_counts_raw <- odeout[,'I'] + rnorm(T+1, 0, sigma)
infec_counts <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

data <- data.frame(time = 0:T, y = infec_counts)



sir_step <- Csnippet("
	double dS = - (R0*r/500.0)*S*I;
	double dI = (R0*r/500.0)*S*I - r*I;
	double dR = r*I;
	S += dt*dS;
	I += dt*dI;
	R += dt*dR;
")

sir_init <- Csnippet("
	S = 500.0 - I0;
	I = I0;
	R = 0;
")

dmeas <- Csnippet("
  	lik = dnorm(y, I, sigma, give_log);
")

rmeas <- Csnippet("
  y = rnorm(I, sigma);
")

sir <- pomp(data = data,
            time = "time",
            t0 = 0,
            initializer = sir_init,
            rprocess = euler.sim(step.fun = sir_step,
                                 delta.t = 1.0/steps),
            dmeasure = dmeas,
            rmeasure = rmeas,
            statenames = c("S","I","R"),
            paramnames = c("R0","r","I0","sigma"))

simStates <- simulate(sir,nsim=1,params=c(R0 = 3.0, r = 0.1, I0 = 5.0, sigma = 10.0),
                      states = TRUE, obvs = TRUE, as = TRUE)

m1 <- mif2(sir,
	      Nmif = 50,
	      start = c(R0 = 3.0, r = 0.1, I0 = 5.0, sigma = 10.0),
	      rw.sd = rw.sd(R0 = 0.3, r = 0.01, I0 = 0.5, sigma = 1.0,
	      cooling.fraction.50 = 0.95,
	      Np = 1000
	      )

plot(m1)
