library(deSolve)
library(rstan)
library(shinystan)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(gridExtra)

# NOTE: to save plots, uncomment the ggsave lines

SIR <- function(Time, State, Pars) {
    
    with(as.list(c(State, Pars)), {
        
        B   <- R0*r/N
        BSI <- B*S*I
        rI  <- r*I
        
        dS = -BSI
        dI = BSI - rI
        dR = rI
        
        return(list(c(dS, dI, dR)))
        
    })
    
}

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

		B <- exp( log(B0) + eta*(log(B) - log(B0)) + rnorm(1, 0, berr) )

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
Tlim 	<- 25

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

infec_counts_raw <- sdeout[,'I'] + rnorm(T+1, 0, sigma)
infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

datapart <- c(infec_counts[1:(Tlim+1)],rep(NA,T-Tlim))

plotdata <- data.frame(times = 0:T, true = sdeout[,'I'], data = datapart)

g <- ggplot(plotdata, aes(times)) +
        geom_line(aes(y = true, colour = "True")) + 
        geom_point(aes(y = data, colour = "Data")) +
        labs(x = "Time", y = "Infection count", color = "") +
        scale_color_brewer(palette="Paired") +
        theme(panel.background = element_rect(fill = "#F0F0F0"))

quartz()
print(g)
#ggsave(g, filename="dataplot.pdf", height=4, width=6.5)

datlen <- T*7 + 1

data <- matrix(data = -1, nrow = T+1, ncol = steps)
data[,1] <- infec_counts
standata <- as.vector(t(data))[1:datlen]

sir_data <- list( T = datlen,   	# simulation time
                  y = standata, 	# infection count data
                  N = 500,      	# population size
                  h = 1/steps )   	# step size per day 
                    
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stan_options <- list(   chains = 1,    		# number of chains
                        iter   = 5000, 		# iterations per chain
                        warmup = 1000, 		# warmup interations
                        thin   = 10,   		# thinning number
                        verbose = TRUE,
                        refresh = 50)

fit <- with(stan_options,
            stan(file  	= "sirode_euler.stan",
	            data    = sir_data,
	            chains  = chains,
	            iter    = iter,
	            warmup  = warmup,
	            thin    = thin,
	            verbose = verbose,
	            refresh = refresh )
        	)

# traceplots
if (stan_options$chains > 1) {

	exfit <- extract(fit, permuted = FALSE, inc_warmup = FALSE)
	plotdata <- melt(exfit[,,'R0'])
	tracefitR0 <- ggplot() +
	              geom_line(data = plotdata,
	                        aes(x = iterations,
	                        y = value,
	                        color = factor(chains, labels = 1:stan_options$chains))) +
	              labs(x = "Sample", y = expression(R[0]), color = "Chain") +
	              scale_color_brewer(palette="Paired") +
	              theme(panel.background = element_rect(fill = "#F0F0F0"))

	quartz()
	print(tracefitR0)
	#ggsave(tracefitR0, filename="traceplotR0.pdf", height=4, width=6.5)

	exfit <- extract(fit, permuted = FALSE, inc_warmup = TRUE)
	plotdata <- melt(exfit[,,2])
	tracefitR0 <- ggplot() +
	              geom_line(data = plotdata,
	                        aes(x = iterations,
	                        y = value,
	                        color = factor(chains, labels = 1:stan_options$chains))) +
	              labs(x = "Sample", y = expression(R[0]), color = "Chain") +
	              scale_color_brewer(palette="Paired") +
	              theme(panel.background = element_rect(fill = "#F0F0F0"))

	quartz()
	print(tracefitR0)
	#ggsave(tracefitR0, filename="traceplotR0_inc.pdf", height=4, width=6.5)

}

exfit <- extract(fit, permuted = FALSE, inc_warmup = FALSE)
paramdata <- data.frame(R0 = melt(exfit[,,'R0'])$value,
               			r = melt(exfit[,,'r'])$value,
               			sigma = melt(exfit[,,'sigma'])$value,
               			eta = melt(exfit[,,'eta'])$value,
               			berr = melt(exfit[,,'berr'])$value,
               			Sinit = melt(exfit[,,4])$value,
               			Iinit = melt(exfit[,,5])$value,
               			Rinit = melt(exfit[,,6])$value )

# sample from parameter distributions
sdeout_true <- sdeout

nTraj 	<- 100
datlen 	<- dim(paramdata)[1]
inds 	<- sample.int(datlen,nTraj,replace = TRUE)
params 	<- paramdata[inds,]

bootstrapdata <- matrix(NA, nrow = nTraj, ncol = T+1)

for (i in 1:nTraj) {

	init_cond <- c(S = params$Sinit[i],
	               I = params$Iinit[i],
	               R = params$Rinit[i])
	pars <- c(R0 = params$R0[i],
	          r = params$r[i],
	          N = 500.0,
	          eta = params$eta[i],
	          berr = params$berr[i])

	sdeout <- StocSIR(init_cond, pars, T, steps)
	colnames(sdeout) <- c('S','I','R','B')

	bootstrapdata[i,] <- sdeout[,'I']

}

meanTraj 	<- colMeans(bootstrapdata)
quantTraj 	<- apply(bootstrapdata, 2, quantile, probs = c(0.025,0.975))

plotdata <- data.frame(times=0:T,true=sdeout_true[,'I'],est=meanTraj,quants=t(quantTraj),datapart=datapart)

g <- ggplot(plotdata, aes(times)) +
		geom_ribbon(aes(ymin = quants.2.5., ymax=quants.97.5.), alpha=0.1) +
        geom_line(aes(y = true, colour = "True")) + 
        geom_line(aes(y = est, colour = "IF2")) +
        geom_point(aes(y = datapart, color = "Data")) +
        labs(x = "Time", y = "Infection count", color = "") +
        scale_color_brewer(palette="Paired") +
        theme(panel.background = element_rect(fill = "#F0F0F0"))

quartz()
print(g)

# kernel plots

f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
kcolours <- f("Paired")


trueval.R0 		<- 3.0
trueval.r 		<- 0.1
trueval.sigma 	<- 10.0
trueval.Iinit 	<- 5
trueval.eta 	<- 0.5
trueval.berr 	<- 0.5

meanval.R0 		<- mean(paramdata$R0)
meanval.r 		<- mean(paramdata$r)
meanval.sigma 	<- mean(paramdata$sigma)
meanval.Iinit 	<- mean(paramdata$Iinit)
meanval.eta 	<- mean(paramdata$eta)
meanval.berr 	<- mean(paramdata$berr)


linecolour <- "grey50"
lineweight <- 0.5

kerdataR0 <- data.frame(R0points = paramdata$R0)
R0kernel <- ggplot(kerdataR0, aes(x = R0points, y = ..scaled..)) +
                geom_density(color = kcolours[1], fill = kcolours[1]) +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = expression(R[0]), y = "Density", color = "") +
                geom_vline(aes(xintercept=trueval.R0), linetype="solid", size=lineweight, color=linecolour) +
                geom_vline(aes(xintercept=meanval.R0), linetype="dashed", size=lineweight, color=linecolour)

kerdatar <- data.frame(rpoints = paramdata$r)
rkernel <- ggplot(kerdatar, aes(x = rpoints, y = ..scaled..)) +
                geom_density(color = kcolours[2], fill = kcolours[2]) +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = "r", y = "", color = "") +
                geom_vline(aes(xintercept=trueval.r), linetype="solid", size=lineweight, color=linecolour) +
                geom_vline(aes(xintercept=meanval.r), linetype="dashed", size=lineweight, color=linecolour)

kerdatasigma <- data.frame(sigmapoints = paramdata$sigma)
sigmakernel <- ggplot(kerdatasigma, aes(x = sigmapoints, y = ..scaled..)) +
                geom_density(color = kcolours[3], fill = kcolours[3]) +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = expression(sigma), y = "Density", color = "") +
                geom_vline(aes(xintercept=trueval.sigma), linetype="solid", size=lineweight, color=linecolour) +
                geom_vline(aes(xintercept=meanval.sigma), linetype="dashed", size=lineweight, color=linecolour)

kerdatainfec <- data.frame(infecpoints = paramdata$Iinit)
infeckernel <- ggplot(kerdatainfec, aes(x = infecpoints, y = ..scaled..)) +
                geom_density(color = kcolours[4], fill = kcolours[4]) +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = "Initial infected", y = "", color = "") +
                geom_vline(aes(xintercept=trueval.Iinit), linetype="solid", size=lineweight, color=linecolour) +
                geom_vline(aes(xintercept=meanval.Iinit), linetype="dashed", size=lineweight, color=linecolour)

kerdataeta <- data.frame(etapoints = paramdata$eta)
etakernel <- ggplot(kerdataeta, aes(x = etapoints, y = ..scaled..)) +
                geom_density(color = kcolours[5], fill = kcolours[5]) +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = expression(eta), y = "", color = "") +
                geom_vline(aes(xintercept=trueval.eta), linetype="solid", size=lineweight, color=linecolour) +
                geom_vline(aes(xintercept=meanval.eta), linetype="dashed", size=lineweight, color=linecolour)

kerdataberr <- data.frame(berrpoints = paramdata$berr)
berrkernel <- ggplot(kerdataberr, aes(x = berrpoints, y = ..scaled..)) +
                geom_density(color = kcolours[6], fill = kcolours[6]) +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = expression(beta[err]), y = "", color = "") +
                geom_vline(aes(xintercept=trueval.berr), linetype="solid", size=lineweight, color=linecolour) +
                geom_vline(aes(xintercept=meanval.berr), linetype="dashed", size=lineweight, color=linecolour)

# show grid
quartz()
grid.arrange(R0kernel, rkernel, sigmakernel, infeckernel, etakernel, berrkernel, ncol = 3, nrow = 2)

# Shiny Stan!!!

sso <- as.shinystan(fit)
sso <- launch_shinystan(sso)