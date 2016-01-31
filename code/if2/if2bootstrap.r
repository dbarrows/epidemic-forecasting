library(deSolve)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(Rcpp)
library(RColorBrewer)

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

set.seed(50721456)

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
		  berr = 0.3) 	# Beta geometric walk noise

true_init_cond <- c(S = N - i_infec,
					I = i_infec,
					R = 0)

sdeout_true <- StocSIR(true_init_cond, pars_true, T, steps)
colnames(sdeout_true) <- c('S','I','R','B')

#quartz()
#plot(1:(T+1), sdeout[,'B'])

infec_counts_raw <- sdeout_true[,'I'] + rnorm(T+1, 0, sigma)
infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

plotdata <- data.frame(times = 0:T, true = sdeout_true[,'I'], data = infec_counts)

g <- ggplot(plotdata, aes(times)) +
        geom_line(aes(y = true, colour = "True")) + 
        geom_point(aes(y = data, colour = "Data")) +
        labs(x = "Time", y = "Infection count", color = "") +
        scale_color_brewer(palette="Paired") +
        theme(panel.background = element_rect(fill = "#F0F0F0"))

if ( Sys.info()['sysname'] == 'Darwin' )
	quartz()
print(g)

#ggsave(g, filename="dataplot.pdf", height=4, width=6.5)

## Rcpp stuff
##

NP <- 30000
nPasses <- 10
coolrate <- 0.9

sourceCpp(paste(getwd(),"if2.cpp",sep="/"))

if2data <- if2(infec_counts[1:(Tlim+1)], Tlim+1, N, NP, nPasses, coolrate)

paramdata <- data.frame( if2data$paramdata )
names(paramdata) <- c("R0", "r", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")

convergedata <- data.frame( if2data$means )
names(convergedata) <- c("R0", "r", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")

statedata <- data.frame( if2data$statemeans )
names(statedata) <- c("S","I","R")

## last particle filter means plot
if ( Sys.info()['sysname'] == 'Darwin' )
	quartz()
qplot(1:(T+1), sdeout_true[,'I'], geom = "line", xlab = "Time", ylab = "Infection counts") +
	geom_point(aes(y = infec_counts)) +
	geom_line(aes(y = statedata[,'I']), linetype = "dashed") +
	theme_bw()

## sample from parameter distributions
##

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
quantTraj 	<- t(apply(bootstrapdata, 2, quantile, probs = c(0.025,0.975)))
colnames(quantTraj) <- c("025","975")

datapart <- c(infec_counts[1:(Tlim+1)],rep(NA,T-Tlim))

g <- qplot(1:(T+1), sdeout_true[,'I'], geom = "line", xlab = "Time", ylab = "Infection count") +
        geom_ribbon(aes(ymin = quantTraj[,'025'], ymax=quantTraj[,'975']), alpha=0.1) +
        geom_line(aes(y = meanTraj), linetype = "dashed") +
        geom_point(aes(y = datapart)) +
        theme_bw()

if ( Sys.info()['sysname'] == 'Darwin' )
	quartz()
print(g)

#ggsave(g, filename="if2plot.pdf", height=4, width=6.5)

f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
kcolours <- f("Paired")


trueval.R0 		<- 3.0
trueval.r 		<- 0.1
trueval.sigma 	<- 10.0
trueval.Iinit 	<- 5
trueval.eta 	<- 0.5
trueval.berr 	<- 0.3

meanval.R0 		<- mean(convergedata$R0)
meanval.r 		<- mean(convergedata$r)
meanval.sigma 	<- mean(convergedata$sigma)
meanval.Iinit 	<- mean(convergedata$Iinit)
meanval.eta 	<- mean(convergedata$eta)
meanval.berr 	<- mean(convergedata$berr)


linecolour <- "grey50"
lineweight <- 0.5

## Convergence 'kernels'
if (FALSE) {
	R0kernel <- qplot(convergedata$R0, geom = "density", xlab = expression(R[0]), ylab = "") +
	    			geom_vline(aes(xintercept=pars_true['R0']), linetype="solid", size=lineweight, color=linecolour) +
	    			geom_vline(aes(xintercept=meanval.R0), linetype="dashed", size=lineweight, color=linecolour) +
					theme_bw()

	rkernel <- qplot(convergedata$r, geom = "density", xlab = "r", ylab = "") +
	    			geom_vline(aes(xintercept=pars_true['r']), linetype="solid", size=lineweight, color=linecolour) +
	    			geom_vline(aes(xintercept=meanval.r), linetype="dashed", size=lineweight, color=linecolour) +
					theme_bw()

	sigmakernel <- qplot(convergedata$sigma, geom = "density", xlab = expression(sigma), ylab = "") +
	    			geom_vline(aes(xintercept=sigma), linetype="solid", size=lineweight, color=linecolour) +
	    			geom_vline(aes(xintercept=meanval.sigma), linetype="dashed", size=lineweight, color=linecolour) +
					theme_bw()

	infeckernel <- qplot(convergedata$Iinit, geom = "density", xlab = "Initial Infected", ylab = "") +
	    			geom_vline(aes(xintercept=i_infec), linetype="solid", size=lineweight, color=linecolour) +
	    			geom_vline(aes(xintercept=meanval.Iinit), linetype="dashed", size=lineweight, color=linecolour) +
					theme_bw()

	etakernel <- qplot(convergedata$eta, geom = "density", xlab = expression(eta), ylab = "") +
	    			geom_vline(aes(xintercept=pars_true['eta']), linetype="solid", size=lineweight, color=linecolour) +
	    			geom_vline(aes(xintercept=meanval.eta), linetype="dashed", size=lineweight, color=linecolour) +
					theme_bw()

	berrkernel <- qplot(convergedata$berr, geom = "density", xlab = expression(epsilon[proc]), ylab = "") +
	    			geom_vline(aes(xintercept=pars_true['berr']), linetype="solid", size=lineweight, color=linecolour) +
	    			geom_vline(aes(xintercept=meanval.berr), linetype="dashed", size=lineweight, color=linecolour) +
					theme_bw()

	# show grid
	if ( Sys.info()['sysname'] == 'Darwin' )
		quartz()
	grid.arrange(R0kernel, rkernel, sigmakernel, infeckernel, etakernel, berrkernel,
	             ncol = 3, nrow = 2, top = "Convergence kernels")
}

# Convergence plots

R0converge <- qplot(1:dim(convergedata)[1], convergedata$R0, geom = "line", xlab = expression(R[0]), ylab = "frequency") +
    			geom_hline(aes(yintercept=pars_true['R0']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_hline(aes(yintercept=meanval.R0), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

rconverge <- qplot(1:dim(convergedata)[1], convergedata$r, geom = "line", xlab = "r", ylab = "") +
    			geom_hline(aes(yintercept=pars_true['r']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_hline(aes(yintercept=meanval.r), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

sigmaconverge <- qplot(1:dim(convergedata)[1], convergedata$sigma, geom = "line", xlab = expression(sigma), ylab = "") +
    			geom_hline(aes(yintercept=sigma), linetype="solid", size=lineweight, color=linecolour) +
    			geom_hline(aes(yintercept=meanval.sigma), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

infecconverge <- qplot(1:dim(convergedata)[1], convergedata$Iinit, geom = "line", xlab = "Initial Infected", ylab = "frequency") +
    			geom_hline(aes(yintercept=i_infec), linetype="solid", size=lineweight, color=linecolour) +
    			geom_hline(aes(yintercept=meanval.Iinit), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

etaconverge <- qplot(1:dim(convergedata)[1], convergedata$eta, geom = "line", xlab = expression(eta), ylab = "") +
    			geom_hline(aes(yintercept=pars_true['eta']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_hline(aes(yintercept=meanval.eta), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

berrconverge <- qplot(1:dim(convergedata)[1], convergedata$berr, geom = "line", xlab = expression(epsilon[proc]), ylab = "") +
    			geom_hline(aes(yintercept=pars_true['berr']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_hline(aes(yintercept=meanval.berr), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

# show grid
if ( Sys.info()['sysname'] == 'Darwin' )
	quartz()
grid.arrange(R0converge, rconverge, sigmaconverge, infecconverge, etaconverge, berrconverge,
             ncol = 3, nrow = 2, top = "Convergence plots")


# Final swarm kernel densities

meanval.R0 		<- mean(paramdata$R0)
meanval.r 		<- mean(paramdata$r)
meanval.sigma 	<- mean(paramdata$sigma)
meanval.Iinit 	<- mean(paramdata$Iinit)
meanval.eta 	<- mean(paramdata$eta)
meanval.berr 	<- mean(paramdata$berr)

linecolour <- "grey50"
lineweight <- 0.5

R0kernel <- qplot(paramdata$R0, geom = "density", xlab = expression(R[0]), ylab = "frequency") +
    			geom_vline(aes(xintercept=pars_true['R0']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_vline(aes(xintercept=meanval.R0), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

rkernel <- qplot(paramdata$r, geom = "density", xlab = "r") +
    			geom_vline(aes(xintercept=pars_true['r']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_vline(aes(xintercept=meanval.r), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

sigmakernel <- qplot(paramdata$sigma, geom = "density", xlab = expression(sigma)) +
    			geom_vline(aes(xintercept=sigma), linetype="solid", size=lineweight, color=linecolour) +
    			geom_vline(aes(xintercept=meanval.sigma), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

infeckernel <- qplot(paramdata$Iinit, geom = "density", xlab = "Initial Infected", ylab = "frequency") +
    			geom_vline(aes(xintercept=i_infec), linetype="solid", size=lineweight, color=linecolour) +
    			geom_vline(aes(xintercept=meanval.Iinit), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

etakernel <- qplot(paramdata$eta, geom = "density", xlab = expression(eta)) +
    			geom_vline(aes(xintercept=pars_true['eta']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_vline(aes(xintercept=meanval.eta), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

berrkernel <- qplot(paramdata$berr, geom = "density", xlab = expression(epsilon[proc])) +
    			geom_vline(aes(xintercept=pars_true['berr']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_vline(aes(xintercept=meanval.berr), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()


# show grid
grid.arrange(R0kernel, rkernel, sigmakernel, infeckernel, etakernel, berrkernel,
             ncol = 3, nrow = 2, top = "Swarm final density kernels")

#pdf("if2kernels.pdf", height = 6.5, width = 6.5)
#grid.arrange(R0kernel, rkernel, sigmakernel, infeckernel, ncol = 2, nrow = 2)
#dev.off()