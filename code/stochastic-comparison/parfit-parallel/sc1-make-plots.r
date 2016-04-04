library(ggplot2)
library(Rcpp)
library(gridExtra)
library(rstan)
library(reshape2)

set.seed(112358)

printvar <- function(v) {
  	name <- deparse(substitute(v))
  	print(paste(name, ":", v))
}


## External files
##

stocsirfile 	<- paste(getwd(), "../../sir-functions", "StocSIR.r", sep = "/")
stocsirstanfile <- paste(getwd(), "../../sir-functions", "StocSIRstan.r", sep = "/")
if2file 		<- paste(getwd(),"../../if2", "if2.cpp", sep="/")
hmcfile 		<- paste(getwd(), "../../hmc", "sirode_euler.stan", sep="/")
source(stocsirfile)
source(stocsirstanfile)


## Beta process plot
###########################################################################################

T <- 7*60
berr <- 0.5
B <- numeric(T)
B0 <- 3.0 * 0.1 / 500;
B[1] <- B0
eta <- 0.5

for (t in 2:T) {
	B[t] <- exp( log(B[t-1]) + eta*(log(B0) - log(B[t-1])) + rnorm(1,0,berr) )
}

bplot <- qplot(1:T, B, geom = "line", xlab = "Time", ylab = expression(beta)) +
			scale_y_continuous(limits=c(0, 0.002)) +
			theme_bw()

ggsave(bplot, filename = "betaplot.pdf", height = 4, width = 6.5)

bdensity <- qplot(B, geom = "density", xlab = expression(beta), ylab = "frequency") +
				scale_x_continuous(limits=c(0, 0.002)) +
				theme_bw()

ggsave(bdensity, filename = "betadensity.pdf", height = 4, width = 6.5)

###########################################################################################

## print Beta theoretical stats

mu <- log(3*0.1/500)
sigma <- 0.5
mean <- exp(mu + sigma^2/2)
printvar(mean)

var <- sqrt( sigma^2 / (1-(1-eta)^2) )
sd <- sqrt( (exp(var^2) - 1) * exp(2*mu + var^2) )
printvar(sd)

## print Beta computed stats

printvar(mean(B))
printvar(sd(B))


## Model behaviour plot
###########################################################################################

set.seed(1004)

T 		<- 60
i_infec <- 5
steps 	<- 7
N 		<- 500
sigma 	<- 10
Tlim 	<- T

## Generate true trajecory and synthetic data

pars_true <- c(R0 = 3.0, 	# new infected people per infected person
          r = 0.1, 		# recovery rate
		  N = 500, 		# population size
		  eta = 0.5, 	# geometric random walk
		  berr = 0.5) 	# Beta geometric walk noise

true_init_cond <- c(S = N - i_infec,
					I = i_infec,
					R = 0)

sdeout_true <- StocSIR(true_init_cond, pars_true, T, steps)
colnames(sdeout_true) <- c('S','I','R','B')

infec_counts_raw <- sdeout_true[,'I'] + rnorm(T+1, 0, sigma)
infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

datapart <- c(infec_counts[1:(Tlim+1)],rep(NA,T-Tlim))

numtraj <- 100

trajData <- matrix(NA, nrow = numtraj, ncol = (T+1))
for(i in 1:numtraj) {
    sdeout <- StocSIR(true_init_cond, pars_true, T, steps)
    colnames(sdeout) <- c('S','I','R','B')
    trajData[i,] <- sdeout[,'I']
}

meanTraj    <- colMeans(trajData)
avTraj      <- meanTraj
quantTraj   <- t(apply(trajData, 2, quantile, probs = c(0.025,0.975)))
colnames(quantTraj) <- c("025","975")

sirmean <- qplot(0:T, sdeout_true[,'I'], geom = "line", xlab = "Time", ylab = "Infection count") +
		    geom_line(aes(x = 0:T, y = meanTraj), linetype = "dotted") +
		    geom_ribbon(aes(ymin = quantTraj[,'025'], ymax=quantTraj[,'975']), alpha=0.1) +
		    geom_point(aes(y = datapart)) +
		    theme_bw()

ggsave(sirmean, filename = "sirmean.pdf", width = 6.5, height = 4)



## IF2 fit
###########################################################################################

## Setup

NP          <- 2500
nPasses     <- 50
coolrate    <- 0.975

## fit

sourceCpp(if2file)
if2time <- system.time( if2data <- if2(infec_counts[1:(Tlim+1)], Tlim+1, N, NP, nPasses, coolrate) )

## extract

paramdata <- data.frame( if2data$paramdata )
names(paramdata) <- c("R0", "r", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")

convergedata <- data.frame( if2data$means )
names(convergedata) <- c("R0", "r", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")

statedata <- data.frame( if2data$statemeans )
names(statedata) <- c("S","I","R")

sdsdata <- data.frame( if2data$sds )
names(sdsdata) <- c("R0", "r", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")

## compute stats

if2means <- colMeans(paramdata[c("R0","r","Iinit","sigma","eta","berr")])
names(if2means) <- c("R0","r","I0","sigma","eta","berr")
printvar(if2means)

truevals <- c( pars_true[c("R0","r")], I0 = i_infec, sigma = sigma, pars_true[c("eta","berr")])
relerrs <- (if2means - truevals) / truevals
printvar(relerrs)

## fitted state plot

if2state <- qplot(0:T, sdeout_true[,'I'], geom = "line", xlab = "Time", ylab = "Infection counts") +
				geom_point(aes(y = infec_counts)) +
				geom_line(aes(y = statedata[,'I']), linetype = "dashed") +
				theme_bw()

ggsave(if2state, filename = "if2state.pdf", width = 6.5, height = 4)


## IF2 Convergence plots
#############################################################################################################

## Parameter values
##

meanval.R0 		<- mean(paramdata$R0)
meanval.r 		<- mean(paramdata$r)
meanval.sigma 	<- mean(paramdata$sigma)
meanval.Iinit 	<- mean(paramdata$Iinit)
meanval.eta 	<- mean(paramdata$eta)
meanval.berr 	<- mean(paramdata$berr)

linecolour <- "grey50"
lineweight <- 0.5

R0converge <- qplot(1:dim(convergedata)[1], convergedata$R0, geom = "line", xlab = "", ylab = expression(R[0])) +
    			geom_hline(aes(yintercept=pars_true['R0']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_hline(aes(yintercept=meanval.R0), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

rconverge <- qplot(1:dim(convergedata)[1], convergedata$r, geom = "line", xlab = "", ylab = "r") +
    			geom_hline(aes(yintercept=pars_true['r']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_hline(aes(yintercept=meanval.r), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

sigmaconverge <- qplot(1:dim(convergedata)[1], convergedata$sigma, geom = "line", xlab = "", ylab = expression(sigma)) +
    			geom_hline(aes(yintercept=sigma), linetype="solid", size=lineweight, color=linecolour) +
    			geom_hline(aes(yintercept=meanval.sigma), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

infecconverge <- qplot(1:dim(convergedata)[1], convergedata$Iinit, geom = "line", xlab = "", ylab = "Initial Infected") +
    			geom_hline(aes(yintercept=i_infec), linetype="solid", size=lineweight, color=linecolour) +
    			geom_hline(aes(yintercept=meanval.Iinit), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

etaconverge <- qplot(1:dim(convergedata)[1], convergedata$eta, geom = "line", xlab = "", ylab = expression(eta)) +
    			geom_hline(aes(yintercept=pars_true['eta']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_hline(aes(yintercept=meanval.eta), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

berrconverge <- qplot(1:dim(convergedata)[1], convergedata$berr, geom = "line", xlab = "", ylab = expression(epsilon[proc])) +
    			geom_hline(aes(yintercept=pars_true['berr']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_hline(aes(yintercept=meanval.berr), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

pdf("if2convergence.pdf", width = 6.5, height = 4)
grid.arrange(R0converge, rconverge, sigmaconverge, infecconverge, etaconverge, berrconverge, ncol = 3, nrow = 2)
dev.off()

## Standard devs
##

meanval.R0 		<- mean(paramdata$R0)
meanval.r 		<- mean(paramdata$r)
meanval.sigma 	<- mean(paramdata$sigma)
meanval.Iinit 	<- mean(paramdata$Iinit)
meanval.eta 	<- mean(paramdata$eta)
meanval.berr 	<- mean(paramdata$berr)

linecolour <- "grey50"
lineweight <- 0.5

R0sds <- qplot(1:dim(sdsdata)[1], sdsdata$R0, geom = "line", xlab = "", ylab = expression(R[0])) +
				theme_bw()

rsds <- qplot(1:dim(sdsdata)[1], sdsdata$r, geom = "line", xlab = "", ylab = "r") +
				theme_bw()

sigmasds <- qplot(1:dim(sdsdata)[1], sdsdata$sigma, geom = "line", xlab = "", ylab = expression(sigma)) +
				theme_bw()

infecsds <- qplot(1:dim(sdsdata)[1], sdsdata$Iinit, geom = "line", xlab = "", ylab = "Initial Infected") +
				theme_bw()

etasds <- qplot(1:dim(sdsdata)[1], sdsdata$eta, geom = "line", xlab = "", ylab = expression(eta)) +
				theme_bw()

berrsds <- qplot(1:dim(sdsdata)[1], sdsdata$berr, geom = "line", xlab = "", ylab = expression(epsilon[proc])) +
				theme_bw()

pdf("if2sdconvergence.pdf", width = 6.5, height = 4)
grid.arrange(R0sds, rsds, sigmasds, infecsds, etasds, berrsds, ncol = 3, nrow = 2)
dev.off()

## IF2 densities

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

rkernel <- qplot(paramdata$r, geom = "density", xlab = "r", ylab = "") +
    			geom_vline(aes(xintercept=pars_true['r']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_vline(aes(xintercept=meanval.r), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

sigmakernel <- qplot(paramdata$sigma, geom = "density", xlab = expression(sigma), ylab = "") +
    			geom_vline(aes(xintercept=sigma), linetype="solid", size=lineweight, color=linecolour) +
    			geom_vline(aes(xintercept=meanval.sigma), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

infeckernel <- qplot(paramdata$Iinit, geom = "density", xlab = "Initial Infected", ylab = "frequency") +
    			geom_vline(aes(xintercept=i_infec), linetype="solid", size=lineweight, color=linecolour) +
    			geom_vline(aes(xintercept=meanval.Iinit), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

etakernel <- qplot(paramdata$eta, geom = "density", xlab = expression(eta), ylab = "") +
    			geom_vline(aes(xintercept=pars_true['eta']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_vline(aes(xintercept=meanval.eta), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

berrkernel <- qplot(paramdata$berr, geom = "density", xlab = expression(epsilon[proc]), ylab = "") +
    			geom_vline(aes(xintercept=pars_true['berr']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_vline(aes(xintercept=meanval.berr), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()


# show grid
pdf("if2kernels.pdf", width = 6.5, height = 4)
grid.arrange(R0kernel, rkernel, sigmakernel, infeckernel, etakernel, berrkernel, ncol = 3, nrow = 2)
dev.off()



## HMCMC Fitting
#############################################################################################################

datlen <- T*7 + 1

data <- matrix(data = -1, nrow = T+1, ncol = steps)
data[,1] <- infec_counts
standata <- as.vector(t(data))[1:datlen]

sir_data <- list( T = datlen,   	# simulation time
                  y = standata, 	# infection count data
                  N = 500,      	# population size
                  h = 1/steps )   	# step size per day 
                    
rstan_options(auto_write = TRUE)

stan_options <- list(   chains = 1,    		# number of chains
                        iter   = 2000, 		# iterations per chain
                        warmup = 1000, 		# warmup interations
                        thin   = 1)   		# thinning number

#hmcfile <- paste(getwd(), "../../code/hmc", "sirode_euler.stan", sep="/")

hmctime <- system.time(fit <- with(stan_options,
		            	stan(file  	= hmcfile,
				            data    = sir_data,
				            chains  = chains,
				            iter    = iter,
				            warmup  = warmup,
				            thin    = thin)
		        		)
			)

exfit <- extract(fit, permuted = FALSE, inc_warmup = FALSE)
paramdata <- data.frame(R0 = melt(exfit[,,'R0'])$value,
               			r = melt(exfit[,,'r'])$value,
               			sigma = melt(exfit[,,'sigma'])$value,
               			eta = melt(exfit[,,'eta'])$value,
               			berr = melt(exfit[,,'berr'])$value,
               			Sinit = melt(exfit[,,'y0[1]'])$value,
               			Iinit = melt(exfit[,,'y0[2]'])$value,
               			Rinit = melt(exfit[,,'y0[3]'])$value )

for (j in 1:datlen) {
	varname <- paste('Bnoise[', j, ']', sep = '')
	paramdata[[varname]] <- melt( exfit[,,varname] )$value
}

## Compute stats

hmcmeans <- colMeans(paramdata[c("R0","r","Iinit","sigma","eta","berr")])
names(hmcmeans) <- c("R0","r","I0","sigma","eta","berr")
printvar(hmcmeans)

truevals <- c( pars_true[c("R0","r")], I0 = i_infec, sigma = sigma, pars_true[c("eta","berr")])
relerrs <- (hmcmeans - truevals) / truevals
printvar(relerrs)


## HMCMC Kernels

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

rkernel <- qplot(paramdata$r, geom = "density", xlab = "r", ylab = "") +
    			geom_vline(aes(xintercept=pars_true['r']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_vline(aes(xintercept=meanval.r), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

sigmakernel <- qplot(paramdata$sigma, geom = "density", xlab = expression(sigma), ylab = "") +
    			geom_vline(aes(xintercept=sigma), linetype="solid", size=lineweight, color=linecolour) +
    			geom_vline(aes(xintercept=meanval.sigma), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

infeckernel <- qplot(paramdata$Iinit, geom = "density", xlab = "Initial Infected", ylab = "frequency") +
    			geom_vline(aes(xintercept=i_infec), linetype="solid", size=lineweight, color=linecolour) +
    			geom_vline(aes(xintercept=meanval.Iinit), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

etakernel <- qplot(paramdata$eta, geom = "density", xlab = expression(eta), ylab = "") +
    			geom_vline(aes(xintercept=pars_true['eta']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_vline(aes(xintercept=meanval.eta), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

berrkernel <- qplot(paramdata$berr, geom = "density", xlab = expression(epsilon[proc]), ylab = "") +
    			geom_vline(aes(xintercept=pars_true['berr']), linetype="solid", size=lineweight, color=linecolour) +
    			geom_vline(aes(xintercept=meanval.berr), linetype="dashed", size=lineweight, color=linecolour) +
				theme_bw()

pdf("hmckernels.pdf", width = 6.5, height = 4)
grid.arrange(R0kernel, rkernel, sigmakernel, infeckernel, etakernel, berrkernel, ncol = 3, nrow = 2)
dev.off()


## HMCMC bootstrapping

nTraj <- 100
# sample from parameter distributions

pardatlen 	<- dim(paramdata)[1]
inds 	    <- sample.int(pardatlen,nTraj,replace = TRUE)
params 	    <- paramdata[inds,]

bootstrapdata <- matrix(NA, nrow = nTraj, ncol = T+1)

for (i in 1:nTraj) {

	paramset <- params[i,]

	init_cond <- c(S = paramset$Sinit,
	               I = paramset$Iinit,
	               R = paramset$Rinit)
	pars <- c(R0 = paramset$R0,
	          r = paramset$r,
	          N = 500.0,
	          eta = paramset$eta,
	          berr = paramset$berr)

	berrvec <- numeric(datlen)
	for (j in 1:datlen) {
		varname <- paste("Bnoise[", j, "]", sep = "")
		berrvec[j] <- paramset[[varname]]
	}

	sdeout <- StocSIRstan(init_cond, pars, T, steps, berrvec, datlen)
	colnames(sdeout) <- c('S','I','R','B')

	bootstrapdata[i,] <- sdeout[,'I']

}

# remove NaN rows (why in the hell are these showing up????)
bootstrapdata <- bootstrapdata[complete.cases(bootstrapdata),]

meanTraj 	<- colMeans(bootstrapdata, na.rm = FALSE)
quantTraj 	<- t( apply(bootstrapdata, 2, quantile, probs = c(0.025,0.975), na.rm = FALSE) )

colnames(quantTraj) <- c("025","975")

hmcboot <- qplot(0:T, sdeout_true[,'I'], geom = "line", xlab = "Time", ylab = "Infection count") +
			geom_ribbon(aes(ymin = quantTraj[,'025'], ymax=quantTraj[,'975']), alpha=0.1) +
			geom_line(aes(y = meanTraj), linetype = "dashed") +
			geom_point(aes(y = datapart)) +
		    geom_line(aes(y = avTraj), linetype = "dotted") +
			theme_bw()

ggsave(hmcboot, filename = "hmcboot.pdf", width = 6.5, height = 4)




## Multi-fit densities
######################################################################################################

load("sc1-multi-time.RData")

# extract if2 data

if2names <- c("if2.R0","if2.r","if2.sigma","if2.eta","if2.berr","if2.Iinit")
if2data <- estmat[,if2names]
colnames(if2data) <- c("R0","r","sigma","eta","berr","Iinit")
if2times <- estmat[,7]

# extract hmc data

hmcnames <- c("hmc.R0","hmc.r","hmc.sigma","hmc.eta","hmc.berr","hmc.Iinit")
hmcdata <- estmat[,hmcnames]
colnames(hmcdata) <- c("R0","r","sigma","eta","berr","Iinit")
hmctimes <- estmat[,14]


# average times

avif2time <- mean(if2times)
avhmctime <- mean(hmctimes)

# times plot

timedf <- data.frame(if2 = if2times, hmc = hmctimes)
#timedf$variable <- factor(timedf$variable, levels = c("if2", "hmc"), ordered = TRUE)
timedata <- melt(timedf)
timeplot <- ggplot(timedata, aes(factor(variable, ordered = TRUE), value)) +
				geom_boxplot() +
				scale_x_discrete(labels = c("IF2", "HMCMC")  ) +
				labs(x = "", y = "Time (seconds)") +
				coord_flip() +
				theme_bw()

ggsave(timeplot, filename = "timeplot.pdf", width = 6.5, height = 4)


# sort results

printvar(avif2time)
printvar(avhmctime)

if2sorted <- apply(if2data, 2, sort)
hmcsorted <- apply(hmcdata, 2, sort)

## take centre 95%

nTrials <- dim(estmat)[1]

cinds <- (0.025*nTrials+1):(0.975*nTrials)
if295 <- if2sorted[cinds,]
hmc95 <- hmcsorted[cinds,]

if2plotdata <- data.frame(if295)
hmcplotdata <- data.frame(hmc95)

linecolour <- "grey50"
lineweight <- 0.5


## IF2 multi densities

R0kernel <- qplot(if2plotdata$R0, geom = "density", xlab = expression(R[0]), ylab = "frequency") +
    geom_vline(aes(xintercept=pars_true['R0']), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=median(if2plotdata$R0)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

rkernel <- qplot(if2plotdata$r, geom = "density", xlab = "r", ylab = "") +
    geom_vline(aes(xintercept=pars_true['r']), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=median(if2plotdata$r)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

sigmakernel <- qplot(if2plotdata$sigma, geom = "density", xlab = expression(sigma), ylab = "") +
    geom_vline(aes(xintercept=sigma), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=median(if2plotdata$sigma)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

infeckernel <- qplot(if2plotdata$Iinit, geom = "density", xlab = "Initial Infected", ylab = "frequency") +
    geom_vline(aes(xintercept=i_infec), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=median(if2plotdata$Iinit)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

etakernel <- qplot(if2plotdata$eta, geom = "density", xlab = expression(eta), ylab = "") +
    geom_vline(aes(xintercept=pars_true['eta']), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=median(if2plotdata$eta)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

berrkernel <- qplot(if2plotdata$berr, geom = "density", xlab = expression(epsilon[proc]), ylab = "") +
    geom_vline(aes(xintercept=pars_true['berr']), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=median(if2plotdata$berr)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()


# show grid
pdf("if2-multi.pdf", width = 6.5, height = 4)
grid.arrange(R0kernel, rkernel, sigmakernel, infeckernel, etakernel, berrkernel, ncol = 3, nrow = 2)
dev.off()

## HMCMC multi densities

R0kernel <- qplot(hmcplotdata$R0, geom = "density", xlab = expression(R[0]), ylab = "frequency") +
    geom_vline(aes(xintercept=pars_true['R0']), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=median(hmcplotdata$R0)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

rkernel <- qplot(hmcplotdata$r, geom = "density", xlab = "r", ylab = "") +
    geom_vline(aes(xintercept=pars_true['r']), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=median(hmcplotdata$r)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

sigmakernel <- qplot(hmcplotdata$sigma, geom = "density", xlab = expression(sigma), ylab = "") +
    geom_vline(aes(xintercept=sigma), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=median(hmcplotdata$sigma)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

infeckernel <- qplot(hmcplotdata$Iinit, geom = "density", xlab = "Initial Infected", ylab = "frequency") +
    geom_vline(aes(xintercept=i_infec), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=median(hmcplotdata$Iinit)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

etakernel <- qplot(hmcplotdata$eta, geom = "density", xlab = expression(eta), ylab = "") +
    geom_vline(aes(xintercept=pars_true['eta']), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=median(hmcplotdata$eta)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

berrkernel <- qplot(hmcplotdata$berr, geom = "density", xlab = expression(epsilon[proc]), ylab = "") +
    geom_vline(aes(xintercept=pars_true['berr']), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=median(hmcplotdata$berr)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()


# show grid
pdf("hmc-multi.pdf", width = 6.5, height = 4)
grid.arrange(R0kernel, rkernel, sigmakernel, infeckernel, etakernel, berrkernel, ncol = 3, nrow = 2)
dev.off()


## Combined plots

R0kernel <- qplot(if2plotdata$R0, geom = "density", xlab = expression(R[0]), ylab = "frequency") +
    geom_density(aes(x = hmcplotdata$R0), color = "grey") +
    geom_vline(aes(xintercept=pars_true['R0']), linetype="solid", size=lineweight) +
    theme_bw()

rkernel <- qplot(if2plotdata$r, geom = "density", xlab = "r", ylab = "") +
    geom_density(aes(x = hmcplotdata$r), color = "grey") +
    geom_vline(aes(xintercept=pars_true['r']), linetype="solid", size=lineweight, color=linecolour) +
    theme_bw()

sigmakernel <- qplot(if2plotdata$sigma, geom = "density", xlab = expression(sigma), ylab = "") +
    geom_density(aes(x = hmcplotdata$sigma), color = "grey") +
    geom_vline(aes(xintercept=sigma), linetype="solid", size=lineweight, color=linecolour) +
    theme_bw()

infeckernel <- qplot(if2plotdata$Iinit, geom = "density", xlab = "Initial Infected", ylab = "frequency") +
    geom_density(aes(x = hmcplotdata$Iinit), color = "grey") +
    geom_vline(aes(xintercept=i_infec), linetype="solid", size=lineweight, color=linecolour) +
    theme_bw()

etakernel <- qplot(if2plotdata$eta, geom = "density", xlab = expression(eta), ylab = "") +
    geom_density(aes(x = hmcplotdata$eta), color = "grey") +
    geom_vline(aes(xintercept=pars_true['eta']), linetype="solid", size=lineweight, color=linecolour) +
    theme_bw()

berrkernel <- qplot(if2plotdata$berr, geom = "density", xlab = expression(epsilon[proc]), ylab = "") +
    geom_density(aes(x = hmcplotdata$berr), color = "grey") +
    geom_vline(aes(xintercept=pars_true['berr']), linetype="solid", size=lineweight, color=linecolour) +
    theme_bw()

pdf("combined-multi.pdf", width = 6.5, height = 4)
grid.arrange(R0kernel, rkernel, sigmakernel, infeckernel, etakernel, berrkernel, ncol = 3, nrow = 2)
dev.off()