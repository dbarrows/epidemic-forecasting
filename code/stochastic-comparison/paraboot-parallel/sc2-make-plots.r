library(ggplot2)
library(reshape2)
library(rstan)

printvar <- function(v) {
  	name <- deparse(substitute(v))
  	print(paste(name, ":", v))
}

load(file = "sc2-script.RData")

q <- qplot(0:T, sdeout_true[,'I'], geom = "line", xlab = "Time", ylab = "Infection count") +
		geom_point(aes(y = c(infec_counts[1:(Tlim+1)], rep(NA,T-Tlim)))) +
		theme_bw()

ggsave(q, filename="dataplot.pdf", height=4, width=6.5)


statedata <- data.frame( if2data$statemeans )
names(statedata) <- c("S","I","R")

statepart <- c(statedata[,'I'],rep(NA,T-Tlim))
qif2 <- qplot(0:T, sdeout_true[,'I'], geom = "line", xlab = "Time", ylab = "Infection counts") +
			geom_point(aes(y = datapart)) +
			geom_line(aes(y = statepart), linetype = "dashed") +
			theme_bw()

ggsave(qif2, filename="if2fit.pdf", height=4, width=6.5)

## Without observation noise
##

pdata <- data.frame(parabootdata)
pforecast <- pdata[,paste("counts", 2:(T-Tlim+1), sep = "")]

meanTraj 	<- colMeans(pforecast)
quantTraj 	<- t(apply(pforecast, 2, quantile, probs = c(0.025,0.975)))
colnames(quantTraj) <- c("025","975")

meanTrajpart 	<- c(rep(NA, Tlim+1), meanTraj)
quantTrajpart 	<- rbind(matrix(NA, nrow = Tlim+1, ncol = 2), quantTraj)

qif2forecast <- qplot(0:T, sdeout_true[,'I'], geom = "line", xlab = "Time", ylab = "Infection count") +
					geom_ribbon(aes(ymin = quantTrajpart[,'025'], ymax=quantTrajpart[,'975']), alpha=0.1) +
					geom_line(aes(y = meanTrajpart), linetype = "dotted") +
    				geom_point(aes(y = datapart)) +
					geom_line(aes(y = statepart), linetype = "dashed") +
					theme_bw()

ggsave(qif2forecast, filename="if2forecast.pdf", height=4, width=6.5)

truefuture  <- sdeout_true[(Tlim+2):(T+1),'I']
estfuture   <- meanTraj
err <- estfuture - truefuture
if2sse <- sum(err^2)
printvar(if2sse)

## With observation noise
##

pdata <- data.frame(parabootdata)
pforecast <- pdata[,paste("counts", 2:(T-Tlim+1), sep = "")]
psigma <- pdata[,"sigma"]

#pcounts <- sapply(1:nTraj, function(i) pforecast[i,1:(T-Tlim)] + rnorm(T-Tlim, 0, psigma[i]))
#pcounts <- ifelse(pcounts < 0, 0, pcounts)
noise <- t(as.matrix(sapply(psigma, function(sigma) rnorm(T-Tlim, 0, sigma))))
pcounts <- pforecast + noise
pcounts[as.matrix(pcounts) < 0] <- 0

meanTraj 	<- colMeans(pcounts)
quantTraj 	<- t(apply(pcounts, 2, quantile, probs = c(0.025,0.975)))
colnames(quantTraj) <- c("025","975")

meanTrajpart 	<- c(rep(NA, Tlim+1), meanTraj)
quantTrajpart 	<- rbind(matrix(NA, nrow = Tlim+1, ncol = 2), quantTraj)

qif2forecast_c <- qplot(0:T, sdeout_true[,'I'], geom = "line", xlab = "Time", ylab = "Infection count") +
					geom_ribbon(aes(ymin = quantTrajpart[,'025'], ymax=quantTrajpart[,'975']), alpha=0.1) +
					geom_line(aes(y = meanTrajpart), linetype = "dotted") +
    				geom_point(aes(y = datapart)) +
					geom_line(aes(y = statepart), linetype = "dashed") +
					theme_bw()

ggsave(qif2forecast_c, filename="if2forecast_c.pdf", height=4, width=6.5)

## Combined

pdata <- data.frame(parabootdata)
pforecast <- pdata[,paste("counts", 2:(T-Tlim+1), sep = "")]

meanTraj 	<- colMeans(pforecast)
quantTraj 	<- t(apply(pforecast, 2, quantile, probs = c(0.025,0.975)))
colnames(quantTraj) <- c("025","975")

meanTrajpart_n 	<- c(rep(NA, Tlim+1), meanTraj)
quantTrajpart_n 	<- rbind(matrix(NA, nrow = Tlim+1, ncol = 2), quantTraj)

qif2combined <- qplot(0:T, sdeout_true[,'I'], geom = "line", xlab = "Time", ylab = "Infection count") +
					geom_ribbon(aes(ymin = quantTrajpart[,'025'], ymax=quantTrajpart[,'975']), alpha=0.1) +
					geom_ribbon(aes(ymin = quantTrajpart_n[,'025'], ymax=quantTrajpart_n[,'975']), alpha=0.2) +
					geom_line(aes(y = meanTrajpart_n), linetype = "dotted") +
    				geom_point(aes(y = datapart)) +
					geom_line(aes(y = statepart), linetype = "dashed") +
					theme_bw()

ggsave(qif2combined, filename="if2combined.pdf", height=4, width=6.5)



## HMCMC
##

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

# remove NaN rows
bootstrapdata <- bootstrapdata[complete.cases(bootstrapdata),]

meanTraj 	<- colMeans(bootstrapdata, na.rm = FALSE)
quantTraj 	<- t( apply(bootstrapdata, 2, quantile, probs = c(0.025,0.975), na.rm = FALSE) )
colnames(quantTraj) <- c("025","975")


qhmcforecast <- qplot(0:T, sdeout_true[,'I'], geom = "line", xlab = "Time", ylab = "Infection count") +
					geom_ribbon(aes(ymin = quantTraj[,'025'], ymax=quantTraj[,'975']), alpha=0.2) +
					geom_line(aes(y = meanTraj), linetype = "dashed") +
					geom_point(aes(y = datapart)) +
					theme_bw()

ggsave(qhmcforecast, filename="hmcforecast.pdf", height=4, width=6.5)

truefuture  <- sdeout_true[(Tlim+2):(T+1),'I']
estfuture   <- meanTraj[(Tlim+2):(T+1)]

err <- estfuture - truefuture
hmcsse <- sum(err^2)
printvar(hmcsse)