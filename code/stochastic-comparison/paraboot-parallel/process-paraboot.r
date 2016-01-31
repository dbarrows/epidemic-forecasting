library(ggplot2)

## load laster data file
load("if2-paraboot.RData")

## this has pretty much everything
parabootdata <- data.frame(trajectories)

## get param star data
paramnames <- c("R0", "r", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")
paramstar <- parabootdata[,paramnames]

## get new infection count data, mean, and 95 quantiles
countnames <- paste("counts",1:(dim(parabootdata)[2]-9), sep = "")
countdata <- parabootdata[,countnames]
countmeans <- colMeans(countdata)
countquants <- t(apply(countdata, 2, quantile, probs = c(0.025,0.975)))
colnames(countquants) <- c("025","975")

## get expected observed data, means, and 95 quantiles
obsdata <- matrix(NA, nrow = dim(countdata)[1], ncol = dim(countdata)[2])
for (i in 1:dim(countdata)[1]) {
	obsdataraw <- countdata[i,] + rnorm(dim(countdata)[2], 0, paramstar[i,'sigma'])
	obsdata[i,] <- unlist(ifelse(obsdataraw < 0, 0, obsdataraw))
}
obsquants <- t(apply(obsdata, 2, quantile, probs = c(0.025,0.975)))
colnames(obsquants) <- c("025","975")


qplot(0:T, sdeout_true[,'I'], geom = "line", xlab = "Time", ylab = "Infection count") +
	geom_point( aes( y = c(infec_counts[1:(Tlim+1)], rep(NA,T-Tlim)) ) ) +
	geom_line(aes(y = c(rep(NA,Tlim+1), countmeans)), linetype = "dotted") +
	geom_ribbon(aes(ymin = c(rep(NA,Tlim+1), countquants[,'025']),
        	    	ymax = c(rep(NA,Tlim+1), countquants[,'975'])),
					alpha=0.15) +
	geom_ribbon(aes(ymin = c(rep(NA,Tlim+1), obsquants[,'025']),
	            	ymax = c(rep(NA,Tlim+1), obsquants[,'975'])),
					alpha=0.1) +
	theme_bw()