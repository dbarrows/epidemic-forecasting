library(ggplot2)
library(reshape2)
library(RColorBrewer)

printvar <- function(v) {
  	name <- deparse(substitute(v))
  	print(paste(name, ":", v))
}

dir <- paste(getwd(), "../../../../Thesis_data","sirs-varimages", sep = "/")
filelist <- list.files(dir)

source("stanbootmean.r")

nTrials <- length(filelist)
first <- TRUE
fctr <- 0

for (filenum in 1:length(filelist)) {

	filename <- filelist[filenum]
	filepath <- paste(dir, filename, sep = "/")

	if (filename != "thread-0.RData" & filename != "varimages.old") {

		fctr <- fctr + 1
		print(paste(fctr, filename))

		e <- new.env()
		load(filepath, e)

		## get sizing info from first file
		if (first) {

			stepsAhead 	<- get("stepsAhead", e)
			T 		   	<- get("T", e)
			Tlim 	   	<- get("Tlim", e)
			nTraj 	   	<- get("nTraj", e)
			datlen 	   	<- get("datlen", e)
			steps 		<- get("steps", e)	

			truetraj <- matrix(NA, nrow = nTrials, ncol = stepsAhead)
			if2traj <- matrix(NA, nrow = nTrials, ncol = stepsAhead)
			hmctraj <- matrix(NA, nrow = nTrials, ncol = stepsAhead)
			smaptraj <- matrix(NA, nrow = nTrials, ncol = stepsAhead)

			hmctraj2 <- matrix(NA, nrow = nTrials, ncol = stepsAhead)

			if2times 	<- numeric(nTrials)
			hmctimes 	<- numeric(nTrials)
			smaptimes 	<- numeric(nTrials)

			first <- FALSE

		}

		## True trajectory
		trueproj <- get("sdeout_true", e)[(Tlim+2):(T+1),'I']
		truetraj[fctr,] <- trueproj

		## IF2
		parabootdata <- get("if2_paraboot_data", e)[,paste("counts", 2:(T-Tlim+1), sep = "")]
		parabootdata <- parabootdata[complete.cases(parabootdata),]
		if2proj <- apply(parabootdata, 2, mean)
		if2traj[fctr,] <- if2proj
		if2times[fctr] <- get("if2time", e)[['user.self']]

		## HMCMC
		hmcbootdata <- get("bootstrapdata", e)[,(Tlim+2):(T+1)]
		hmcproj <- apply(hmcbootdata, 2, mean)
		hmctraj[fctr,] <- hmcproj
		hmctimes[fctr] <- get("hmctime", e)[['user.self']]
		# test
		stanfit <- get("fit", e)
		hmcbootdata2 <- stanbootmean(stanfit, nTraj, 500.0, Tlim, steps)[,(Tlim+2):(T+1)]
		hmcproj2 <- apply(hmcbootdata2, 2, mean)
		hmctraj2[fctr,] <- hmcproj2

		## S-map
		smapproj <- get("predictions", e)
		smaptraj[fctr,] <- smapproj
		smaptimes[fctr] <- get("smaptime", e)[['user.self']]
		
	}
}

## trim
if2traj   <- if2traj[1:fctr,]
hmctraj   <- hmctraj[1:fctr,]
smaptraj  <- smaptraj[1:fctr,]
truetraj  <- truetraj[1:fctr,]
if2times  <- if2times[1:fctr]
hmctimes  <- hmctimes[1:fctr]
smaptimes <- smaptimes[1:fctr]

#temp
hmctraj2  <- hmctraj2[1:fctr,]

## average times
if2meantime <- mean(if2times)
hmcmeantime <- mean(hmctimes)
smapmeantime <- mean(smaptimes)

printvar( if2meantime )
printvar( hmcmeantime )
printvar( smapmeantime )

if2sse <- colMeans((abs(if2traj - truetraj)^2))
hmcinds <- complete.cases(hmctraj)
hmcsse <- colMeans((abs(hmctraj[hmcinds,] - truetraj[hmcinds,])^2))
smapsse <- colMeans((abs(smaptraj - truetraj)^2))

#temp
hmcinds2 <- complete.cases(hmctraj2)
hmcsse2 <- colMeans((abs(hmctraj2[hmcinds2,] - truetraj[hmcinds2,])^2))


if2sseraw 	<- (abs(if2traj - truetraj)^2)
if2quant  	<- t(apply(if2sseraw, 2, quantile, probs = c(0.025,0.975)))
hmcsseraw 	<- (abs(hmctraj[hmcinds,] - truetraj[hmcinds,])^2)
hmcquant  	<- t(apply(hmcsseraw, 2, quantile, probs = c(0.025,0.975)))
smapsseraw 	<- (abs(smaptraj - truetraj)^2)
smapquant 	<- t(apply(smapsseraw, 2, quantile, probs = c(0.025,0.975)))

quantdf <- data.frame(time = 1:dim(if2quant)[1], if2 = log(if2quant), hmc = log(hmcquant), smap = log(smapquant))
#quantdf[quantdf < 0] <- 0
quantplot <- ggplot(quantdf, aes(x = time)) +
				geom_ribbon(aes(ymin = quantdf[,"if2.2.5."], ymax = quantdf[,"if2.97.5."]), alpha = 0.01) +
				geom_ribbon(aes(ymin = quantdf[,"hmc.2.5."], ymax = quantdf[,"hmc.97.5."]), alpha = 0.02) +
				geom_ribbon(aes(ymin = quantdf[,"smap.2.5."], ymax = quantdf[,"smap.97.5."]), alpha = 0.03) +
				theme_bw()

## SSE plot
##########################################################################################

df <- data.frame(time = 1:(T-Tlim), if2 = log(if2sse), hmc = log(hmcsse), smap = log(smapsse))
plotdata <- melt(df, id = "time")

q <- qplot(data = plotdata, x = time, y = value, geom = "line", color = variable, xlab = "Weeks ahead", ylab = "Average error (log)") +
		scale_color_grey(labels = c("IF2","HMCMC","S-map")) +
		theme_bw() +
		theme(legend.title=element_blank())

ggsave(q, filename = "sseplot.pdf", width = 6.5, height = 4)

df2 <- data.frame(time = 1:(T-Tlim), if2 = log(if2sse), hmc = log(hmcsse2), smap = log(smapsse))
plotdata2 <- melt(df2, id = "time")

q2 <- qplot(data = plotdata2, x = time, y = value, geom = "line", color = variable, xlab = "Weeks ahead", ylab = "Average error (log)") +
		scale_color_grey(labels = c("IF2","HMCMC","S-map")) +
		theme_bw() +
		theme(legend.title=element_blank())

ggsave(q2, filename = "sseplot2.pdf", width = 6.5, height = 4)


# times plot

timedf <- data.frame(if2 = if2times, hmc = hmctimes, smap = smaptimes)
timedata <- melt(timedf)
timeplot <- ggplot(timedata, aes(factor(variable, ordered = TRUE), value)) +
				geom_boxplot() +
				scale_x_discrete(labels = c("IF2", "HMCMC", "S-map")  ) +
				labs(x = "", y = "Time (seconds)") +
				coord_flip() +
				theme_bw()

ggsave(timeplot, filename = "timeplot.pdf", width = 6.5, height = 4)



#save.image("fsim.RData")