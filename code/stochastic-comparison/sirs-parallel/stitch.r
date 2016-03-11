library(ggplot2)
library(reshape2)
library(RColorBrewer)

dir <- paste(getwd(), "varimages", sep = "/")
filelist <- list.files(dir)

nTrials <- length(filelist)
first <- TRUE
fctr <- 0

for (filenum in 1:length(filelist)) {

	filename <- filelist[filenum]
	filepath <- paste(getwd(), "varimages", filename, sep = "/")

	if (filename != "thread-0.RData") {

		fctr <- fctr + 1
		print(paste(fctr, filename))

		e <- new.env()
		load(filepath, e)

		## get sizing info from first file
		if (first) {
			stepsAhead <- get("stepsAhead", e)
			T 		   <- get("T", e)
			Tlim 	   <- get("Tlim", e)
			truetraj <- matrix(NA, nrow = nTrials, ncol = stepsAhead)
			if2traj <- matrix(NA, nrow = nTrials, ncol = stepsAhead)
			hmctraj <- matrix(NA, nrow = nTrials, ncol = stepsAhead)
			smaptraj <- matrix(NA, nrow = nTrials, ncol = stepsAhead)
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

		## HMCMC
		hmcbootdata <- get("bootstrapdata", e)[,(Tlim+2):(T+1)]
		hmcproj <- apply(hmcbootdata, 2, mean)
		hmctraj[fctr,] <- hmcproj

		## S-map
		smapproj <- get("predictions", e)
		smaptraj[fctr,] <- smapproj
		
	}
}

## trim
if2traj <- if2traj[1:fctr,]
hmctraj <- hmctraj[1:fctr,]
smaptraj <- smaptraj[1:fctr,]
truetraj <- truetraj[1:fctr,]

if2sse <- colMeans((abs(if2traj - truetraj)^2))
hmcinds <- complete.cases(hmctraj)
hmcsse <- colMeans((abs(hmctraj[hmcinds,] - truetraj[hmcinds,])^2))
smapsse <- colMeans((abs(smaptraj - truetraj)^2))

q <- qplot(1:(T-Tlim), log(if2sse), geom = "line") +
		geom_line(aes(y = log(hmcsse)), linetype = "dashed") +
		geom_line(aes(y = log(smapsse)), linetype = "dotted") +
		theme_bw()

ggsave(q, filename = "sseplot.pdf", width = 6.5, height = 3)

save.image("fsim.RData")