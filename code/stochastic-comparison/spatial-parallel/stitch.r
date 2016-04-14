library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)

printvar <- function(v) {
  	name <- deparse(substitute(v))
  	print(paste(name, ":", v))
}

mine <- function(L){
    if (length(L) == 2) {
        ctr <<- ctr + 1
        masterlist[[ctr]] <<- L[[2]]
        mine(L[[1]])
    } else {
        ctr <<- ctr + 1
        masterlist[[ctr]] <<- L
    }
}

## to write out
writefile <- source(paste(getwd(), "../../cuIF2", "writedata.r", sep = "/"))
datdir <- paste("../../cuIF2/data")

dir <- paste(getwd(), "../../../../Thesis_data/spatial-varimages", sep = "/")
filelist <- list.files(dir)

nTrials <- length(filelist)
first <- TRUE
fctr <- 0

for (filenum in 1:length(filelist)) {

    filename <- filelist[filenum]
    filepath <- paste(dir, filename, sep = "/")

    if (filename != "thread-0.RData") {

        fctr <- fctr + 1
        print(paste(fctr, filename))

        e <- new.env()
        load(filepath, e)

        ## get sizing info from first file
        if (first) {

            stepsAhead  <- get("stepsAhead", e)
            T           <- get("T", e)
            Tlim        <- get("Tlim", e)
            nloc        <- get("nloc", e)
            nTraj       <- get("nTraj", e)

            truetraj    <- array(NA, c(length(filelist), nloc, stepsAhead) )
            if2traj     <- array(NA, c(length(filelist), nloc, stepsAhead) )
            hmctraj     <- array(NA, c(length(filelist), nloc, stepsAhead) )
            smaptraj    <- array(NA, c(length(filelist), nloc, stepsAhead) )

            countmeans   <- matrix(NA, nloc, T-Tlim+1)
            parabootdata <- array(NA, c(nTraj, nloc, T-Tlim+1))

            if2times    <- numeric(nTrials)
            hmctimes    <- numeric(nTrials)
            smaptimes   <- numeric(nTrials)

            if2fittimes <- numeric(nTrials)
            hmcfittimes <- numeric(nTrials)
            cuIF2times  <- numeric(nTrials)

            if2covdata <- data.frame(matrix(NA, nTrials, 6))
			names(if2covdata) <- c("mean2", "lower2", "upper2", "mean10", "lower10", "upper10")
			hmccovdata <- data.frame(matrix(NA, nTrials, 6))
			names(hmccovdata) <- c("mean2", "lower2", "upper2", "mean10", "lower10", "upper10")

            first <- FALSE

        }

        ## True trajectory
        trueproj <- get("ssdeout_true", e)[,'I',(Tlim+2):(T+1)]
        truetraj[fctr,,] <- trueproj

        # write data out for cuIF2
        data <- get("ssdeout_true", e)[,'I',1:(Tlim+1)]
        neinum <- get("neinum", e)
        neibmat <- get("neibmat", e)
        trial <- get("trial", e)
        writedata(data, neinum, neibmat, datdir, trial)

        ## IF2 extract data
        masterlist <- list()
        ctr <- 0
        if2_paraboot_data <- get("if2_paraboot_data", e)
        mine(if2_paraboot_data)
        # get mean of trajectories
        for (i in 1:nTraj) {
            datalist <- masterlist[[i]]
            trajectories <- datalist$counts
            parabootdata[i,,] <- trajectories
        }
        for (i in 1:nloc) {
            locslice <- parabootdata[,i,]
            locslice <- locslice[complete.cases(locslice),]
            locslice[locslice < 0] <- 0
            countmeans[i,] <- colMeans(locslice)
        }
        if2traj[fctr,,] <- countmeans[,-1]
        if2times[fctr] <- get("if2time", e)[['user.self']]
        if2fittimes[fctr] <- get("if2time1", e)[['user.self']]
        # coverage
        # location 8, times 2 and 10
        real2 	 <- trueproj[8,2]
		real10 	 <- trueproj[8,10]
		quants2 <- quantile(parabootdata[,8,2], probs = c(0.05, 0.95), na.rm = TRUE)
		quants10 <- quantile(parabootdata[,8,10], probs = c(0.05, 0.95), na.rm = TRUE)
		mean2   <- countmeans[8,3] - real2
		diff2 	 <- quants2 - real2
		mean10 	 <- countmeans[8,11] - real10
		diff10 	 <- quants10 - real10
		if2covdata[fctr,] <- c(mean2, diff2, mean10, diff10)


        ## HMCMC
        hmcbootdata <- get("bootstrapdata", e)[,,(Tlim+2):(T+1)]
        meanTraj <- matrix(NA, nloc, T-Tlim)
        # in case of explosion
        for (loc in 1:nloc) {
            locdata <- hmcbootdata[,loc,]
            locdata <- locdata[complete.cases(locdata),]
            locdata[locdata < 0] <- 0
            meanTraj[loc,] <- colMeans(locdata)
        }
        hmctraj[fctr,,] <- meanTraj
        hmctimes[fctr] <- get("hmctime", e)[['user.self']]
        hmcfittimes[fctr] <- get("hmctime1", e)[['user.self']]
        # coverage
        # location 8, times 2 and 10
        real2 	 <- trueproj[8,2]
		real10 	 <- trueproj[8,10]
		quants2 <- quantile(hmcbootdata[,8,2], probs = c(0.05, 0.95), na.rm = TRUE)
		quants10 <- quantile(hmcbootdata[,8,10], probs = c(0.05, 0.95), na.rm = TRUE)
		mean2   <- meanTraj[8,2] - real2
		diff2 	 <- quants2 - real2
		mean10 	 <- meanTraj[8,10] - real10
		diff10 	 <- quants10 - real10
		hmccovdata[fctr,] <- c(mean2, diff2, mean10, diff10)


        ## S-map
        smapproj <- get("predictions", e)
        smapproj[smapproj < 0] <- 0
        smaptraj[fctr,,] <- smapproj
        smaptimes[fctr] <- get("smaptime", e)[['user.self']]
        
    }
}

## cuIF2 get data
cudir 	<- paste(getwd(), "../../cuIF2/log", sep = "/")
cufiles <- list.files(cudir)

for (filenum in 1:length(cufiles)) {
	filepath <- paste(cudir, cufiles[filenum], sep = "/")
	sysout <- system( paste("cat ", filepath, " | grep Rawtime  | sed -e s/[^0-9.]//g" , sep = ""), intern = TRUE)
	cuIF2times[filenum] <- as.numeric(sysout)
}


## trim
if2traj   <- if2traj[1:fctr,,]
hmctraj   <- hmctraj[1:fctr,,]
smaptraj  <- smaptraj[1:fctr,,]
truetraj  <- truetraj[1:fctr,,]
if2times  <- if2times[1:fctr]
hmctimes  <- hmctimes[1:fctr]
smaptimes <- smaptimes[1:fctr]

## average times
if2meantime <- mean(if2times)
hmcmeantime <- mean(hmctimes)
smapmeantime <- mean(smaptimes)

printvar( if2meantime )
printvar( hmcmeantime )
printvar( smapmeantime )

## average fitting times
if2fitmeantime <- mean(if2fittimes)
hmcfitmeantime <- mean(hmcfittimes)
cuIF2meantime <- mean(cuIF2times)

printvar( if2fitmeantime )
printvar( hmcfitmeantime )
printvar( cuIF2meantime )

if2sses <- matrix(NA, nloc, stepsAhead)
hmcsses <- matrix(NA, nloc, stepsAhead)
smapsses <- matrix(NA, nloc, stepsAhead)

for (loc in 1:nloc) {
    if2sses[loc,] <- colMeans((abs(if2traj[,loc,] - truetraj[,loc,])^2))
    hmcinds <- complete.cases(hmctraj[,loc,])
    hmcsses[loc,] <- colMeans((abs(hmctraj[hmcinds,loc,] - truetraj[hmcinds,loc,])^2))
    smapsses[loc,] <- colMeans((abs(smaptraj[,loc,] - truetraj[,loc,])^2))
}

if2ssecollapse <- colMeans(if2sses)
hmcssecollapse <- colMeans(hmcsses)
smapssecollapse <- colMeans(smapsses)

## SSE plot
##########################################################################################

df <- data.frame(time = 1:(T-Tlim), if2 = log10(if2ssecollapse), hmc = log10(hmcssecollapse), smap = log10(smapssecollapse))
plotdata <- melt(df, id = "time")

q <- qplot(data = plotdata, x = time, y = value, geom = "line", color = variable, xlab = "Weeks ahead", ylab = expression(SSE~(log[10]))) +
		scale_color_grey(labels = c("IF2","HMCMC","S-map")) +
		theme_bw() +
		theme(legend.title=element_blank())

ggsave(q, filename = "sseplot.pdf", width = 6.5, height = 4)

## times plot
##########################################################################################

timedf <- data.frame(if2 = if2times, hmc = hmctimes, smap = smaptimes)
timedata <- melt(timedf)
timeplot <- ggplot(timedata, aes(factor(variable, ordered = TRUE), value)) +
				geom_boxplot() +
				scale_x_discrete(labels = c("IF2", "HMCMC", "S-map")  ) +
				labs(x = "", y = "Time (seconds)") +
				coord_flip() +
				theme_bw()

ggsave(timeplot, filename = "timeplot.pdf", width = 6.5, height = 4)

## second times plot with cuIF2

timedf2 <- data.frame(if2 = log10(if2fittimes), hmc = log10(hmcfittimes), cuif2 = log10(cuIF2times))
timedata2 <- melt(timedf2)
timeplot2 <- ggplot(timedata2, aes(factor(variable, ordered = TRUE), value)) +
				geom_boxplot() +
				scale_x_discrete(labels = c("IF2", "HMCMC", "cuIF2")  ) +
				labs(x = "", y = expression(Time~(log[10]~seconds)) )  +
				coord_flip() +
				theme_bw()

ggsave(timeplot2, filename = "timeplot2.pdf", width = 6.5, height = 4)

#save.image("fsim.RData")


## Let's try this coverage plot thing
##########################################################################################


inds <- complete.cases(hmccovdata)
nCovSets <- sum(inds)

if2pdraw <- if2covdata[inds,][1:nCovSets,]
hmcpdraw <- hmccovdata[inds,][1:nCovSets,]

## 10 weeks ahead
##

pd2if2 <- data.frame(setnum = 1:nCovSets, method = rep("if2",nCovSets), if2pdraw[order(if2pdraw[,"mean2"]),])
pd2hmc <- data.frame(setnum = 1:nCovSets, method = rep("hmc",nCovSets), hmcpdraw[order(if2pdraw[,"mean2"]),])
cnames <- names(pd2if2)
pd2if2melted <- melt(pd2if2, id = cnames)
pd2hmcmelted <- melt(pd2hmc, id = cnames)
pd2 <- rbind(pd2if2melted, pd2hmcmelted)

pd <- position_dodge(0.5)
cov2plot <- ggplot(pd2, aes(x = setnum, color = method)) +
					geom_hline(aes(yintercept = 0)) +
					geom_errorbar(aes(ymin = lower2, ymax = upper2), width = 0, position = pd) +
					geom_point(aes(y = mean2), position = pd) +
					labs(x = "", y = "Estimate - True") +
					theme_bw() +
					scale_color_grey() +
					theme(axis.ticks.x=element_blank(),
					      axis.text.x=element_blank(),
					      legend.position='none')

## 45 weeks ahead
##

pd10if2 <- data.frame(setnum = 1:nCovSets, method = rep("if2",nCovSets), if2pdraw[order(if2pdraw[,"mean10"]),])
pd10hmc <- data.frame(setnum = 1:nCovSets, method = rep("hmc",nCovSets), hmcpdraw[order(if2pdraw[,"mean10"]),])
cnames <- names(pd10if2)
pd10if2melted <- melt(pd10if2, id = cnames)
pd10hmcmelted <- melt(pd10hmc, id = cnames)
pd10 <- rbind(pd10if2melted, pd10hmcmelted)

pd <- position_dodge(0.5)
cov10plot <- ggplot(pd10, aes(x = setnum, color = method)) +
					geom_hline(aes(yintercept = 0)) +
					geom_errorbar(aes(ymin = lower10, ymax = upper10), width = 0, position = pd) +
					geom_point(aes(y = mean10), position = pd) +
					labs(x = "Data set", y = "Estimate - True") +
					theme_bw() +
					scale_color_grey() +
					theme(axis.ticks.x=element_blank(),
					      axis.text.x=element_blank(),
					      legend.position='none')


## save combined plot
##
pdf("coverage.pdf", width = 6.5, height = 4)
grid.arrange(cov2plot, cov10plot, ncol = 1, nrow = 2)
dev.off()