library(ggplot2)
library(reshape2)
library(RColorBrewer)

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

        ## S-map
        smapproj <- get("predictions", e)
        smapproj[smapproj < 0] <- 0
        smaptraj[fctr,,] <- smapproj
        smaptimes[fctr] <- get("smaptime", e)[['user.self']]
        
    }
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

df <- data.frame(time = 1:(T-Tlim), if2 = log(if2ssecollapse), hmc = log(hmcssecollapse), smap = log(smapssecollapse))
plotdata <- melt(df, id = "time")

q <- qplot(data = plotdata, x = time, y = value, geom = "line", color = variable, xlab = "Weeks ahead", ylab = "Average error (log)") +
		scale_color_grey(labels = c("IF2","HMCMC","S-map")) +
		theme_bw() +
		theme(legend.title=element_blank())

ggsave(q, filename = "sseplot.pdf", width = 6.5, height = 4)

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