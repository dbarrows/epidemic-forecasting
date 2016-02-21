library(ggplot2)
library(gridExtra)
library(reshape2)

source(paste(getwd(), "../../", "sir-functions/StocSSIR.r", sep="/"))
source("multiplot.r")

T 		<- 60
N 		<- 500
i_infec <- 5
nloc  	<- 10
steps   <- 7
sigma   <- 10

pars <- c(	R0 	 = 3.0, # new infected people per infected person
          	r 	 = 0.1, # recovery rate
          	N 	 = 500, # population size
          	eta  = 0.5, # geometric random walk
          	berr = 0.5, # Beta geometric walk noise
          	phi  = 0.1) # degree of interconnectivity

initcond <- matrix(NA, nloc, 3)
colnames(initcond) <- c("S","I","R")

## default (everyone susceptible)
initcond[,"S"] <- N
initcond[,"I"] <- 0
initcond[,"R"] <- 0

## one region infected
initcond[2,] <- c(N - i_infec, i_infec, 0);

## set up number of neighbors vector and neighbor index matrix

## ring topology:
neinum <- rep(2,nloc)
neibmat <- matrix(NA, nloc, 3)
neibmat[,1] <- c( ((0:(nloc-1)-1)%%nloc)+1 )
neibmat[,2] <- c( ((0:(nloc-1)+1)%%nloc)+1 )

ssdeout <- StocSSIR(initcond, pars, T, steps, neinum, neibmat) 

df <- data.frame(time = 0:T, loc = t(ssdeout[,'I',]))
plotdata <- melt(df, id = "time")

## full
qplot(data = plotdata, x = time, y = value, geom = "line", color = variable) + theme_bw()

## 1 v. 6, no bridge
quartz()
qplot(0:T, ssdeout[1,'I',], geom = "line") +
	geom_line(aes(y = ssdeout[6,'I',])) +
	theme_bw()


## 1 v. 6 bridge

## add some extra connections
neinum[1] <- neinum[1] + 1
neibmat[1,3] <- 6

ssdeout <- StocSSIR(initcond, pars, T, steps, neinum, neibmat) 

quartz()
qplot(0:T, ssdeout[1,'I',], geom = "line") +
	geom_line(aes(y = ssdeout[6,'I',])) +
	theme_bw()