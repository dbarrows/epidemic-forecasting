library(deSolve)
library(ggplot2)
library(RColorBrewer)
library(pracma)
library(reshape2)

set.seed(106745)

## external files
#################################################################################

stoc_spa_sir_file 	<- paste(getwd(), "../../sir-functions", "StocSSIR.r", sep = "/")
dewdrop_file 		<- paste(getwd(), "dewdrop.r", sep = "/")
source(stoc_spa_sir_file)
source(dewdrop_file)

## parameters
#################################################################################

T 		<- 100
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

## Simulate data
#################################################################################

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

infec_counts <- ssdeout[,'I',] + matrix(rnorm(nloc*(T+1), 0, sigma), nloc, T+1)
infec_counts <- ifelse(infec_counts < 0, 0, infec_counts)

## Plot and save data
#################################################################################

df <- data.frame(time = 0:T, loc = t(ssdeout[,'I',]))
plotdata <- melt(df, id = "time")

## full
q <- qplot(data = plotdata, x = time, y = value, geom = "line", color = variable, xlab = "Time", ylab = "Infection Count") +
		scale_color_grey() +
		theme_bw()

ggsave(q, filename="dataplot.pdf", height=4, width=6.5)

## S-map
#################################################################################

## parameters
##

E <- 10
theta <- 3
stepsAhead <- 10

## comform data

Tlim <- 50
smapdata <- ssdeout[,'I',1:(Tlim + 1)]

## run S-map
##

predictions <- dewdrop(smapdata, E, theta, stepsAhead)

qp <- qplot(0:(Tlim+stepsAhead), c(smapdata[2,], rep(NA, stepsAhead)), geom = "line") +
		geom_line(aes(y = c(rep(NA, Tlim+1),predictions[2,])), linetype = "dashed") +
		theme_bw()

#p <- qplot(0:(T+stepsAhead), c(infec_counts, rep(NA, stepsAhead)), geom = "line", xlab = "Time", ylab = "Infection count") +
#		geom_line(aes(y = c(rep(NA, T+1), predictions)), linetype = "dotted") +
#		theme_bw()

#print(p)
#ggsave(p, filename="smap-project.pdf", height=4, width=6.5)