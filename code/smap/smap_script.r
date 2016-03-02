library(deSolve)
library(ggplot2)
library(RColorBrewer)
library(pracma)

set.seed(106745)

## external files
##
stoc_sirs_file 	<- paste(getwd(), "../sir-functions", "StocSIRS.r", sep = "/")
smap_file 		<- paste(getwd(), "smap.r", sep = "/")
source(stoc_sirs_file)
source(smap_file)

## parameters
##
T 		<- 5*52
i_infec <- 10
steps 	<- 7
N 		<- 500
sigma 	<- 5

true_pars <- c(	R0 = 3.0, 	# new infected people per infected person
          		r = 0.1, 	# recovery rate
		    	N = 500,    # population size
		    	eta = 0.5, 	# geometric random walk
		    	berr = 0.5, # Beta geometric walk noise
          		re = 1)  	# resuceptibility rate

true_init_cond <- c(S = N - i_infec,
                    I = i_infec,
                    R = 0)

## get true trajectory
##
sdeout <- StocSIRS(true_init_cond, true_pars, T, steps)

## perturb to get data
##
infec_counts_raw <- sdeout[,'I'] + rnorm(T+1,0,sigma)
infec_counts 	 <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)


# plot data and true states
g <- qplot(0:T, sdeout[,'I'], geom = "line", xlab = "Time", ylab = "Infection count") +
		geom_point(aes(y = infec_counts)) +
		theme_bw()

print(g)
ggsave(g, filename="dataplot.pdf", height=4, width=6.5)

## S-map parameters
##
E <- 14
theta <- 3
stepsAhead <- 52

## run S-map
##
predictions <- smap(infec_counts, E, theta, stepsAhead)

p <- qplot(0:(T+stepsAhead), c(infec_counts, rep(NA, stepsAhead)), geom = "line", xlab = "Time", ylab = "Infection count") +
		geom_line(aes(y = c(rep(NA, T+1), predictions)), linetype = "dotted") +
		theme_bw()

print(p)
ggsave(p, filename="smap-project.pdf", height=4, width=6.5)