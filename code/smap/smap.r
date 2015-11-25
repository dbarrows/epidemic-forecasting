library(deSolve)
library(ggplot2)
library(RColorBrewer)
library(pracma)

# base simulation ODE RHS function

SIR <- function(Time, State, Pars) {
    
    with(as.list(c(State, Pars)), {
        
        B   <- R0*r/N
        BSI <- B*S*I
        rI  <- r*I
        
        dS = -BSI
        dI = BSI - rI
        dR = rI
        
        return(list(c(dS, dI, dR)))
        
    })
    
}

# parameters

T 		<- 100
N 		<- 500
sigma 	<- 10
i_infec <- 1
replicates <- 5

# get true underlying states

true_init_cond <- c(S = N - i_infec,
                    I = i_infec,
                    R = 0)

true_pars <- c(R0 = 3.0,
               r = 0.1,
               N = 500.0)

odeout <- ode(true_init_cond, 0:(T-1), SIR, true_pars)
trueTraj <- odeout[,3]

# replicate to mimic cycles

multiTrueTraj <- repmat(trueTraj, 1, replicates)

# perturb with normally distributed noise

set.seed(1000)

infec_counts_raw <- multiTrueTraj + rnorm(replicates*T, 0, sigma)
infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

# plot data and true states

plotdata <- data.frame(times = 1:(replicates*T),
                       true = t(multiTrueTraj),
                       data = t(infec_counts) )

g <- ggplot(plotdata, aes(times)) +
        geom_line(aes(y = true, colour = "True")) +
        geom_point(aes(y = data, color = "Data")) +
        labs(x = "Time", y = "Infection count", color = "") +
        scale_color_brewer(palette="Paired") +
        theme(panel.background = element_rect(fill = "#F0F0F0"))

print(g)


## SMAP

# trim off a bit of the time series - TEMPORARY
lim <- length(as.vector(infec_counts)) - 24
infec_counts <- as.vector(infec_counts)[1:lim]

E <- 10
theta <- 10
stepsAhead <- 200

# construct library
tseries <- as.vector(infec_counts)
liblen 	<- length(tseries) - E + 1 - stepsAhead
lib 	<- matrix(NA, liblen, E)

for (i in 1:E) {
	lib[,i] <- tseries[(E-i+1):(liblen+E-i)]
}

# predict from the last index
tslen <- length(tseries)
predictee <- rev(t(as.matrix(tseries[(tslen-E+1):tslen])))
predictions <- numeric(stepsAhead)

#allPredictees <- matrix(NA, stepsAhead, E)

# for each prediction index (number of steps ahead)
for(i in 1:stepsAhead) {

	# set up weight calculation
	predmat <- repmat(predictee, liblen, 1)
	distances <- sqrt( rowSums( abs(lib - predmat)^2 ) )
	meanDist <- mean(distances)

	# calculate weights
	weights <- exp( - (theta * distances) / meanDist )

	# construct A, B

	preds <- tseries[(E+i):(liblen+E+i-1)]

	A <- cbind( rep(1.0, liblen), lib ) * repmat(as.matrix(weights), 1, E+1)
	B <- as.matrix(preds * weights)

	# solve system for C

	Asvd <- svd(A)
	C <- Asvd$v %*% diag(1/Asvd$d) %*% t(Asvd$u) %*% B

	# get prediction

	predsum <- sum(C * c(1,predictee))

	# save

	predictions[i] <- predsum

	# next predictee

	#predictee <- c( predsum, predictee[-E] )
	#allPredictees[i,] <- predictee

}

plotdata <- data.frame(times = 1:(length(tseries)+stepsAhead),
                       tseries = c(tseries, rep(NA,stepsAhead) ),
                       predictions = c( rep(NA,length(tseries)), predictions),
                       predictee = c(rep(NA,length(tseries)-E),rev(predictee),rep(NA,stepsAhead)) )


p <- ggplot(plotdata, aes(times)) +
        geom_line(aes(y = tseries, colour = "Data")) +
        geom_line(aes(y = predictions, color = "Predictions")) +
        geom_point(aes(y = predictee, color = "Predictee")) +
        labs(x = "Time", y = "Infection count", color = "") +
        scale_color_brewer(palette="Paired") +
        theme(panel.background = element_rect(fill = "#F0F0F0"))

print(p)