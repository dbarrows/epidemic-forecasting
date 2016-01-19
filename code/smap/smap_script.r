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

#print(g)


## SMAP

# trim off a bit of the time series - TEMPORARY
lim <- length(as.vector(infec_counts)) - 24
infec_counts <- as.vector(infec_counts)[1:lim]

## smap function form external file
source("smap.r")

E <- 10
theta <- 10
stepsAhead <- 200

predictions <- smap(infec_counts, E, theta, stepsAhead)

dlen <- length(infec_counts)
dlen
predictee <- infec_counts[(dlen-E+1):dlen]
predictee

plotdata <- data.frame(times = 1:(length(infec_counts)+stepsAhead),                                         # master index of times
                       tseries = c(infec_counts, rep(NA,stepsAhead) ),                                      # time series data, padded
                       predictions = c( rep(NA,length(infec_counts)), predictions),                         # predictions
                       predictee = c(rep(NA,length(infec_counts)-E),predictee,rep(NA,stepsAhead)) )    # last predictee vector


p <- ggplot(plotdata, aes(times)) +
        geom_line(aes(y = tseries), linetype = "solid") +
        geom_line(aes(y = predictions), linetype = "twodash") +
        geom_point(aes(y = predictee)) +
        labs(x = "Time", y = "Infection count", color = "") +
        theme_bw()

print(p)