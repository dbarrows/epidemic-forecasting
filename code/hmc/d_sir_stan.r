## Dexter Barrows
## McMaster University
## 2016

library(deSolve)
library(rstan)
library(shinystan)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

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

pars  <- c(R0  <- 3.0,    # average number of new infected individuals per infectious person
           r   <- 0.1,    # recovery rate
           N   <- 500)    # population size

T <- 100
y_ini <- c(S = 495, I = 5, R = 0)
times <- seq(0, T, by = 1)

odeout <- ode(y_ini, times, SIR, pars)

set.seed(1001)
sigma <- 10
infec_counts_raw <- odeout[,3] + rnorm(101, 0, sigma)
infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

g <- qplot(0:T, odeout[,3], geom = "line", xlab = "Time (weeks)", ylab = "Infection Count") +
		geom_point(aes(y = infec_counts)) +
		theme_bw()

print(g)
ggsave(g, filename="dataplot.pdf", height=4, width=6.5)

sPw <- 7
datlen <- (T-1)*7 + 1

data <- matrix(data = -1, nrow = T+1, ncol = sPw)
data[,1] <- infec_counts
standata <- as.vector(t(data))[1:datlen]

sir_data <- list( T = datlen,   # simulation time
                  y = standata, # infection count data
                  N = 500,      # population size
                  h = 1/sPw )   # step size per day 
                    
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stan_options <- list(   chains = 4,    # number of chains
                        iter   = 2000, # iterations per chain
                        warmup = 1000, # warmup interations
                        thin   = 2)   # thinning number
fit <- stan(file    = "d_sirode_euler.stan",
            data    = sir_data,
            chains  = stan_options$chains,
            iter    = stan_options$iter,
            warmup  = stan_options$warmup,
            thin    = stan_options$thin )

exfit <- extract(fit, permuted = TRUE, inc_warmup = FALSE)

R0points <- exfit$R0
R0kernel <- qplot(R0points, geom = "density", xlab = expression(R[0]), ylab = "frequency") +
				geom_vline(aes(xintercept=R0), linetype="dashed", size=1, color="grey50") +
				theme_bw()

print(R0kernel)
ggsave(R0kernel, filename="kernelR0.pdf", height=3, width=3.25)

rpoints <- exfit$r
rkernel <- qplot(rpoints, geom = "density", xlab = "r", ylab = "frequency") +
				geom_vline(aes(xintercept=r), linetype="dashed", size=1, color="grey50") +
				theme_bw()

print(rkernel)
ggsave(rkernel, filename="kernelr.pdf", height=3, width=3.25)

sigmapoints <- exfit$sigma
sigmakernel <- qplot(sigmapoints, geom = "density", xlab = expression(sigma), ylab = "frequency") +
				geom_vline(aes(xintercept=sigma), linetype="dashed", size=1, color="grey50") +
				theme_bw()

print(sigmakernel)
ggsave(sigmakernel, filename="kernelsigma.pdf", height=3, width=3.25)

infecpoints <- exfit$y0[,2]
infeckernel <- qplot(infecpoints, geom = "density", xlab = "Initial Infected", ylab = "frequency") +
				geom_vline(aes(xintercept=y_ini[['I']]), linetype="dashed", size=1, color="grey50") +
				theme_bw()

print(infeckernel)
ggsave(infeckernel, filename="kernelinfec.pdf", height=3, width=3.25)

exfit <- extract(fit, permuted = FALSE, inc_warmup = FALSE)
plotdata <- melt(exfit[,,"R0"])
tracefitR0 <- ggplot() +
              geom_line(data = plotdata,
                        aes(x = iterations,
                        y = value,
                        color = factor(chains, labels = 1:stan_options$chains))) +
              labs(x = "Sample", y = expression(R[0]), color = "Chain") +
              scale_color_brewer(palette="Greys") +
              theme_bw()

print(tracefitR0)
ggsave(tracefitR0, filename="traceplotR0.pdf", height=4, width=6.5)

exfit <- extract(fit, permuted = FALSE, inc_warmup = TRUE)
plotdata <- melt(exfit[,,"R0"])
tracefitR0 <- ggplot() +
              geom_line(data = plotdata,
                        aes(x = iterations,
                        y = value,
                        color = factor(chains, labels = 1:stan_options$chains))) +
              labs(x = "Sample", y = expression(R[0]), color = "Chain") +
              scale_color_brewer(palette="Greys") +
              theme_bw()

print(tracefitR0)
ggsave(tracefitR0, filename="traceplotR0_inc.pdf", height=4, width=6.5)

sso <- as.shinystan(fit)
sso <- launch_shinystan(sso)