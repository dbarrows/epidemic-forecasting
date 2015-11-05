library(deSolve)
library(rstan)
library(shinystan)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

# NOTE: to save plots, uncomment the ggsave lines

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
times <- 0:100

odeout <- ode(y_ini, times, SIR, pars)

set.seed(1001)
sigma <- 10
infec_counts_raw <- odeout[,3] + rnorm(101, 0, sigma)
infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

plotdata <- data.frame(times=1:(T+1),true=odeout[,3],data=infec_counts)

g <- ggplot(plotdata, aes(times)) +
        geom_line(aes(y = true, colour = "True")) + 
        geom_point(aes(y = data, colour = "Data")) +
        labs(x = "Time", y = "Infection count", color = "") +
        scale_color_brewer(palette="Paired") +
        theme(panel.background = element_rect(fill = "#F0F0F0"))

print(g)
#ggsave(g, filename="dataplot.pdf", height=4, width=6.5)

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
stan_options <- list(   chains = 10,    # number of chains
                        iter   = 5000, # iterations per chain
                        warmup = 1000, # warmup interations
                        thin   = 10)   # thinning number
fit <- stan(file    = "sirode_euler.stan",
            data    = sir_data,
            chains  = stan_options$chains,
            iter    = stan_options$iter,
            warmup  = stan_options$warmup,
            thin    = stan_options$thin )

exfit <- extract(fit, permuted = TRUE, inc_warmup = FALSE)

R0points <- exfit$R0
kerdataR0 <- data.frame(R0points)
R0kernel <- ggplot(kerdataR0, aes(x = R0points, y = ..scaled..)) +
                geom_density(colour="#8dd3c7", fill="#8dd3c7") +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = expression(R[0]), y = "Density", color = "") +
                geom_vline(aes(xintercept=R0), linetype="dashed", size=1, color="grey50")

print(R0kernel)
#ggsave(R0kernel, filename="kernelR0.pdf", height=3, width=3.25)

rpoints <- exfit$r
kerdatar <- data.frame(rpoints)
rkernel <- ggplot(kerdatar, aes(x = rpoints, y = ..scaled..)) +
                geom_density(colour="#ffffb3", fill="#ffffb3") +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = "r", y = "Density", color = "") +
                geom_vline(aes(xintercept=r), linetype="dashed", size=1, color="grey50")

print(rkernel)
#ggsave(rkernel, filename="kernelr.pdf", height=3, width=3.25)

sigmapoints <- exfit$sigma
kerdatasigma <- data.frame(sigmapoints)
sigmakernel <- ggplot(kerdatasigma, aes(x = sigmapoints, y = ..scaled..)) +
                geom_density(colour="#bebada", fill="#bebada") +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = expression(sigma), y = "Density", color = "") +
                geom_vline(aes(xintercept=sigma), linetype="dashed", size=1, color="grey50")

print(sigmakernel)
#ggsave(sigmakernel, filename="kernelsigma.pdf", height=3, width=3.25)

infecpoints <- exfit$y0[,2]
kerdatainfec <- data.frame(infecpoints)
infeckernel <- ggplot(kerdatainfec, aes(x = infecpoints, y = ..scaled..)) +
                geom_density(colour="#fb8072", fill="#fb8072") +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = "Initial infected", y = "Density", color = "") +
                geom_vline(aes(xintercept=5), linetype="dashed", size=1, color="grey50")

print(infeckernel)
#ggsave(infeckernel, filename="kernelinfec.pdf", height=3, width=3.25)

exfit <- extract(fit, permuted = FALSE, inc_warmup = FALSE)
plotdata <- melt(exfit[,,2])
tracefitR0 <- ggplot() +
              geom_line(data = plotdata,
                        aes(x = iterations,
                        y = value,
                        color = factor(chains, labels = 1:stan_options$chains))) +
              labs(x = "Sample", y = expression(R[0]), color = "Chain") +
              scale_color_brewer(palette="Paired") +
              theme(panel.background = element_rect(fill = "#F0F0F0"))

print(tracefitR0)
#ggsave(tracefitR0, filename="traceplotR0.pdf", height=4, width=6.5)

exfit <- extract(fit, permuted = FALSE, inc_warmup = TRUE)
plotdata <- melt(exfit[,,2])
tracefitR0 <- ggplot() +
              geom_line(data = plotdata,
                        aes(x = iterations,
                        y = value,
                        color = factor(chains, labels = 1:stan_options$chains))) +
              labs(x = "Sample", y = expression(R[0]), color = "Chain") +
              scale_color_brewer(palette="Paired") +
              theme(panel.background = element_rect(fill = "#F0F0F0"))

print(tracefitR0)
#ggsave(tracefitR0, filename="traceplotR0_inc.pdf", height=4, width=6.5)

exfit <- extract(fit, permuted = FALSE, inc_warmup = FALSE)
paramdata <- data.frame(R0 = melt(exfit[,,1])$value,
               			r = melt(exfit[,,2])$value,
               			sigma = melt(exfit[,,3])$value,
               			Sinit = melt(exfit[,,4])$value,
               			Iinit = melt(exfit[,,5])$value,
               			Rinit = melt(exfit[,,6])$value )

# sample from parameter distributions
nTraj <- 1000
datlen <- dim(paramdata)[1]
inds <- sample.int(datlen,nTraj,replace = TRUE)
params <- paramdata[inds,]

bootstrapdata <- matrix(NA, nrow = nTraj, ncol = (T + 1))

for (i in 1:nTraj) {
	init_cond <- c(S = params$Sinit[i],
	               I = params$Iinit[i],
	               R = params$Rinit[i])
	pars <- c(R0 = params$R0[i],
	          r = params$r[i],
	          N = 500.0)
	odeout <- ode(init_cond, 0:T, SIR, pars)
	bootstrapdata[i,] <- odeout[,3]
}

meanTraj <- colMeans(bootstrapdata)
sdTraj <- apply(bootstrapdata, 2, sd)

true_init_cond <- c(S = 495,
                    I = 5,
                    R = 0)
true_pars <- c(R0 = 3.0,
               r = 0.1,
               N = 500.0)
odeout <- ode(true_init_cond, 0:T, SIR, true_pars)
trueTraj <- odeout[,3]

plotdata <- data.frame(times=1:(T+1),true=trueTraj,data=meanTraj,sds=sdTraj)

g <- ggplot(plotdata, aes(times)) +
		geom_ribbon(aes(ymin = data-sds, ymax=data+sds), alpha=0.1) +
        geom_line(aes(y = true, colour = "True")) + 
        geom_line(aes(y = data, colour = "HMCMC")) +
        labs(x = "Time", y = "Infection count", color = "") +
        scale_color_brewer(palette="Paired") +
        theme(panel.background = element_rect(fill = "#F0F0F0"))

print(g)

sso <- as.shinystan(fit)
sso <- launch_shinystan(sso)