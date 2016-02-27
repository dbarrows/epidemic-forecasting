##	Author: Dexter Barrows
##	Github: dbarrows.github.io

library(deSolve)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(Rcpp)

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

T 		<- 100
N 		<- 500
sigma 	<- 10
i_infec <- 5

## Generate true trajecory and synthetic data
##

true_init_cond <- c(S = N - i_infec,
                    I = i_infec,
                    R = 0)

true_pars <- c(R0 = 3.0,
               r = 0.1,
               N = 500.0)

odeout <- ode(true_init_cond, 0:T, SIR, true_pars)
trueTraj <- odeout[,3]

set.seed(1001)

infec_counts_raw <- odeout[,3] + rnorm(T+1, 0, sigma)
infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

g <- qplot(0:T, odeout[,3], geom = "line", xlab = "Time (weeks)", ylab = "Infection Count") +
		geom_point(aes(y = infec_counts)) +
		theme_bw()

print(g)
ggsave(g, filename="dataplot.pdf", height=4, width=6.5)

## Rcpp stuff
##

sourceCpp(paste(getwd(),"d_if2.cpp",sep="/"))

paramdata <- data.frame(if2(infec_counts, T+1, N))
colnames(paramdata) <- c("R0", "r", "sigma", "Sinit", "Iinit", "Rinit")

## Parameter density kernels
##

R0points <- paramdata$R0
R0kernel <- qplot(R0points, geom = "density", xlab = expression(R[0]), ylab = "frequency") +
				geom_vline(aes(xintercept=true_pars[["R0"]]), linetype="dashed", size=1, color="grey50") +
				theme_bw()

print(R0kernel)
ggsave(R0kernel, filename="kernelR0.pdf", height=3, width=3.25)

rpoints <- paramdata$r
rkernel <- qplot(rpoints, geom = "density", xlab = "r", ylab = "frequency") +
				geom_vline(aes(xintercept=true_pars[["r"]]), linetype="dashed", size=1, color="grey50") +
				theme_bw()

print(rkernel)
ggsave(rkernel, filename="kernelr.pdf", height=3, width=3.25)

sigmapoints <- paramdata$sigma
sigmakernel <- qplot(sigmapoints, geom = "density", xlab = expression(sigma), ylab = "frequency") +
				geom_vline(aes(xintercept=sigma), linetype="dashed", size=1, color="grey50") +
				theme_bw()

print(sigmakernel)
ggsave(sigmakernel, filename="kernelsigma.pdf", height=3, width=3.25)

infecpoints <- paramdata$Iinit
infeckernel <- qplot(infecpoints, geom = "density", xlab = "Initial Infected", ylab = "frequency") +
				geom_vline(aes(xintercept=true_init_cond[['I']]), linetype="dashed", size=1, color="grey50") +
				theme_bw()

print(infeckernel)
ggsave(infeckernel, filename="kernelinfec.pdf", height=3, width=3.25)

# show grid
grid.arrange(R0kernel, rkernel, sigmakernel, infeckernel, ncol = 2, nrow = 2)

pdf("if2kernels.pdf", height = 6.5, width = 6.5)
grid.arrange(R0kernel, rkernel, sigmakernel, infeckernel, ncol = 2, nrow = 2)
dev.off()
#ggsave(filename="if2kernels.pdf", g2,  height=6.5, width=6.5)