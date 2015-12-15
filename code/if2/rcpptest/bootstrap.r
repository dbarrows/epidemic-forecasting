library(deSolve)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(Rcpp)
library(RColorBrewer)

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
Tlim 	<- T

## Generate true trajecory and synthetic data
##

true_init_cond <- c(S = N - i_infec,
                    I = i_infec,
                    R = 0)

true_pars <- c(R0 = 3.0,
               r = 0.1,
               N = 500.0)

odeout <- ode(true_init_cond, 0:(T-1), SIR, true_pars)
trueTraj <- odeout[,3]

set.seed(1000)

infec_counts_raw <- odeout[,3] + rnorm(T, 0, sigma)
infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

plotdata <- data.frame(times=1:T,true=trueTraj,data=infec_counts)

g <- ggplot(plotdata, aes(times)) +
        geom_line(aes(y = true, colour = "True")) +
        geom_point(aes(y = data, color = "Data")) +
        labs(x = "Time", y = "Infection count", color = "") +
        scale_color_brewer(palette="Paired") +
        theme(panel.background = element_rect(fill = "#F0F0F0"))

print(g)
ggsave(g, filename="dataplot.pdf", height=4, width=6.5)

## Rcpp stuff
##

sourceCpp(paste(getwd(),"if2.cpp",sep="/"))

paramdata <- data.frame(if2(infec_counts[1:Tlim], Tlim, N))
colnames(paramdata) <- c("R0", "r", "sigma", "Sinit", "Iinit", "Rinit")

## sample from parameter distributions
##

nTraj 	<- 100
datlen 	<- dim(paramdata)[1]
inds 	<- sample.int(datlen,nTraj,replace = TRUE)
params 	<- paramdata[inds,]

bootstrapdata <- matrix(NA, nrow = nTraj, ncol = T)

for (i in 1:nTraj) {

	init_cond <- c(S = params$Sinit[i],
	               I = params$Iinit[i],
	               R = params$Rinit[i])
	pars <- c(R0 = params$R0[i],
	          r = params$r[i],
	          N = 500.0)

	odeout <- ode(init_cond, 0:(T-1), SIR, pars)

	bootstrapdata[i,] <- odeout[,3]

}


meanTraj 	<- colMeans(bootstrapdata)
quantTraj 	<- apply(bootstrapdata, 2, quantile, probs = c(0.025,0.975))

datapart <- c(infec_counts[1:Tlim],rep(NA,T-Tlim))

plotdata <- data.frame(times=1:T,true=trueTraj,est=meanTraj,quants=t(quantTraj),datapart=datapart)

g <- ggplot(plotdata, aes(times)) +
		geom_ribbon(aes(ymin = quants.2.5., ymax=quants.97.5.), alpha=0.1) +
        geom_line(aes(y = true, colour = "True")) + 
        geom_line(aes(y = est, colour = "IF2")) +
        geom_point(aes(y = datapart, color = "Data")) +
        labs(x = "Time", y = "Infection count", color = "") +
        scale_color_brewer(palette="Paired") +
        theme(panel.background = element_rect(fill = "#F0F0F0"))

print(g)
ggsave(g, filename="if2plot.pdf", height=4, width=6.5)

f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
kcolours <- f("Paired")


trueval.R0 		<- 3.0
trueval.r 		<- 0.1
trueval.sigma 	<- 10.0
trueval.Iinit 	<- 5

meanval.R0 		<- mean(paramdata$R0)
meanval.r 		<- mean(paramdata$r)
meanval.sigma 	<- mean(paramdata$sigma)
meanval.Iinit 	<- mean(paramdata$Iinit)

linecolour <- "grey50"
lineweight <- 0.5

kerdataR0 <- data.frame(R0points = paramdata$R0)
R0kernel <- ggplot(kerdataR0, aes(x = R0points, y = ..scaled..)) +
                geom_density(color = kcolours[1], fill = kcolours[1]) +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = expression(R[0]), y = "Density", color = "") +
                geom_vline(aes(xintercept=trueval.R0), linetype="solid", size=lineweight, color=linecolour) +
                geom_vline(aes(xintercept=meanval.R0), linetype="dashed", size=lineweight, color=linecolour)

kerdatar <- data.frame(rpoints = paramdata$r)
rkernel <- ggplot(kerdatar, aes(x = rpoints, y = ..scaled..)) +
                geom_density(color = kcolours[2], fill = kcolours[2]) +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = "r", y = "", color = "") +
                geom_vline(aes(xintercept=trueval.r), linetype="solid", size=lineweight, color=linecolour) +
                geom_vline(aes(xintercept=meanval.r), linetype="dashed", size=lineweight, color=linecolour)

kerdatasigma <- data.frame(sigmapoints = paramdata$sigma)
sigmakernel <- ggplot(kerdatasigma, aes(x = sigmapoints, y = ..scaled..)) +
                geom_density(color = kcolours[3], fill = kcolours[3]) +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = expression(sigma), y = "Density", color = "") +
                geom_vline(aes(xintercept=trueval.sigma), linetype="solid", size=lineweight, color=linecolour) +
                geom_vline(aes(xintercept=meanval.sigma), linetype="dashed", size=lineweight, color=linecolour)

kerdatainfec <- data.frame(infecpoints = paramdata$Iinit)
infeckernel <- ggplot(kerdatainfec, aes(x = infecpoints, y = ..scaled..)) +
                geom_density(color = kcolours[4], fill = kcolours[4]) +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = "Initial infected", y = "", color = "") +
                geom_vline(aes(xintercept=trueval.Iinit), linetype="solid", size=lineweight, color=linecolour) +
                geom_vline(aes(xintercept=meanval.Iinit), linetype="dashed", size=lineweight, color=linecolour)

# show grid
grid.arrange(R0kernel, rkernel, sigmakernel, infeckernel, ncol = 2, nrow = 2)

pdf("if2kernels.pdf", height = 6.5, width = 6.5)
grid.arrange(R0kernel, rkernel, sigmakernel, infeckernel, ncol = 2, nrow = 2)
dev.off()