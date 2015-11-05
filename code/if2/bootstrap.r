library(deSolve)
library(ggplot2)
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

data <- read.table("plotdata.dat", col.names = c("times","true", "data","estimate") )
paramdata <- read.table("pfdata.dat", col.names = c("R0","r", "sigma","Sinit","Iinit","Rinit") )

datapart <- data$data[data$data >= 0]
Tlim <- length(datapart)



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

plotdata <- data.frame(times=1:(T+1),true=trueTraj,est=meanTraj,sds=sdTraj,datapart=c(datapart,rep(NA,T-Tlim+1)))

g <- ggplot(plotdata, aes(times)) +
		geom_ribbon(aes(ymin = est-sds, ymax=est+sds), alpha=0.1) +
        geom_line(aes(y = true, colour = "True")) + 
        geom_line(aes(y = est, colour = "IF2")) +
        geom_point(aes(y = datapart, color = "Data")) +
        labs(x = "Time", y = "Infection count", color = "") +
        scale_color_brewer(palette="Paired") +
        theme(panel.background = element_rect(fill = "#F0F0F0"))

print(g)
ggsave(g, filename="if2plot.pdf", height=4, width=6.5)


trueval.R0 <- 3.0
trueval.r <- 0.1
trueval.sigma <- 10.0
trueval.Iinit <- 5

kerdataR0 <- data.frame(R0points = paramdata$R0)
R0kernel <- ggplot(kerdataR0, aes(x = R0points, y = ..scaled..)) +
                geom_density(colour="#8dd3c7", fill="#8dd3c7") +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = expression(R[0]), y = "Density", color = "") +
                geom_vline(aes(xintercept=trueval.R0), linetype="dashed", size=1, color="grey50")

print(R0kernel)

kerdatar <- data.frame(rpoints = paramdata$r)
rkernel <- ggplot(kerdatar, aes(x = rpoints, y = ..scaled..)) +
                geom_density(colour="#ffffb3", fill="#ffffb3") +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = "r", y = "Density", color = "") +
                geom_vline(aes(xintercept=trueval.r), linetype="dashed", size=1, color="grey50")

print(rkernel)

kerdatasigma <- data.frame(sigmapoints = paramdata$sigma)
sigmakernel <- ggplot(kerdatasigma, aes(x = sigmapoints, y = ..scaled..)) +
                geom_density(colour="#bebada", fill="#bebada") +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = expression(sigma), y = "Density", color = "") +
                geom_vline(aes(xintercept=trueval.sigma), linetype="dashed", size=1, color="grey50")

print(sigmakernel)

kerdatainfec <- data.frame(infecpoints = paramdata$Iinit)
infeckernel <- ggplot(kerdatainfec, aes(x = infecpoints, y = ..scaled..)) +
                geom_density(colour="#fb8072", fill="#fb8072") +
                theme(panel.background = element_rect(fill = "#F0F0F0")) +
                scale_color_brewer(palette="Paired") +
                labs(x = "Initial infected", y = "Density", color = "") +
                geom_vline(aes(xintercept=trueval.Iinit), linetype="dashed", size=1, color="grey50")

print(infeckernel)