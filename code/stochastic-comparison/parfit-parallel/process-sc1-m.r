library(ggplot2)
library(gridExtra)

load("sc1-multi-time.RData")

# extract if2 data

if2names <- c("if2.R0","if2.r","if2.sigma","if2.eta","if2.berr","if2.Iinit")
if2data <- estmat[,if2names]
colnames(if2data) <- c("R0","r","sigma","eta","berr","Iinit")
if2times <- estmat[,7]

# extract hmc data

hmcnames <- c("hmc.R0","hmc.r","hmc.sigma","hmc.eta","hmc.berr","hmc.Iinit")
hmcdata <- estmat[,hmcnames]
colnames(hmcdata) <- c("R0","r","sigma","eta","berr","Iinit")
hmctimes <- estmat[,14]

# average times

avif2time <- mean(if2times)
avhmctime <- mean(hmctimes)

avif2time
avhmctime

# sort results

if2sorted <- apply(if2data, 2, sort)
hmcsorted <- apply(hmcdata, 2, sort)

## take centre 95%

nTrials <- dim(estmat)[1]

cinds <- (0.025*nTrials+1):(0.975*nTrials)
if295 <- if2sorted[cinds,]
hmc95 <- hmcsorted[cinds,]

if2plotdata <- data.frame(if295)
hmcplotdata <- data.frame(hmc95)

linecolour <- "grey50"
lineweight <- 0.5


## if2 density plots

R0kernel <- qplot(if2plotdata$R0, geom = "density", xlab = expression(R[0]), ylab = "frequency") +
    geom_vline(aes(xintercept=pars_true['R0']), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=mean(if2plotdata$R0)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

rkernel <- qplot(if2plotdata$r, geom = "density", xlab = "r", ylab = "") +
    geom_vline(aes(xintercept=pars_true['r']), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=mean(if2plotdata$r)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

sigmakernel <- qplot(if2plotdata$sigma, geom = "density", xlab = expression(sigma), ylab = "") +
    geom_vline(aes(xintercept=sigma), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=mean(if2plotdata$sigma)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

infeckernel <- qplot(if2plotdata$Iinit, geom = "density", xlab = "Initial Infected", ylab = "frequency") +
    geom_vline(aes(xintercept=i_infec), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=mean(if2plotdata$Iinit)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

etakernel <- qplot(if2plotdata$eta, geom = "density", xlab = expression(eta), ylab = "") +
    geom_vline(aes(xintercept=pars_true['eta']), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=mean(if2plotdata$eta)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

berrkernel <- qplot(if2plotdata$berr, geom = "density", xlab = expression(epsilon[proc]), ylab = "") +
    geom_vline(aes(xintercept=pars_true['berr']), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=mean(if2plotdata$berr)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()


# show grid
grid.arrange(R0kernel, rkernel, sigmakernel, infeckernel, etakernel, berrkernel, ncol = 3, nrow = 2)





## hmc density plots

R0kernel <- qplot(hmcplotdata$R0, geom = "density", xlab = expression(R[0]), ylab = "frequency") +
    geom_vline(aes(xintercept=pars_true['R0']), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=mean(hmcplotdata$R0)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

rkernel <- qplot(hmcplotdata$r, geom = "density", xlab = "r", ylab = "") +
    geom_vline(aes(xintercept=pars_true['r']), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=mean(hmcplotdata$r)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

sigmakernel <- qplot(hmcplotdata$sigma, geom = "density", xlab = expression(sigma), ylab = "") +
    geom_vline(aes(xintercept=sigma), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=mean(hmcplotdata$sigma)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

infeckernel <- qplot(hmcplotdata$Iinit, geom = "density", xlab = "Initial Infected", ylab = "frequency") +
    geom_vline(aes(xintercept=i_infec), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=mean(hmcplotdata$Iinit)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

etakernel <- qplot(hmcplotdata$eta, geom = "density", xlab = expression(eta), ylab = "") +
    geom_vline(aes(xintercept=pars_true['eta']), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=mean(hmcplotdata$eta)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()

berrkernel <- qplot(hmcplotdata$berr, geom = "density", xlab = expression(epsilon[proc]), ylab = "") +
    geom_vline(aes(xintercept=pars_true['berr']), linetype="solid", size=lineweight, color=linecolour) +
    geom_vline(aes(xintercept=mean(hmcplotdata$berr)), linetype="dashed", size=lineweight, color=linecolour) +
    theme_bw()


# show grid
grid.arrange(R0kernel, rkernel, sigmakernel, infeckernel, etakernel, berrkernel, ncol = 3, nrow = 2)



## combined plots

R0kernel <- qplot(hmcplotdata$R0, geom = "density", xlab = expression(R[0]), ylab = "frequency") +
    geom_density(aes(x = if2plotdata$R0), linetype = "dashed") +
    geom_vline(aes(xintercept=pars_true['R0']), linetype="solid", size=lineweight, color=linecolour) +
    theme_bw()

rkernel <- qplot(hmcplotdata$r, geom = "density", xlab = "r", ylab = "") +
    geom_density(aes(x = if2plotdata$r), linetype = "dashed") +
    geom_vline(aes(xintercept=pars_true['r']), linetype="solid", size=lineweight, color=linecolour) +
    theme_bw()

sigmakernel <- qplot(hmcplotdata$sigma, geom = "density", xlab = expression(sigma)), ylab = "") +
    geom_density(aes(x = if2plotdata$sigma), linetype = "dashed") +
    geom_vline(aes(xintercept=sigma), linetype="solid", size=lineweight, color=linecolour) +
    theme_bw()

infeckernel <- qplot(hmcplotdata$Iinit, geom = "density", xlab = "Initial Infected", ylab = "frequency") +
    geom_density(aes(x = if2plotdata$Iinit), linetype = "dashed") +
    geom_vline(aes(xintercept=i_infec), linetype="solid", size=lineweight, color=linecolour) +
    theme_bw()

etakernel <- qplot(hmcplotdata$eta, geom = "density", xlab = expression(eta), ylab = "") +
    geom_density(aes(x = if2plotdata$eta), linetype = "dashed") +
    geom_vline(aes(xintercept=pars_true['eta']), linetype="solid", size=lineweight, color=linecolour) +
    theme_bw()

berrkernel <- qplot(hmcplotdata$berr, geom = "density", xlab = expression(epsilon[proc]), ylab = "") +
    geom_density(aes(x = if2plotdata$berr), linetype = "dashed") +
    geom_vline(aes(xintercept=pars_true['berr']), linetype="solid", size=lineweight, color=linecolour) +
    theme_bw()

grid.arrange(R0kernel, rkernel, sigmakernel, infeckernel, etakernel, berrkernel, ncol = 3, nrow = 2)
