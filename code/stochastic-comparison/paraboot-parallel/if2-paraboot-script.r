# Dexter Barrows
#
# run: Rscript if2-paraboot-script.r > if2-paraboot-script.log 2>&1

library(foreach)
library(parallel)
library(doParallel)
library(ggplot2)
library(Rcpp)

set.seed(8265813)

# R source files
stoc_sir_file 		<- paste(getwd(),"../../sir-functions", "StocSIR.r", sep = "/")
if2_paraboot_file 	<- paste(getwd(), "if2-paraboot.r", sep = "/")
source(stoc_sir_file)
source(if2_paraboot_file)

# C++ source files
if2file <- paste(getwd(),"../../if2", "if2.cpp", sep="/")
if2_s_file <- paste(getwd(),"../../if2", "if2-s.cpp", sep="/")

## data params
T       <- 60
i_infec <- 5
steps   <- 7
N       <- 500
sigma   <- 10
Tlim    <- round(T * 1/2)

## Generate true trajecory and synthetic data
##

pars_true <- c(R0 = 3.0,    # new infected people per infected person
          r = 0.1,      # recovery rate
          N = 500,      # population size
          eta = 0.5,    # geometric random walk
          berr = 0.5)   # Beta geometric walk noise

true_init_cond <- c(S = N - i_infec,
                    I = i_infec,
                    R = 0)

# get true trajectory
sdeout_true <- StocSIR(true_init_cond, pars_true, T, steps)
colnames(sdeout_true) <- c('S','I','R','B')

# get initial data
infec_counts_raw <- sdeout_true[,'I'] + rnorm(T+1, 0, sigma)
infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

## parallel setup
##

nTrials <- 10

numCores <- detectCores()
cl <- makeCluster(min(numCores, nTrials))
registerDoParallel(cl)

## if2 setup
##

NP          <- 3000
nPasses     <- 50
coolrate    <- 0.975

## initial fit
##



#if (FALSE) {
sourceCpp(if2file)
if2time <- system.time( if2data <- if2(infec_counts[1:(Tlim+1)], Tlim+1, N, NP, nPasses, coolrate) )

# get parameter estimates

paramdata <- data.frame( if2data$paramdata )
names(paramdata) <- c("R0", "r", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")
parmeans <- colMeans(paramdata)
names(parmeans) <- c("R0", "r", "sigma", "eta", "berr", "Sinit", "Iinit", "Rinit")

if2_paraboot_data <- if2_paraboot(if2data,
                                  T, Tlim, steps, N, nTrials,
                                  if2file, if2_s_file, stoc_sir_file,
                                  NP, nPasses, coolrate)

stopCluster(cl)

save.image(file = "if2-paraboot.RData")