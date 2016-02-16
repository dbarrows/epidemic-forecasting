library(Rcpp)
library(ggplot2)
library(rstan)
library(reshape2)
library(beepr)

## External files
source(paste(getwd(), "../../sir-functions", "StocSIRSstan.r", sep = "/"))
stanfile <- "sirs-euler.stan"

#set.seed(6287502)

T 		<- 200
i_infec <- 5
steps 	<- 7
N 		<- 500
sigma 	<- 3

true_pars <- c(R0 = 3.0,  	# new infected people per infected person
          	   r = 0.1, 	# recovery rate
		  	   N = 500,    	# population size
		  	   eta = 0.5, 	# geometric random walk
		  	   berr = 0.5, 	# Beta geometric walk noise
          	   re = 1)    # resuceptibility rate

true_init_cond <- c(S = N - i_infec,
					I = i_infec,
					R = 0)

sdeout <- StocSIRS(true_init_cond, true_pars, T, steps)
colnames(sdeout) <- c('S','I','R','B')

infec_counts_raw <- sdeout[,'I'] + rnorm(T+1, 0, sigma)
infec_counts     <- ifelse(infec_counts_raw < 0, 0, infec_counts_raw)

#quartz()
qplot(0:T, sdeout[,'I'], geom = "line", xlab = "week", ylab = "Infection Count") +
    geom_point(aes(y = infec_counts)) +
    theme_bw()

## Stan setup
##

datlen <- T*steps + 1

data <- matrix(data = -1, nrow = T+1, ncol = steps)
data[,1] <- infec_counts
standata <- as.vector(t(data))[1:datlen]

sir_data <- list( T = datlen,   	# simulation time
                  y = standata, 	# infection count data
                  N = 500,      	# population size
                  h = 1/steps )   	# step size per day 
                    
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stan_options <- list(   chains = 1,    		# number of chains
                        iter   = 2000, 		# iterations per chain
                        warmup = 1000, 		# warmup interations
                        thin   = 1 )   		# thinning number


## fit
##
hmctime <- system.time(fit <- with(stan_options,
                        stan( file    = stanfile,
                              data    = sir_data,
                              chains  = chains,
                              iter    = iter,
                              warmup  = warmup,
                              thin    = thin)
                        )
            )

# Done beep
beep()


## extract fit
##
exfit <- extract(fit, permuted = FALSE, inc_warmup = FALSE)

paramdata <- data.frame(R0 = melt(exfit[,,'R0'])$value,
                        r = melt(exfit[,,'r'])$value,
                        re = melt(exfit[,,'re'])$value,
                        sigma = melt(exfit[,,'sigma'])$value,
                        eta = melt(exfit[,,'eta'])$value,
                        berr = melt(exfit[,,'berr'])$value,
                        Iinit = melt(exfit[,,'Iinit'])$value )

for (j in 1:datlen) {
	varname <- paste('Bnoise[', j, ']', sep = '')
	paramdata[[varname]] <- melt( exfit[,,varname] )$value
}



#parnames <- c("R0","r","re","sigma","eta","berr","Iinit")
#parvars <- var(paramdata[parnames])
parmeans <- colMeans(paramdata)
#names(parmeans) <- parnames

berrvec <- numeric(datlen)
for (j in 1:datlen) {
	varname <- paste("Bnoise[", j, "]", sep = "")
	berrvec[j] <- parmeans[[varname]]
}

## Setup for generated trajectory
pars_fit <- with( as.list(parmeans),
              c(R0 = R0,
      			r = r,
      			N = N,
      			eta = eta,
      			berr = berr,
      			re = re) )

init_cond_fit <- with( as.list(parmeans),
                   c(S = N - Iinit,
               		 I = Iinit,
               		 R = 0.0) )

sdeout_fit <- StocSIRSstan(init_cond_fit, pars_fit, T, steps, berrvec, datlen)
colnames(sdeout_fit) <- c('S','I','R','B')

#quartz()
qplot(0:T, sdeout[,'I'], geom = "line", xlab = "week", ylab = "Infection Count") +
	geom_line(aes(y = sdeout_fit[,'I']), linetype = "dashed") +
	#geom_line(aes(y = if2data$statemeans[,2]), linetype = "dotted") +
    geom_point(aes(y = infec_counts)) +
    theme_bw()