library(Rcpp)
library(ggplot2)

source(paste(getwd(), "../../sir-functions", "StocSIRS.r", sep = "/"))

#set.seed(6287502)

T 		<- 200
i_infec <- 5
steps 	<- 7
N 		<- 500
sigma 	<- 5

true_pars <- c(R0 = 3.0,  	# new infected people per infected person
          	   r = 0.1, 	# recovery rate
		  	   N = 500,    	# population size
		  	   eta = 0.5, 	# geometric random walk
		  	   berr = 0.5, 	# Beta geometric walk noise
          	   re = 0.05)    # resuceptibility rate

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

NP          <- 2500
nPasses     <- 50
coolrate    <- 0.975

if2_sirs_file = paste(getwd(), "if2-sirs.cpp", sep = "/")
sourceCpp(if2_sirs_file)

if2time <- system.time( if2data <- if2_sirs(infec_counts, T+1, N, NP, nPasses, coolrate) )

paramdata <- if2data$paramdata
names(paramdata) <- c("R0","r","re","sigma","eta","berr","Sinit","Iinit","Rinit")
parmeans <- colMeans(paramdata)
names(parmeans) <- c("R0","r","re","sigma","eta","berr","Sinit","Iinit","Rinit")

pars_fit <- with( as.list(parmeans),
              c(R0 = R0,
      			r = r,
      			N = N,
      			eta = eta,
      			berr = berr,
      			re = re) )

init_cond_fit <- with( as.list(parmeans),
                   c(S = Sinit,
               		 I = Iinit,
               		 R = Rinit) )

sdeout_fit <- StocSIRS(init_cond_fit, pars_fit, T, steps)
colnames(sdeout_fit) <- c('S','I','R','B')

#quartz()
qplot(0:T, sdeout[,'I'], geom = "line", xlab = "week", ylab = "Infection Count") +
	geom_line(aes(y = sdeout_fit[,'I']), linetype = "dashed") +
	geom_line(aes(y = if2data$statemeans[,2]), linetype = "dotted") +
    geom_point(aes(y = infec_counts)) +
    theme_bw()