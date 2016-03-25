library(rstan)

stanbootmean <- function(fit, nTraj, N, Tlim, steps) {

	source("../../sir-functions/StocSIRSstan.r")

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

	parnames <- c("R0","r","re","sigma","eta","berr","Iinit")
	parvars <- var(paramdata[parnames])
	hmcparmeans <- colMeans(paramdata[parnames])
	names(hmcparmeans) <- paste("hmc.", parnames,sep="")

	bootstrapdata <- matrix(NA, nrow = nTraj, ncol = (T+1) )

	pardatlen 	<- dim(paramdata)[1]
	inds 	    <- sample.int(pardatlen,nTraj,replace = TRUE)
	params 	    <- paramdata[inds,]

	# get mean final states
	final_states <- matrix(NA, nTraj, 4)

	for (i in 1:nTraj) {

		paramset <- params[i,]

		init_cond <- c(S = N - paramset$Iinit,
		               I = paramset$Iinit,
		               R = 0.0 )
		pars <- c(R0 = paramset$R0,
		          r = paramset$r,
		          re = paramset$re,
		          N = 500.0,
		          eta = paramset$eta,
		          berr = paramset$berr)

		berrvec <- numeric(datlen)
		for (j in 1:datlen) {
			varname <- paste("Bnoise[", j, "]", sep = "")
			berrvec[j] <- paramset[[varname]]
		}

		sdeout <- StocSIRSstan(init_cond, pars, Tlim, steps, berrvec, datlen)
		colnames(sdeout) <- c('S','I','R','B')

		final_states[i,] <- sdeout[Tlim+1,]

		bootstrapdata[i,1:(Tlim+1)] <- sdeout[,'I']

	}

	final_state_means <- colMeans(final_states)
	names(final_state_means) <- c('S','I','R','B')

	# bootstrap into future

	for (i in 1:nTraj) {

		paramset <- params[i,]

		init_cond <- c(	S = final_state_means[['S']],
		               	I = final_state_means[['I']],
		               	R = final_state_means[['R']])
		pars <- c(R0 = paramset$R0,
		          r = paramset$r,
		          re = paramset$re,
		          N = 500.0,
		          eta = paramset$eta,
		          berr = paramset$berr)

		berrvec <- numeric(datlen)
		for (j in 1:datlen) {
			varname <- paste("Bnoise[", j, "]", sep = "")
			berrvec[j] <- paramset[[varname]]
		}

		sdeout <- StocSIRSstan(init_cond, pars, T-Tlim, steps, berrvec, datlen)
		colnames(sdeout) <- c('S','I','R','B')

		#print( length( bootstrapdata[i,(Tlim+2):(T+1)] ) )
		#print( length( sdeout[-1,'I'] ) )
		bootstrapdata[i,(Tlim+2):(T+1)] <- sdeout[-1,'I']

	}

	# in case of explosion
	bootstrapdata <- bootstrapdata[complete.cases(bootstrapdata),]
	goodbootsize <- dim(bootstrapdata)[1]

	return( bootstrapdata )

}