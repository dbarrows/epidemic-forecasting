library(pracma)

# data will be a mtrix where the rows are individual time series
# (assumes all time series same length, this could be relaxed in the future)

dewdrop <- function(data, E, theta, stepsAhead) {

    # construct library
    #tseries <- as.matrix(data)
    liblen  <- dim(data)[2] - E + 1 - stepsAhead
    libnum  <- dim(data)[1]
    lib     <- matrix(NA, libnum*liblen, E)
    preds 	<- matrix(NA, libnum*liblen, stepsAhead)

    ## make library and predictions array (flattened)

    libctr <- 1
    for (lidx in 1:libnum) {
    	tseries <- data[lidx,]
	    for (i in 1:liblen) {
	    	lstart <- i
	    	lstop  <- lstart + E - 1
	    	lib[libctr,] <- tseries[lstart:lstop]
	    	pstart <- lstop + 1
	    	pstop  <- pstart + stepsAhead - 1
	    	preds[libctr,] <- tseries[pstart:pstop]
	    	libctr <- libctr + 1
	    }
	}

	predictees  <- matrix(NA, libnum, E)

    # predict from the last index
	tslen <- dim(data)[2]
    for (lidx in 1:libnum) {
    	tseries <- data[lidx,]
	    predictees[lidx,] <- data[lidx,(tslen-E+1):(tslen)]
	}

	## predictions matrix

	predictions <- matrix(NA, libnum, stepsAhead)

    # for each prediction index (number of steps ahead)
	for (preidx in 1:libnum) {

		predictee <- predictees[preidx,]

	    for(i in 1:stepsAhead) {

	        # set up weight calculation
	        predmat <- repmat(predictee, liblen*libnum, 1)
	        #print(dim(predmat))
	        #print(dim(lib))
	        distances <- sqrt( rowSums( abs(lib - predmat)^2 ) )
	        meanDist <- mean(distances)

	        # calculate weights
	        weights <- exp( - (theta * distances) / meanDist )

	        # construct A, B

	        predvec <- preds[,i]

	        A <- cbind( rep(1.0, liblen*libnum), lib ) * repmat(as.matrix(weights), 1, E+1)
	        B <- as.matrix(predvec * weights)

	        # solve system for C

	        Asvd <- svd(A)
	        C <- Asvd$v %*% diag(1/Asvd$d) %*% t(Asvd$u) %*% B

	        # get prediction

	        predsum <- sum(C * c(1,predictee))

	        # save

	        predictions[preidx,i] <- predsum

	    }

	}

    return(predictions)

}