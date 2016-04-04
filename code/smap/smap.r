## Dexter Barrows
## dbarrows.github.io
## McMaster University
## 2016

library(pracma) # needed for tiling function

smap <- function(data, E, theta, stepsAhead) {

    # construct library
    tseries <- as.vector(data)
    liblen  <- length(tseries) - E + 1 - stepsAhead
    lib     <- matrix(NA, liblen, E)

    for (i in 1:E) {
        lib[,i] <- tseries[(E-i+1):(liblen+E-i)]
    }

    # predict from the last index
    tslen <- length(tseries)
    predictee <- rev(t(as.matrix(tseries[(tslen-E+1):tslen])))
    predictions <- numeric(stepsAhead)

    # for each prediction index (number of steps ahead)
    for(i in 1:stepsAhead) {

        # set up weight calculation
        predmat <- repmat(predictee, liblen, 1)
        distances <- sqrt( rowSums( abs(lib - predmat)^2 ) )
        meanDist <- mean(distances)

        # calculate weights
        weights <- exp( - (theta * distances) / meanDist )

        # construct A, B

        preds <- tseries[(E+i):(liblen+E+i-1)]

        A <- cbind( rep(1.0, liblen), lib ) * repmat(as.matrix(weights), 1, E+1)
        B <- as.matrix(preds * weights)

        # solve system for C

        Asvd <- svd(A)
        C <- Asvd$v %*% diag(1/Asvd$d) %*% t(Asvd$u) %*% B

        # get prediction

        predsum <- sum(C * c(1,predictee))

        # save

        predictions[i] <- predsum

    }

    return(predictions)

}