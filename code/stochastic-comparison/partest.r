library(foreach)
library(parallel)
library(doParallel)

# setup cluster
numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)

matrix <- foreach(i = 1:10, .combine = rbind) %dopar% {

	retvec <- i*(1:10)

	return(retvec)

}

matrix