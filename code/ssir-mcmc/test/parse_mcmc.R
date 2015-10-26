execute_C <- function(exec=file.path(".","sir-mcmc"),
                        tmpfile=tempfile())
{
    command=paste(exec,"data",">",tmpfile)
    system(command)
    read.table(tmpfile)
}

data <- execute_C()

plot( density(data$V1) )
plot( density(data$V2) )
plot( density(data$V3) )
plot( density(data$V4) )
plot( density(data$V5) )
plot( density(data$V6) )

R0mean = mean(data$V1); R0mean
rmean = mean(data$V2); rmean
sigmean = mean(data$V3); sigmean
Smean = mean(data$V4); Smean
Imean = mean(data$V5); Imean
Rmean = mean(data$V6); Rmean