library(ggplot2)

load("sc2-multi.RData")

ssedata <- data.frame(SSEmat)
if2mat <- ssedata[,1:nTrials]
if2means <- rowMeans(if2mat, na.rm = TRUE)
hmcmat <- ssedata[,(nTrials+1):(2*nTrials)]
hmcmeans <- rowMeans(hmcmat, na.rm = TRUE)


times <- 1:dim(ssedata)[1]

qplot(times, log(if2means), geom = "line", xlab = "Truncation", ylab = "log(SSE)") +
	geom_line(aes(y = log(hmcmeans)), linetype = "dashed") +
	theme_bw()