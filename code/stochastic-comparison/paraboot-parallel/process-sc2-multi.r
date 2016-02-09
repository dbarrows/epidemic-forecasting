library(ggplot2)

load("sc2-multi.RData")

ssedata <- data.frame(SSEmat)
names(ssedata) <- c("IF2","HMCMC")

times <- 1:dim(ssedata)[1]

qplot(times, log(ssedata$IF2), geom = "line", xlab = "Truncation", ylab = "log(SSE)") +
	geom_line(aes(y = log(ssedata$HMCMC)), linetype = "dashed") +
	theme_bw()