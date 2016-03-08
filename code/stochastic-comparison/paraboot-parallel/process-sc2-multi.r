library(ggplot2)
library(reshape2)
library(RColorBrewer)

load("sc2-multi.RData")

ssedata <- data.frame(SSEmat)
if2mat <- ssedata[,1:nTrials]
if2means <- rowMeans(if2mat, na.rm = TRUE)
hmcmat <- ssedata[,(nTrials+1):(2*nTrials)]
hmcmeans <- rowMeans(hmcmat, na.rm = TRUE)

times <- 1:dim(ssedata)[1]

df <- data.frame(time = times, if2 = log(if2means), hmc = log(hmcmeans))
plotdata <- melt(df, id = "time")

## full
q <- qplot(data = plotdata, x = time, y = value, geom = "line", color = variable, xlab = "Truncation", ylab = "log(SSE)") +
		scale_colour_manual(values = c("black","grey"),
		                    labels = c("IF2", "HMCMC"),
		                    name = "Method") +
		theme_bw()

ggsave(q, filename="truncation.pdf", height=4, width=6.5)
