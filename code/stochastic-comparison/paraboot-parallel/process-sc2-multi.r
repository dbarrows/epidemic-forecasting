library(ggplot2)
library(reshape2)
library(RColorBrewer)

load("sc2-multi.RData")

ssedata <- data.frame(SSEmat)
if2mat <- ssedata[,1:nTrials]
if2mat[if2mat == Inf] <- NaN
if2means <- rowMeans(if2mat, na.rm = TRUE)
hmcmat <- ssedata[,(nTrials+1):(2*nTrials)]
hmcmat[hmcmat == Inf] <- NaN
hmcmeans <- rowMeans(hmcmat, na.rm = TRUE)

if2meds <- apply(if2mat, 1, median)
hmcmeds <- apply(hmcmat, 1, median)

times <- 1:dim(ssedata)[1]

df <- data.frame(time = times, if2 = log10(if2means), hmc = log10(hmcmeans))
plotdata <- melt(df, id = "time")

## full
q <- qplot(data = plotdata, x = time, y = value, geom = "line", color = variable, xlab = "Truncation", ylab = expression(SSE~(log[10]))) +
		scale_colour_manual(values = c("black","grey"),
		                    labels = c("IF2", "HMCMC"),
		                    name = "Method") +
		theme_bw()

ggsave(q, filename="truncation.pdf", height=4, width=6.5)


df2 <- data.frame(time = times, if2 = log10(if2meds), hmc = log10(hmcmeds))
plotdata2 <- melt(df2, id = "time")

q2 <- qplot(data = plotdata2, x = time, y = value, geom = "line", color = variable, xlab = "Truncation", ylab = "log(SSE)") +
		scale_colour_manual(values = c("black","grey"),
		                    labels = c("IF2", "HMCMC"),
		                    name = "Method") +
		theme_bw()

ggsave(q2, filename="truncation2.pdf", height=4, width=6.5)

