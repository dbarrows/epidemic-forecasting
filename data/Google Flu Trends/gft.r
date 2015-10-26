require(ggplot2)
library(scales)

print.eval  = TRUE

data <- read.csv(file="data.csv",head=TRUE,sep=",")

col_names <- names(data)

dates <- data[,1]
print(dates)

datlen <- length(dates)

mock_data = runif(datlen)

cdat1 = data[,2]
cdat9 = data[,9]

plot(dates,cdat9)

cdat9[is.na(cdat9)] <- 0

write(cdat9, ncolumns = 1, file = "Canada_gft_ili.dat")

#dates <- data[,1]
#print(dates)








# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------
if (FALSE) {
    cpu_means <- rep(0, 10);
    cublas_means <- rep(0,10);
    cuda_means <- rep(0,10);
    cudaopt_means <- rep(0,10);

    for (i in 1:10) {
    	start <- (i-1)*10;
    	end <- i*10;
    	cpu_means[i] = mean(run_data$CPU[start:end]);
    	cublas_means[i] = mean(run_data$CUBLAS[start:end]);
    	cuda_means[i] = mean(run_data$CUDA[start:end]);
    	cudaopt_means[i] = mean(run_data$CUDA.opt[start:end]);
    }

    xx <- 2^c(1:10);

    df <- data.frame(xx,cpu_means,cublas_means, cuda_means, cudaopt_means);

    ymax <-  max(cpu_means, cublas_means, cudaopt_means, cuda_means);
    ymin <-  min(cpu_means, cublas_means, cudaopt_means, cuda_means);
    xmax <- max(xx);
    xmin <- min(xx);
    ydiff <- ymax - ymin;
    xdiff <- xmax - xmin;
    yrat <- ymax/ymin;
    xrat <- xmax/xmin;

    g1 <- ggplot(df, aes(xx)) +                  			
        geom_line(aes(y=cpu_means, col="CPU")) +
        geom_line(aes(y=cublas_means, col="CUBLAS")) +
        geom_line(aes(y=cuda_means, col="CUDA")) +
        geom_line(aes(y=cudaopt_means, col="CUDA opt")) +
        xlab("Number of tiles (32x32 float)") +
        ylab("Run time (s)") +
        ggtitle("Number of tiles vs run times - log scale") +
        theme(legend.title=element_blank()) +
        scale_x_continuous(trans = log2_trans(),
            breaks = trans_breaks("log2", function(x) 2^x),
            labels = trans_format("log2", math_format(2^.x))) +
        scale_y_continuous(trans = log10_trans(),
            breaks = trans_breaks("log10", function(x) 10^x),
            labels = trans_format("log10", math_format(10^.x))) +
        theme(aspect.ratio = 1)

    g2 <- ggplot(df, aes(xx)) +                  			
        geom_line(aes(y=cpu_means, col="CPU")) +
        geom_line(aes(y=cublas_means, col="CUBLAS")) +
        geom_line(aes(y=cuda_means, col="CUDA")) +
        geom_line(aes(y=cudaopt_means, col="CUDA opt")) +
        ggtitle("Number of tiles vs run times") +
        xlab("Number of tiles (32x32 float)") +
        ylab("Run time (s)") +
        ggtitle("Tile size vs run times") +
        theme(legend.title=element_blank()) +
        theme(aspect.ratio = 1)

    ggsave(g1, filename="times1.pdf")
    ggsave(g2, filename="times2.pdf")
}



