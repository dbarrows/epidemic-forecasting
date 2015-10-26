# load SSR functions
source("SSR_functions.R")

# CCM analysis for Didinium-Paramecium
if(FALSE)
{
 # load data
  model_data <- read.table("vr.raw.txt")
  para <- model_data[11:71,2]
  didi <- model_data[11:71,3]
  lib <- c(1, 60)
  pred <- c(1, 60)
  lib_sizes <- seq(5, 60, 5)
  E <- 3
  
  x_xmap_y <- ccm(para, didi, lib_sizes, lib, pred, E)
  y_xmap_x <- ccm(didi, para, lib_sizes, lib, pred, E)
  
  # compute mean rhos at each L
  x_xmap_y$L <- as.factor(x_xmap_y$L)
  x_xmap_y_means <- do.call(rbind, lapply(split(x_xmap_y, x_xmap_y$L), function(x){max(0, mean(x$rho))}))
  y_xmap_x$L <- as.factor(y_xmap_x$L)
  y_xmap_x_means <- do.call(rbind, lapply(split(y_xmap_x, y_xmap_x$L), function(x){max(0, mean(x$rho))}))
  
  # produce plot (Figure 5D from Sugihara et al. 2012)
  plot(lib_sizes, x_xmap_y_means, type = "l", col = "royalblue", lwd = 2, 
       xlim = c(min(lib_sizes), max(lib_sizes)), ylim = c(0.4, 1.0), xlab = "L", ylab = expression(rho))
  lines(lib_sizes, y_xmap_x_means, col = "red3", lwd = 2)
  legend(max(lib_sizes), 0.4, legend = c("para xmap didi", "didi xmap para"), 
         xjust = 1, yjust = 0, lty = 1, lwd = 2, col = c("royalblue", "red3"))
}

# CCM analysis for 2sp model time
if(FALSE)
{
  # load data
  model_data <- read.table("2sp.txt")
  x <- model_data[,2]
  y <- model_data[,3]
  lib <- c(1, 1000)
  pred <- c(1, 1000)
  lib_sizes <- c(50, 100, 200, 350, 500, 1000)
  E <- 2
  
  x_xmap_y <- ccm(x, y, lib_sizes, lib, pred, E)
  y_xmap_x <- ccm(y, x, lib_sizes, lib, pred, E)
  
  # compute mean rhos at each L
  x_xmap_y$L <- as.factor(x_xmap_y$L)
  x_xmap_y_means <- do.call(rbind, lapply(split(x_xmap_y, x_xmap_y$L), function(x){max(0, mean(x$rho))}))
  y_xmap_x$L <- as.factor(y_xmap_x$L)
  y_xmap_x_means <- do.call(rbind, lapply(split(y_xmap_x, y_xmap_x$L), function(x){max(0, mean(x$rho))}))
  
  # produce plot (Figure 5D from Sugihara et al. 2012)
  plot(lib_sizes, x_xmap_y_means, type = "l", col = "royalblue", lwd = 2, 
       xlim = c(min(lib_sizes), max(lib_sizes)), ylim = c(0.0, 1.0), xlab = "L", ylab = expression(rho))
  lines(lib_sizes, y_xmap_x_means, col = "red3", lwd = 2)
  legend(max(lib_sizes), 0.0, legend = c("x xmap y", "y xmap x"), 
         xjust = 1, yjust = 0, lty = 1, lwd = 2, col = c("royalblue", "red3"))
}

# CCM analysis for 5sp model
if(TRUE)
{
  # load data
  model_data <- read.table("5sp.txt")
  lib_sizes <- c(50, 100, 200, 300)
  lib <- c(1, 300)
  pred <- c(301, 600)
  E <- 3
  
  par(mfrow = c(5, 5), mar = c(1,1,1,1), oma = c(2, 4, 2, 0))
  for(i in 1:5)
  {
    for(j in 1:5)
    {
      if(i == j)
      {
        plot(NA, NA, type = "n", 
             xlim = c(min(lib_sizes), max(lib_sizes)), ylim = c(0, 1), 
             xaxt = if(i<5) "n" else "s", yaxt = if(j==1) "s" else "n")
        text(min(lib_sizes)+max(lib_sizes)/2, 0.5, paste("species ", i, sep = ""), adj = c(0.5, 0.5)) 
      }
      else
      {
        x <- model_data[1:600, i]
        y <- model_data[1:600, j]
        
        x_xmap_y <- ccm(x, y, lib_sizes, lib, pred, E)
        
        # compute mean rhos at each L
        x_xmap_y$L <- as.factor(x_xmap_y$L)
        x_xmap_y_means <- do.call(rbind, lapply(split(x_xmap_y, x_xmap_y$L), function(x){max(0, mean(x$rho))}))
        
        # produce plot (Figure 5D from Sugihara et al. 2012)
        plot(lib_sizes, x_xmap_y_means, type = "l", col = "royalblue", lwd = 2, 
                     xlim = c(min(lib_sizes), max(lib_sizes)), ylim = c(0.0, 1.0), 
             xlab = "L", ylab = expression(rho), 
             xaxt = if(i<5) "n" else "s", yaxt = if(j==1) "s" else "n")
      }
    }
  }
  
  mtext(c("cross mapping from", "cross mapping to"), c(2, 3), line = c(2, 1), outer = TRUE)
}