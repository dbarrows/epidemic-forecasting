# install.packages("rstan")

Sys.setenv(MAKEFLAGS = "-j4") 
source('http://mc-stan.org/rstan/install.R', echo = TRUE, max.deparse.length = 2000)
install_rstan()

# alternative code to install
options(repos = c(getOption("repos"), rstan = "http://rstan.org/repo/"))
install.packages("rstan", type = 'source')

