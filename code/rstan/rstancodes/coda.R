library(coda)

coda.obj = mcmc.list(lapply (1:ncol(fit), function(x) {mcmc(as.array(fit)[,x,])}))

summary(coda.obj)
summary(coda.obj)$statistics
summary(coda.obj)$statistics[,'Mean']
summary(coda.obj)$statistics[,'SD']
summary(coda.obj)$statistics[,'Naive SE']

summary(coda.obj)$quantiles

summary(coda.obj)$start
summary(coda.obj)$end

summary(coda.obj)$thin

summary(coda.obj)$nchain

traceplot(coda.obj)
