# TODO
# add sourcing
# rewrite
# create figures

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXAMINE OUTPUT AND ASSESS CONVERGENCE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OUT <- list(`out-1`, `out-2`, `out-3`) %>% as.mcmc.list()
library(coda)
OUT <- out
rhat <- gelman.diag(OUT, multivariate = FALSE) 
rhat.proc <- rhat$psrf %>% as.data.frame() %>% filter(!is.nan(`Point est.`))


## CONVERT SUMMARY TO A DATA FRAME
OUT<-as.data.frame(IPMrun$summary$all.chains)
OUT$parameter<-row.names(OUT)


## CONVERT SAMPLES TO MCMC LIST FOR USE IN CODA TO CALCULATE R-HAT
IPMdiag<-as.mcmc.list(IPMrun$samples)
gelman.plot(IPMdiag, confidence = 0.95, transform=FALSE, autoburnin=TRUE,multivariate=TRUE)
traceplot(IPMdiag, smooth = FALSE, col = 1:6, type = "l", xlab = "Iterations", ylab = "")
#coda::plot.mcmc(EVdiag)
Rhat<-gelman.diag(IPMdiag, confidence = 0.95, transform=FALSE, autoburnin=TRUE,multivariate=TRUE)
OUT$R_hat<-Rhat$psrf[,1]


## EXPORT OUTPUT
write.table(OUT,"YNAL_estimates_nimble.csv", sep=",")

