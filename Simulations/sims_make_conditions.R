##########################################
## CONDITIONS
##########################################

PATH <- "path where your source files are located"

noisevalues <- c(3, 5.5, 7)
alphavalues <- c(0.001, 0.05)
betavalues <- c(0.1, 0.2, 0.3)

noise <- rep(noisevalues, each=length(alphavalues)*length(betavalues))
alpha <- rep(alphavalues, each=length(betavalues))
beta <- rep(betavalues, times = length(noisevalues)*length(alphavalues))


conditions <- as.data.frame(cbind(noise,alpha,beta))
colnames(conditions) <- c("noise", "alpha","beta")
conditions


write.table(conditions,file=paste(PATH,"conditions_results1.txt",sep=""),row.names=FALSE,col.names=TRUE)


## Conditions for effect size computation in sims_results_anova.R

scanvalues <- c(50, 100, 150)
noisevalues <- c(3, 5.5, 7)
alphavalues <- c(0.001, 0.05)
betavalues <- c(0.1, 0.2, 0.3)

scans <- rep(scanvalues, each=length(noisevalues)*length(alphavalues)*length(betavalues))
noise <- rep(noisevalues, each=length(alphavalues)*length(betavalues))
#alpha <- rep(alphavalues, each=length(alphavalues)*length(betavalues))
alpha <- rep(alphavalues, each=length(betavalues))
beta <- rep(betavalues, times = length(noisevalues)*length(alphavalues))


conditions <- as.data.frame(cbind(scans,noise,alpha,beta))
colnames(conditions) <- c("scans","noise", "alpha","beta")
conditions


write.table(conditions,file=paste(PATH,"conditions_results2.txt",sep=""),row.names=FALSE,col.names=TRUE)
