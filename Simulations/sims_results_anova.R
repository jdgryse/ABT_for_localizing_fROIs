######################################################
## REVIEW SIMULATIONS: RESULTS WITH 50 SCANS ##
######################################################

results <- "PATH1" #where the storage arrays with the results of your simulations are stored per simulation condition
setwd(results)


## LIBRARIES
library(neuRosim)
library(AnalyzeFMRI)
#library(methods)
library(gplots)
library(lattice)
library(fmri)
library(matrixStats)
library(oro.nifti)
library(lsr)


## SOURCE FILES

source("sims_stimfunction2.R")
source("sims_functions.GLM.R")
conditions_results <- read.table("sims_conditions_results1.txt", header = TRUE)


##Overall result arrays

active <- array(NA, dim=c(0,dim(conditions_results)[2]+3))
colnames(active) <- c("scans", "noise", "alpha", "beta", "ES1", "ES2", "ES0")

inactive <- array(NA, dim=c(0,dim(conditions_results)[2]+3))
colnames(inactive) <- c("scans", "noise", "alpha", "beta", "ES1", "ES2", "ES0")

uncertain <- array(NA, dim=c(0,dim(conditions_results)[2]+3))
colnames(uncertain) <- c("scans", "noise", "alpha", "beta", "ES1", "ES2", "ES0")

practins <- array(NA, dim=c(0,dim(conditions_results)[2]+3))
colnames(practins) <- c("scans", "noise", "alpha", "beta", "ES1", "ES2", "ES0")



for (i in 1:dim(conditions_results)[1]) {

  print(i)

  scans <- conditions_results[i,1]
  noise <- conditions_results[i,2]
  alpha <- conditions_results[i,3]
  beta <- conditions_results[i,4]

  tmp_act <- read.csv(paste("ALL_COUNTS_ACTIVE_",scans,"_noise_",noise,"_alpha_",alpha,"_beta_",beta,sep=""))
  tmp_act[,1] <- rep(scans, times=dim(tmp_act)[1])
  active <- rbind(active, tmp_act)

  tmp_inact <- read.csv(paste("ALL_COUNTS_INACTIVE_",scans,"_noise_",noise,"_alpha_",alpha,"_beta_",beta,sep=""))
  tmp_inact[,1] <- rep(scans, times=dim(tmp_inact)[1])
  inactive <- rbind(inactive, tmp_inact)

  tmp_uncertain <- read.csv(paste("ALL_COUNTS_UNCERTAIN_",scans,"_noise_",noise,"_alpha_",alpha,"_beta_",beta,sep=""))
  tmp_uncertain[,1] <- rep(scans, times=dim(tmp_uncertain)[1])
  uncertain <- rbind(uncertain, tmp_uncertain)

  tmp_practins <- read.csv(paste("ALL_COUNTS_PRACTINS_",scans,"_noise_",noise,"_alpha_",alpha,"_beta_",beta,sep=""))
  tmp_practins[,1] <- rep(scans, times=dim(tmp_practins)[1])
  practins <- rbind(practins, tmp_practins)

}


for (i in 1:4) {

	active[,i] <- factor(active[,i])
  inactive[,i] <- factor(inactive[,i])
  uncertain[,i] <- factor(uncertain[,i])
  practins[,i] <- factor(practins[,i])

}


active_anova1 <- aov(active[,"ES.1"]~active[,"X"]+active[,"noise"]+active[,"alpha"]+active[,"beta"]+(active[,"X"]*active[,"noise"])+(active[,"X"]*active[,"alpha"])+(active[,"X"]*active[,"beta"])+
	(active[,"noise"]*active[,"alpha"])+(active[,"noise"]*active[,"beta"])+(active[,"alpha"]*active[,"beta"]))
active_eta1 <- etaSquared(x=active_anova1, type=2, anova=TRUE)
active_anova2 <- aov(active[,"ES.2"]~active[,"X"]+active[,"noise"]+active[,"alpha"]+active[,"beta"]+(active[,"X"]*active[,"noise"])+(active[,"X"]*active[,"alpha"])+(active[,"X"]*active[,"beta"])+
                       (active[,"noise"]*active[,"alpha"])+(active[,"noise"]*active[,"beta"])+(active[,"alpha"]*active[,"beta"]))
active_eta2 <- etaSquared(x=active_anova2, type=2, anova=TRUE)
active_anova3 <- aov(active[,"ES.0"]~active[,"X"]+active[,"noise"]+active[,"alpha"]+active[,"beta"]+(active[,"X"]*active[,"noise"])+(active[,"X"]*active[,"alpha"])+(active[,"X"]*active[,"beta"])+
                       (active[,"noise"]*active[,"alpha"])+(active[,"noise"]*active[,"beta"])+(active[,"alpha"]*active[,"beta"]))
active_eta3 <- etaSquared(x=active_anova3, type=2, anova=TRUE)

inactive_anova1 <- aov(inactive[,"ES.1"]~inactive[,"X"]+inactive[,"noise"]+inactive[,"alpha"]+inactive[,"beta"]+(inactive[,"X"]*inactive[,"noise"])+(inactive[,"X"]*inactive[,"alpha"])+(inactive[,"X"]*inactive[,"beta"])+
                         (inactive[,"noise"]*inactive[,"alpha"])+(inactive[,"noise"]*inactive[,"beta"])+(inactive[,"alpha"]*inactive[,"beta"]))
inactive_eta1 <- etaSquared(x=inactive_anova1, type=2, anova=TRUE)
inactive_anova2 <- aov(inactive[,"ES.2"]~inactive[,"X"]+inactive[,"noise"]+inactive[,"alpha"]+inactive[,"beta"]+(inactive[,"X"]*inactive[,"noise"])+(inactive[,"X"]*inactive[,"alpha"])+(inactive[,"X"]*inactive[,"beta"])+
                         (inactive[,"noise"]*inactive[,"alpha"])+(inactive[,"noise"]*inactive[,"beta"])+(inactive[,"alpha"]*inactive[,"beta"]))
inactive_eta2 <- etaSquared(x=inactive_anova2, type=2, anova=TRUE)
inactive_anova3 <- aov(inactive[,"ES.0"]~inactive[,"X"]+inactive[,"noise"]+inactive[,"alpha"]+inactive[,"beta"]+(inactive[,"X"]*inactive[,"noise"])+(inactive[,"X"]*inactive[,"alpha"])+(inactive[,"X"]*inactive[,"beta"])+
                         (inactive[,"noise"]*inactive[,"alpha"])+(inactive[,"noise"]*inactive[,"beta"])+(inactive[,"alpha"]*inactive[,"beta"]))
inactive_eta3 <- etaSquared(x=inactive_anova3, type=2, anova=TRUE)

uncertain_anova1 <- aov(uncertain[,"ES.1"]~uncertain[,"X"]+uncertain[,"noise"]+uncertain[,"alpha"]+uncertain[,"beta"]+(uncertain[,"X"]*uncertain[,"noise"])+(uncertain[,"X"]*uncertain[,"alpha"])+(uncertain[,"X"]*uncertain[,"beta"])+
                          (uncertain[,"noise"]*uncertain[,"alpha"])+(uncertain[,"noise"]*uncertain[,"beta"])+(uncertain[,"alpha"]*uncertain[,"beta"]))
uncertain_eta1 <- etaSquared(x=uncertain_anova1, type=2, anova=TRUE)
uncertain_anova2 <- aov(uncertain[,"ES.2"]~uncertain[,"X"]+uncertain[,"noise"]+uncertain[,"alpha"]+uncertain[,"beta"]+(uncertain[,"X"]*uncertain[,"noise"])+(uncertain[,"X"]*uncertain[,"alpha"])+(uncertain[,"X"]*uncertain[,"beta"])+
                          (uncertain[,"noise"]*uncertain[,"alpha"])+(uncertain[,"noise"]*uncertain[,"beta"])+(uncertain[,"alpha"]*uncertain[,"beta"]))
uncertain_eta2 <- etaSquared(x=uncertain_anova2, type=2, anova=TRUE)
uncertain_anova3 <- aov(uncertain[,"ES.0"]~uncertain[,"X"]+uncertain[,"noise"]+uncertain[,"alpha"]+uncertain[,"beta"]+(uncertain[,"X"]*uncertain[,"noise"])+(uncertain[,"X"]*uncertain[,"alpha"])+(uncertain[,"X"]*uncertain[,"beta"])+
                          (uncertain[,"noise"]*uncertain[,"alpha"])+(uncertain[,"noise"]*uncertain[,"beta"])+(uncertain[,"alpha"]*uncertain[,"beta"]))
uncertain_eta3 <- etaSquared(x=uncertain_anova3, type=2, anova=TRUE)

practins_anova1 <- aov(practins[,"ES.1"]~practins[,"X"]+practins[,"noise"]+practins[,"alpha"]+practins[,"beta"]+(practins[,"X"]*practins[,"noise"])+(practins[,"X"]*practins[,"alpha"])+(practins[,"X"]*practins[,"beta"])+
                         (practins[,"noise"]*practins[,"alpha"])+(practins[,"noise"]*practins[,"beta"])+(practins[,"alpha"]*practins[,"beta"]))
practins_eta1 <- etaSquared(x=practins_anova1, type=2, anova=TRUE)
practins_anova2 <- aov(practins[,"ES.2"]~practins[,"X"]+practins[,"noise"]+practins[,"alpha"]+practins[,"beta"]+(practins[,"X"]*practins[,"noise"])+(practins[,"X"]*practins[,"alpha"])+(practins[,"X"]*practins[,"beta"])+
                         (practins[,"noise"]*practins[,"alpha"])+(practins[,"noise"]*practins[,"beta"])+(practins[,"alpha"]*practins[,"beta"]))
practins_eta2 <- etaSquared(x=practins_anova2, type=2, anova=TRUE)
practins_anova3 <- aov(practins[,"ES.0"]~practins[,"X"]+practins[,"noise"]+practins[,"alpha"]+practins[,"beta"]+(practins[,"X"]*practins[,"noise"])+(practins[,"X"]*practins[,"alpha"])+(practins[,"X"]*practins[,"beta"])+
                         (practins[,"noise"]*practins[,"alpha"])+(practins[,"noise"]*practins[,"beta"])+(practins[,"alpha"]*practins[,"beta"]))
practins_eta3 <- etaSquared(x=practins_anova3, type=2, anova=TRUE)

interaction.plot(active[,"alpha"], active[,"noise"], active[,"ES.0"])
interaction.plot(active[,"noise"], active[,"beta"], active[,"ES.1"])
interaction.plot(active[,"alpha"], active[,"noise"], active[,"ES.2"])

interaction.plot(inactive[,"X"], inactive[,"noise"], inactive[,"ES.1"])
interaction.plot(inactive[,"alpha"], inactive[,"noise"], inactive[,"ES.1"])
interaction.plot(inactive[,"beta"], inactive[,"noise"], inactive[,"ES.2"])
interaction.plot(inactive[,"noise"], inactive[,"X"], inactive[,"ES.0"])
interaction.plot(inactive[,"noise"], inactive[,"beta"], inactive[,"ES.0"])

interaction.plot(practins[,"X"], practins[,"noise"], practins[,"ES.1"])
interaction.plot(practins[,"alpha"], practins[,"noise"], practins[,"ES.1"])
interaction.plot(practins[,"beta"], practins[,"noise"], practins[,"ES.1"])
interaction.plot(practins[,"beta"], practins[,"noise"], practins[,"ES.2"])
interaction.plot(practins[,"beta"], practins[,"alpha"], practins[,"ES.2"])
interaction.plot(practins[,"alpha"], practins[,"noise"], practins[,"ES.0"])

interaction.plot(uncertain[,"alpha"], uncertain[,"noise"], uncertain[,"ES.2"])
interaction.plot(uncertain[,"X"], uncertain[,"noise"], uncertain[,"ES.0"])
interaction.plot(uncertain[,"beta"], uncertain[,"noise"], uncertain[,"ES.0"])




