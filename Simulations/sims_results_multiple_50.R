######################################################
## REVIEW SIMULATIONS: RESULTS WITH 50 SCANS ##
######################################################

simulations <- "PATH1" #where simulated brain images are located
results <- "PATH2" #where you want the results to be stored
masks <- "PAHT3" #where a general grey matter mask is stored, sims_rmask.txt
setwd(simulations)


## LIBRARIES
library(neuRosim)
library(AnalyzeFMRI)
#library(methods)
library(gplots)
library(lattice)
library(fmri)
library(matrixStats)
library(oro.nifti)


## SOURCE FILES
source("sims_stimfunction2.R")
source("sims_functions.GLM.R")
conditions_results <- read.table("sims_conditions_results.txt", header = TRUE)



##MASKS
setwd(masks)
noisemask <- fmri::read.NIFTI("rmask.nii")
noisemask <- fmri::extract.data(noisemask, what = "data")
noisemask <- noisemask[,,,1]
nobrain <- which(noisemask == 0)
nvoxels <- sum(noisemask==1)
noisemask[noisemask==0] <- NA

## PRELIMINARY

base <- 100
#qval <- 0.05
mudelta1 <- 0.015
tau <- 0.005
nsim <- 500

fwhm <- 8
sigma <- fwhm/sqrt(8*log(2))
width <- 5


#Image characteristics
dim1 <- 64
dim2 <- 64
dim3 <- 40
imdim <- c(dim1,dim2,dim3) #image dimensions
voxdim <- c(3.5,3.5,3.51) #voxelsize
ext <- 4
nscan <- 50
TR <- 2
nregio <- 2
coord1 <- c(22,17,11)
coord2 <- c(40,18,11)
coord <- list(coord1,coord2)



#CREATING GROUND TRUTH

regions1 <- simprepSpatial(regions=1, coord=coord1, radius=ext, form="sphere", fading=0)
truthdesign <- simprepTemporal(1,1,onsets=1,effectsize = 1, durations=1,TR=1,acc=0.1)
truth1 <- simVOLfmri(design=truthdesign,image=regions1,dim=imdim,SNR=1,noise="none", template=noisemask)[,,,1]
truth1 <- ifelse(truth1 > 0, 1, 0)

regions2 <- simprepSpatial(regions=1, coord=coord2, radius=ext, form="sphere", fading=0)
truthdesign <- simprepTemporal(1,1,onsets=1,effectsize = 1, durations=1,TR=1,acc=0.1)
truth2 <- simVOLfmri(design=truthdesign,image=regions2,dim=imdim,SNR=1,noise="none", template=noisemask)[,,,1]
truth2 <- ifelse(truth2 > 0, 1, 0)


#Result arrats
results.counts.act <- array(NA, c(nrow(conditions_results), 9))
results.counts.inact <- array(NA, c(nrow(conditions_results), 9))
results.counts.unc <- array(NA, c(nrow(conditions_results), 9))
results.counts.practins <- array(NA, c(nrow(conditions_results),9))

colnames(results.counts.act) <- c("noise","alpha","beta","ES=1","ES=2","ES=0","SD=1", "SD=2", "SD=3")
colnames(results.counts.inact) <- c("noise","alpha","beta","ES=1","ES=2","ES=0","SD=1", "SD=2", "SD=3")
colnames(results.counts.unc) <- c("noise","alpha","beta","ES=1","ES=2","ES=0","SD=1", "SD=2", "SD=3")
colnames(results.counts.practins) <- c("noise","alpha","beta","ES=1","ES=2","ES=0","SD=1", "SD=2", "SD=3")
 

#Temporary result arrays within one condition
tmp.counts.act <- array(NA, c(nsim, 6))
tmp.counts.inact <- array(NA, c(nsim, 6))
tmp.counts.unc <- array(NA, c(nsim, 6))
tmp.counts.practins <- array(NA, c(nsim,6))

colnames(tmp.counts.act) <- c("noise","alpha","beta","ES=1","ES=2","ES=0")
colnames(tmp.counts.inact) <- c("noise","alpha","beta","ES=1","ES=2","ES=0")
colnames(tmp.counts.unc) <- c("noise","alpha","beta","ES=1","ES=2","ES=0")
colnames(tmp.counts.practins) <- c("noise","alpha","beta","ES=1","ES=2","ES=0")


#arrays for each simulation
emap.act <- array(NA, c(dim1, dim2, dim3))
emap.inact <- array(NA, c(dim1, dim2, dim3))
emap.unc <- array(NA, c(dim1, dim2, dim3))
emap.practins <- array(NA, c(dim1, dim2, dim3))



## COMPUTING RESULTS
setwd(simulations)
for(i in 1:nrow(conditions_results)) {

	print(i)
	#i<-2
	sigmnoise <- conditions_results[i,1]
	alpha <- conditions_results[i,2]
	beta <- conditions_results[i,3]

	for(s in 1:nsim) {

		#Read in simulated images
		print(s)
		#s<- 1
		setwd(simulations)
		tmap <- f.read.nifti.volume(paste("TMAP_Multiple_null_simulation_",s,"_noise_",sigmnoise,"_scans_50.nii", sep=""))[,,,1]
		tmap <- tmap*noisemask
		sb1 <- f.read.nifti.volume(paste("SB1_Multiple_simulation_",s,"_noise_",sigmnoise,"_scans_50.nii", sep=""))[,,,1]
		sb1 <- sb1*noisemask

		#Null
		pmap.null <- pt(tmap, nscan-2,lower.tail=FALSE)
		active.null <- ifelse(pmap.null<=alpha, 1, 0)

		#Alternative
		pmap.alt <- pnorm(tmap, (mudelta1/sb1), sqrt(((sb1^2 + tau^2)/sb1^2)))
		active.alt <- ifelse(pmap.alt < beta, 0, 1)
		

		## LAYERED STATISTICAL PARAMETRIC MAP: which information in which layer?

		#Make layers
		active <- ifelse(pmap.null <= alpha & pmap.alt >= beta, 1, NA)
		inactive <- ifelse(pmap.null > alpha & pmap.alt < beta, 1, NA)
		uncertain <- ifelse(pmap.null > alpha & pmap.alt >= beta, 1, NA)
		practins <- ifelse(pmap.null <= alpha & pmap.alt < beta, 1, NA)

		#Number of voxels with ES=1, ES=2 en ES=0 in layer
		emap.act <- ifelse(active==1 & truth1==1, 1, NA) #ES=1, practically irrelevant
		emap.act[(active==1 & truth2==1)] <- 2 #ES=2, practically relevant
		emap.act[(active==1 & truth1==0 & truth2==0)] <- 0 #ES=0, inactive		

		emap.inact <- ifelse(inactive==1 & truth1==1, 1, NA) #ES=1, practically irrelevant
		emap.inact[(inactive==1 & truth2==1)] <- 2 #ES=2, practically relevant
		emap.inact[(inactive==1 & truth1==0 & truth2==0)] <- 0 #ES=0, inactive	

		emap.unc <- ifelse(uncertain==1 & truth1==1, 1, NA) #ES=1, practically irrelevant
		emap.unc[(uncertain==1 & truth2==1)] <- 2 #ES=2, practically relevant
		emap.unc[(uncertain==1 & truth1==0 & truth2==0)] <- 0 #ES=0, inactive	

		emap.practins <- ifelse(practins==1 & truth1==1, 1, NA) #ES=1, practically irrelevant
		emap.practins[(practins==1 & truth2==1)] <- 2 #ES=2, practically relevant
		emap.practins[(practins==1 & truth1==0 & truth2==0)] <- 0 #ES=0, inactive	

		

		## Write away data

		tmp.counts.act[s,4] <- sum(emap.act[!is.na(emap.act)]==1)
		tmp.counts.act[s,5] <- sum(emap.act[!is.na(emap.act)]==2)
		tmp.counts.act[s,6] <- sum(emap.act[!is.na(emap.act)]==0)
		tmp.counts.act[s,1:3] <- c(sigmnoise,alpha,beta)

		tmp.counts.inact[s,4] <- sum(emap.inact[!is.na(emap.inact)]==1)
		tmp.counts.inact[s,5] <- sum(emap.inact[!is.na(emap.inact)]==2)
		tmp.counts.inact[s,6] <- sum(emap.inact[!is.na(emap.inact)]==0)
		tmp.counts.inact[s,1:3] <- c(sigmnoise,alpha,beta)

		tmp.counts.unc[s,4] <- sum(emap.unc[!is.na(emap.unc)]==1)
		tmp.counts.unc[s,5] <- sum(emap.unc[!is.na(emap.unc)]==2)
		tmp.counts.unc[s,6] <- sum(emap.unc[!is.na(emap.unc)]==0)
		tmp.counts.unc[s,1:3] <- c(sigmnoise,alpha,beta)

		tmp.counts.practins[s,4] <- sum(emap.practins[!is.na(emap.practins)]==1)
		tmp.counts.practins[s,5] <- sum(emap.practins[!is.na(emap.practins)]==2)
		tmp.counts.practins[s,6] <- sum(emap.practins[!is.na(emap.practins)]==0)
		tmp.counts.practins[s,1:3] <- c(sigmnoise,alpha,beta)

	
	}

	#Write away temporary arrays
	setwd(results)
	write.csv(tmp.counts.act, paste("ALL_COUNTS_ACTIVE_50_noise_",sigmnoise,"_alpha_",alpha,"_beta_",beta,sep=""))
	write.csv(tmp.counts.inact, paste("ALL_COUNTS_INACTIVE_50_noise_",sigmnoise,"_alpha_",alpha,"_beta_",beta,sep=""))
	write.csv(tmp.counts.unc, paste("ALL_COUNTS_UNCERTAIN_50_noise_",sigmnoise,"_alpha_",alpha,"_beta_",beta,sep=""))
	write.csv(tmp.counts.practins, paste("ALL_COUNTS_PRACTINS_50_noise_",sigmnoise,"_alpha_",alpha,"_beta_",beta,sep=""))

	#Summary measures
	results.counts.act[i,4:6] <- round(colMeans(tmp.counts.act[,4:6]))
	results.counts.act[i,7:9] <- colSds(tmp.counts.act[,4:6])
	results.counts.act[i,1:3] <- tmp.counts.act[1,1:3]

	results.counts.inact[i,4:6] <- round(colMeans(tmp.counts.inact[,4:6]))
	results.counts.inact[i,7:9] <- colSds(tmp.counts.inact[,4:6])
	results.counts.inact[i,1:3] <- tmp.counts.inact[1,1:3]

	results.counts.unc[i,4:6] <- round(colMeans(tmp.counts.unc[,4:6]))
	results.counts.unc[i,7:9] <- colSds(tmp.counts.unc[,4:6])
	results.counts.unc[i,1:3] <- tmp.counts.unc[1,1:3]

	results.counts.practins[i,4:6] <- round(colMeans(tmp.counts.practins[,4:6]))
	results.counts.practins[i,7:9] <- colSds(tmp.counts.practins[,4:6])
	results.counts.practins[i,1:3] <- tmp.counts.practins[1,1:3]
}



## WRITE THE RESULTING TABLES AWAY
setwd(results)
write.csv(results.counts.act, "COUNTS_ACTIVE_50")
write.csv(results.counts.inact, "COUNTS_INACTIVE_50")
write.csv(results.counts.unc, "COUNTS_UNCERTAIN_50")
write.csv(results.counts.practins, "COUNTS_PRACTINS_50")

