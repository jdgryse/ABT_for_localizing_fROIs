#########################################################################################
##### Evaluate spatial accuracy of thresholded run with ground truth of effect sizes #####
#########################################################################################

## PRELIMINARY

#Libraries
library(oro.nifti)
library(lattice)
library(xtable)

#Working directories
data1 <- "PATH1" #results path from real_abt.R, so where the ABT maps are stored
data2 <- "PATH2" #where ground truths of effect sizes are stored
results <- "PATH3" #path where you want to store results concerning voxel counts etc.

#Parameter values
df <- 165-2
mudelta1 <- 75 #manually adjust to 25, 50 or 75

conditions <- read.table("real_conditions_results.txt", header = TRUE)
mask <- readNIfTI("real_OVERALLMASK.nii")
mask[mask==0] <- NA


#Results matrices

results.act <- array(NA, dim=c(6,8))
colnames(results.act) <- c("alpha","beta","meanES","medES","sdES",">0.5/truth",">0.5/layer","dice")
results.inact <- array(NA, dim=c(6,8))
colnames(results.inact) <- c("alpha","beta","meanES","medES","sdES",">0.5/truth",">0.5/layer","dice")
results.unc <- array(NA, dim=c(6,8))
colnames(results.unc) <- c("alpha","beta","meanES","medES","sdES",">0.5/truth",">0.5/layer","dice")
results.pract <- array(NA, dim=c(6,8))
colnames(results.pract) <- c("alpha","beta","meanES","medES","sdES",">0.5/truth",">0.5/layer","dice")
results.actnull <- array(NA, dim=c(6,8))
colnames(results.actnull) <- c("alpha","beta","meanES","medES","sdES",">0.5/truth",">0.5/layer","dice")
results.inactnull <- array(NA, dim=c(6,8))
colnames(results.inactnull) <- c("alpha","beta","meanES","medES","sdES",">0.5/truth",">0.5/layer","dice")

tmp.act <- array(NA, dim=c(100,6))
colnames(tmp.act) <- c("meanES","medES","sdES",">0.5/truth",">0.5/layer","dice")
tmp.inact <- array(NA, dim=c(100,6))
colnames(tmp.inact) <- c("meanES","medES","sdES",">0.5/truth",">0.5/layer","dice")
tmp.unc <- array(NA, dim=c(100,6))
colnames(tmp.unc) <- c("meanES","medES","sdES",">0.5/truth",">0.5/layer","dice")
tmp.pract <- array(NA, dim=c(100,6))
colnames(tmp.pract) <- c("meanES","medES","sdES",">0.5/truth",">0.5/layer","dice")
tmp.actnull <- array(NA, dim=c(100,6))
colnames(tmp.actnull) <- c("meanES","medES","sdES",">0.5/truth",">0.5/layer","dice")
tmp.inactnull <- array(NA, dim=c(100,6))
colnames(tmp.inactnull) <- c("meanES","medES","sdES",">0.5/truth",">0.5/layer","dice")

for(r in 1:nrow(conditions)) {

	print(r)

	alpha <- conditions[r,1]
	beta <- conditions[r,2]

	for (i in 1:100) {

		print(i)

		#Ground truth of ESs
		truthcope <- readNIfTI(paste("PATH2/truth_ES_",i,".nii.gz",sep=""))
		truthcope <- mask*truthcope

		#Counts of number of voxels > mudelta1 in ground truth
		func.rel.truth <- ifelse(truthcope>=mudelta1, 1, 0)
		truthvox <- sum(func.rel.truth,na.rm=TRUE)
		func.irrel.truth <- ifelse(truthcope<mudelta1,1,0)
		intruthvox <- sum(func.irrel.truth,na.rm=TRUE)

		#Distribution ES in NHST
		cope <- readNIfTI(paste("PATH2/cope_",i,".nii",sep=""),reorient=FALSE)
		runmask <- ifelse(cope!=0, 1, NA)
		varcope <- readNIfTI(paste("PATH2/varcope_",i,".nii",sep=""),reorient=FALSE)
		se <- sqrt(varcope)
		tmap <- cope/se
		tmap <- tmap*runmask
		pmap.null <- 1-pt(tmap,df=df)
		active.null <- ifelse(pmap.null <= alpha, 1, NA)
		inactive.null <- ifelse(pmap.null > alpha, 1, NA)

		gem.actnull <- mean(truthcope[active.null==1],na.rm=TRUE)
		med.actnull <- median(truthcope[active.null==1],na.rm=TRUE)
		std.actnull <- sd(truthcope[active.null==1],na.rm=TRUE)

		gem.inactnull <- mean(truthcope[inactive.null==1],na.rm=TRUE)
		med.inactnull <- median(truthcope[inactive.null==1],na.rm=TRUE)
		std.inactnull <- sd(truthcope[inactive.null==1],na.rm=TRUE)

		#Counts of active and inactive layer
		rel.active.null <- ifelse(truthcope[active.null==1]>=mudelta1,1,NA)
		relvox.act.null <- sum(rel.active.null, na.rm=TRUE)
		relvox.actnulltruth <- relvox.act.null/truthvox
		relvox.actnulllay <- relvox.act.null/sum(active.null, na.rm=TRUE)

		rel.inactive.null <- ifelse(truthcope[inactive.null==1]<mudelta1,1,NA)
		relvox.inact.null <- sum(rel.inactive.null, na.rm=TRUE)
		relvox.inactnulltruth <- relvox.inact.null/intruthvox 
		relvox.inactnulllay <- relvox.inact.null/sum(inactive.null, na.rm=TRUE)


		#Read in different layers
		setwd(data)
		active <- readNIfTI(paste("ACTIVE_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""))
		active <- active * runmask
		inactive <- readNIfTI(paste("INACTIVE_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""))
		inactive <- inactive * runmask
		uncertain <- readNIfTI(paste("UNCERTAIN_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""))
		uncertain <- uncertain * runmask
		practins <- readNIfTI(paste("PRACTINS_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""))
		practins <- practins * runmask


		##COMPARE LAYERS WITH GROUND TRUTH

		#Distributions ES LSPM
		gem.act <- round(mean(truthcope[active==1],na.rm=TRUE)/100, digits=3)
		med.act <- round(median(truthcope[active==1],na.rm=TRUE)/100, digits=3)
		std.act <- round(sd(truthcope[active==1],na.rm=TRUE)/100, digits=3)

		gem.inact <- round(mean(truthcope[inactive==1],na.rm=TRUE)/100, digits=3)
		med.inact <- round(median(truthcope[inactive==1],na.rm=TRUE)/100, digits=3)
		std.inact <- round(sd(truthcope[inactive==1],na.rm=TRUE)/100, digits=3)

		gem.unc <- round(mean(truthcope[uncertain==1],na.rm=TRUE)/100, digits=3)
		med.unc <- round(median(truthcope[uncertain==1],na.rm=TRUE)/100, digits=3)
		std.unc <- round(sd(truthcope[uncertain==1],na.rm=TRUE)/100, digits=3)
		
		gem.pract <- round(mean(truthcope[practins==1],na.rm=TRUE)/100, digits=3)
		med.pract <- round(median(truthcope[practins==1],na.rm=TRUE)/100, digits=3)
		std.pract <- round(sd(truthcope[practins==1],na.rm=TRUE)/100, digits=3)


		#Counts of number of voxels > 50 in ground truth found in layer and how much of the ground truth is found in layer
		rel.active <- ifelse(truthcope[active==1]>=mudelta1,1,NA)
		relvox.act <- sum(rel.active, na.rm=TRUE)
		relvox.acttruth <- relvox.act/truthvox
		relvox.actlay <- relvox.act/sum(active, na.rm=TRUE)

		rel.inactive <- ifelse(truthcope[inactive==1]<mudelta1,1,NA)
		relvox.inact <- sum(rel.inactive, na.rm=TRUE)
		relvox.inacttruth <- relvox.inact/intruthvox
		relvox.inactlay <- relvox.inact/sum(inactive, na.rm=TRUE)

		rel.uncertain <- ifelse(truthcope[uncertain==1]>=mudelta1,1,NA)
		relvox.unc <- sum(rel.uncertain, na.rm=TRUE)
		relvox.unctruth <- relvox.unc/truthvox
		relvox.unclay <- relvox.unc/sum(uncertain, na.rm=TRUE)

		rel.practins <- ifelse(truthcope[practins==1]>=mudelta1,1,NA)
		relvox.pract <- sum(rel.practins, na.rm=TRUE)
		relvox.practtruth <- relvox.pract/truthvox
		relvox.practlay <- relvox.pract/sum(practins, na.rm=TRUE)

		#Results storage in temporary arrays
		tmp.act[i,1] <- gem.act
		tmp.act[i,2] <- med.act
		tmp.act[i,3] <- std.act
		tmp.act[i,4] <- relvox.acttruth
		tmp.act[i,5] <- relvox.actlay
		tmp.act[i,6] <- (2*relvox.act)/(sum(active, na.rm=TRUE)+ truthvox)

		tmp.inact[i,1] <- gem.inact
		tmp.inact[i,2] <- med.inact
		tmp.inact[i,3] <- std.inact
		tmp.inact[i,4] <- relvox.inacttruth
		tmp.inact[i,5] <- relvox.inactlay
		tmp.inact[i,6] <- (2*relvox.inact)/(sum(inactive,na.rm=TRUE) + intruthvox)

		tmp.unc[i,1] <- gem.unc
		tmp.unc[i,2] <- med.unc
		tmp.unc[i,3] <- std.unc
		tmp.unc[i,4] <- relvox.unctruth
		tmp.unc[i,5] <- relvox.unclay
		tmp.unc[i,6] <- (2*relvox.unc)/(sum(uncertain, na.rm=TRUE)+ truthvox)

		tmp.pract[i,1] <- gem.pract
		tmp.pract[i,2] <- med.pract
		tmp.pract[i,3] <- std.pract
		tmp.pract[i,4] <- relvox.practtruth
		tmp.pract[i,5] <- relvox.practlay		
		tmp.pract[i,6] <- (2*relvox.pract)/(sum(practins, na.rm=TRUE)+ truthvox)

		tmp.actnull[i,1] <- gem.actnull
		tmp.actnull[i,2] <- med.actnull
		tmp.actnull[i,3] <- std.actnull
		tmp.actnull[i,4] <- relvox.actnulltruth
		tmp.actnull[i,5] <- relvox.actnulllay
		tmp.actnull[i,6] <- (2*relvox.act.null)/(sum(active.null, na.rm=TRUE)+ truthvox)

		tmp.inactnull[i,1] <- gem.inactnull
		tmp.inactnull[i,2] <- med.inactnull
		tmp.inactnull[i,3] <- std.inactnull
		tmp.inactnull[i,4] <- relvox.inactnulltruth
		tmp.inactnull[i,5] <- relvox.inactnulllay
		tmp.inactnull[i,6] <- (2*relvox.inact.null)/(sum(inactive.null, na.rm=TRUE)+ intruthvox)


	}


		results.act[r,1] <- alpha
		results.act[r,2] <- beta
		results.act[r,3] <- mean(tmp.act[,1],na.rm=T)
		results.act[r,4] <- mean(tmp.act[,2],na.rm=T)
		results.act[r,5] <- mean(tmp.act[,3],na.rm=T)
		results.act[r,6] <- mean(tmp.act[,4],na.rm=T)
		results.act[r,7] <- mean(tmp.act[,5],na.rm=T)
		results.act[r,8] <- mean(tmp.act[,6],na.rm=T)

		results.inact[r,1] <- alpha
		results.inact[r,2] <- beta
		results.inact[r,3] <- mean(tmp.inact[,1],na.rm=T)
		results.inact[r,4] <- mean(tmp.inact[,2],na.rm=T)
		results.inact[r,5] <- mean(tmp.inact[,3],na.rm=T)
		results.inact[r,6] <- mean(tmp.inact[,4],na.rm=T)
		results.inact[r,7] <- mean(tmp.inact[,5],na.rm=T)
		results.inact[r,8] <- mean(tmp.inact[,6],na.rm=T)

		results.unc[r,1] <- alpha
		results.unc[r,2] <- beta
		results.unc[r,3] <- mean(tmp.unc[,1],na.rm=T)
		results.unc[r,4] <- mean(tmp.unc[,2],na.rm=T)
		results.unc[r,5] <- mean(tmp.unc[,3],na.rm=T)
		results.unc[r,6] <- mean(tmp.unc[,4],na.rm=T)
		results.unc[r,7] <- mean(tmp.unc[,5],na.rm=T)
		results.unc[r,8] <- mean(tmp.unc[,6],na.rm=T)

		results.pract[r,1] <- alpha
		results.pract[r,2] <- beta
		results.pract[r,3] <- mean(tmp.pract[,1],na.rm=T)
		results.pract[r,4] <- mean(tmp.pract[,2],na.rm=T)
		results.pract[r,5] <- mean(tmp.pract[,3],na.rm=T)
		results.pract[r,6] <- mean(tmp.pract[,4],na.rm=T)
		results.pract[r,7] <- mean(tmp.pract[,5],na.rm=T)
		results.pract[r,8] <- mean(tmp.pract[,6],na.rm=T)

		results.actnull[r,1] <- alpha
		results.actnull[r,2] <- beta
		results.actnull[r,3] <- mean(tmp.actnull[,1],na.rm=T)
		results.actnull[r,4] <- mean(tmp.actnull[,2],na.rm=T)
		results.actnull[r,5] <- mean(tmp.actnull[,3],na.rm=T)
		results.actnull[r,6] <- mean(tmp.actnull[,4],na.rm=T)
		results.actnull[r,7] <- mean(tmp.actnull[,5],na.rm=T)
		results.actnull[r,8] <- mean(tmp.actnull[,6],na.rm=T)

		results.inactnull[r,1] <- alpha
		results.inactnull[r,2] <- beta
		results.inactnull[r,3] <- mean(tmp.inactnull[,1],na.rm=T)
		results.inactnull[r,4] <- mean(tmp.inactnull[,2],na.rm=T)
		results.inactnull[r,5] <- mean(tmp.inactnull[,3],na.rm=T)
		results.inactnull[r,6] <- mean(tmp.inactnull[,4],na.rm=T)
		results.inactnull[r,7] <- mean(tmp.inactnull[,5],na.rm=T)
		results.inactnull[r,8] <- mean(tmp.inactnull[,6],na.rm=T)

}

setwd(results)
write.csv(results.act, paste("Active_mudelta_",mudelta1,sep=""))
write.csv(results.inact, paste("Inactive_mudelta_",mudelta1,sep=""))
write.csv(results.unc, paste("Uncertain_mudelta_",mudelta1,sep=""))
write.csv(results.pract, paste("Practins_mudelta_",mudelta1,sep=""))
write.csv(results.actnull,paste("Active_null_mudelta_",mudelta1,sep=""))
write.csv(results.inactnull, paste("Inactive_null_mudelta_",mudelta1,sep=""))




