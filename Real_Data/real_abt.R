##################################################################################
##### Alternative Hypothesis Testing with uncorrected and FDR-corrected NHST #####
##################################################################################

## ADAPT BLOCKS, ADAPT DF ##

## PRELIMINARY

#Libraries
library(oro.nifti)
library(lattice)
library(qvalue)
library(locfdr)


#Working directories
data <- "PATH1" #Where your copes and varcopes for each run are located as well as your ground truths of effect sizes
results <- "PATH2" #Where you want to store your results

conditions <- read.table("real_conditions_results.txt", header = TRUE)

#Parameter values
mudelta1 <- 25  #10000 is the grand mean intensity in FSL for a voxel; this has to be changed manually to 25, 50 or 75
tau <- 21
df <- 165-2



for(r in 1:nrow(conditions)) { #loop over conditions consisting of a configuration of parameters

   	print(r)

   	alpha <- conditions[r,1]
   	beta <- conditions[r,2]


	for(i in 1:100) { #loop over runs

		setwd(data)

		#Read in cope and varcope of a run + uncorrected NHST under null distribution
		print(i)
		cope <- readNIfTI(paste("cope_",i,".nii",sep=""),reorient=FALSE)
		varcope <- readNIfTI(paste("varcope_",i,".nii",sep=""),reorient=FALSE)
		se <- sqrt(varcope)
		tmap <- cope/se
   		pmap.null <- 1 - pt(tmap,df)

	    #FDR NHST	
		if(alpha == 0.05) 
			qval <- alpha
		else 
			qval <- 0.01

	    lfdr <- locfdr(tmap[!is.na(tmap)])
		tmap2 <- (tmap-lfdr$fp0[3,1])/lfdr$fp0[3,2]
		pmap2 <- 1-pt(tmap2,df=74)
		pvalms <- sort(pmap2)
		orderpvalms <- rank(pmap2)
		FDRqval <- (orderpvalms/length(pmap2))*qval
		pr <- ifelse(pvalms[orderpvalms] < FDRqval,1,0)
		pthres.FDR <-ifelse(sum(pr,na.rm=TRUE)==0,0,max(FDRqval[pr==1],na.rm=TRUE))

		#Alternative hypothesis testing + uncorrected NHST
 		pmap.alt <- pnorm(tmap, (mudelta1/se), sqrt(((varcope + tau^2)/varcope)))

		#Creating LSPM
		active <- ifelse(pmap.null <= alpha & pmap.alt >= beta, 1, 0)
		active[is.na(active)] <- 0

		inactive <- ifelse(pmap.null > alpha & pmap.alt < beta, 1, 0)
		inactive[is.na(inactive)] <- 0

		uncertain <- ifelse(pmap.null > alpha & pmap.alt >= beta, 1, 0)
		uncertain[is.na(uncertain)] <- 0

		practins <- ifelse(pmap.null <= alpha & pmap.alt < beta, 1, 0)
		practins[is.na(practins)] <- 0

		es.active <- ifelse(active==1, cope, NA)				
		es.inactive <- ifelse(inactive==1, cope, NA)
		es.uncertain <- ifelse(uncertain==1, cope, NA)
		es.practins <- ifelse(practins==1, cope, NA)

		
		setwd(results)
		writeNIfTI(active, paste("ACTIVE_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""),gzipped=TRUE)
		writeNIfTI(inactive, paste("INACTIVE_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""),gzipped=TRUE)
		writeNIfTI(uncertain, paste("UNCERTAIN_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""),gzipped=TRUE)
		writeNIfTI(practins, paste("PRACTINS_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""),gzipped=TRUE)

		writeNIfTI(es.active, paste("ES_ACTIVE_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""),gzipped=TRUE)
		writeNIfTI(es.inactive, paste("ES_INACTIVE_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""),gzipped=TRUE)
		writeNIfTI(es.uncertain, paste("ES_UNCERTAIN_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""),gzipped=TRUE)
		writeNIfTI(es.practins, paste("ES_PRACTINS_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""),gzipped=TRUE)

		
		##Alternative hypothesis testing + FDR-corrected NHST under null distribution

		#Creating LSPM
		active <- ifelse(pmap.null <= pthres.FDR & pmap.alt >= beta, 1, 0)
		active[is.na(active)] <- 0

		inactive <- ifelse(pmap.null > pthres.FDR & pmap.alt < beta, 1, 0)
		inactive[is.na(inactive)] <- 0

		uncertain <- ifelse(pmap.null > pthres.FDR & pmap.alt >= beta, 1, 0)
		uncertain[is.na(uncertain)] <- 0

		practins <- ifelse(pmap.null <= pthres.FDR & pmap.alt < beta, 1, 0)
		practins[is.na(practins)] <- 0

		es.active <- ifelse(active==1, cope, NA)				
		es.inactive <- ifelse(inactive==1, cope, NA)
		es.uncertain <- ifelse(uncertain==1, cope, NA)
		es.practins <- ifelse(practins==1, cope, NA)


		setwd(results)

		writeNIfTI(active, paste("ACTIVE_FDR_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""),gzipped=TRUE)
		writeNIfTI(inactive, paste("INACTIVE_FDR_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""),gzipped=TRUE)
		writeNIfTI(uncertain, paste("UNCERTAIN_FDR_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""),gzipped=TRUE)
		writeNIfTI(practins, paste("PRACTINS_FDR_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""),gzipped=TRUE)

		writeNIfTI(es.active, paste("ES_ACTIVE_FDR_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""),gzipped=TRUE)
		writeNIfTI(es.inactive, paste("ES_INACTIVE_FDR_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""),gzipped=TRUE)
		writeNIfTI(es.uncertain, paste("ES_UNCERTAIN_FDR_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""),gzipped=TRUE)
		writeNIfTI(es.practins, paste("ES_PRACTINS_FDR_subject_1_run_",i,"_alpha_",alpha,"_beta_",beta,"_",mudelta1,sep=""),gzipped=TRUE)

		
	}
}
