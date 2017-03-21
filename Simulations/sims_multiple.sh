COMMENT: PATH is where your simulation scripts are stored

#!/bin/sh

#PBS -o outputSamen/outputMultiple.file

#PBS -e errorSamen/errorMultiple.file

#PBS -l walltime=71:00:00

#PBS -l nodes=1:ppn=1

#PBS -l vmem=10GB

module load R/3.0.2-ictce-5.5.0

NOISEVALUES=(3 5.5 7)
SCANVALUES=(50 100 150)

for NSIM in $(seq 1 100);

	do

	for COND1 in $(seq 1 3);
       
   		do

		HLP1=$(($COND1 - 1))
   		NOISE=${NOISEVALUES[$HLP1]}
   			
   			for COND2 in $(seq 1 3);
   		       
   		   		do

   				HLP2=$(($COND2 - 1))
   		   		SCANS=${SCANVALUES[$HLP2]}
  
   				Rscript ~PATH/Sim_Multiple_Regions.R  $PBS_ARRAYID $NSIM $NOISE $SCANS

        	done
	done

done	

echo "job finished"
