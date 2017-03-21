# Simulations of alternative-based thresholding for the definition of fROIs

## sims_multiple.R

In this file, the subject maps are simulated.

## sims_multiple.sh

This file is used to submit to the linux-based supercomputer in order to simulate the brain maps.

## sims_results.R

These files are used to analyze the simulated brain maps with both ABT and NHST. Voxels in different layers are counted.

## sims_make_conditions.R

This file makes a table with all parameter configurations needed in sims_results. The resulting table is conditions_results.txt

## sims_conditions_results….txt

These files represent the arrays needed for the loops in sims_results_multiple_….R and sims_results_anova.R respectively.

## sims_rmask.txt

The mask with which the brain maps are simulated

## sims_function.GLM.R

This file contains the function to estimate the parameters of the GLM, their standard error and the tmap

## sims_results_multiple_….R

In these files, the information within each layer of the LSPM after ABT is evaluated. this is done for the different number of scans separately.

## sims_results_anova.R

In this file anovas were performed in order to compute effect sizes for the simulation parameters contributing to the simulation design.


# Real data example and evaluation of alternative-based thresholding for the definition of fROIs

## Data

A toy data set in order to perform the cross-validation performed in the original study, as well as an overall brain mask for this “simulated individual”.

## design.

These files contain information concerning the design of the fixed effects analysis in FSL, needed for the flameo command in real_es_groundtruth.sh

## real_abt.R

In this file, NHST and ABT is performed on the left out run under different parameter configutations.

## real_designCrossVal.R

Here, the design files mentioned above are made.

## real_es_accuracy.R

In this file, the left out run is compared with the functionally relevant ground truth.

# real_es_groundtruth.sh

Here, the ground truth of effect sizes for each step of the cross-validation is constructed.

## real_writefiles.R

In this file, the overall brain mask for the “simulated individual” is constructed.
