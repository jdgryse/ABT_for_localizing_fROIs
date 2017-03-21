# Simulations of alternative-based thresholding for the definition of fROIs

## sims_multiple.R

In this file, the subject maps are simulated.

## sims_multiple.sh

This file is used to submit to the linux-based supercomputer in order to simulate the brain maps.

## sims_results.R

These files are used to analyze the simulated brain maps with both ABT and NHST. Voxels in different layers are counted.

## make_conditions.R

This file makes a table with all parameter configurations needed in sims_results. The resulting table is conditions_results.txt

## rmask.txt

The mask with which the brain maps are simulated

## function.GLM.R

This file contains the function to estimate the parameters of the GLM, their standard error and the tmap