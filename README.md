exHet_Supplement
================

Example script and simulation data to use the Truong and McCormick et al. (2014) variable heterozygosity model.

File descriptions:

   genotype_probabilities_with_excess_heterozygosity.m - Matlab M file to derive the general solution for the probability of a marker in a genotype class given a selfing population, a finite generation interval, and heterozygosity maintained per generation other than 0.5.

   exampleForExHetEstimation.R - R script to obtain and compile the modified R/qtl code base that contains the est.rf.exHet() function and examples to run it on the provided simulated dataset.

   simulated_genotypes_manuscript.csv - Error-free simulated genotypes with the recombination fractions and map distances (Haldane) that were simulated to generate the genotypes.

   simulated_genotypes_manuscript_rqtl.csv - Error-free simulated genotypes in R/qtl format.

   simulated_genotypes_maniscrupt_errors_and_missing_rqtl.csv - Simulated genotypes with 1% error rate and 5% missingness rate in R/qtl format.

   simulated_genotypes_manuscript_errors_and_missing_rqtl_mapEstimated.csv - Simulated genotypes with 1% error rate and 5% missingness rate in R/qtl format output after an initial map estimation prior to removal of tight double recombinations (aka short double crossovers, or SDCOs).

   simulated_genotypes_manuscript_errors_and_missing_rqtl_mapEstimated.csv-SDCOsToMissingV02 - Simulated genotypes with 1% error rate and 5% missingness rate in R/qtl format after SDCOs less than 2cM have been removed.

   simulateHetGenosV02.py - Python script to generate a simulated dataset. Not a deterministic simulation.

   RqtlSDCOsToMissing.py - Python script to take R/qtl csvr formatted data and set short double crossovers to missing.

   simulated_genotypes_manuscript_rqtl_DEV.csv - A small subset of simulated_genotypes_manuscript.rqtl.csv for testing.
