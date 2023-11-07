# Permutation testing for 15s data

1. Combine all data into one matrix (e.g. for PSD, generate a matrix that is trials x channel x freqband). 
2. Generate a list of labels corresponding to entries in matrix
3. Load this dataset in `run_permutation.m`
    1. this script generates binary labels depending on the conditions of interest, shuffles the data, and generates a matrix of t-statistic values 
    2. fxn `permutations.m` can perform permutation testing for ROIs or channels, and `fxn_handle` is currently configured for `Tfunc` for t-tests but this can be changed for ANOVA, etc. 
4. first entry of t-statistic matrix actual data, rest is shuffled data  
5. p_value is calculated in `hypothesis_testing.m`

Other notes: 

- pipeline was validated using simulated data with known probability distribution
    1. `run_sim_data.m` runs simulated data to get t-statistic values, p-values using steps listed in 1-5.