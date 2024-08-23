1. Ensure that the SU_Waveform_Output_Extraction.mat, stochastic responsiveness .mat, FilterCluster_SUMU.mat, TimeStamps.mat is in the kilosort folder.
2. Run PE_PI.m with desired animal ID and threshold values specified. (for HPC, you need to change the SU_Waveform_Output_Extraction file name to _HPC_SU_Waveform_Output_Extraction which the script loads and output file names to include _HPC as well).
3a. For first site, run stochasticFiring.m with the correct subject ID and change:
    cortexThresh = XXXX; 
    hippThreshLower = XXXX; 
    hippThreshUpper = XXXX; 
3b. For second site (HPC), run stochasticFiringHPC.m with the correct subject ID and change:
    hippThresh_lower = XXXX;
    hippThresh_upper = XXXX;
4. Put the correct stimResponse_gui script (first, second, third cohort) into the kilosort folder and run it by typing the file name of the script into the command window.