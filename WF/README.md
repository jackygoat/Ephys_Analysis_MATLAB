REQUIRES COMPUTECANADA!

1. Run FilterCluster_SUMU and put the resulting file in kilosort_output
2. Copy load_open_ephys_binary.m, readNPY.m, readNPYheader.m, loadKSdir.m, loadParamsPy.m, readClusterGroupsCSV.m into kilosort_output folder
3. Upload animal folder with kilosort_output to SCRATCH on narval
4. Modify getWaveformExtraction.m, change paths and animalID for each specific animal
5. Create a getWF.sl file using nano and copy the content into nano
6. sbatch the getWF.sl file and queue the job
7. Resulting WF output file will be in SCRATCH folder