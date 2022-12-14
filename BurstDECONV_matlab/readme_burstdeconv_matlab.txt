
The steps to perform the analysis of a short movie are

1) run Read_data.m  
	-reads the excel data files (resulting from short movies)
	-all data files should be in the folder data, the code reads them all and treats them as a single phenotype
	-plots of data traces go to output/images 
2) run parameters_droso.m
	-loads data specific hyperparameters
	-change values in this script if needed
3)  run Deconv.m  
	-depending on your Matlab version you may want to use sumSignal1_par.m (older Matlab versions), or sum Signalv_par.m (vectorized version using repelem)
	-output goes to output/matfiles/result_*.mat
	-before next step make sure that the output/matfiles folder contains only the result files corresponding to the same phenotype and not other 
	result files from a previous run
4) run AllFit 
	-loads all result files from output/matfiles
	-output in output/FitResults 
