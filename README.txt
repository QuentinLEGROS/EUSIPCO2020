This project provides the code corresponding to the reconstruction method 
presented in the paper 

Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen, "Stochastic EM
      algorithm for fast analysis of single waveform multi-spectral Lidar 
      data," 2020 28th European Signal Processing Conference (EUSIPCO), 
      Amsterdam, 2021, pp. 2413-2417, doi: 10.23919/Eusipco47968.2020.9287414


The main function 'EM_ranging.m' is the mail algorithm performing the estimation
of the reflectivity depth profile.

A demo is available in 'main_synthetic_multi.m'.
This demo uses pre generated data. However, two codes are also provided to generate
synthetic data:

   - 'gener_data4git.m' is the code used to generate the ground truth depth and 
	reflectivity profiles

   - 'initdata.m' generates histograms of photon counts from the ground truth 
	generated in 'gener_data4git.m'.