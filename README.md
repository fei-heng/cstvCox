# cstvCox
Codes for the paper:
Heng, F., Sun, Y., Hyun, S. et al. Analysis of the time-varying Cox model for the cause-specific hazard functions with missing causes. Lifetime Data Anal 26, 731–760 (2020). https://doi.org/10.1007/s10985-020-09497-y

This folder contains the MATLAB code and a simulated data set:
	* sim_data.mat: a simulated dataset
	* MATLAB functions:
		+ main.m: main function
		+ esta.m, vara.m: compute estimates and estimated standard errors using two-stage AIPW method
		+ esti.m, vari.m: compute estimates and estimated standard errors using IPW method
		+ estf.m, varf.m: compute estimates and estimated standard errors for full data
		+ epanker.m: (1/h)*K((t_k-t)/h), K is Epanechnikov kernel

In main.m, we use csvread() to load the simulated dataset: csvread('simdata.csv').

It includes input variables as shown below. 
	
INPUT:
------
	time: X=min{T,C}, where T is the failure time and C is the censoring time
	delta: I(T<=C), censoring indicator, 1 for the death/failure is observed, 0 for censored
	z1: covariate matrix
	R: missing indicator, 0 for that a death/failure is observed but cause is missing, otherwise 1.
	A: a binary auxialiry covariate
	cause: cause of failure V
	
PARAMETER SETTINGS:
-------------------
	missingmodel: working model for the missing probability
	h: bandwidth
	tau: administration time
	tstep: steps of grid points
	
OUTPUT:
-------
* Estimates of covariate coefficients: sbeta_acc
* Estimated standard errors: sstd_acc
* p-values: 
    + pv_test1_a1: for the alternative H1a using a supremum type test statistic
	+ pv_test1_a2: for the alternative H1a using a integrated type test statistic
	+ pv_test1_m1: for the alternative H1m using a supremum type test statistic
	+ pv_test1_m2: for the alternative H1m using a integrated type test statistic
	+ pv_test2_a1: for the alternative H2a using a supremum type test statistic
	+ pv_test2_a2: for the alternative H2a using a integrated type test statistic
	+ pv_test2_m1: for the alternative H2m using a supremum type test statistic
	+ pv_test2_m2: for the alternative H2m using a integrated type test statistic
	
COMPUTATION TIME:
-----------------
The required computation time for sample size 1200 is about 1 minute running on the High Performance Computing Cluster at UNC Charlotte.****
