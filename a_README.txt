This README describes the data and code files used in:
E. G. Ozaktas et al. "Aperiodicity and Disorder as Systematic Spectral Tuning Mechanisms for Plasmonic Nanostructures", 2025.

There are several categories, each corresponding to approximately one figure of the paper. Each is accompanied by one plotting code and a folder containing the relevant data files. There is also a function in a separate file calc_smoothness.m that calculates local variance.

BasicPhysics:
BasicPhysics_Plot.m: This plots the extinction spectra in Fig. 1.
BasicPhysics_a_1200nm_r100_h200.mat: This file contains the frequency and transmittance data for an array of Au plasmonic nanodisks with period 1200 nm, radius 100 nm, height 200 nm.

Bernoulli_2D:
Bernoulli_2D_Plots.m: This plots extinction profiles, calculates local variance, and visualizes geometry for plasmonic nanodisks distributed according to a Bernoulli point process (Fig. 3).
The data in the folder "Bernoulli_2D" are named as follows:
"Bernoulli_2D_divX_Y.mat" where X = [4,5,6,7,8,10,20] is the factor corresponding to how many grid points are in one period length of the unit cell, and Y corresponds to the trial number (Y = 1 does not exist, the first trial filename has _1 omitted).

Bernoulli_L:
Bernoulli_L_Plots.m: This plots extinction profiles, calculates local variance, and visualizes geometry for plasmonic nanodisks distributed according to a Bernoulli point process (Fig. 4).
The data in the folder "Bernoulli_L" are named as follows:
"Bernoulli_L_ZxZ_nX_Y.mat" where Z = [3,4,5,6] is the factor corresponding to how many grid points are in one period length of the unit cell for p = 1. X = [15, 20, 25, 30] corresponds to how many points are in one period length for p = 0.04, and Y corresponds to the trial number (Y = 1 does not exist, the first trial filename has _1 omitted).

FP:
FP_Plots.m: This plots extinction profiles, calculates local variance, and visualizes geometry for plasmonic nanodisks distributed according to a frozen phonon process (Fig. 5).
The data in the folder "FP" are named as follows:
"FP_Xa_Y.mat" where X = [00, 05, 10, 20, 30, 40, 50] is b parametrizing the variance (see manuscript), and Y corresponds to the trial number (Y = 1 does not exist, the first trial filename has _1 omitted).

LGC:
LGC_Plots.m: This plots extinction profiles, calculates local variance, and visualizes geometry for plasmonic nanodisks distributed according to a log Gaussian Cox process (Fig. 8).
The data in the folder "LGC" are named as follows:
"LGC_sX_l02_N30_Y.mat" where X = [05, 1, 15, 2, 25, 3, 4, 5] (two digit numbers should be interpreted with a decimal point in between) is s parametrizing the amplitude of the Gaussian (see manuscript), and Y corresponds to the trial number (Y = 1 does not exist, the first trial filename has _1 omitted).

LR:
LR_Plots.m: This plots extinction profiles, calculates local variance, and visualizes geometry for plasmonic nanodisks distributed according to a long range frozen phonon process (Fig. 6).
The data in the folder "LR" are named as follows:
"LR_Xa_Y.mat" where X = [00, 05, 10, 20, 30, 40, 50] is b parametrizing the variance (see manuscript), and Y corresponds to the trial number (Y = 1 does not exist, the first trial filename has _1 omitted).

QCoct:
QCoct_Plots.m: This plots extinction profiles, calculates local variance, and visualizes geometry for plasmonic nanodisks distributed according to an octogonal quasicrystal (Fig. 2).
The data in the folder "LR" are named as follows:
"QCoct_nX.mat" where X = [0, 1, 2, 3] corresponds to the approximant number of the quasicrystal.

Strauss:
Strauss_Plots.m: This plots extinction profiles, calculates local variance, and visualizes geometry for plasmonic nanodisks distributed according to a Strauss point process (Fig. 7).
The data in the folder "LR" are named as follows:
"Strauss_gX_r01_mu30_Y.mat" where X = [00, 02, 04, 06, 08, 10] is gamma parametrizing the process (see manuscript), and Y corresponds to the trial number (Y = 1 does not exist, the first trial filename has _1 omitted).

