# PHYS 499
The code used in my PHYS 499. It's all in Python.
The different names for the files will be explained below in what they do.


bayes_linear.py - Linear MCMC code written by Erik Rosolowsky to be used for Bayesian analysis of linear data

Emcee_AgeMetBinLxM.py - MCMC analysis of all three dependent parameters (age, metallicity, and binary fraction) versus x-ray emissivitiy
Emcee_BinMetLxM.py - MCMC analysis of only binary fraction and metallicity (no age) versus x-ray emissivity

ODR_AgeLxM.py - Plot of age versus x-ray emissivity. Also included are the plots of binary fraction and metallicity with age being kept constant versus x-ray emissivity. Orthogonal Distance Regression (ODR) was used to create the linear fit.
ODR_BinLxM.py - Plot of binary fraction versus x-ray emissivity. Also included are the plots of age and metallicity with binary fraction being kept constant versus x-ray emissivity. Orthogonal Distance Regression (ODR) was used to create the linear fit.
ODR_MetLxM.py - Plot of metallicity versus x-ray emissivity. Orthogonal Distance Regression (ODR) was used to create the linear fit. 

ODR_BinRo.py - Plot of binary fraction versus density of low-density stellar environments. Orthogonal Distance Regression (ODR) was used to create the linear fit.
