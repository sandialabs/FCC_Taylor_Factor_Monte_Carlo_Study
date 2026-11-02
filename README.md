# FCC Taylor Factor Monte Carlo Study

Compute mean Taylor factors and Lankford ratios (R-values) for FCC materials by sampling
Euler orientations and minimizing the average Taylor factor over a lateral
shrinking ratio *q*.

## Repository structure
├── LICENSE  
├── THIRD_PARTY_NOTICES.md  
├── Monte_Carlo_Taylor_Generator_Parallelized.m    # Driver script  
└── FCC_Taylor_Factor.m              # GPL-2.0-licensed Taylor-factor function  

## Prerequisites
- MATLAB R2014a or newer  
- Parallel Computing Toolbox  
- Optimization Toolbox  
- Statistics and Machine Learning Toolbox for `datasample`  

## Usage
1. Place your Euler‐angle data file (3 columns: phi1, Phi, phi2 in radians) in the same folder, or edit the `orientationFile` variable in `Monte_Carlo_Taylor_Generator_Parallelized.m` to point to your file.
2. Adjust the `numWorkers` variable based on your system.
3. Adjust the `desiredIterations` and `grains_parameter` variables based on the number of individual texture realizations and number of grains in each texture realization.

## Copyright
Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
