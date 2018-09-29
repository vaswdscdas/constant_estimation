# Rigorous interpolation error constant estimation (alpha version)

Xuefeng LIU

First version: 2018/09/26

## Running environment

- The codes run on latest MATLAB environment. 
- Since Octave does not supply **ldl** function for indefinite symmetric matrix, the rigorous evaluation is not available for Octave.

## Configuration

Edit **my_intlab_mode_config.m** to configure the computing environment.

```matlab
%The path of the codes for switch between verified computing and approximate computing.
addpath('./mode_swith_interface/');  

%The path of the library of verified eigenvalue estimation for matrix.
addpath('./verified_eig_estimation/'); 

%The path of INTLAB toolbox and initialization.
addpath('~/app/Intlab_V9/');
startintlab;

global INTERVAL_MODE;

%INTERVAL_MODE=1; for rigorous computing based on interval arithmetic.
%INTERVAL_MODE=0; for approximate computing with rounding error inside.
INTERVAL_MODE=1;
```

Remeber to run "my_intlab_mode_config" before other codes.

## How to use it?

For the purpose of rigorous computation, please run the code in MATLAB along with the INTLAB toolbox.

The code has two running mode: approximate evaluation and verified evaluation.

To swith between each other, please set the value of  **INTERVAL_MODE**.

- INTERVAL_MODE=0: approximate computation mode.
- INTERVAL_MODE=1: verified computation mode. It takes longer time that approximate computation mode. So be careful with the mesh size, which should not be too small. (To start around h=1/16 is recommended.)

A sample usage can be found in **main.m**.
