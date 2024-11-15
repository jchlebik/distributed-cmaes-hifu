# PDD CMAES arcipelego

## Introduction

This repository contains the implementation of two distributed computation models for HIFU treatment optimizations using CMAES. The models in question are Farmer-workers model and Island evolution model.


## Repository structure

    .
    +--GECCO                            - Paper created from results of this project for GECCO2021
    +--Intel Tools Analysis Screenshots - Screenshot of performance metrics from varios intel toolchain applications
    +--Sources                          - Root folder for the sources.
    Readme.md     - Read me file


## Build instruction
Build on Salomon supercomputer using the provided Makefile. WIP: parametrization of optimization arguments. Currently, to change parameters used during optimization (i.e., number of immigrants, max runtime, ...) one must edit the constant variables near to top of the main source file.

Build using module intel/2019a.
Performance metrics collected using modules VTune/2019_update4 and itac/9.1.2.024

make hifu_pp
make hifu_Bcas

## Usage instruction
The 8 hour benchmark was run using the following series of commands. JobScript using these commands procided in Sources/Jobscripts folder. Jobscripts need to be updated for your use.

qsub -A projectID -q qprod -l select=N,walltime=08:04:00 -I  
ml intel/2019a  
mpirun -np \[N\] ./hifuPp \[seed\] \[OutputDir\]  

OR  

mpirun -np \[N\] ./hifuBcast \[seed\] \[OutputDir\]  


## Author information

 * Name: Jakub Chlebík
 * Email: jakub.chlebik@gmail.com, ichlebik@fit.vutbr.cz
 * Data: 2020/2021
