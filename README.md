# 🌟 Highlights
- Implementation of [PECS](https://ieeexplore.ieee.org/document/9928573/), [FISL](https://www.sciencedirect.com/science/article/abs/pii/S1051200421000300), and CDMM. 
- Runtime and performance comparison of CDMM with PECS and FISL. 
- Sensitivity analysis of polynomial phase model

# ℹ️ Overview
Welcome to the simulation-base for the paper: "Accelerated ISL Reduction of Polynomial Phase Sequences: A NUFFT-Driven Majorize-Minimize Coordinate Approach". 
All simulative figures in it can be reproduced using this repo.

# ✍️ Author
Hi, I'm Marcel Follmann, an early career researcher at the SPARC group in University Luxembourg. My research interests include cognitive radar, waveform optimization, control theory and signal processing as a whole.

# 🚀 Usage
1) At the start the startup.m file should be run once such that all necessary folders for figure and data storage are created. 

2) The repo is structured into different folders. Specifically the `+Algorithm` folder serves as a package to all algorithms compared in the paper. 

3) The folders `04_Experiments` and `05_Plotting` can be used to re-run the experiments in the paper. Simply run an experiment in the experiment folder first and then start the respective plot script in the plot folder to visualize the generated data. 