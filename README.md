# NUFFT based ISL minimization of polynomial phase sequences
This repositroy is the simulation-base for the paper: "Accelerated ISL Reduction of Polynomial Phase Sequences: A NUFFT-Driven Majorize-Minimize Coordinate Approach"

## Folder Structure

**01_Figures**: Storage for generated figures.

**02_Data**: Storage for generated experiment data.

**03_Tests**: Internal function tests.

**04_Experiments**: Pure data generation.
+ **cd_improvement_cdmm**: How much does ISL improve in each coordinate?
+ **mi_cdmm**: Simulate mutual interference of polynomial phase models.
+ **perf_cdmm_pecs**: ISL results after exemplary CDMM / PECS optimizations.
+ **q2t_cdmm_pecs**: Quality budget to runtime simulation
    1) CDMM gives quality budget and respective runtime
    2) Analyze runtime of PECS with given quality goal
+ **rt_cdmm_fisl**: Plain runtime simulation of CDMM, CDMM-OOD, and FISL.
+ **rt_cdmm_pecs**: Plain runtime simulation of CDMM and PECS.
+ **sen_cdmm**: Sensitivity analysis of polynomial phase model.
+ **t2q_cdmm_pecs**: Runtime budget to isl performance simulation
    1) CDMM gives runtime budget and respective reached ISL.
    2) Analyze ISL performance of PECS after given time budget.

**05_Plotting**: Figure generation.

**99_Utilities**: Storage for Submodules.

## Relation between scripts and figures in the paper: 
