# PSM-Work-Mapping
Patient-specific modeling assessment of myocardial work with simple and robust clinical datasets

This repo is split into two directories, a Patient-Specific Modeling directory and a Simplified Work Estimates directory. Documentation for each directory is below.


## Patient-Specific Modeling
Adarsh: put your documentation here



## Simplified Work Estimates
This directory contains data and scripts for generating simplified myocardial work measurements from finite element meshes of patient-specific LV endocardium. These geometric models for ech patient are stored in the *LV Geometric Models* folder and include a time-series of patient-specific geometries throughout the cardiac cycle. Scripts for generating and anayzing data are stores in the *scripts* folder. Datasets are stored in the *data* folder.  

These myocardial work measurements approximate the regional stress-strain loop area along the LV endocardial surface.

# Calculating regional strain
Regional strain is computed in *CalculateTriStrain.m*. This script inputs the meshes from the *LV Geometric Models* folder and computes strain as the change in area of each finite element over the course of the cardiac cycle. Regional stran is stored as *RS_CT.mat* This script also computes and stores AHA-based segmental strain for each patient as *seg_strain_allpats.mat*.

# Calculating regional work
In this study, we evaluated five different myocardial work estimates that utilize the same regional strain measurement, but explore different approximation of regional stress.

1. Left heart catheterization-based LV pressure-strain area
  
