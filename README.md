# PSM-Work-Mapping
Patient-specific modeling assessment of myocardial work with simple and robust clinical datasets

This repo is split into two directories, a Patient-Specific Modeling directory and a Simplified Work Estimates directory. Documentation for each directory is below.


## Patient-Specific Modeling
Adarsh: put your documentation here



## Simplified Work Estimates
This directory contains data and scripts for generating simplified myocardial work measurements from finite element meshes of patient-specific LV endocardium. These geometric models for ech patient are stored in the *LV Geometric Models* folder and include a time-series of patient-specific geometries throughout the cardiac cycle. Scripts for generating and anayzing data are stores in the *scripts* folder. Datasets are stored in the *data* folder.  

These myocardial work measurements approximate the regional stress-strain loop area along the LV endocardial surface.

### Calculating regional strain
Regional strain is computed in *CalculateTriStrain.m*. This script inputs the meshes from the *LV Geometric Models* folder and computes strain as the change in area of each finite element over the course of the cardiac cycle. Regional stran is stored as *RS_CT.mat* This script also computes and stores AHA segment-based segmental strain for each patient in data file *seg_strain_allpats.mat*.

### Calculating regional work
In this study, we evaluated five different myocardial work estimates that utilize the same regional strain measurement, but explore different approximation of regional stress.

1. Left heart catheterization-based LV pressure-strain area (P_L_H_CSA)
   P_L_H_CSA is computed in *CalculateTriMW.m*. This script inputs the regional strain results *RS_CT.mat* and the patient-specific LV pressure waveforms recorded from pre-CRT left-heart catheterization. These waveforms are stored in *WorkPSMLBBBCRT.mat*, which stores the pre-CRT PSM results. P_L_H_CSA is computed as the area of the LV pressure-regional strain loop area for each element on the mesh, and the result is stores in *MWCT_tris.mat*.
   This script also computes AHA segment-based segmental P_L_H_CSA and PSM-derived work. These results are stored as variables *segMWCT_allpats* and *segWork_allpats* respectively in data file *all_seg_work_all_pats.mat*.
   
  
