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

***1. Left heart catheterization-based LV pressure-strain area (P<sub>LHC</sub>SA)***

   Scripts needed to compute work estimates:  *CalculateTriMW.m*

   Data files needed: *RS_CT.mat*, *WorkPSMLBBBCRT.mat*

   Scripts where necessary data is generated: *PrincipalCurvatureAnalysis2.mat*
   
   File where results are stored: *MWCT_tris.mat*

   ***Description:*** P<sub>LHC</sub>SA is computed in *CalculateTriMW.m*. This script inputs the regional strain results *RS_CT.mat* and the patient-specific LV pressure waveforms recorded from pre-CRT left-heart catheterization. These waveforms are stored in *WorkPSMLBBBCRT.mat*, which stores the pre-CRT PSM results. P<sub>LHC</sub>SA is computed as the area of the LV pressure-regional strain loop area for each element on the mesh, and the result is stored in *MWCT_tris.mat*.

   This script also computes AHA segment-based segmental P<sub>LHC</sub>SA and PSM-derived work. These results are stored as variables *segMWCT_allpats* and *segWork_allpats* respectively in data file *all_seg_work_all_pats.mat*.

***2. End-diastolic wall stress-strain area (WS<sub>ED</sub>SA)***
***3. Time-varying wall stress-strain area (WS<sub>TV</sub>SA)***
   
   Scripts needed to compute work estimates: *calculateSimpleWork.m*
   
   Data files needed: *Laplace_measurements.mat*
   
   Scripts where necessary data is generated: *PrincipalCurvatureAnalysis2.mat*
   
   File where results are stored: *Laplace_work_all_patches.mat*
   
   ***Description:*** These two work estimates leverage the law of Laplace to compute regional wall stress, which states that wall stress is proportional to cavity pressure, wall thickness, and cavity radius. WS<sub>ED</sub>SA applied  regional shape information at end-diastole, while WS<sub>TV</sub>SA applied regional shape information for the whole cardiac cycle. These wall stresses are computed in *calculateSimpleWork.m* using the function *calculatePatchWorkEstimates*. (WS<sub>ED</sub>SA) is stored as variable *MWCTLP* and (WS<sub>TV</sub>SA) is stored as variable *MWCTLPTV* in the data file *Laplace_work_all_patches.mat*.

   To calculate these two approximations, regional radius and wall thickness information is needed. Both variables are generated in script *PrincipalCurvatureAnalysis2.mat*. Regional radius is computed as the regional effective radius, which is calculated for each AHA segment by fitting an ellipsoid to the patient-specific geometric model. The effective radius is stored as variable *r_eff_allpats*. Regional wall thickness information was deried from the patient-specific geometric model. Wall thickness data is stored in the *Segmental Wall Thickness Measurements* folder and saved as variable *thickness_data*. Both variables are saved in data file *Laplace_measurements.mat*
  
   
  
