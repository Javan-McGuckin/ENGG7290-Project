# ENGG7290-Project
## Overview
This project was part of the authors placement project in the course ENGG7290 at the University of Queensland which looks at using an existing implemenation of the ADM1 in MATLAB to build a model to fit with data supplied from a pilot anaerobic digestion plant. The origional code this is built from is avaliable at https://github.com/wwtmodels/Anaerobic-Digestion-Models (note the base code used here was retrieved from 1. Default Matlab/Simulink ADM1)

The changes made to this work included the addition of a fixed SRT in the model (defined in the initilisation file) which was included as an attempted measure to better model systems with low HRTs. Additionall MSS states were added into the model. Some additional code for calculating the influent conditions from measured field data (equations and work based on the methods outlined in https://www.sciencedirect.com/science/article/pii/S0043135409000475 ), error calculations for model against measured data, post processing to calcluate TSS, VSS, and sCOD, and plots is also included.

## The following list describes each file included

### DATA
* **AD_dynamicconc_pilot.mat** - This file holds the influet data which was calculated based on the measuered values collected through running of the pilot plant and is a direct input into the model
* **AD_steppedflow_pilot.mat** - This file holds influnet data with constant parameters for all values except flow rate, these parameters were calculated using averages of the supplied data values
* **PilotData_Effluent.mat** - This file holds the effluent data that was used for error calculations and for plotting charts to compare model preformance
* **PilotData_measuredInfluientProperties.mat** - Raw data used to calculate the influent model values

### Initilisation files
* **adm1init_bsm2.m** - Initilisation parameters and conditions
* **adm1init_bsm2_sensitivity.m** - Same as above but SRT and K_hydrolysis are defined in run file for sensititivity analysis
* **initilisationFlow.m** - Takes the data from PilotData_measuredInfluientProperties.mat and calculates the model inputs
* **initilisationFlow_sensitivity.m** - same as above but defined degradability fraction in run file so it can be adjusted eaiser and used in loops

### Run Files
* **Running_Model.m** - This script runs the model 
* **Running_Model+Error_Analysis.m** - Runs model and compares error for TSS sCOD Gasflow and gas quality
* **Sensitivity_Analysis** - Creates surface plots of error for different input parameters


### Post Processing Files
* **post_processing.m** - Calculates TSS, VSS, sCOD, gasflow, and corrects for constant SRT
* **plotting_chartss.m** - Generates charts of model preformance over time


### Other
* **adm1_ODE_bsm2.c** - C code for the system (ADM1)
* **adm1_ss.mdl** - Simulink model
* **mexall_adm1.m** - Code to convert C code to matlab readable format

