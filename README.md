# ENGG7290-Project
This project looks at using the ADM1 in MATLAB to build a model for a pilot anaerobic digestion facility. The files are described below.

AD_dynamicconc_pilot.mat 

This is a matrix which stores the influent paramerts which were measured in the pilot facility through its operation which are used in modelling the system.

AD_steppedflow_pilot.mat

This is a matrix which has constant paremeter for all values except the flow rate into the unit which was used for modelling just the effect of flow with constant parameters and concentrations in the influent.

Error_Analysis.m

This m file is used to test the error in sCOD, VSS/TSS ratio, gas flow, and gas composition, for different degradability coeffieicents, SRT's, and Hydrolysis Coefficients.

Running_Model.m

This m file is used to run the model.

Sensitivity_Analysis.m

This file is used to compare the SSE for a range of combinations of different degradability coeffieicents, SRT's, and Hydrolysis Coefficients.

adm1_ODE_bsm2.c

This is the C code which is used in the ODE solver.

adm1init_Pilot.m and adm1init_BSM2.

These files define the initial and system parameters for the model.

adm1init_bsm_sensitivity.m

This is the same as adm1init_BSM2 however it passes the Hydrolysis coefficient, SRT, and degradibility fraction from the running codes rather then as a constent.

initilisationFlow.m

This file takes the measured flow parameters and converts them to the input paremeters required in the modelling based on the set degradibility fraction.
