clear all; clc;
% DESCRIPTION: This file is used to run the model for the pilot system.
% Note well this is the model of ADM1 for simulink/MATLAB from the work ()
% with some adjustments for this project.

%% SET UP 

adm1init_bsm2;  % Defines the initial conditions for the model and solver
load AD_dynamicconc_pilot.mat; % Loads the storage matrix for flow conditions
adm1_ss % Loads Simulink File

%% RUNNING & PLOTTING

disp('Simulating ADM1 with constant influent, Solver = ode15s ')
disp('********************************************************')
disp(' ')
start=clock; 
disp(['Start time for simulation (hour:min:sec) = ', num2str(round(start(4:6)))]); %Display simulation start time 

sim('adm1_ss'); %Simulate the Pilot under dynamic influent

%% DATA PRINTING AND PROCESSING

% File which prints key final values and balances into the command window
plotting_results

% File which plots desired charts

plotting_charts

