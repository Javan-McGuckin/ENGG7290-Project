%% DESCRIPTION
% In order to change the degradability of the influent data for the trials
% need to do the data processing in MATLAB. This code takes the raw data
% from the effleunet (which has had some preprocessing including linearly
% interpolating any blank data, and correction of units) and then saves
% this as the input to the system AD_dynamicconc_pilot.mat

%% SET UP

load pilotData_measuredInfluentProperties.mat

% Column 1 is the time values, 2 is s_COD, 3 is p_COD, 4 is acetic acid, 5
% is propanoic acid, 6 is MSS, and 7 is TKN. Units for column 1 is days,
% Units for column 2 through 6 are kg COD/m^3, and column 7 is in kmol
% N/m^3.

% Define storage matrix to fill 
AD_dynamicconc_pilot = zeros(length(pilotData_measuredInfluentProperties)+1,94);

Degradability_Fraction = 0.7;


%% FILLING MATRIX

% TIME VALUES
AD_dynamicconc_pilot([2:end],1)= pilotData_measuredInfluentProperties(:,1);

% CONSTANTS
% Defining Variables which remain constant through the process.
AD_dynamicconc_pilot([2:end],11)  = 0.008;    % IC (kM C / M^3)
AD_dynamicconc_pilot([2:end],13) = 0.03;      % S_I (kg COD / M^3)
AD_dynamicconc_pilot([2:end],27) = 20;        % Temperature (C)
AD_dynamicconc_pilot([2:end],48) = 0.014;     % Sodium (Na)
AD_dynamicconc_pilot([2:end],50) = 0.0002;    % Chlorine (Cl)

% DYNAMIC
   
% S_pro (kg COD / M^3)
AD_dynamicconc_pilot([2:end],7) = pilotData_measuredInfluentProperties(:,5);    

% S_ac (kg COD / M^3)
AD_dynamicconc_pilot([2:end],8) = pilotData_measuredInfluentProperties(:,4);   

% S_su (kg COD / M^3)
AD_dynamicconc_pilot([2:end],2) = (pilotData_measuredInfluentProperties(:,2)...
    - AD_dynamicconc_pilot([2:end],7) - AD_dynamicconc_pilot([2:end],8)...
    - AD_dynamicconc_pilot([2:end],13))/2;  

% S_aa (kg COD / M^3)
AD_dynamicconc_pilot([2:end],3) = AD_dynamicconc_pilot([2:end],2); 

% X_ch (kg COD / M^3)
AD_dynamicconc_pilot([2:end],15) = Degradability_Fraction*0.4*pilotData_measuredInfluentProperties(:,3);  

% X_pr (kg COD / M^3)
AD_dynamicconc_pilot([2:end],16) = Degradability_Fraction*0.4*pilotData_measuredInfluentProperties(:,3);   

% X_li (kg COD / M^3)
AD_dynamicconc_pilot([2:end],17) = Degradability_Fraction*0.2*pilotData_measuredInfluentProperties(:,3);     

% X_I (kg COD / M^3)
AD_dynamicconc_pilot([2:end],25) = (1-Degradability_Fraction)*pilotData_measuredInfluentProperties(:,3);     

% MSS (kg COD / M^3)
AD_dynamicconc_pilot([2:end],28) = pilotData_measuredInfluentProperties(:,6);     

% IN (kM N / M^3)
AD_dynamicconc_pilot([2:end],12) = pilotData_measuredInfluentProperties(:,7)...
    -(0.007*(AD_dynamicconc_pilot([2:end],3)+...
    AD_dynamicconc_pilot([2:end],16)));     

% FLOW RATE
% for this work the flow rate was stepped up from 8 day HRT, to 4 days and
% then 2 days through the running of the pilot plant.
AD_dynamicconc_pilot([2:151],26)   = 0.09;     % Q_in (M^3 / Day)
AD_dynamicconc_pilot([152:276],26) = 0.18;     
AD_dynamicconc_pilot([277:end],26) = 0.36;     

% FIRST ROW
AD_dynamicconc_pilot(1,[2:end])= AD_dynamicconc_pilot(2,[2:end]);
% This sets the first row values to the same as the seccond allowing the
% model to run at a constant influent and reach steady state before using
% the dynamic flow.