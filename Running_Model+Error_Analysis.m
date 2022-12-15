% This code is a run File which includes the error analysis based on the
% measured data stored in PilotDATA_Eddluent. This data was daily so the
% model takes the first predicted value from each day to calculate the
% error.

%% Variables
Degradability_Fraction = 0.5;
K_Hydrolysis = 0.1;
SRT = 50;

%% RUNNING SIMULATION
initilisationFlow_sensitivity; % Loads the storage matrix for flow conditions
adm1init_bsm2_sensitivity;  % Defines the initial conditions for the model and solver

%adm1_ss % Loads Simulink File

sim('adm1_ss'); %Simulate the Pilot under dynamic influent

post_processing;

% Load known effluent data from storage matrix (C1 = Time, C2 = sCOD, C3 =
% pCOD, C4 = TSS, C5 = VSS, C6 = Gas Flow, C7 = Gas Quality)
clf
load PilotData_Effluent.mat

% Calculate values from Model
Total_gas_flow              = digesterout(:,54);
Methane_flow_volume         = (digesterout(:,43)./digesterout(:,53)*P_atm*16/(R_cte*T_op)).*digesterout(:,54)*1.4;
fraction_CH4                = Methane_flow_volume/Total_gas_flow;


%% CALCULATING ERROR
% Because simulink is using a different time interval to the collected data
% cant subtract the values directly. As such am rounding all of the time
% values in the time output down and taking the position of the first cell
% with the time value (using find) and using this to get the modelled 
% outputs which are stored in Modelled_Output.

FLOOR_TIME = floor(time);

Modelled_Output = zeros(440,7);

for j = 127:566;
k = find(FLOOR_TIME==j,1);
Modelled_Output(j-126,1) = j ;
Modelled_Output(j-126,2) = sCOD_out(k,1);
Modelled_Output(j-126,3) = pCOD_out(k,1);
Modelled_Output(j-126,4) = TSS_out(k,1);
Modelled_Output(j-126,5) = VSS_out(k,1);
Modelled_Output(j-126,6) = Total_gas_flow(k,1);
Modelled_Output(j-126,7) = (Methane_flow_volume(k,1)./Total_gas_flow(k,1));
end

% Difference between observed and modelled outputs (C1=sCOD, C2=VSS/TSS, C3
% = GAS FLOW, C4 = CH4 fraction)
Difference = zeros(440,4);
Difference(:,1) = PilotData_Effluent(:,2) - Modelled_Output(:,2);
Difference(:,2) = PilotData_Effluent(:,5)./PilotData_Effluent(:,4) - ...
                    Modelled_Output(:,5)./Modelled_Output(:,4);
Difference(:,3) = PilotData_Effluent(:,6) - Modelled_Output(:,6);
Difference(:,4) = PilotData_Effluent(:,7)/100 - Modelled_Output(:,7);

% Square the difference
Difference_Squared = Difference.^2;
% Sum the squared errors to get SSE
SSE = sum(Difference_Squared,'omitnan');

%% PRINT RESULTS
clc
disp([' The Sum Squared Errors for Parameters'])
disp([' '])
disp(['Degradability          = ',num2str(Degradability_Fraction)])
disp(['Hydrolysis Coefficient = ',num2str(K_Hydrolysis)])
disp(['SRT                    = ',num2str(SRT)])
disp([' '])
disp(['SSE (sCOD)             = ',num2str(SSE(1,1))])
disp(['SSE (VSS/TSS)          = ',num2str(SSE(1,2))])
disp(['SSE (Gas Flow)         = ',num2str(SSE(1,3))])
disp(['SSE (CH_4 Fraction)    = ',num2str(SSE(1,4))])


%% MAKING PLOTS
% This plots the modelled values against the observed values and also
% includes plots of the residuals  below each respective plot.

clf
subplot(2,4,1)
hold on
scatter(PilotData_Effluent(:,1),PilotData_Effluent(:,2))
plot(time,sCOD_out)
axis([127 566 0 0.5])

title("sCOD")
ylabel("kg m^-^3")
xlabel("Days")
grid on
legend("Observed Data","Model")

subplot(2,4,2)
hold on
scatter(PilotData_Effluent(:,1),PilotData_Effluent(:,5)./PilotData_Effluent(:,4))
plot(time, VSS_out./TSS_out)
axis([127 566 0 1])
title("VSS/TSS")
ylabel("Unitless")
xlabel("Days")
grid on
legend("Observed Data","Model")

subplot(2,4,3)
hold on
scatter(PilotData_Effluent(:,1),PilotData_Effluent(:,6))
plot(time,Total_gas_flow)
axis([127 566 0 0.1])
title("Gas Flow")
ylabel("m^3 Day^-^1")
xlabel("Days")
grid on
legend("Observed Data","Model")

subplot(2,4,4)
hold on
scatter(PilotData_Effluent(:,1),PilotData_Effluent(:,7)/100)
plot(time,Methane_flow_volume./Total_gas_flow)

title("Gas Quality")
ylabel("% CH_4 ")
xlabel("Days")
grid on
legend("Observed Data","Model")

subplot(2,4,5)

scatter(Modelled_Output(:,1), Difference(:,1))
yline(0)
xlabel("Time")
ylabel("Residual")
grid on
title("Residual Plot")

subplot(2,4,6)

scatter(Modelled_Output(:,1), Difference(:,2))
yline(0)
xlabel("Time")
ylabel("Residual")
grid on
title("Residual Plot")

subplot(2,4,7)

scatter(Modelled_Output(:,1), Difference(:,3))
yline(0)
xlabel("Time")
ylabel("Residual")
grid on
title("Residual Plot")

subplot(2,4,8)

scatter(Modelled_Output(:,1), Difference(:,4))
yline(0)
xlabel("Time")
ylabel("Residual")
grid on
title("Residual Plot")


