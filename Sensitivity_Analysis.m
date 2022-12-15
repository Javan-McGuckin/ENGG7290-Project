%% Variables

K_Hydrolysis                 = 0.1;
SRT_Range                    = [4:2:50];
Degradability_Fraction_Range = [0.2:0.1:0.90];

SSE_SCOD = zeros(length(Degradability_Fraction_Range),length(SRT_Range));
SSE_TSVS = zeros(length(Degradability_Fraction_Range),length(SRT_Range));
SSE_FLOW = zeros(length(Degradability_Fraction_Range),length(SRT_Range));
SSE_QUALITY = zeros(length(Degradability_Fraction_Range),length(SRT_Range));
count1=0;
count2=0;


for Degradability_Fraction = Degradability_Fraction_Range
count1 = count1+1
count2=0;

for SRT = SRT_Range
count2 = count2+1
%% RUNNING SIMULATION
initilisationFlow_sensitivity; % Loads the storage matrix for flow conditions
adm1init_bsm2_sensitivity;  % Defines the initial conditions for the model and solver

adm1_ss % Loads Simulink File

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

%% STORING VALUES
SSE_SCOD(count1,count2) = SSE(1,1);
SSE_TSVS(count1,count2) = SSE(1,2);
SSE_FLOW(count1,count2) = SSE(1,3);
SSE_QUALITY(count1,count2) = SSE(1,4);

end

end
%% 
subplot(2,2,1)
surf(SRT_Range,Degradability_Fraction_Range,SSE_SCOD)
xlabel("SRT")
ylabel("Degradability Fraction")
zlabel("SSE in sCOD values")
title("sCOD Error")

subplot(2,2,2)
surf(SRT_Range,Degradability_Fraction_Range,SSE_TSVS)
xlabel("SRT")
ylabel("Degradability Fraction")
zlabel("SSE in TSS/VSS values")
title("TSS/VSS Error")

subplot(2,2,3)
surf(SRT_Range,Degradability_Fraction_Range,SSE_FLOW)
xlabel("SRT")
ylabel("Degradability Fraction")
zlabel("SSE in Flow values")
title("Flow Error")

subplot(2,2,4)
surf(SRT_Range,Degradability_Fraction_Range,SSE_QUALITY)
xlabel("SRT")
ylabel("Degradability Fraction")
zlabel("SSE in CH_4 Composition values")
title("Quality Error")