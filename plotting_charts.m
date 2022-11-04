%% USER: Javan McGuckin - November 2022
% This code is used to plot outputs of the model into charts. Each chart
% has been allocated to a section in the code and includes plotting gas
% flows, plotting the COD and SS values and comparing the model before and
% after correction for SRT in the outlet.

close all;

%% FIGURE 1 - PLOTTING GAS FLOWS
% VARIABLE DEFINITION
pH  = digesterout(:,28);

Total_gas_flow              = digesterout(:,54);
Methane_flow_volume         = (digesterout(:,43)./digesterout(:,53)*P_atm*16/(R_cte*T_op)).*digesterout(:,54)*1.4;
Hydrogen_flow_Volume        = (digesterout(:,42)./digesterout(:,53)*P_atm*16/(R_cte*T_op)).*digesterout(:,54)*11.126;
CarbonDioxide_flow_Volume   = (digesterout(:,44)./digesterout(:,53)*P_atm*16/(R_cte*T_op)).*digesterout(:,54)*1.836;

Methane_flow_mass         = (digesterout(:,43)./digesterout(:,53)*P_atm*16/(R_cte*T_op)).*digesterout(:,54);
Hydrogen_flow_mass        = (digesterout(:,42)./digesterout(:,53)*P_atm*16/(R_cte*T_op)).*digesterout(:,54);
CarbonDioxide_flow_mass   = (digesterout(:,44)./digesterout(:,53)*P_atm*16/(R_cte*T_op)).*digesterout(:,54);
% Plot of the composition and volume of gas flow over the whole time period
% of the data
figure
load GAS_DATA.mat
subplot(3,3,1:2)

area(time, [Methane_flow_volume Hydrogen_flow_Volume CarbonDioxide_flow_Volume])

title("Modelled Gas flow"," composition chart by Volume")
grid on
xlabel("Time (Days)")
ylabel("Total Gas Flow (m^3/Day)")
legend("Methane","Hydrogen","Carbon Dioxide")
axis([127 566 0 0.125])

subplot(3,3,4:5)

area(GAS_DATA(:,1), [GAS_DATA(:,3) GAS_DATA(:,4)])

title("Data Gas flow"," composition chart by Volume")
grid on
xlabel("Time (Days)")
ylabel("Total Gas Flow (m^3/Day)")
legend("Methane","Other")
axis([127 566 0 0.125])

% Plot of the total gas flow against the measured gas flow for a HRT of 8
% days
subplot(3,3,7:9)
hold on
scatter(GAS_DATA(:,1),GAS_DATA(:,2))
plot(time,Total_gas_flow,'k')
grid on
axis([127 566 0 0.1])

title("Modelled Gas FLow vs Gas Flow Data")
xlabel("Time (Days)")
ylabel("Total Gas Flow (m^3/Day)")
legend("Data","Modelled Preformance")

% Plotting pH
subplot(3,3,[3 6])

plot(time,pH)

xline(127)
xline(277)
xline(401)
axis([0 600 4 9])
title("pH over time")
grid on
xlabel("Time (Days)")
ylabel("pH")
legend("pH")

%% FIGURE 2 - COD and SS Charts influent and Effluent
load COD_SS_DATA.mat
figure
% COD
subplot(5,2,1)

scatter(COD_SS_DATA(:,1),COD_SS_DATA(:,2))

title("COD","COD in influent (Modelled")
xlabel("Time (Days)")
ylabel("COD (kg m^-^3)")
grid on
axis([127 566 0 2])

subplot(5,2,2)
hold on
scatter(COD_SS_DATA(:,1),COD_SS_DATA(:,7))
plot(time,tCOD_out)

title("COD","COD in effluent (Modelled and Data)")
xlabel("Time (Days)")
ylabel("COD (kg m^-^3)")
grid on
axis([127 566 0 2])
legend("Data","Model")

% sCOD
subplot(5,2,3)

scatter(COD_SS_DATA(:,1),COD_SS_DATA(:,3))

title("sCOD","sCOD in influent (Modelled")
xlabel("Time (Days)")
ylabel("sCOD (kg m^-^3)")
grid on
axis([127 566 0 1.5])

subplot(5,2,4)
hold on
scatter(COD_SS_DATA(:,1),COD_SS_DATA(:,8))
plot(time,sCOD_out)

title("sCOD","sCOD in effluent (Modelled and Data)")
xlabel("Time (Days)")
ylabel("sCOD (kg m^-^3)")
grid on
axis([127 566 0 1.5])
legend("Data","Model")

% pCOD
subplot(5,2,5)

scatter(COD_SS_DATA(:,1),COD_SS_DATA(:,4))

title("pCOD","pCOD in influent (Modelled")
xlabel("Time (Days)")
ylabel("pCOD (kg m^-^3)")
grid on
axis([127 566 0 1.5])

subplot(5,2,6)
hold on
scatter(COD_SS_DATA(:,1),COD_SS_DATA(:,9))
plot(time,pCOD_out)

title("pCOD","pCOD in effluent (Modelled and Data)")
xlabel("Time (Days)")
ylabel("pCOD (kg m^-^3)")
grid on
axis([127 566 0 1.5])
legend("Data","Model")

% TSS
subplot(5,2,7)

scatter(COD_SS_DATA(:,1),COD_SS_DATA(:,5))

title("TSS","TSS in influent (Modelled")
xlabel("Time (Days)")
ylabel("TSS (kg m^-^3)")
grid on
axis([127 566 0 1])

subplot(5,2,8)
hold on
scatter(COD_SS_DATA(:,1),COD_SS_DATA(:,10))
plot(time,TSS_out)

title("TSS","TSS in effluent (Modelled and Data)")
xlabel("Time (Days)")
ylabel("TSS (kg m^-^3)")
grid on
axis([127 566 0 1])
legend("Data","Model")

% VSS
subplot(5,2,9)

scatter(COD_SS_DATA(:,1),COD_SS_DATA(:,6))

title("VSS","VSS in influent (Modelled")
xlabel("Time (Days)")
ylabel("VSS (kg m^-^3)")
grid on
axis([127 566 0 1])

subplot(5,2,10)
hold on
scatter(COD_SS_DATA(:,1),COD_SS_DATA(:,11))
plot(time,VSS_out)

title("VSS","VSS in effluent (Modelled and Data)")
xlabel("Time (Days)")
ylabel("VSS (kg m^-^3)")
grid on
axis([127 566 0 1])
legend("Data","Model")

%% FIGURE 3 - Comparing corrected and not corrected values
figure
subplot(3,2,1)

hold on
scatter(COD_SS_DATA(:,1),COD_SS_DATA(:,7))
plot(time,tCOD_out)
plot(time,sum(digesterout(:,1:9),2)+sum(digesterout(:,14:24),2))
title("COD","COD in effluent (Modelled and Data)")
xlabel("Time (Days)")
ylabel("COD (kg m^-^3)")
grid on
axis([127 566 0 5])
legend("Data","Model Corrected","Model not Corrected for SRT")

% sCOD


subplot(3,2,3)
hold on
scatter(COD_SS_DATA(:,1),COD_SS_DATA(:,8))
plot(time,sCOD_out)
plot(time,sum(digesterout(:,1:9),2))

title("sCOD","sCOD in effluent (Modelled and Data)")
xlabel("Time (Days)")
ylabel("sCOD (kg m^-^3)")
grid on
axis([127 566 0 3.5])
legend("Data","Model Corrected","Model not Corrected for SRT")

% pCOD

subplot(3,2,5)
hold on
scatter(COD_SS_DATA(:,1),COD_SS_DATA(:,9))
plot(time,pCOD_out)
plot(time,sum(digesterout(:,14:24),2))

title("pCOD","pCOD in effluent (Modelled and Data)")
xlabel("Time (Days)")
ylabel("pCOD (kg m^-^3)")
grid on
axis([127 566 0 1.5])
legend("Data","Model Corrected","Model not Corrected for SRT")

% TSS


subplot(3,2,2)
hold on
scatter(COD_SS_DATA(:,1),COD_SS_DATA(:,10))
plot(time,TSS_out)
plot(time,sum(digesterout(:,14:16),2)/1.5686 + ...
          sum(digesterout(:,17:23),2)/1.3072 + ...
          digesterout(:,24)/1.5686 + digesterout(:,33))

title("TSS","TSS in effluent (Modelled and Data)")
xlabel("Time (Days)")
ylabel("TSS (kg m^-^3)")
grid on
axis([127 566 0 1])
legend("Data","Model Corrected","Model not Corrected for SRT")

% VSS


subplot(3,2,4)
hold on
scatter(COD_SS_DATA(:,1),COD_SS_DATA(:,11))
plot(time,VSS_out)
plot(time,sum(digesterout(:,14:16),2)/1.5686 + ...
          sum(digesterout(:,17:23),2)/1.3072 + ...
          digesterout(:,24)/1.5686)
title("VSS","VSS in effluent (Modelled and Data)")
xlabel("Time (Days)")
ylabel("VSS (kg m^-^3)")
grid on
axis([127 566 0 1])
legend("Data","Model Corrected","Model not Corrected for SRT")

subplot(3,2,6)
hold on
scatter(COD_SS_DATA(:,1),COD_SS_DATA(:,10)-COD_SS_DATA(:,11))
plot(time,MSS_out)
plot(time,digesterout(:,33))
title("VSS","VSS in effluent (Modelled and Data)")
xlabel("Time (Days)")
ylabel("VSS (kg m^-^3)")
grid on
axis([127 566 0 0.3])
legend("Data","Model Corrected","Model not Corrected for SRT")
