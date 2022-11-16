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

load PilotData_Effluent.mat

subplot(5,1,1)
hold on
scatter(PilotData_Effluent(:,1),PilotData_Effluent(:,2))
plot(time,sCOD_out)
axis([127 566 0 0.5])

title("sCOD")
ylabel("kg m^-^3")
xlabel("Days")
grid on
legend("Observed Data","Model")

subplot(5,1,2)
hold on
scatter(PilotData_Effluent(:,1),PilotData_Effluent(:,5)./PilotData_Effluent(:,4))
plot(time, VSS_out./TSS_out)
axis([127 566 0 1])
title("VSS/TSS")
ylabel("Unitless")
xlabel("Days")
grid on
legend("Observed Data","Model")

subplot(5,1,3)
hold on
scatter(PilotData_Effluent(:,1),PilotData_Effluent(:,6))
plot(time,Total_gas_flow)
axis([127 566 0 0.1])
title("Gas Flow")
ylabel("m^3 Day^-^1")
xlabel("Days")
grid on
legend("Observed Data","Model")

subplot(5,1,4)
hold on
scatter(PilotData_Effluent(:,1),PilotData_Effluent(:,7)/100)
plot(time,Methane_flow_volume./Total_gas_flow)

title("Gas Quality")
ylabel("% CH_4 ")
xlabel("Days")
grid on
legend("Observed Data","Model")

subplot(5,1,5)
hold on
plot(time,pH)


title("pH over time")
ylabel("pH")
xlabel("Days")
grid on
