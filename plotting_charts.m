% USER: Javan McGuckin - October 2022

% DESCRIPTION: This code takes the data from the digestor out matrix from
% the simulink simulation and plots it as a time series for visualisation

close all;
%% VARIABLE DEFINITION
load PILOT_DATA_ALL.mat

pH  = digesterout(:,28);

Total_gas_flow              = digesterout(:,54);
Methane_flow_volume         = (digesterout(:,43)./digesterout(:,53)*P_atm*16/(R_cte*T_op)).*digesterout(:,54)*1.4;
Hydrogen_flow_Volume        = (digesterout(:,42)./digesterout(:,53)*P_atm*16/(R_cte*T_op)).*digesterout(:,54)*11.126;
CarbonDioxide_flow_Volume   = (digesterout(:,44)./digesterout(:,53)*P_atm*16/(R_cte*T_op)).*digesterout(:,54)*1.836;

Methane_flow_mass         = (digesterout(:,43)./digesterout(:,53)*P_atm*16/(R_cte*T_op)).*digesterout(:,54);
Hydrogen_flow_mass        = (digesterout(:,42)./digesterout(:,53)*P_atm*16/(R_cte*T_op)).*digesterout(:,54);
CarbonDioxide_flow_mass   = (digesterout(:,44)./digesterout(:,53)*P_atm*16/(R_cte*T_op)).*digesterout(:,54);

Acetate     = digesterout(:,7);
Propionate  = digesterout(:,6);
Butyrate    = digesterout(:,5);
Valerate    = digesterout(:,4);
%% PLOTTING GAS FLOWS
% Plot of the composition and volume of gas flow over the whole time period
% of the data
figure
subplot(2,3,1:2)

area(time, [Methane_flow_volume Hydrogen_flow_Volume CarbonDioxide_flow_Volume])

title("Modelled Gas flow composition chart by Volume")
grid on
xlabel("Time (Days)")
ylabel("Total Gas Flow (m^3/Day)")
legend("Methane","Hydrogen","Carbon Dioxide")
axis([0 600 0 0.125])

% Plot of the total gas flow against the measured gas flow for a HRT of 8
% days
subplot(2,3,4)
hold on
scatter(PILOT_DATA_ALL(:,1)+127,PILOT_DATA_ALL(:,4))
plot(time,Total_gas_flow,'k')
grid on
axis([127 277 0 0.1])

title("Modelled Gas FLow vs Gas Flow Data","8 Day HRT")
xlabel("Time (Days)")
ylabel("Total Gas Flow (m^3/Day)")
legend("Data","Modelled Preformance")

% Plot of the total gas flow against the measured gas flow for a HRT of 4
% days
subplot(2,3,5)
hold on
scatter(PILOT_DATA_ALL(:,1)+127,PILOT_DATA_ALL(:,4))
plot(time,Total_gas_flow,'k')

grid on
axis([277 401 0 0.1])

title("Modelled Gas FLow vs Gas Flow Data","4 Day HRT")
xlabel("Time (Days)")
ylabel("Total Gas Flow (m^3/Day)")
legend("Data","Modelled Preformance")

% Plot of the total gas flow against the measured gas flow for a HRT of 2
% days
subplot(2,3,6)
hold on
scatter(PILOT_DATA_ALL(:,1)+127,PILOT_DATA_ALL(:,4))
plot(time,Total_gas_flow,'k')

grid on
axis([401 566 0 0.1])

title("Modelled Gas FLow vs Gas Flow Data","2 Day HRT")
xlabel("Time (Days)")
ylabel("Total Gas Flow (m^3/Day)")
legend("Data","Modelled Preformance")

%% Plotting pH
subplot(2,3,3)

plot(time,pH)

xline(127)
xline(277)
xline(401)
%axis([0 600 0 8])
title("pH over time")
grid on
xlabel("Time (Days)")
ylabel("pH")
legend("pH")

%% Solids
figure

% CARBOHYDRATES IN
subplot(3,2,1)
plot(time,digesterin(:,14),'r')
title("Carbohydrate concentration in","From Dynamic Input Data")

grid on
xlabel("Time (Days)")
ylabel("Concentration (COD / m^3)")
axis([127 566 0 0.25])

% CARBOHYDRATES OUT
subplot(3,2,2)
plot(time,digesterout(:,14),'r')
title("Carbohydrate concentration out","From Model")

grid on
xlabel("Time (Days)")
ylabel("Concentration (COD / m^3)")
axis([127 566 0 0.25])

% PROTEINES IN
subplot(3,2,3)
plot(time,digesterin(:,15),'g')
title("Proteins concentration in","From Dynamic Input Data")
grid on
xlabel("Time (Days)")
ylabel("Concentration (COD / m^3)")
axis([127 566 0 0.25])

% PROTEINS OUT
subplot(3,2,4)
plot(time,digesterout(:,15),'g')
title("Proteins concentration out","From Model")

grid on
xlabel("Time (Days)")
ylabel("Concentration (COD / m^3)")
axis([127 566 0 0.25])

% LIPIDS IN
subplot(3,2,5)
plot(time,digesterin(:,16),'b')
title("Lipids concentration in","From Dynamic Input Data")

grid on
xlabel("Time (Days)")
ylabel("Concentration (COD / m^3)")
axis([127 566 0 0.125])

% LIPIDS OUT
subplot(3,2,6)
plot(time,digesterout(:,16),'b')
title("Lipids concentration out","From Model")

grid on
xlabel("Time (Days)")
ylabel("Concentration (COD / m^3)")
axis([127 566 0 0.125])



