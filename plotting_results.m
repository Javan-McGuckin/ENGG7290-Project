%% ORIGIONAL CODE
[m n]=size(digesterout);

disp(' ')
disp('This script will provide the final values of the ADM1 ')
disp(' ')

disp('ADM1 influent conditions ')
disp('*********************')
disp(['   S_su =  ', num2str(digesterin(m,1)), ' kg COD/m3']);
disp(['   S_aa =  ', num2str(digesterin(m,2)), ' kg COD/m3']);
disp(['   S_fa =  ', num2str(digesterin(m,3)), ' kg COD/m3']);
disp(['   S_va =  ', num2str(digesterin(m,4)), ' kg COD/m3']);
disp(['   S_bu =  ', num2str(digesterin(m,5)), ' kg COD/m3']);
disp(['   S_pro =  ', num2str(digesterin(m,6)), ' kg COD/m3']);
disp(['   S_ac =  ', num2str(digesterin(m,7)), ' kg COD/m3']);
disp(['   S_h2 =  ', num2str(digesterin(m,8)), ' kg COD/m3']);
disp(['   S_ch4 =  ', num2str(digesterin(m,9)), ' kg COD/m3']);
disp(['   S_IC =  ', num2str(digesterin(m,10)), ' kmol C/m3']);
disp(['   S_IN =  ', num2str(digesterin(m,11)), ' kmol N/m3']);
disp(['   S_I =  ', num2str(digesterin(m,12)), ' kg COD/m3']);


disp(['   X_ch =  ', num2str(digesterin(m,14)), ' kg COD/m3']);
disp(['   X_pr =  ', num2str(digesterin(m,15)), ' kg COD/m3']);
disp(['   X_li =  ', num2str(digesterin(m,16)), ' kg COD/m3']);
disp(['   X_su =  ', num2str(digesterin(m,17)), ' kg COD/m3']);
disp(['   X_aa =  ', num2str(digesterin(m,18)), ' kg COD/m3']);
disp(['   X_fa =  ', num2str(digesterin(m,19)), ' kg COD/m3']);
disp(['   X_c4 =  ', num2str(digesterin(m,20)), ' kg COD/m3']);
disp(['   X_pro =  ', num2str(digesterin(m,21)), ' kg COD/m3']);
disp(['   X_ac =  ', num2str(digesterin(m,22)), ' kg COD/m3']);
disp(['   X_h2 =  ', num2str(digesterin(m,23)), ' kg COD/m3']);
disp(['   X_I =  ', num2str(digesterin(m,24)), ' kg COD/m3']);

disp(' ')
disp(['   VFA =  ', num2str(sum(digesterin(m,4:7))), ' kg COD/m3']);
disp(['   CODsol =  ', num2str(sum(digesterin(m,1:7))), ' kg COD/m3']);
disp(['   CODpart =  ', num2str(sum(digesterin(m,12:24))), ' kg COD/m3']);
disp(['   Qin =  ', num2str(digesterin(m,25)), ' m3/d']);
disp(' ')
disp('ADM1 state variables ')
disp('*********************')
disp(['   S_su =  ', num2str(digesterout(m,1)), ' kg COD/m3']);
disp(['   S_aa =  ', num2str(digesterout(m,2)), ' kg COD/m3']);
disp(['   S_fa =  ', num2str(digesterout(m,3)), ' kg COD/m3']);
disp(['   S_va =  ', num2str(digesterout(m,4)), ' kg COD/m3']);
disp(['   S_bu =  ', num2str(digesterout(m,5)), ' kg COD/m3']);
disp(['   S_pro =  ', num2str(digesterout(m,6)), ' kg COD/m3']);
disp(['   S_ac =  ', num2str(digesterout(m,7)), ' kg COD/m3']);
disp(['   S_h2 =  ', num2str(digesterout(m,8)), ' kg COD/m3']);
disp(['   S_ch4 =  ', num2str(digesterout(m,9)), ' kg COD/m3']);
disp(['   S_IC =  ', num2str(digesterout(m,10)), ' kmol C/m3']);
disp(['   S_IN =  ', num2str(digesterout(m,11)), ' kmol N/m3']);
disp(['   S_I =  ', num2str(digesterout(m,12)), ' kg COD/m3']);


disp(['   X_ch =  ', num2str(digesterout(m,14)), ' kg COD/m3']);
disp(['   X_pr =  ', num2str(digesterout(m,15)), ' kg COD/m3']);
disp(['   X_li =  ', num2str(digesterout(m,16)), ' kg COD/m3']);
disp(['   X_su =  ', num2str(digesterout(m,17)), ' kg COD/m3']);
disp(['   X_aa =  ', num2str(digesterout(m,18)), ' kg COD/m3']);
disp(['   X_fa =  ', num2str(digesterout(m,19)), ' kg COD/m3']);
disp(['   X_c4 =  ', num2str(digesterout(m,20)), ' kg COD/m3']);
disp(['   X_pro =  ', num2str(digesterout(m,21)), ' kg COD/m3']);
disp(['   X_ac =  ', num2str(digesterout(m,22)), ' kg COD/m3']);
disp(['   X_h2 =  ', num2str(digesterout(m,23)), ' kg COD/m3']);
disp(['   X_I =  ', num2str(digesterout(m,24)), ' kg COD/m3']);
disp(' ')
disp(['   VFA =  ', num2str(sum(digesterout(m,4:7))), ' kg COD/m3']);
disp(['   CODsol =  ', num2str(sum(digesterout(m,1:7))), ' kg COD/m3']);
disp(['   CODpart =  ', num2str(sum(digesterout(m,12:24))), ' kg COD/m3']);
disp(['   Qin =  ', num2str(digesterout(m,25)), ' m3/d']);
disp(' ')

disp(' ')
disp('ADM1 Scat(+)/San(-) extension state variables ')
disp('**********************************************')
disp(' ')
disp(['   S_Na =  ', num2str(digesterout(m,75)), ' kmol /m3']);
disp(['   S_K =  ', num2str(digesterout(m,76)),  ' kmol /m3']);
disp(['   S_Cl =  ', num2str(digesterout(m,77)), ' kmol /m3']);


disp(' ')
disp('Gas transpher variables ')
disp('************************')


Methanevec = digesterout(m,43)./digesterout(m,53)*P_atm*16/(R_cte*T_op); %kg CH4/m3
Methaneflowvec = Methanevec.*digesterout(m,54);


Hydrogenvec = digesterout(m,42)./digesterout(m,53)*P_atm*2/(R_cte*T_op); %kg H2/m3
Hydrogenflowvec = Hydrogenvec.*digesterout(m,54);

Carbondioxidevec = digesterout(m,44)./digesterout(m,53)*P_atm*44/(R_cte*T_op); %kg CO2/m3
Carbondioxideflowvec = Carbondioxidevec.*digesterout(m,54);

Qgasvec= digesterout(m,54);

disp(['Average methane production = ',num2str(Methaneflowvec),' kg CH4/d = ', num2str(Methaneflowvec*50.014/3.6),' kWh/d'])
disp(['Average hydrogen gas production (kg H2/d) = ', num2str(Hydrogenflowvec),' kg H2/d']);
disp(['Average carbon dioxide gas production (kg CO2/d) = ', num2str(Carbondioxideflowvec),' kg CO2/d']);
disp(['Average total gas flow rate (AD, normalized to P_atm) = ', num2str(Qgasvec), ' m3/d']);

disp(' ')
disp('Display the H+ and pH ')
disp('**************************************')
disp(' ')

disp(['   S_H+ =  ', num2str(digesterout(m,29)), ' mol/L']);
disp(['   pH =  ', num2str(digesterout(m,28)), ' ']);
disp(' ')
disp(['   SVa- =  ', num2str(digesterout(m,78)), ' kgCOD/m3 ']);
disp(['   SBu- =  ', num2str(digesterout(m,79)), ' kgCOD/m3 ']);
disp(['   Spro- =  ', num2str(digesterout(m,80)), ' kgCOD/m3 ']);
disp(['   SAc- =  ', num2str(digesterout(m,81)), ' kgCOD/m3 ']);
disp(' ')
disp(['   ShVa =  ', num2str(digesterout(m,4) - digesterout(m,78)), ' kgCOD/m3 ']);
disp(['   ShBu =  ', num2str(digesterout(m,5) - digesterout(m,79)), ' kgCOD/m3 ']);
disp(['   Shpro =  ', num2str(digesterout(m,6) -digesterout(m,80)), ' kgCOD/m3 ']);
disp(['   ShAc =  ', num2str(digesterout(m,7) -digesterout(m,81)), ' kgCOD/m3 ']);
disp(' ')
disp(['   Shco3- =  ', num2str(digesterout(m,82)), ' kmol/m3 ']);
disp(['   Snh3 =  ', num2str(digesterout(m,83)), ' kmol/m3 ']);
disp(['   Sco2 =  ', num2str(digesterout(m,10)-digesterout(m,82)), ' kmol/m3 ']);
disp(['   Snh4 =  ', num2str(digesterout(m,11)-digesterout(m,83)), ' kmol/m3 ']);
disp(' ')
disp('COD balance ')
disp('*********************')
disp(['   COD load in =  ',  num2str((sum(digesterin(m,1:9)) +   sum(digesterin(m,12:24))  + digesterin(m,28) +  digesterin(m,30) +   sum(digesterin(m,38:43 )))*digesterin(m,25)  + AD_dynamicconc_pilot(end,91) + AD_dynamicconc_pilot(end,92)), ' kg COD/d']);
disp(['   COD load out = ',  num2str((sum(digesterout(m,1:9)) + sum(digesterout(m,12:24)) +  digesterout(m,56) + digesterout(m,58) +  sum(digesterout(m,66:71)))*digesterout(m,25) + Methaneflowvec.*64/16 + Hydrogenflowvec.*16/2), ' kg COD/d']);
disp(['   COD conversion (gas) = ',  num2str(Methaneflowvec.*64/16/((sum(digesterout(m,1:9)) + sum(digesterout(m,12:24)) +  digesterout(m,56) + digesterout(m,58) +  sum(digesterout(m,66:71)))*digesterout(m,25) + Methaneflowvec.*64/16 + Hydrogenflowvec.*16/2)), ' kg COD/d']);
disp(['   COD conversion (liq) = ',  num2str(((sum(digesterout(m,1:9)) + sum(digesterout(m,12:24)) +  digesterout(m,56) + digesterout(m,58) +  sum(digesterout(m,66:71)))*digesterout(m,25) + Hydrogenflowvec.*16/2)/((sum(digesterout(m,1:9)) + sum(digesterout(m,12:24)) +  digesterout(m,56) + digesterout(m,58) +  sum(digesterout(m,66:71)))*digesterout(m,25) + Methaneflowvec.*64/16 + Hydrogenflowvec.*16/2)), ' kg COD/d']);
disp(' ')
disp('C balance ')
disp('*********************')
totCinfAD = (digesterin(m,1)*C_su + digesterin(m,2)*C_aa + digesterin(m,3)*C_fa + digesterin(m,4)*C_va + digesterin(m,5)*C_bu + digesterin(m,6)*C_pro + digesterin(m,7)*C_ac + digesterin(m,9)*C_ch4 + digesterin(m,10) + digesterin(m,12)*C_xI  + (sum(digesterin(m,17:23)) + digesterin(m,30)+ sum(digesterin(m,39:43)))*C_bac + digesterin(m,14)*C_ch + digesterin(m,15)*C_pr +  digesterin(m,16)*C_li + digesterin(m,24)*C_xI  )*digesterin(m,25);
totCeffAD = (digesterout(m,1)*C_su + digesterout(m,2)*C_aa + digesterout(m,3)*C_fa + digesterout(m,4)*C_va + digesterout(m,5)*C_bu + digesterout(m,6)*C_pro + digesterout(m,7)*C_ac + digesterout(m,9)*C_ch4 + digesterout(m,10) + digesterout(m,12)*C_xI  + (sum(digesterout(m,17:23)) + digesterout(m,58)+ sum(digesterout(m,67:71)))*C_bac + digesterout(m,14)*C_ch + digesterout(m,15)*C_pr +  digesterout(m,16)*C_li + digesterout(m,24)*C_xI  )*digesterout(m,25) + Methaneflowvec.*1/16 + Carbondioxideflowvec*1/44;

disp(['   C load in =  ',  num2str(totCinfAD), ' kmol C/d']);
disp(['   C load out = ',  num2str(totCeffAD), ' kmol C/d']);
disp(' ')
disp('N balance ')
disp('*********************')

totNinfAD = (digesterin(m,11) + digesterin(m,2)*N_aa + digesterin(m,15)*N_aa + (sum(digesterin(m,17:23)) + digesterin(m,30)+ sum(digesterin(m,39:43)))*N_bac + digesterin(m,12)*N_I + digesterin(m,24)*N_I)*digesterin(m,25);
totNeffAD = (digesterout(m,11) + digesterout(m,2)*N_aa + digesterout(m,15)*N_aa + (sum(digesterout(m,17:23)) + digesterout(m,58)+ sum(digesterout(m,67:71)))*N_bac + digesterout(m,12)*N_I + digesterout(m,24)*N_I + digesterout(m,101))*digesterout(m,25);

disp(['   N load in =  ',  num2str(totNinfAD), ' kmol N/d']);
disp(['   N load out = ',  num2str(totNeffAD), ' kmol N/d']);

%% CORRECTING FOR SRT and ADDITIONAL OUTPUTS

% This section of the code is set up to correct the outflows of the model
% for the SRT. In the C code currently the outputs are mapped directly to
% the state which is impacted by the seperatre SRT and doesnt represent the
% actual outflow from the model;

% Define Storage matrix to populate with corrected values
Corrected_Output_Values = ones(length(digesterout),24);

Corrected_Output_Values(:,1:24) = (digesterout(:,1:24))*V_liq./(digesterout(:,25)*DIGESTERPAR(101));

% Calculate pCOD out of the model by summing all of the Xi terms in the
% model outputs, note well does not include MSS. Similarly can calculate
% the sCOD and TCOD

pCOD_out = sum(Corrected_Output_Values(:,14:24),2);
sCOD_out = sum(Corrected_Output_Values(:,1:9),2);
tCOD_out = pCOD_out + sCOD_out;

% Calculate the corrected MSS value using the same equation used to correct
% other values 
MSS_out = (digesterout(:,33))*V_liq./(digesterout(:,25)*DIGESTERPAR(101));

% Calculating VSS Note well here the values being used were taken from a 
% later version of the model which was used to account for PHA.

VSS_out = sum(Corrected_Output_Values(:,14:16),2)/1.5686 + ...
          sum(Corrected_Output_Values(:,17:23),2)/1.3072 + ...
          Corrected_Output_Values(:,24)/1.5686;

TSS_out = MSS_out + VSS_out;



