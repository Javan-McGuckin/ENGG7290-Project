/*
 * adm1_ODE_bsm2.c is a C-file S-function for IAWQ AD Model No 1. 
 * In addition to the ADM1, temperature dependency and dummy states are added.
 * Some details are adjusted for BSM2 (pH inhibition, gas flow output etc).
 * In this the model derivatives of ion states anf hydrogen gas are set to zero.
 * Instead they are calculated using algebraic equations and a pH solver
 * using the Newton-Raphson method nd an equivalent H2-solver). This way the main stiffness is removed.
 * In addition to the ADM1, temperature dependency and dummy states are added.
 * Some details are adjusted for BSM2 (pH inhibition, gas flow output etc).
 *
 * Copyright (2006):
 * Dr Christian Rosen, Dr Darko Vrecko and Dr Ulf Jeppsson
 * Dept. Industrial Electrical Engineering and Automation (IEA)
 * Lund University, Sweden
 * http://www.iea.lth.se/ 
 *
 * Modified by Xavier Flores-Alsina in order to include PC corrections (and other extensions)
 * Xc is removed from the list of state variables. Decay process is automatically mapped into Xch, Xli, Xpro, Xi and Si
 * May 2015
 *
 */

#define S_FUNCTION_NAME adm1_ODE_bsm2

#include "simstruc.h"
#include <math.h>

#define XINIT   ssGetArg(S,0)
#define PAR  	ssGetArg(S,1)
#define V	    ssGetArg(S,2)


/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 105);  /* number of continuous states           */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states             */
    ssSetNumInputs(        S, 93);  /* number of inputs                      */
    ssSetNumOutputs(       S, 150);  /* number of outputs                     */
    ssSetDirectFeedThrough(S, 1);   /* direct feedthrough flag               */
    ssSetNumSampleTimes(   S, 1);   /* number of sample times                */
    ssSetNumSFcnParams(    S, 3);   /* number of input arguments             */
    ssSetNumRWork(         S, 0);   /* number of real work vector elements   */
    ssSetNumIWork(         S, 0);   /* number of integer work vector elements*/
    ssSetNumPWork(         S, 0);   /* number of pointer work vector elements*/
}


/*
 * mdlInitializeSampleTimes - initialize the sample times array
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}


/*
 * mdlInitializeConditions - initialize the states
 */
static void mdlInitializeConditions(double *x0, SimStruct *S)
{
int i;

for (i = 0; i < 105; i++) {
   x0[i] = mxGetPr(XINIT)[i];
}
}


/*
 * mdlOutputs - compute the outputs
 */
static void mdlOutputs(double *y, double *x, double *u, SimStruct *S, int tid)
{
  double R, T_base, T_op, P_atm, p_gas_h2o, P_gas, k_P, q_gas, V_liq, procT8, procT9, procT10, p_gas_h2, p_gas_ch4, p_gas_co2, t_res;
  double kLa, K_H_h2_base, K_H_ch4_base, K_H_co2_base, phi, S_H_ion, pK_w_base, K_H_h2o_base, K_H_h2, K_H_ch4, K_H_co2, K_w, factor;
  double S_h2, S_ch4, S_co2;             
  double S_gas_h2, S_gas_ch4, S_gas_co2; 
  double pH_op;
  
  int i;
  
  R = mxGetPr(PAR)[77];
  T_base = mxGetPr(PAR)[78];
  T_op = mxGetPr(PAR)[79];
  P_atm = mxGetPr(PAR)[93];
  p_gas_h2o = mxGetPr(PAR)[94];
  V_liq = mxGetPr(V)[0];
  kLa = mxGetPr(PAR)[95];
  K_H_h2_base = mxGetPr(PAR)[98];
  K_H_ch4_base = mxGetPr(PAR)[97];
  K_H_co2_base = mxGetPr(PAR)[96];
  K_H_h2o_base = mxGetPr(PAR)[95];
  pK_w_base = mxGetPr(PAR)[80];
  k_P = mxGetPr(PAR)[99];
  t_res = mxGetPr(PAR)[100];

  factor = (1.0/T_base - 1.0/T_op)/(100.0*R);
  K_H_h2 = K_H_h2_base*exp(-4180.0*factor);     /* T adjustment for K_H_h2 */
  K_H_ch4 = K_H_ch4_base*exp(-14240.0*factor);  /* T adjustment for K_H_ch4 */
  K_H_co2 = K_H_co2_base*exp(-19410.0*factor);  /* T adjustment for K_H_co2 */
  K_w = pow(10,-pK_w_base)*exp(55900.0*factor); /* T adjustment for K_w */
  p_gas_h2o = K_H_h2o_base*exp(5290.0*(1.0/T_base - 1.0/T_op));  /* T adjustement for water vapour saturation pressure */
  
  for (i = 0; i < 7; i++) {
      y[i] = x[i];
  }
  
  y[7] = x[7];   /* Sh2 */
  
  for (i = 0; i < 16; i++) {
      y[i+8] = x[i+8];
  }
   
  y[24] = u[24];   /* flow */
  
  y[25] = T_op - 273.15;      /* Temp = 35 degC */
  
  y[26] = 0.0;             /* H2 not in use in this implementation  */
  
  
  phi = x[56]+(x[10]-x[64])-x[63]-x[62]/64.0-x[61]/112.0-x[60]/160.0-x[59]/208.0-x[58]; // charge balance (positive and negative charges)
  S_H_ion = -phi*0.5+0.5*sqrt(phi*phi+4.0*K_w); /* SH+ */
  pH_op = -log10(S_H_ion); /* pH */
  
  y[27] = pH_op;    /* pH */
  y[28] = S_H_ion;            /* SH+ */  
  
  /* MASS TRANSFER EQUATIONS */  
  
  y[29] = x[24];   /* Sh2 */   
  y[30] = x[25];   /* Sch4 */
  y[31] = x[26];   /* Sco2 */
  y[32] = x[27];   /* D2 - Now using for mineral solids*/
  y[33] = x[28];   /* D3 */
  y[34] = x[29];   /* D4 */   
  y[35] = x[30];   /* D5 */
  y[36] = x[31];   /* D6 */
  y[37] = x[32];   /* D7 */
  y[38] = x[33];   /* D8 */
  y[39] = x[34];   /* D7 */
  y[40] = x[35];   /* D8 */
  
 
  S_gas_h2 =  x[24];
  S_gas_ch4 = x[25];
  S_gas_co2 = x[26];       
         
  p_gas_h2 =  S_gas_h2*R*T_op/16.0;
  p_gas_ch4 = S_gas_ch4*R*T_op/64.0;
  p_gas_co2 = S_gas_co2*R*T_op;
  P_gas = p_gas_h2 + p_gas_ch4 + p_gas_co2 + p_gas_h2o;
  
  q_gas = k_P*(P_gas - P_atm);
  if (q_gas < 0)
    q_gas = 0.0;

  
  y[41] = p_gas_h2;               
  y[42] = p_gas_ch4;               
  y[43] = p_gas_co2;              
  y[44] = 0.0;
  y[45] = 0.0;
  y[46] = 0.0;
  y[47] = 0.0;
  y[48] = 0.0;
  y[49] = 0.0;
  y[50] = 0.0;
  y[51] = 0.0;
  
  y[52] = P_gas;                /* total head space pressure from H2, CH4, CO2 and H2O */
  y[53] = q_gas * P_gas/P_atm;  /* The output gas flow is recalculated to atmospheric pressure (normalization) */
  
  
  /* Bio P equations (future use)*/  
 
  for (i = 0; i < 10; i++) {
      y[i+54] = x[36+i];
  }
  
  /* Bio S equations (future use)*/  
  
  for (i = 0; i < 10; i++) {
      y[i+64] = x[46+i];
  }
  
  /* Scat+ and San- (future use)*/  
  
  for (i = 0; i < 20; i++) {
      y[i+74] = x[56+i];
  }
  
  /* Multiple minerals(future use) */  
  
  for (i = 0; i < 28; i++) {
      y[i+94] = 0.0;
  }
  
  /* Additional outputs */  
  
  y[122] = 0.0;
  y[123] = 0.0;
  y[124] = 0.0;
  y[125] = 0.0;
  y[126] = 0.0;
  y[127] = 0.0;
  y[128] = 0.0;        
  y[129] = 0.0;
  y[130] = 0.0;    
  y[131] = 0.0;
  y[132] = 0.0;      
  y[133] = 0.0;
  y[134] = 0.0;
  y[135] = 0.0;
  y[136] = 0.0;
  y[137] = 0.0;
  y[138] = 0.0;
  y[139] = 0.0; 
  y[140] = 0.0;
  y[141] = 0.0;
  y[142] = 0.0;
  y[143] = 0.0;
  y[144] = 0.0;
  y[145] = 0.0;
  y[146] = 0.0;
  y[147] = 0.0;
  y[148] = 0.0;
  y[149] = 0.0;        
}
 

  
  



/*
 * mdlUpdate - perform action at major integration time step
 */
static void mdlUpdate(double *x, double *u, SimStruct *S, int tid)
{
}


/*
 * mdlDerivatives - compute the derivatives
 */
static void mdlDerivatives(double *dx, double *x, double *u, SimStruct *S, int tid)
{

double f_sI_xc, f_xI_xc, f_ch_xc, f_pr_xc, f_li_xc, N_xc, N_I, N_aa, C_xc, C_sI, C_ch;
double C_pr, C_li, C_xI, C_su, C_aa, f_fa_li, C_fa, f_h2_su, f_bu_su, f_pro_su, f_ac_su;
double N_bac, C_bu, C_pro, C_ac, C_bac, Y_su, f_h2_aa, f_va_aa, f_bu_aa, f_pro_aa, f_ac_aa;
double C_va, Y_aa, Y_fa, Y_c4, Y_pro, C_ch4, Y_ac, Y_h2;
double k_dis, k_hyd_ch, k_hyd_pr, k_hyd_li, K_S_IN, k_m_su, K_S_su, pH_UL_aa, pH_LL_aa;
double k_m_aa, K_S_aa, k_m_fa, K_S_fa, K_Ih2_fa, k_m_c4, K_S_c4, K_Ih2_c4, k_m_pro, K_S_pro;
double K_Ih2_pro, k_m_ac, K_S_ac, K_I_nh3, pH_UL_ac, pH_LL_ac, k_m_h2, K_S_h2, pH_UL_h2, pH_LL_h2;
double k_dec_Xsu, k_dec_Xaa, k_dec_Xfa, k_dec_Xc4, k_dec_Xpro, k_dec_Xac, k_dec_Xh2;
double R, T_base, T_op, pK_w_base, pK_a_va_base, pK_a_bu_base, pK_a_pro_base, pK_a_ac_base, pK_a_co2_base, pK_a_IN_base;
double K_w, K_a_va, K_a_bu, K_a_pro, K_a_ac, K_a_co2, K_a_IN, K_H_co2, K_H_ch4, K_H_h2;
double K_A_Bva, K_A_Bbu, K_A_Bpro, K_A_Bac, K_A_Bco2, K_A_BIN;
double P_atm, p_gas_h2o, P_gas, k_P, kLa, K_H_h2o_base, K_H_co2_base, K_H_ch4_base, K_H_h2_base, factor;
double V_liq, V_gas, t_res;
double eps, pH_op, phi, S_H_ion;

double proc1, proc2, proc3, proc4, proc5, proc6, proc7, proc8, proc9, proc10, proc11, proc12, proc13;
double proc14, proc15, proc16, proc17, proc18, proc19, procA4, procA5, procA6, procA7, procA10, procA11;
double procT8, procT9, procT10;

double I_pH_aa, I_pH_ac, I_pH_h2, I_IN_lim, I_h2_fa, I_h2_c4, I_h2_pro, I_nh3;

double reac1, reac2, reac3, reac4, reac5, reac6, reac7, reac8, reac9, reac10, reac11, reac12, reac13;
double reac14, reac15, reac16, reac17, reac18, reac19, reac20, reac21, reac22, reac23, reac24;

double N_pr;

double stoich1_c, stoich2_c, stoich3_c, stoich4_c, stoich5_c, stoich6_c, stoich7_c, stoich8_c, stoich9_c, stoich10_c, stoich11_c, stoich12_c, stoich13_c; // NEW
double stoich1_n, stoich2_n, stoich3_n, stoich4_n, stoich5_n, stoich6_n, stoich7_n, stoich8_n, stoich9_n, stoich10_n, stoich11_n, stoich12_n, stoich13_n; // NEW

double xtemp[105], inhib[6];

double S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_I, X_xc, X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, S_Na, S_K, S_Cl, S_co2, S_nh3; //NEW
double S_gas_h2, S_gas_ch4, S_gas_co2; //NEW
double p_gas_h2, p_gas_ch4, p_gas_co2, q_gas;
double pHLim_aa, pHLim_ac, pHLim_h2, a_aa, a_ac, a_h2, n_aa, n_ac, n_h2;

int i;

eps = 0.000001;

f_sI_xc = mxGetPr(PAR)[0];
f_xI_xc = mxGetPr(PAR)[1];
f_ch_xc = mxGetPr(PAR)[2];
f_pr_xc = mxGetPr(PAR)[3];
f_li_xc = mxGetPr(PAR)[4];
N_xc = mxGetPr(PAR)[5];
N_I = mxGetPr(PAR)[6];
N_aa = mxGetPr(PAR)[7];
C_xc = mxGetPr(PAR)[8];
C_sI = mxGetPr(PAR)[9];
C_ch = mxGetPr(PAR)[10];
C_pr = mxGetPr(PAR)[11];
C_li = mxGetPr(PAR)[12];
C_xI = mxGetPr(PAR)[13];
C_su = mxGetPr(PAR)[14];
C_aa = mxGetPr(PAR)[15];
f_fa_li = mxGetPr(PAR)[16];
C_fa = mxGetPr(PAR)[17];
f_h2_su = mxGetPr(PAR)[18];
f_bu_su = mxGetPr(PAR)[19];
f_pro_su = mxGetPr(PAR)[20];
f_ac_su = mxGetPr(PAR)[21];
N_bac = mxGetPr(PAR)[22];
C_bu = mxGetPr(PAR)[23];
C_pro = mxGetPr(PAR)[24];
C_ac = mxGetPr(PAR)[25];
C_bac = mxGetPr(PAR)[26];
Y_su = mxGetPr(PAR)[27];
f_h2_aa = mxGetPr(PAR)[28];
f_va_aa = mxGetPr(PAR)[29];
f_bu_aa = mxGetPr(PAR)[30];
f_pro_aa = mxGetPr(PAR)[31];
f_ac_aa = mxGetPr(PAR)[32];
C_va = mxGetPr(PAR)[33];
Y_aa = mxGetPr(PAR)[34];
Y_fa = mxGetPr(PAR)[35];
Y_c4 = mxGetPr(PAR)[36];
Y_pro = mxGetPr(PAR)[37];
C_ch4 = mxGetPr(PAR)[38];
Y_ac = mxGetPr(PAR)[39];
Y_h2 = mxGetPr(PAR)[40];
k_dis = mxGetPr(PAR)[41];
k_hyd_ch = mxGetPr(PAR)[42];
k_hyd_pr = mxGetPr(PAR)[43];
k_hyd_li = mxGetPr(PAR)[44];
K_S_IN = mxGetPr(PAR)[45];
k_m_su = mxGetPr(PAR)[46];
K_S_su = mxGetPr(PAR)[47];
pH_UL_aa = mxGetPr(PAR)[48];
pH_LL_aa = mxGetPr(PAR)[49];
k_m_aa = mxGetPr(PAR)[50];
K_S_aa = mxGetPr(PAR)[51];
k_m_fa = mxGetPr(PAR)[52];
K_S_fa = mxGetPr(PAR)[53];
K_Ih2_fa = mxGetPr(PAR)[54];
k_m_c4 = mxGetPr(PAR)[55];
K_S_c4 = mxGetPr(PAR)[56];
K_Ih2_c4 = mxGetPr(PAR)[57];
k_m_pro = mxGetPr(PAR)[58];
K_S_pro = mxGetPr(PAR)[59];
K_Ih2_pro = mxGetPr(PAR)[60];
k_m_ac = mxGetPr(PAR)[61];
K_S_ac = mxGetPr(PAR)[62];
K_I_nh3 = mxGetPr(PAR)[63];
pH_UL_ac = mxGetPr(PAR)[64];
pH_LL_ac = mxGetPr(PAR)[65];
k_m_h2 = mxGetPr(PAR)[66];
K_S_h2 = mxGetPr(PAR)[67];
pH_UL_h2 = mxGetPr(PAR)[68];
pH_LL_h2 = mxGetPr(PAR)[69];
k_dec_Xsu = mxGetPr(PAR)[70];
k_dec_Xaa = mxGetPr(PAR)[71];
k_dec_Xfa = mxGetPr(PAR)[72];
k_dec_Xc4 = mxGetPr(PAR)[73];
k_dec_Xpro = mxGetPr(PAR)[74];
k_dec_Xac = mxGetPr(PAR)[75];
k_dec_Xh2 = mxGetPr(PAR)[76];
R = mxGetPr(PAR)[77];
T_base = mxGetPr(PAR)[78];
T_op = mxGetPr(PAR)[79];
pK_w_base = mxGetPr(PAR)[80];
pK_a_va_base = mxGetPr(PAR)[81];
pK_a_bu_base = mxGetPr(PAR)[82];
pK_a_pro_base = mxGetPr(PAR)[83];
pK_a_ac_base = mxGetPr(PAR)[84];
pK_a_co2_base = mxGetPr(PAR)[85];
pK_a_IN_base = mxGetPr(PAR)[86];
K_A_Bva = mxGetPr(PAR)[87];
K_A_Bbu = mxGetPr(PAR)[88];
K_A_Bpro = mxGetPr(PAR)[89];
K_A_Bac = mxGetPr(PAR)[90];
K_A_Bco2 = mxGetPr(PAR)[91];
K_A_BIN = mxGetPr(PAR)[92];
P_atm = mxGetPr(PAR)[93];
kLa = mxGetPr(PAR)[94];
K_H_h2o_base = mxGetPr(PAR)[95];
K_H_co2_base = mxGetPr(PAR)[96];
K_H_ch4_base = mxGetPr(PAR)[97];
K_H_h2_base = mxGetPr(PAR)[98];
k_P = mxGetPr(PAR)[99];
t_res = mxGetPr(PAR)[100];

V_liq = mxGetPr(V)[0];
V_gas = mxGetPr(V)[1];


// Nitrogen extension

N_pr = 0.007;         // NEW



for (i = 0; i < 105; i++) {
   if (x[i] < 0)
     xtemp[i] = 0;
   else
     xtemp[i] = x[i];
}

S_su =  xtemp[0];
S_aa =  xtemp[1];
S_fa =  xtemp[2];
S_va =  xtemp[3];
S_bu = xtemp[4];
S_pro = xtemp[5];
S_ac =  xtemp[6];
S_h2 =  xtemp[7];
S_ch4 = xtemp[8];
S_IC = xtemp[9];
S_IN = xtemp[10];
S_I = xtemp[11];
X_xc = xtemp[12];
X_ch = xtemp[13];
X_pr = xtemp[14];
X_li = xtemp[15];
X_su = xtemp[16];
X_aa = xtemp[17];
X_fa = xtemp[18];
X_c4 = xtemp[19];
X_pro = xtemp[20];
X_ac = xtemp[21];
X_h2 =  xtemp[22];
X_I =  xtemp[23];

S_gas_h2 =  xtemp[24];
S_gas_ch4 = xtemp[25];
S_gas_co2 = xtemp[26];
M_S = xtemp[27];

S_co2 = (xtemp[9]-xtemp[63]);
S_nh3 = xtemp[64];

factor = (1.0/T_base - 1.0/T_op)/(100.0*R);
K_w = pow(10,-pK_w_base)*exp(55900.0*factor); /* T adjustment for K_w */
K_a_va = pow(10,-pK_a_va_base);
K_a_bu = pow(10,-pK_a_bu_base);
K_a_pro = pow(10,-pK_a_pro_base);
K_a_ac = pow(10,-pK_a_ac_base);
K_a_co2 = pow(10,-pK_a_co2_base)*exp(7646.0*factor); /* T adjustment for K_a_co2 */
K_a_IN = pow(10,-pK_a_IN_base)*exp(51965.0*factor); /* T adjustment for K_a_IN */

K_H_h2 = K_H_h2_base*exp(-4180.0*factor);     /* T adjustment for K_H_h2 */
K_H_ch4 = K_H_ch4_base*exp(-14240.0*factor);  /* T adjustment for K_H_ch4 */
K_H_co2 = K_H_co2_base*exp(-19410.0*factor);  /* T adjustment for K_H_co2 */
p_gas_h2o = K_H_h2o_base*exp(5290.0*(1.0/T_base - 1.0/T_op));  /* T adjustement for water vapour saturation pressure */
  
// S_H_ion = u[94];
// pH_op = -log10(u[94]);   /* pH */

p_gas_h2 = S_gas_h2*R*T_op/16.0;
p_gas_ch4 = S_gas_ch4*R*T_op/64.0;
p_gas_co2 = S_gas_co2*R*T_op;

P_gas = p_gas_h2 + p_gas_ch4 + p_gas_co2 + p_gas_h2o;

phi = xtemp[56]+(xtemp[10]-xtemp[64])-xtemp[63]-xtemp[62]/64.0-xtemp[61]/112.0-xtemp[60]/160.0-xtemp[59]/208.0-xtemp[58]; // charge balance (positive and negative charges)
S_H_ion = -phi*0.5+0.5*sqrt(phi*phi+4.0*K_w); /* SH+ */
pH_op = -log10(S_H_ion); /* pH */

/* Hill function on SH+ used within BSM2, ADM1 Workshop, Copenhagen 2005. */
pHLim_aa = pow(10,(-(pH_UL_aa + pH_LL_aa)/2.0));
pHLim_ac = pow(10,(-(pH_UL_ac + pH_LL_ac)/2.0));
pHLim_h2 = pow(10,(-(pH_UL_h2 + pH_LL_h2)/2.0));
n_aa = 3.0/(pH_UL_aa-pH_LL_aa);
n_ac = 3.0/(pH_UL_ac-pH_LL_ac);
n_h2 = 3.0/(pH_UL_h2-pH_LL_h2);
I_pH_aa = pow(pHLim_aa,n_aa)/(pow(S_H_ion,n_aa)+pow(pHLim_aa ,n_aa));
I_pH_ac = pow(pHLim_ac,n_ac)/(pow(S_H_ion,n_ac)+pow(pHLim_ac ,n_ac));
I_pH_h2 = pow(pHLim_h2,n_h2)/(pow(S_H_ion,n_h2)+pow(pHLim_h2 ,n_h2));

I_IN_lim = 1.0/(1.0+K_S_IN/S_IN);
I_h2_fa = 1.0/(1.0+S_h2/K_Ih2_fa);
I_h2_c4 = 1.0/(1.0+S_h2/K_Ih2_c4);
I_h2_pro = 1.0/(1.0+S_h2/K_Ih2_pro);
I_nh3 = 1.0/(1.0+ S_nh3/K_I_nh3);

inhib[0] = I_pH_aa*I_IN_lim;
inhib[1] = inhib[0]*I_h2_fa;
inhib[2] = inhib[0]*I_h2_c4;
inhib[3] = inhib[0]*I_h2_pro;
inhib[4] = I_pH_ac*I_IN_lim*I_nh3;
inhib[5] = I_pH_h2*I_IN_lim;

proc1 = 0.0;                                                                            // this kinetics has been modified : Xc is omitted
proc2 = k_hyd_ch*X_ch;
proc3 = k_hyd_pr*X_pr;
proc4 = k_hyd_li*X_li;
proc5 = k_m_su*S_su /(K_S_su+S_su )*X_su*inhib[0];
proc6 = k_m_aa*S_aa/(K_S_aa+S_aa)*X_aa*inhib[0];
proc7 = k_m_fa*S_fa/(K_S_fa+S_fa)*X_fa*inhib[1];
proc8 = k_m_c4*S_va/(K_S_c4+S_va)*X_c4*S_va/(S_va+S_bu+eps)*inhib[2];
proc9 = k_m_c4*S_bu/(K_S_c4+S_bu)*X_c4*S_bu/(S_va+S_bu+eps)*inhib[2];
proc10 = k_m_pro*S_pro/(K_S_pro+S_pro)*X_pro*inhib[3];
proc11 = k_m_ac*S_ac/(K_S_ac+S_ac)*X_ac*inhib[4];
proc12 = k_m_h2*(S_h2)/(K_S_h2+(S_h2))*X_h2*inhib[5];
proc13 = k_dec_Xsu*X_su;
proc14 = k_dec_Xaa*X_aa;
proc15 = k_dec_Xfa*X_fa;
proc16 = k_dec_Xc4*X_c4;
proc17 = k_dec_Xpro*X_pro;
proc18 = k_dec_Xac*X_ac;
proc19 = k_dec_Xh2*X_h2;

procA4 = K_A_Bva*(xtemp[59]*(K_a_va+S_H_ion)-K_a_va*xtemp[3]);
procA5 = K_A_Bbu*(xtemp[60]*(K_a_bu+S_H_ion)-K_a_bu*xtemp[4]);
procA6 = K_A_Bpro*(xtemp[61]*(K_a_pro+S_H_ion)-K_a_pro*xtemp[5]);
procA7 = K_A_Bac*(xtemp[62]*(K_a_ac+S_H_ion)-K_a_ac*xtemp[6]);
procA10 = K_A_Bco2*(xtemp[63]*(K_a_co2+S_H_ion)-K_a_co2*xtemp[9]);
procA11 = K_A_BIN*(xtemp[64]*(K_a_IN+S_H_ion)-K_a_IN*xtemp[10]);

procT8 =  kLa*(S_h2- 16.0*K_H_h2*p_gas_h2);
procT9 =  kLa*(S_ch4-64.0*K_H_ch4*p_gas_ch4);
procT10 = kLa*(S_co2-     K_H_co2*p_gas_co2);

stoich1_c = 0.0;																		// this stoich has been modified : Xc is omitted
stoich2_c = -C_ch+C_su;
stoich3_c = -C_pr+C_aa;
stoich4_c = -C_li+(1.0-f_fa_li)*C_su+f_fa_li*C_fa;
stoich5_c = -C_su+(1.0-Y_su)*(f_bu_su*C_bu+f_pro_su*C_pro+f_ac_su*C_ac)+Y_su*C_bac;
stoich6_c = -C_aa+(1.0-Y_aa)*(f_va_aa*C_va+f_bu_aa*C_bu+f_pro_aa*C_pro+f_ac_aa*C_ac)+Y_aa*C_bac;
stoich7_c = -C_fa+(1.0-Y_fa)*0.7*C_ac+Y_fa*C_bac;
stoich8_c = -C_va+(1.0-Y_c4)*0.54*C_pro+(1.0-Y_c4)*0.31*C_ac+Y_c4*C_bac;
stoich9_c = -C_bu+(1.0-Y_c4)*0.8*C_ac+Y_c4*C_bac;
stoich10_c = -C_pro+(1.0-Y_pro)*0.57*C_ac+Y_pro*C_bac;
stoich11_c = -C_ac+(1.0-Y_ac)*C_ch4+Y_ac*C_bac;
stoich12_c = (1.0-Y_h2)*C_ch4+Y_h2*C_bac;
stoich13_c = -C_bac+f_sI_xc*C_sI+f_ch_xc*C_ch+f_pr_xc*C_pr+f_li_xc*C_li+f_xI_xc*C_xI;	// this stoich has been modified : Xc is omitted


stoich1_n = 0.0;																		// this stoich has been modified : Xc is omitted
stoich2_n = 0.0;
stoich3_n = -N_pr+N_aa;
stoich4_n = 0.0;
stoich5_n =  +Y_su*N_bac;
stoich6_n =  -N_aa +Y_aa*N_bac;
stoich7_n =  Y_fa*N_bac;
stoich8_n =  Y_c4*N_bac;
stoich9_n =  Y_c4*N_bac;
stoich10_n = Y_pro*N_bac;
stoich11_n = Y_ac*N_bac;
stoich12_n = Y_h2*N_bac;
stoich13_n = -N_bac+f_sI_xc*N_I+f_pr_xc*N_pr+f_xI_xc*N_I;                                           // this stoich has been modified : Xc is omitted


/* Ssu */    reac1 = proc2+(1.0-f_fa_li)*proc4-proc5;
/* Saa */    reac2 = proc3-proc6;
/* Sfa */    reac3 = f_fa_li*proc4-proc7;
/* Sva */    reac4 = (1.0-Y_aa)*f_va_aa*proc6-proc8;
/* Sbu */    reac5 = (1.0-Y_su)*f_bu_su*proc5+(1.0-Y_aa)*f_bu_aa*proc6-proc9;
/* Spro */   reac6 = (1.0-Y_su)*f_pro_su*proc5+(1.0-Y_aa)*f_pro_aa*proc6+(1.0-Y_c4)*0.54*proc8-proc10;
/* Sac */    reac7 = (1.0-Y_su)*f_ac_su*proc5+(1.0-Y_aa)*f_ac_aa*proc6+(1.0-Y_fa)*0.7*proc7+(1.0-Y_c4)*0.31*proc8+(1.0-Y_c4)*0.8*proc9+(1.0-Y_pro)*0.57*proc10-proc11;
/* Sh2 */    reac8 = (1.0-Y_su)*f_h2_su*proc5+(1.0-Y_aa)*f_h2_aa*proc6+(1.0-Y_fa)*0.3*proc7+(1.0-Y_c4)*0.15*proc8+(1.0-Y_c4)*0.2*proc9+(1.0-Y_pro)*0.43*proc10-proc12-procT8;
/* Sch4 */   reac9 = (1.0-Y_ac)*proc11+(1.0-Y_h2)*proc12-procT9;
/* SIC */    reac10 = -stoich1_c*proc1-stoich2_c*proc2-stoich3_c*proc3-stoich4_c*proc4-stoich5_c*proc5- stoich6_c*proc6-stoich7_c*proc7-stoich8_c*proc8-stoich9_c*proc9-stoich10_c*proc10-stoich11_c*proc11-stoich12_c*proc12-stoich13_c*proc13-stoich13_c*proc14-stoich13_c*proc15-stoich13_c*proc16-stoich13_c*proc17-stoich13_c*proc18-stoich13_c*proc19-procT10;
/* SIN */    reac11 = -stoich1_n*proc1-stoich2_n*proc2-stoich3_n*proc3-stoich4_n*proc4-stoich5_n*proc5- stoich6_n*proc6-stoich7_n*proc7-stoich8_n*proc8-stoich9_n*proc9-stoich10_n*proc10-stoich11_n*proc11-stoich12_n*proc12-stoich13_n*proc13-stoich13_n*proc14-stoich13_n*proc15-stoich13_n*proc16-stoich13_n*proc17-stoich13_n*proc18-stoich13_n*proc19;
/* SI */     reac12 = f_sI_xc*(proc13 + proc14 + proc15 + proc16 + proc17 + proc18 + proc19);
/* XC */     reac13 = 0.0;                                                                          // this reac has been modified :direct mapping from biomass decay amd not longer accounted
/* Xch */    reac14 = f_ch_xc*(proc13 + proc14 + proc15 + proc16 + proc17 + proc18 + proc19) -proc2; // this reac has been modified: direct mapping from biomass decay
/* Xpro */   reac15 = f_pr_xc*(proc13 + proc14 + proc15 + proc16 + proc17 + proc18 + proc19) -proc3; // thus reac has been modified :direct mapping from biomass decay
/* Xli */    reac16 = f_li_xc*(proc13 + proc14 + proc15 + proc16 + proc17 + proc18 + proc19)- proc4; // this reac has been modified :direct mapping from biomass decay
/* Xsu */    reac17 = Y_su*proc5-proc13;
/* Xaa */    reac18 = Y_aa*proc6-proc14;
/* Xfa */    reac19 = Y_fa*proc7-proc15;
/* Xc4 */    reac20 = Y_c4*proc8+Y_c4*proc9-proc16;
/* Xpro */   reac21 = Y_pro*proc10-proc17;
/* Xac */    reac22 = Y_ac*proc11-proc18;
/* Xh2 */    reac23 = Y_h2*proc12-proc19;
/* XI */     reac24 = f_xI_xc*(proc13 + proc14 + proc15 + proc16 + proc17 + proc18 + proc19);        // this reac has been modified :direct mapping from biomass decay

q_gas = k_P*(P_gas-P_atm);
if (q_gas < 0)
   q_gas = 0.0;
   
dx[0] = 1.0/V_liq*(u[24]*u[0])-(x[0]/(t_res+V_liq/u[24]))+reac1;
dx[1] = 1.0/V_liq*(u[24]*u[1])-(x[1]/(t_res+V_liq/u[24]))+reac2;
dx[2] = 1.0/V_liq*(u[24]*u[2])-(x[2]/(t_res+V_liq/u[24]))+reac3;
dx[3] = 1.0/V_liq*(u[24]*u[3])-(x[3]/(t_res+V_liq/u[24]))+reac4;
dx[4] = 1.0/V_liq*(u[24]*u[4])-(x[4]/(t_res+V_liq/u[24]))+reac5;
dx[5] = 1.0/V_liq*(u[24]*u[5])-(x[5]/(t_res+V_liq/u[24]))+reac6; 
dx[6] = 1.0/V_liq*(u[24]*u[6])-(x[6]/(t_res+V_liq/u[24]))+reac7;

/* calculated in Sh2solv.c */
dx[7] = 1.0/V_liq*(u[24]*u[7])-(x[7]/(t_res+V_liq/u[24]))+reac8;                                 

dx[8]  = 1.0/V_liq*(u[24]*u[8])-(x[8]/(t_res+V_liq/u[24]))+reac9;
dx[9]  = 1.0/V_liq*(u[24]*u[9])-(x[9]/(t_res+V_liq/u[24]))+reac10;    /* SIC */
dx[10] = 1.0/V_liq*(u[24]*u[10])-(x[10]/(t_res+V_liq/u[24]))+reac11; /* SIN */
dx[11] = 1.0/V_liq*(u[24]*u[11])-(x[11]/(t_res+V_liq/u[24]))+reac12;
dx[12] = 1.0/V_liq*(u[24]*u[12])-(x[12]/(t_res+V_liq/u[24]))+reac13;
dx[13] = 1.0/V_liq*(u[24]*u[13])-(x[13]/(t_res+V_liq/u[24]))+reac14;
dx[14] = 1.0/V_liq*(u[24]*u[14])-(x[14]/(t_res+V_liq/u[24]))+reac15;
dx[15] = 1.0/V_liq*(u[24]*u[15])-(x[15]/(t_res+V_liq/u[24]))+reac16;
dx[16] = 1.0/V_liq*(u[24]*u[16])-(x[16]/(t_res+V_liq/u[24]))+reac17;
dx[17] = 1.0/V_liq*(u[24]*u[17])-(x[17]/(t_res+V_liq/u[24]))+reac18;
dx[18] = 1.0/V_liq*(u[24]*u[18])-(x[18]/(t_res+V_liq/u[24]))+reac19;
dx[19] = 1.0/V_liq*(u[24]*u[19])-(x[19]/(t_res+V_liq/u[24]))+reac20;
dx[20] = 1.0/V_liq*(u[24]*u[20])-(x[20]/(t_res+V_liq/u[24]))+reac21;
dx[21] = 1.0/V_liq*(u[24]*u[21])-(x[21]/(t_res+V_liq/u[24]))+reac22;
dx[22] = 1.0/V_liq*(u[24]*u[22])-(x[22]/(t_res+V_liq/u[24]))+reac23;
dx[23] = 1.0/V_liq*(u[24]*u[23])-(x[23]/(t_res+V_liq/u[24]))+reac24;

/* MASS TRANSFER EQUATIONS */

dx[24] = -S_gas_h2* q_gas/V_gas+procT8*V_liq/V_gas;        /* Sh2 */
dx[25] = -S_gas_ch4*q_gas/V_gas+procT9*V_liq/V_gas;        /* Sch4 */
dx[26] = -S_gas_co2*q_gas/V_gas+procT10*V_liq/V_gas;       /* Sco2 */
dx[27] = 1.0/V_liq*(u[24]*u[27])-(x[27]/(t_res+V_liq/u[27]));  /* D2 - Now using for Mineral Solids*/ 
dx[28] = 0.0;  /* D3 */
dx[29] = 0.0;  /* D4 */
dx[30] = 0.0;  /* D5 */
dx[31] = 0.0;  /* D6 */
dx[32] = 0.0;  /* D7 */
dx[33] = 0.0;  /* D8 */
dx[34] = 0.0;  /* D9 */
dx[35] = 0.0;  /* D10 */

/* Bio P reactions (reac25) */ 

dx[36] = 1.0/V_liq*(u[24]*(u[26]-x[36])); /* D1 */
dx[37] = 0.0; /* D2 */
dx[38] = 0.0; /* D3 */
dx[39] = 0.0; /* D4 */
dx[40] = 0.0; /* D5 */
dx[41] = 0.0; /* D6 */
dx[42] = 0.0; /* D7 */
dx[43] = 0.0; /* D8 */
dx[44] = 0.0; /* D9 */
dx[45] = 0.0; /* D10 */

/* Bio S reactions (reac35) */

dx[46] = 1.0/V_liq*(u[24]*(u[36]-x[46]));  /* D1 */
dx[47] = 0.0;  /* D2 */ 
dx[48] = 0.0;  /* D3 */
dx[49] = 0.0;  /* D4 */
dx[50] = 0.0;  /* D5 */
dx[51] = 0.0;  /* D6 */ 
dx[52] = 0.0;  /* D7 */ 
dx[53] = 0.0;  /* D8 */
dx[54] = 0.0;  /* D9 */
dx[55] = 0.0;  /* D10 */

/* Scations(+) and Sanions(-) */

dx[56] = 1.0/V_liq*(u[24]*(u[46]-x[56]));  /* Na */ 
dx[57] = 1.0/V_liq*(u[24]*(u[47]-x[57]));  /* K */ 
dx[58] = 1.0/V_liq*(u[24]*(u[48]-x[58]));  /* Cl */

/* weak acidbase chemistry */

dx[59] = -procA4;  /* Sva- */
dx[60] = -procA5;  /* Sbu- */
dx[61] = -procA6;  /* Spro- */
dx[62] = -procA7;  /* Sac- */
dx[63] = -procA10; /* SHCO3- */
dx[64] = -procA11; /* SNH3 */

dx[65] = 0.0;  /* D10 */
dx[66] = 0.0;  /* D11 */
dx[67] = 0.0;  /* D12 */ 
dx[68] = 0.0;  /* D13 */
dx[69] = 0.0;  /* D14 */
dx[70] = 0.0;  /* D15 */
dx[71] = 0.0;  /* D16 */ 
dx[72] = 0.0;  /* D17 */ 
dx[73] = 0.0;  /* D18 */
dx[74] = 0.0;  /* D19 */
dx[75] = 0.0;  /* D20 */

/* Multiple minerals */

dx[76] = 0.0;  /* D1 */
dx[77] = 0.0;  /* D2 */ 
dx[78] = 0.0;  /* D3 */
dx[79] = 0.0;  /* D4 */
dx[80] = 0.0;  /* D5 */
dx[81] = 0.0;  /* D6 */ 
dx[82] = 0.0;  /* D7 */ 
dx[83] = 0.0;  /* D8 */
dx[84] = 0.0;  /* D9 */
dx[85] = 0.0;  /* D10 */
dx[86] = 0.0;  /* D11 */
dx[87] = 0.0;  /* D12 */ 
dx[88] = 0.0;  /* D13 */
dx[89] = 0.0;  /* D14 */
dx[90] = 0.0;  /* D15 */
dx[91] = 0.0;  /* D16 */ 
dx[92] = 0.0;  /* D17 */ 
dx[93] = 0.0;  /* D18 */
dx[94] = 0.0;  /* D19 */
dx[95] = 0.0;  /* D20 */
dx[96] = 0.0;  /* D21 */
dx[97] = 0.0;  /* D22 */
dx[98] = 0.0;  /* D23 */
dx[99] = 0.0;  /* D24 */
dx[100] = 0.0;  /* D25 */ 
dx[101] = 0.0;  /* D26 */ 
dx[102] = 0.0;  /* D27 */
dx[103] = 0.0;  /* D28 */
dx[104] = 0.0;  /* D29 */

}


/*
 * mdlTerminate - called when the simulation is terminated.
 */
static void mdlTerminate(SimStruct *S)
{
}

#ifdef	MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif

