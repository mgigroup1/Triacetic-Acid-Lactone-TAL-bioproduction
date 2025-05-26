function [t,F,c,v]=TAL_model_v1; 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% TAL model                                                           %%% 
%%%                                                                     %%% 
%%% CF                                                                  %%% 
%%% V1: version of 11/09/2023                                           %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
clear all 
clc 
clf 
warning off
global Mwx Mws Mwg Mwco2 Mwn Mwo Mwph Mwa Mwc Mwt Mwm Mww % global molecular weights
global Cnf Csf Cgf Cmf Caf % global feed percent mass composition
global Frac % global feed carbon source percentage
global qgmax Ks1 qamax Ksa Ki1_g Ki1_s Ki1_m qsmax Ks2 qmmax Ksm qcupmax Ks6 Ki3_g Ki3_a Ki3_s Ki3_m msmax r_evap % global uptake, cell maintenance, evaporation rate kinetics
global qtmax_g Ks4_g Kiqx_g qtmax_a Ks4_a Kiqx_a qtmax_s Ks4_s Kiqx_s qtmax_m Ks4_m Kiqx_m % global TAL production kinetics
global qcmax_g Ks5_g Ki2_g Kiqx_g qcmax_a Ks5_a Ki2_a Kiqx_a qcmax_s Ks5_s Ki2_s Kiqx_s qcmax_m Ks5_m Ki2_m Kiqx_m % global citrate production kinetics
global alpha1 beta1 alpha2 beta2 zeta kd qbm % global cell growth/death kinetics
global dG1 dG2 dG3 dG4 dG5 dG6 dG7 dG8 dG9 dG10 dG11 % global energy deamand
global PS FB SG % global substrate and fed-batch switches

format long
Tend  = 288;                                                               % Duration of fermentation (h)

% =========== INPUTS ====================================================== 

%User-defined inputs
% PS: primary substrate
% FB: Specify batch or fed-batch
PS = 1;                                                                    %Primary substrate: 1 = glucose, 2 = acetic acid, 3 = ethanol, 4 = methanol
FB = 0;
Frac = [1, 0, 0, 0];


% =========== FIXED & CALCULATED CONSTANTS ================================ 
  
% Fixed constants:
Mwx = 24.6;                                                                % Molweight biomass (g/mol) 
Mws = 46;                                                                  % Molweight ethanol (g/mol) 
Mwg = 180;                                                                 % Molweight glucose (g/mol)
Mwco2 = 44;                                                                % Molweight CO2 (g/mol)
Mwn = 17;                                                                  % Molweight NH3 (g/mol) 
Mwo = 32;                                                                  % Molweight O2 (g/mol) 
Mwph = 96;                                                                 % Molweight pyrophosphate ion (g/mol)
Mwa = 60;                                                                  % Molweight acetic acid; (g/mol)
Mwc = 192;                                                                 % Molweight citrate; (g/mol)
Mwt = 126;                                                                 % Molweight TAL; (g/mol)
Mwm = 32;                                                                  % Molweight methanol; (g/mol)
Mww = 18;                                                                  % Molweight water; (g/mol)

Cnf = 0.29;                                                                % Mass fraction of NH3 in ammonium feed
Csf = 1;                                                                   % Mass fraction of ethanol in ethanol feed
Cgf = 0.75;                                                                % Mass fraction of glucose in glucose feed
Cmf = 1;                                                                   % Mass fraction of methanol in methanol feed
Caf = 1;                                                                   % Mass fraction of acetic acid in acetic acid feed

% Kinetics: 

%Glucose uptake
qgmax = 0.11;
Ks1 = 0.1;
Ks1 = 1;


%Acetic acid uptake
qamax = 0.111904350402778;
Ksa = 1.1;
Ki1_g = 1;
Ki1_s = 1;
Ki1_m = 1;


%Ethanol uptake
qsmax = 0.0608402689216667;
Ks2 = 1;


%Methanol uptake
qmmax = 0.121341209873889;
Ksm = 1;


%Citrate uptake
qcupmax = 0.0068;
Ks6 =0.25;
Ki3_g = 10;
Ki3_a = 6.66666666666667;
Ki3_s = 5.11445611111111;
Ki3_m = 10.67075;


%Cell maintenance
msmax = 0.001;


%Citrate production
%from glucose
qcmax_g = 0.06;
Ks5_g = 2.3;                                                               % MM-constant ethanol for citrate formation (g ethanol/kg broth)
Ki2_g = 20;
Kiqx_g = 0.01;


%from acetic acid
qcmax_a = 0.046;
Ks5_a = 1.533333333333333;
Ki2_a = 30;
Kiqx_a = 0.01;


%from ethanol
qcmax_s = 0.030;
Ks5_s = 1.17632490555556; 
Ki2_s = 38.3584208333333;
Kiqx_s = 0.01;


%from methanol
qcmax_m = 0.06;
Ks5_m = 2.4542725;
Ki2_m = 24.542725;
Kiqx_m = 0.01;


%TAL production
%from glucose
qtmax_g = 0.017;
Ks4_g = 1;
Kiqx_g = 0.01;


%from acetic acid
qtmax_a = 0.0168466666666667;
Ks4_a = 1;
Kiqx_a = 0.01;


%from ethanol
qtmax_s = 0.0171;
Ks4_s = 1;
Kiqx_s = 0.01;


%from methanol
qtmax_m = 0.0176;
Ks4_m = 1;
Kiqx_m = 0.01;


%Biomass formation
alpha1 = 0.6373;
beta1 = -15.19;
alpha2 = 0.03543;
beta2 = -1.739;
zeta = 30;
kd = -0.002; %cell death rate


%Evaporation rates
r_evap = 0.4;


% ATP demands
dG1 = 6.16;                                                                % Energy generation for glucose upatke (mol ATP/mol glucose)
dG2 = 2.08;                                                                % Energy generation for ethanol upatke (mol ATP/mol ethanol)
dG3 = 10.2;                                                                % Energy generation by acetic acid catabolism (mol ATP/mol acetic acid)
dG4 = 1.5;                                                                 % Energy consumption by biomass synthesis (mol ATP/mol biomass)
dG5 = 2;                                                                   % Energy consumption by TAL synthesis (mol ATP/mol TAL)
dG6 = 3.5;                                                                 % Energy consumption by citrate synthesis (mol ATP/mol citrate)
dG7 = 1;                                                                   % Energy consumption by acetic acid uptake (mol ATP/mol acetic acid)
dG8 = 15;                                                                  % Energy consumption by methanol uptake (mol ATP/mol methanol)
dG9 = 3.5;                                                                 % Energy generation by citrate uptake (mol ATP/mol citrate)
dG10 = 1;                                                                  % Energy consumption by acetic acid uptake (mol ATP/mol acetic acid)
dG11 = 0.6666666666;                                                       % Energy generation by methanol uptake (mol ATP/mol methanol)


%Initial conditions
SG = 0;
qbm = 0;


WB0 = 0.00045;                                                             % Initial broth weight

Cx0 = 3;                                                                   % Initial biomass concentration (g/kg)
if PS == 1 && FB == 0
    Cg0 = 180;
elseif PS ~= 1 && FB == 0
    Cg0 = 0;
else
    Cg0 = 80;
end
if PS == 2 && FB == 0
    Caint0 = 183.11620975;                                                 %Batch acetic acid for model fitting of acetic acid parameters
elseif PS ~= 2 && FB == 0
    Caint0 = 9.9;                                                          %Batch acetic acid for model fitting of other substrate (i.e. parameters
else
    Caint0 = 0;
end
if PS == 3 && FB == 0
    Cs0 = 99.55680369;
else
    Cs0 = 0;                                                               % Initial ethanol concentration (g/kg) 
end
if PS == 4 && FB == 0
    Cm0 = 192.0735;
else
    Cm0 = 0;
end


Ca0 = 0;
Ct0 = 0;                                                                   % Initial retinoid concentration (g/kg)
Cc0 = 0;                                                                   % Initial citrate concentration (g/kg)
Cco20 = 0;                                                                 % Initial CO2 concentration (g/kg)
Cn0 = 0.67;                                                                % Initial NH3 concentration (g/kg, NBS3) 
Cz0 = 1e1;                                                                 % Initial Z concentration (mmol/kg, arbitrary value)
Cph0 = 16;                                                                 % Initial P concentration (mmol/kg, arbitrary value for now)



Mx0 = Cx0*WB0;                                                             % Initial biomass weight (kg) 
Ms0 = Cs0*WB0;                                                             % Initial ethanol weight (kg)
Mg0 = Cg0*WB0;                                                             % Initial glucose weight (kg)
Ma0 = Ca0*WB0;                                                             % Initial acetate weight (kg)
Mm0 = Cm0*WB0;                                                             % Initial methanol weight (kg)
Mt0 = Ct0*WB0;                                                             % Initial TAL weight (kg)
Mc0 = Cc0*WB0;                                                             % Initial citrate weight (kg)
Mn0 = Cn0*WB0;                                                             % Initial NH3 weight (kg) 
Mo0 = 0;                                                                   % Initial O2 amount (kg) 
Mco20 = Cco20*WB0;                                                         % Initial CO2 amount (kg) 
Mz0 = Cz0*WB0;                                                             % Initial Z amount (mol) 
Mph0 = Cph0*WB0;                                                           % Initial P amount (kg) 
Mph0 = 0;
Mw0 = 0;                                                                   % Initial water lost to evaporation (kg)
Tn0 = 0;                                                                   % Totalized NH3 dosed (kg) 
Ts0 = Ms0;                                                                 % Totalized ethanol dosed (kg) 
To0 = 0;                                                                   % Totalized O2 consumed (kg) 
Tco20 = Cco20*WB0;                                                         % Totalized CO2 produced (kg) 
Maint0 = Caint0*WB0;
M_Afed0 = Maint0;
M_Mfed0 = Mm0;
Tw0 = 0;                                                                   %Totalized water produced (kg)
Tg0 = Mg0;

[t,F] = ode45(@diffeq,[0:0.1:Tend],[WB0;Mx0;Ms0;Mg0;Mt0;Mc0;Mn0;Mo0;Mco20;Mz0;Tn0;Ts0;To0;Tco20;Mph0;Mw0;WB0;Maint0;Mm0;M_Afed0;M_Mfed0;Tw0;Tg0]); 

SG = 0;
qbm = 0;
count=0; 
c=[]; 
v=[]; 
for n=t' 
    count=count+1; 
    [c_,v_]=kineq(t(count),F(count,:)); 
    c = [c;c_]; 
    v = [v;v_]; 
end

Min = WB0*1000 + F(end,13) + F(end,11)/Cnf + ((F(end,23)-Mg0)/Cgf) + (F(end,20) - M_Afed0)/Caf + (F(end,21) - M_Mfed0)/Cmf + (F(end,12) - Ms0)/Csf;
Mout = F(end,1)*1000 + F(end,9) + F(end,16);

Mass_balance = 100*(Min - Mout)/Min;

Cin = 1*Mx0/Mwx + 6*F(end,23)/Mwg + 6*Mt0/Mwt + 6*Mc0/Mwc + 2*F(end,20)/Mwa + 1*F(end,21)/Mwm + 2*F(end,12)/Mws;
Cout = 1*F(end,2)/Mwx + 6*F(end,4)/Mwg + 6*F(end,5)/Mwt + 6*F(end,6)/Mwc + 1*F(end,9)/Mwco2 + 2*F(end,18)/Mwa + 1*F(end,19)/Mwm + 2*F(end,3)/Mws;
Carbon_balance = 100*(Cin - Cout)/Cin;

figure(1);
subplot(3,4,1)
plot(t,c(:,8),'-k','LineWidth',3)
title('glucose concentration (g/L)')
xlabel('Fermentation age (h)')
ylabel('Concentration (g/L)')
set(gca,'FontSize',12)

subplot(3,4,2)
plot(t,c(:,4),'-k','LineWidth',3)
title('citrate concentration (g/L)')
xlabel('Fermentation age (h)')
ylabel('Concentration (g/L)')
set(gca,'FontSize',12)

subplot(3,4,3)
if PS == 2 && FB == 0
    plot(t,c(:,9),'-k','LineWidth',3)
else
    plot(t,c(:,9),'-k','LineWidth',3)
end
title('Acetate concentration (g/L)')
xlabel('Fermentation age (h)')
ylabel('Concentration (g/L)')
set(gca,'FontSize',12)

subplot(3,4,4)
plot(t,c(:,3),'-k','LineWidth',3)
title('TAL concentration (g/L)')
xlabel('Fermentation age (h)')
ylabel('Concentration (g/L)')
set(gca,'FontSize',12)

subplot(3,4,5)
plot(t,c(:,1),'-k','LineWidth',3)
title('Biomass concentration (g/L)')
xlabel('Fermentation age (h)')
ylabel('Concentration (g/L)')
set(gca,'FontSize',12)

subplot(3,4,6)
plot(t,F(:,1)*1000,'-k','LineWidth',3)
title('Broth weight (kg)')
xlabel('Fermentation age (h)')
ylabel('Mass (kg)')
set(gca,'FontSize',12)

subplot(3,4,7)
plot(t,v(:,7)./Mwo.*c(:,1)*1000,'-k','LineWidth',3)
title('OUR (mmol/kg/h)')
xlabel('Fermentation age (h)')
ylabel('OUR (mmol/kg/h)')
set(gca,'FontSize',12)

subplot(3,4,8)
plot(t,v(:,8)./Mwco2.*c(:,1)*1000,'-k','LineWidth',3)
title('CER (mmol/kg/h)')
xlabel('Fermentation age (h)')
ylabel('CER (mmol/kg/h)')
set(gca,'FontSize',12)

subplot(3,4,9)
plot(t,(v(:,8)./Mwco2)./(v(:,7)./Mwo),'-k','LineWidth',2) 
title('Respiration quotient')
xlabel('Fermentation age (h)') 
ylabel('RQ (-)')
set(gca,'FontSize',12)

subplot(3,4,10)
plot(t,c(:,2),'-k','LineWidth',3)
title('Ethanol concentration (g/L)')
xlabel('Fermentation age (h)')
ylabel('Concentration (g/L)')
set(gca,'FontSize',12)

subplot(3,4,11)
plot(t,c(:,10),'-k','LineWidth',3)
title('Methanol concentration (g/L)')
xlabel('Fermentation age (h)')
ylabel('Concentration (g/L)')
set(gca,'FontSize',12)

subplot(3,4,12)
plot(t,c(:,5),'-k','LineWidth',3)
title('Ammonium concentration (g/L)')
xlabel('Fermentation age (h)')
ylabel('Concentration (g/L)')
set(gca,'FontSize',12)


    placeholder = 1;

end

function f=diffeq(time,f);
global Mwx Mws Mwg Mwco2 Mwn Mwo Mwph Mwa Mwc Mwt Mwm %global molecular weights
global Cnf Csf Cgf Cmf Caf %global feed percent mass composition
global Frac_a Frac_s Frac_m %global feed carbon source percentage
global qgmax Ks1 qamax Ksa Ki1_g Ki1_s Ki1_m qsmax Ks2 qmmax Ksm qcupmax Ks6 Ki3_g Ki3_a Ki3_s Ki3_m msmax r_evap %global uptake, cell maintenance, evaporation rate kinetics
global qtmax_g Ks4_g Kiqx_g qtmax_a Ks4_a Kiqx_a qtmax_s Ks4_s Kiqx_s qtmax_m Ks4_m Kiqx_m %global TAL production kinetics
global qcmax_g Ks5_g Ki2_g Kiqx_g qcmax_a Ks5_a Ki2_a Kiqx_a qcmax_s Ks5_s Ki2_s Kiqx_s qcmax_m Ks5_m Ki2_m Kiqx_m %global citrate production kinetics
global alpha1 beta1 alpha2 beta2 zeta kd qbm %global cell growth/death kinetics
global dG1 dG2 dG3 dG4 dG5 dG6 dG7 dG8 dG9 dG10 dG11 %global energy deamand
global FB SG
  
% State variables: 
Mx =   f(2);                                                               % Biomass amount in broth (kg LDW) 

% Reaction rates: 
[dummy,v] = kineq(time,f); 
clear dummy 
qx = v(1);                                                                 % Specific growth rate (1/h) 
qs = v(2);                                                                 % Specific ethanol uptake rate (g ethanol/g LDW/h) 
ms = v(3);                                                                 % Maintenance (g acetic acid/g LDW/h) 
qt = v(4);                                                                 % Specific TAL production rate (g TAL/g LDW/h) 
qc = v(5);                                                                 % Specific citrate production rate (g citrate/g LDW/h) 
qn = v(6);                                                                 % Specific NH3 consumption rate (g NH3/g LDW/h) 
qo = v(7);                                                                 % Specific O2 consumption rate (mol O2/g LDW/h) 
qco2 = v(8);                                                                 % Specific CO2 production rate (mol CO2/g LDW/h) 
qz = v(9);                                                                 % Specific Z consumption rate (mol Z/g LDW/h) 
Fs = v(10);                                                                 % Ethanol feed rate (kg ethanol/h) 
Fn = v(11);                                                                % NH3 feed rate (g NH3/h) 
Fg = v(12);
qph = v(13);                                                                % Specific P consumption rate (mol P/kg LDW/h)
qglc = v(14);                                                              % Specific glucose consumption rate
qaup = v(17);
qcup = v(18);
qm = v(19);
Fa = v(20);
Fm = v(21);
qH2O = v(22);

% Differentials: 
dWBdt = (Fg/Cgf + Fs/Csf + Fn/Cnf + Fa/Caf + Fm/Cmf + qo*Mx - qco2*Mx - r_evap/1000)/1000;     % Change in broth weight (ton broth/h) 
%dWBdt = (qo*Mx - qco2*Mx - r_evap/1000)/1000;
dMxdt = qx*Mx;                                                             % Change in biomass (kg LDW/h) 
dMsdt = Fs-qs*Mx;                                                          % Change in ethanol (kg ethanol/h)
dMgdt = Fg-qglc*Mx;                                                        % Change in glucose (kg glucose/h)
dMtdt = qt*Mx;                                                             % Change in TAL (kg TAL/h) 
dMcdt = qc*Mx;                                                             % Change in citrate (kg citrate/h)
dMndt = Fn-qn*Mx;                                                          % Change in NH3 (kg NH3/h) 
dModt = -qo*Mx;                                                            % Change in O2 (kg O2/h) 
dMco2dt = qco2*Mx;                                                           % Change in CO2 (kg CO2/h) 
dMzdt = -qz*Mx*1000;                                                       % Change in Z (mol Z/h) 
dTndt = Fn;                                                                % Change in totalized NH3 feed (kg NH3/h) 
dTsdt = Fs;                                                                % Change in totalized ethanol feed (kg ethanol/h) 
dTodt = qo*Mx;                                                             % Change in totalized O2 (kg O2/h) 
dTco2dt = qco2*Mx;                                                           % Change in totalized CO2 (kg CO2/h)
dMphdt = qph*Mx;                                                            % Change in P (kg P/h) 
dMwdt = r_evap/1000;                                                       % Water lost to evaporation (kg H2O/h)
dTWBdt = (Fg/Cgf + Fs/Csf + Fn/Cnf + Fa/Csf + Fm/Cmf + qo*Mx - qc*Mx - r_evap/1000)/1000;                  % Change in broth weight (ton broth/h) 
dAAdt = Fa - qaup*Mx;
dMdt = Fm - qm*Mx;
dTAdt = Fa;
dTMdt = Fm;
dTwdt = qH2O*Mx;
dTgdt = Fg;
f=[dWBdt;dMxdt;dMsdt;dMgdt;dMtdt;dMcdt;dMndt;dModt;dMco2dt;dMzdt;dTndt;dTsdt;dTodt;dTco2dt;dMphdt;dMwdt;dTWBdt;dAAdt;dMdt;dTAdt;dTMdt;dTwdt;dTgdt]; 


end



% Kinetic rate expressions: 
function [c,v]=kineq(time,f); 

global Mwx Mws Mwg Mwco2 Mwn Mwo Mwph Mwa Mwc Mwt Mwm Mww %global molecular weights
global Cnf Csf Cgf Cmf Caf %global feed percent mass composition
global Frac %global feed carbon source percentage
global qgmax Ks1 qamax Ksa Ki1_g Ki1_s Ki1_m qsmax Ks2 qmmax Ksm qcupmax Ks6 Ki3_g Ki3_a Ki3_s Ki3_m msmax r_evap %global uptake, cell maintenance, evaporation rate kinetics
global qtmax_g Ks4_g Kiqx_g qtmax_a Ks4_a Kiqx_a qtmax_s Ks4_s Kiqx_s qtmax_m Ks4_m Kiqx_m %global TAL production kinetics
global qcmax_g Ks5_g Ki2_g Kiqx_g qcmax_a Ks5_a Ki2_a Kiqx_a qcmax_s Ks5_s Ki2_s Kiqx_s qcmax_m Ks5_m Ki2_m Kiqx_m %global citrate production kinetics
global alpha1 beta1 alpha2 beta2 zeta kd qbm %global cell growth/death kinetics
global dG1 dG2 dG3 dG4 dG5 dG6 dG7 dG8 dG9 dG10 dG11 %global energy deamand
global PS FB SG

% State variables: 
WB =   f(1);                                                               % Broth weight (ton) 
Mx =   f(2);                                                               % Biomass in broth (kg LDW) 
Ms =   f(3);                                                               % Ethanol in broth (kg glucose) 
Mg = f(4);                                                                 % Glucose in broth (kg glucose), do not allow to decrease below zero
Mt =   f(5);                                                               % TAL in broth (kg retinoids) 
Mc =  f(6);                                                                % Citrate in broth (kg citrate)
Mn =   f(7);                                                               % NH3 in broth (kg NH3) 
Mz =   f(10);                                                              % Z in broth (mol Z) 
Mph =   f(15);                                                             % P in broth (mol P) 
Ma =  f(18);                                                               % Acetic acid in broth (kg acetic acid)
Mm = f(19);                                                                % Methanol in broth (kg methanol)

Cx =   Mx/WB;                                                              % Biomass concentration in broth (g/kg) 
Cs =   Ms/WB;                                                              % Ethanol concentration in broth (g/kg) 
Cg =   Mg/WB;                                                              % Glucose concentration in broth (g/kg) 
Ct =   Mt/WB;                                                              % TAL concentration in broth (g/kg) 
Cc = Mc/WB;                                                                % Citrate concentration in broth (g/kg)
Cn =   Mn/WB;                                                              % NH3 concentration in broth (g/kg) 
Cz =   Mz/WB*1e-3;                                                         % Z concentration in broth (mol/kg) 
Cph =   Mph/WB;                                                            % P concentration in broth (mol/kg) 
Cm = Mm/WB;                                                                % methanol concentration in broth (mol/kg)

if PS ~= 2 && FB == 0
    if time < 36
        Caint = 0;
    else
        Caint = Ma/WB;
    end
else
    Caint = Ma/WB;
end

% Reaction rates:                                   

%Uptake rates
qg = qgmax*Cg/(Cg+Ks1);                                                    % Specific glucose uptake rate (g glucose/gDW/h) 
qaup = qamax*Caint/(Caint + Ksa)*Ki1_g/(Ki1_g + Cg)*Ki1_s/(Ki1_s + Cs)*Ki1_m/(Ki1_m + Cm); %Acetate uptake used in model fitting
qs = qsmax*Cs/(Cs+Ks2);                                                    % Specific ethanol uptake rate (g ethanol/gDW/h) 
qm = qmmax*Cm/(Cm + Ksm);                                                  %Specific methanol uptake rate(g methanol/gDW/h)
qcup = qcupmax*Cc/(Ks6 + Cc)*Ki3_g/(Ki3_g + Cg)*Ki3_s/(Ki3_s + Cs)*Ki3_a/(Ki3_a + Caint)*Ki3_m/(Ki3_m + Cm);

% Maintinence rate kJ/gDW/h, dependent on ethanol only to avoid numerator going to zero when glucose is depleted
ma = msmax; % Maintinence rate mol ATP/gDW/h

%Citrate production rate
qc_g = qcmax_g*Cg/((Ks5_g + Cg))*Ki2_g/(Ki2_g + Cg)*Kiqx_g/(Kiqx_g + qbm);
if PS == 2 && FB == 0
    qc_a = qcmax_a*Caint/((Ks5_a + Caint))*Ki2_a/(Ki2_a + Caint)*Kiqx_a/(Kiqx_a + qbm);
elseif FB == 1
    qc_a = qcmax_a*Caint/((Ks5_a + Caint))*Ki2_a/(Ki2_a + Caint)*Kiqx_a/(Kiqx_a + qbm);
else
    qc_a = 0;
end

qc_s = qcmax_s*Cs/((Ks5_s + Cs))*Ki2_s/(Ki2_s + Cs)*Kiqx_s/(Kiqx_s + qbm);
qc_m = qcmax_m*Cm/((Ks5_m + Cm))*Ki2_m/(Ki2_m + Cm)*Kiqx_m/(Kiqx_m + qbm);
qc = qc_g + qc_a + qc_s + qc_m;

%TAL production rate
qt_g = qtmax_g*Cg/(Ks4_g + Cg)*Kiqx_g/(Kiqx_g + qbm);
if PS == 2 && FB == 0
    qt_a = qtmax_a*Caint/(Ks4_a + Caint)*Kiqx_a/(Kiqx_a + qbm);
elseif FB == 1
    qt_a = qtmax_a*Caint/(Ks4_a + Caint)*Kiqx_a/(Kiqx_a + qbm);
else
    qt_a = 0;
end
qt_s = qtmax_s*Cs/(Ks4_s + Cs)*Kiqx_s/(Kiqx_s + qbm);
qt_m = qtmax_m*Cm/(Ks4_m + Cm)*Kiqx_m/(Kiqx_m + qbm);
qt = qt_g + qt_a + qt_s + qt_m;

% Cell growth/death

if FB == 1
    if Cg < 1 && SG == 0
        SG = 1;
    end
    if SG == 0
        qbm = (alpha1*exp(beta1*qc) + alpha2*exp(beta2*qc))/zeta;
        Fs = 0;
        Fa = 0;
        Fm = 0;
        Fg = 0;
    else
        qbm = 0;
        Fg = Frac(1)*0.000187;                                             % Glucose feed rate (kg glucose/h)
        Fa = Frac(2)*0.000163;                                              % Acetic acid feed rate (kg acetic acid/h)
        Fs = Frac(3)*0.000105;                                               % Ethanol feed rate (kg ethanol/h)
        Fm = Frac(4)*0.000198;                           
    end
end

if FB == 0
    Fs = 0;
    Fa = 0;
    Fm = 0;
    Fg = 0;
    if time < 71
        qbm = (alpha1*exp(beta1*qc) + alpha2*exp(beta2*qc))/zeta;
    else
        qbm = 0;
    end
end
qd = kd;
qx = qbm + qd;

% Maintenance energy and carbon catabolism by mass and energy balance
ma = qg*dG1/Mwg + qs*dG2/Mws + qm*dG11/Mwm + qcup*dG9/Mwc + dG3*(1/3*qm/Mwm + 2*qg/Mwg + 1*qs/Mws + 1*qaup/Mwa - 3*qt/Mwt - 2.25*qc/Mwc) - qt*dG5/Mwt - qc*dG6/Mwc - qaup*dG10/Mwa - qx*((dG4 + 0.525*dG3)/Mwx);
qcat = 1/3*qm/Mwm + 2*qg/Mwg + 1*qs/Mws  + 1*qaup/Mwa - 0.525*qx/Mwx - 3*qt/Mwt - 2.25*qc/Mwc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Network stoichiometry
% Glucose uptake
% 1 C6H12O6 + 2 O2 -> 2 C2H4O2 + 2 CO2 + 2 H2O + dG1
%
% Ethanol uptake
% 1 C2H6O + 1 O2 -> 1 C2H4O2 + 1 H2O + dG2
% 
% Catabolism: 
% 1 C2H4O2 + 2 O2 -> 2 CO2 + 2 H2O + dG3 
%
% Growth reaction:  
% 0.525 C2H4O2 + 0.2 NH3 + 0.001 Z + 0.01 P + dG4 -> 1
% C1H1.8O0.5N0.2Z0.001P0.01 + 0.45 H2O + 0.05 CO2
% 


% Stoichiometric matrix, row = metabolite, column = reaction
S = [-1 0 0 0 0 0 0 0 0; ...                                               % Glucose
    0 -1 0 0 0 0 0 0 0; ...                                                % Ethanol
    0 0 -0.2 0 0 0 0 0 0; ...                                              % Ammonium
    -2 -1 0.02 -2 0 0 -0.5 0 -5/6; ...                           % Oxygen
    2 1 -0.525 -1 -3 -2.25 0 1 1/3; ...                        % Acetic acid
    0 0 1 0 0 0 0 0 0; ...                                                 % Biomass
    2 0 0.05 2 0 -1.5 2 0 1/3; ...                             % Carbon dioxide
    0 0 -0.001 0 0 0 0 0 0; ...                                            % Pseudometabolite Z
    0 0 -0.01 0 0 0 0 0 0; ...                                             % Phosphate
    0 0 0 0 1 0 2/3 0 0; ...                                    % TAL
    0 0 0 0 0 1 -1 0 0; ...                                                % citrate
    0 0 0 0 0 0 0 -1 0; ...                                                %extracellular acetic acid
    0 0 0 0 0 0 0  0 -1; ...                                               %methanol
    2 1 0.45 2 3 0.5 2 0 4/3];                               %H2O


flx = [qg/Mwg;qs/Mws;qx/Mwx;qcat;qt/Mwt;qc/Mwc;qcup/Mwc;qaup/Mwa;qm/Mwm];  % Flux vector converted to mol/gDW/h
qvec = S*flx;                                                              % dC/dt  (mol/gDW/h)

qglc = -qvec(1)*Mwg;                                                       % g glucose0/gDW/h
qetoh = -qvec(2)*Mws;                                                      % g ethanol/gDW/h
qn = -qvec(3)*Mwn;                                                         % g NH3/gDW/h
qo = -qvec(4)*Mwo;                                                         % g O2/gDW/h
qa = qvec(5);                                                              % mole acetic acid/kgDW/h
qco2 = qvec(7)*Mwco2;                                                      % g CO2/gDW/h
qz = -qvec(8);                                                             % mole Z/gDW/h
qph = -qvec(9)*Mwph;                                                       % g phosphate/gDW/h
qt = qvec(10)*Mwt;                                                         % g TAL/gDW/h
qc = qvec(11)*Mwc;                                                         % g citrate/gDW/h
qaup = -qvec(12)*Mwa;                                                      %g acetic acid/gDW/h
qm = -qvec(13)*Mwm;                                                        % g methanol/gDW/h
Fn = qn*Mx;                                                                % kg NH3/h
qH2O = qvec(14)*Mww;                                                       % g H2O/gDW/h



c=[Cx,Cs,Ct,Cc,Cn,Cz,Cph,Cg,Caint,Cm];                                              % Concentration vector 
v=[qx,qetoh,ma,qt,qc,qn,qo,qco2,qz,Fs,Fn,Fg,qph,qglc,qcat,qa,qaup,qcup,qm,Fa,Fm,qH2O];                % Rate vector (qa included as a Cmol balance check)  
  
end
