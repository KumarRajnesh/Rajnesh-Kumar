clear all
clc
% Topic C: Study on Optimization of Gas Turbine

% Data
gamma = 1.4;
R = 8.314; % J/(mol.K)
M_air = (0.79*28+0.21*32)/1000;  % kg/mol
r = R/M_air;  % J/(kg.K)
Cp = gamma*r/(gamma-1);
Cv = Cp-r;
rho = 10;

% Initialize Matrices
T_2_values = [];
T_3_values = [];
W_comp_values = [];
W_turb_values = [];
W_turb_s_values = [];
Ratio_values = [];
rho_values = [];
eta_reg_values = [];



for rho = 1:0.1:20

rho_values = [rho_values, rho];

% Turbine
eta_is = 0.85;   % isentropic efficiency of turbine

% Combustion Chamber
LHV = 42.3e6; % J/kg

% Point 1
T_1 = 25 + 273.15; % K
P_1 = 1; % bar

% Point 2
T_2 = T_1*(rho)^((gamma-1)/gamma);
T_2_values = [T_2_values, T_2];
P_2 = rho*P_1;

% Point 3
T_3 = T_2 + LHV/(17*Cp);
T_3_values = [T_3_values, T_3];
P_3 = P_2;

% Point 4
T_4s = T_3*(rho)^((1-gamma)/gamma);
T_4 = eta_is*(T_4s - T_3) + T_3;

%Compressor Work
W_comp = Cp*(T_2 - T_1);
W_comp_values=[W_comp_values,W_comp];

%Isentropic work of turbine
W_turb_s = abs(Cp*(T_4s - T_3));
W_turb_s_values = [W_turb_s_values, W_turb_s];

%Real Work of turbine
W_turb = abs(Cp*(T_4-T_3));
W_turb_values = [W_turb_values, W_turb];

% Work Ratio
Ratio = W_comp/W_turb;
Ratio_values = [Ratio_values, Ratio];

% Net Work
Net_Work = W_turb_values - W_comp_values;

% Heat Input
q_in = Cp*(T_3-T_2);

% Thermal Efficiency
eta_th_simple = Net_Work/q_in;

end

% Optimum Values
Net_Work_opt = max(Net_Work')
[row, col] = find(Net_Work' == Net_Work_opt);
W_turb_opt = W_turb_values(col,row)
W_comp_opt = W_comp_values(col,row)
rho_opt = rho_values(col,row)
T_2_opt = T_2_values(col,row);
T_3_opt = T_3_values(col,row);
T_4_opt = T_3_opt - W_turb_opt/Cp;
T_4s_opt = (T_4_opt - T_3_opt)/eta_is + T_3_opt;
q_in_opt = Cp*(T_3_opt - T_2_opt);
eta_th_simple_opt = Net_Work_opt/q_in_opt;
Ratio_opt = W_comp_opt/W_turb_opt;

% Cycle with Regenerator
for eta_reg = 0:0.1:1

eta_reg_values = [eta_reg_values, eta_reg];

    T_5 = eta_reg_values*(T_4_opt-T_2_opt) + T_2_opt;

    q_in_reg = Cp*(T_5 - T_2_opt);

    % new heat input after installing regeneator
    q_in_53 = Cp*(T_3_opt-T_5);

    % new Thermal Efficiency
    eta_th_reg = Net_Work_opt./q_in_53;

end

% Characteristics at eta_reg = 0.7
eta_reg = 0.7;
T_5_a = eta_reg*(T_4_opt-T_2_opt) + T_2_opt;

q_in_reg_a = Cp*(T_5_a - T_2_opt);

% new heat input after installing regeneator
q_in_53_a = Cp*(T_3_opt-T_5_a);

% new Thermal Efficiency
eta_th_reg_a = Net_Work_opt./q_in_53_a;


%Inilialize Matices
P_int_values = [];
W_comp_LP_values = [];
W_comp_HP_values = [];
rho_stage_values = [];


% Staged Compression
for P_int = 1:0.1:rho_opt*P_1
    P_int_values = [P_int_values, P_int];

% Point 1
P_1_stage = P_1; %bar
T_1_stage = T_1; %K

% Point 2
P_2_stage = P_int; % bar
T_2_stage = T_1_stage*(P_2_stage./P_1_stage)^((gamma-1)/gamma);  %K

% Point 3
P_3_stage = P_int; % bar
T_3_stage = T_1_stage;

% Point 4
P_4_stage = rho_opt*P_1_stage;
T_4_stage = T_3_stage*(P_3_stage./P_4_stage)^((1-gamma)/gamma);

% Point 5
T_5_stage = T_3_opt;

% Point 6
T_6_stage = T_4_opt;

% Point 7
T_7_stage = eta_reg*(T_6_stage-T_4_stage) + T_4_stage;

%Compressor # 1 Work
W_comp_LP = Cp*(T_2_stage - T_1_stage);
W_comp_LP_values = [W_comp_LP_values, W_comp_LP];

%Compressor # 2 Work
W_comp_HP = Cp*(T_4_stage - T_3_stage);
W_comp_HP_values = [W_comp_HP_values, W_comp_HP];

%Total Compressor Work
W_comp_total = W_comp_LP_values + W_comp_HP_values;

%Net work
Net_Work_stage = W_turb_opt - W_comp_total;

%Ratio
Ratio_stage = W_comp_total./W_turb_opt;

%Heat provided to combustion chamber
q_in_75 = Cp*(T_5_stage - T_7_stage);


end

%Optimum Values
P_int_opt=sqrt(P_1_stage*P_4_stage);

T_2_stage_opt = T_1_stage*(P_int_opt./P_1_stage)^((gamma-1)/gamma);

T_4_stage_opt = T_3_stage*(P_int_opt./P_4_stage)^((1-gamma)/gamma);

W_comp_total_opt = Cp*(T_2_stage_opt - T_1_stage) + Cp*(T_4_stage_opt - T_3_stage);

Net_Work_stage_opt = W_turb_opt - W_comp_total_opt;

Ratio_stage_opt = W_comp_total_opt/W_turb_opt;

q_in_75_opt = Cp*(T_5_stage - T_7_stage);

eta_th_stage_opt = Net_Work_stage_opt./q_in_75_opt;



figure(1)
plot(rho_values, W_turb_s_values);
hold on
plot(rho_values, W_turb_values);
xlabel("Pressure Ratio")
ylabel("Turbine Work Output (J/kg)")
legend("Isentropic", "Non-Isentropic");

figure(2)  %3.3.1
plot(rho_values,Ratio_values);
xlabel("Pressure Ratio")
ylabel("Work Ratio")

figure(3) % 3.3.2
plot(rho_values, W_comp_values);
hold on
plot(rho_values, W_turb_values);
xlabel("Pressure Ratio")
ylabel("Work (J/kg)")
legend("Wcomp", "Wturb");

figure(4)  %3.4
plot(rho_values, Net_Work);
xlabel("Pressure Ratio")
ylabel("Net Work (J/kg)")

figure(5) %3.4
plot(rho_values, eta_th_simple);
xlabel("Pressure Ratio")
ylabel("Thermal Efficiency (%)")

figure(6) % 4.1.2
plot(eta_reg_values, T_5);
xlabel("eta,reg")
ylabel("T_5 (K)")

figure(7) %4.1.3
plot(eta_reg_values,q_in_53);
xlabel("eta,reg")
ylabel("q,in,reg (J/kg)")

figure(8) %4.1.3
plot(eta_reg_values, eta_th_reg);
xlabel("eta,reg")
ylabel("Thermal Efficiency with Regenerator (%)")

figure(9) %2.4.2
plot(P_int_values, W_comp_total);
xlabel("Intermediate Pressure (bars)")
ylabel("Total Compressor Work (J/kg)")

figure(10) %2.4.2
plot(P_int_values, Ratio_stage);
xlabel("Intermediate Pressure (bar)")
ylabel("Work Ratio")




