%Case Based Module: PEM Electrolyzer Parameters (Electrolyzer Pressure, Current Density, Membrane Thickness) Optimization to acheive Optimum Levelized Cost of Hydrogen (Euros/kg).

clear all
clc

T = 60+273.15; %[K] -- Electrolyzer operating temperature

Henry_H2= ((0.255*10^5*exp(170/T))*(101325/1000000))^-1;
Henry_O2= ((1.33*10^6*exp(-666/T))*(101325/1000000))^-1;

D_m_H2 = 13.1e-9;
D_m_O2 = 5.7e-9;  

lambda_m = 20; % [mol of H2O/mol of SO3]

F = 96500; %[C/mol] -- Faraday Constant
R = 8.314; %[J/(mol.K)] -- General Gas Constant
gamma = 1.4; % Specific heat capacity ratio
HHV_H2 = 285800; %[J/mol] -- Higher Heating value of Hydrogen
M_H2 = 2.016e-3; % [kg/mol] -- Molar mass of hydrogen H2
M_O2 = 32.0e-3;  % [kg/mol] -- Molar mass of Oxygen O2

b = 0.0227; % Tafel Slope
I_0 = 0.0062e4; % [A/m2] Standard Current Density                  

R_el = 20e-7; % [Ohms.m2] -- Electrical Resistance
P_0 = 1e5; % [Pa] -- Atmospheric pressure
PO2_a = 2e5; % [Pa] -- Oxygen Pressure at Anode

Efficiency_comp = 0.8; % Compressor efficiency


Capex_el = 1.5; %[€/W] -- Per watt cost of electrolyzer
Celopex = 3/100; % Percentage of Capex that goes to operating cost of electrolyzer per year
LTel = 79000; % [hours] -- life time of electrolyzer

Capex_comp= 3.8; %[€/W] -- Per watt cost of compressor
Ccompopex = 3/100; % Percentage of Capex that goes to operating cost of compressor per year   
LTcomp = 79000;% [hours] -- life time of compressor


%Parametric Study

PH2_del = [500]*1e5; % [Pa] -- Delivery Pressure of Hydrogen
Electricity_price = [0.12];  % [€/kWh]
t_m = [265]*1e-6; % [m]
PH2_c = [1:1:500]*1e5; % [Pa] -- Hydrogen Pressure at Cathode  
I = [0.05:0.1:0.7]*1e4; % [A/m2] -- Current Density


for i_Ep = 1:length(Electricity_price)
for i_tm = 1:length(t_m)
for i_PH2_c = 1:length(PH2_c)
    for i_I = 1:length(I)

        % Difussion
        N_a_H2(i_I,i_PH2_c) = Henry_H2*D_m_H2*PH2_c(i_PH2_c)./(R*T*t_m(i_tm)); %[mol/(m2.s)] -- Diffusion cathode to anode
        N_c_O2(i_I,i_PH2_c) = Henry_O2*D_m_O2*PO2_a./(R*T*t_m(i_tm)); %[mol/(m2.s)] -- Diffusion anode to cathode

        % Hydrogen Produced at cathode
        N_prod_H2(i_I,i_PH2_c) = I(i_I)/(2*F)-N_a_H2(i_I,i_PH2_c)-2*N_c_O2(i_I,i_PH2_c);

        % Anodic Ratio: Percent of hydogen present anode
        ANODICRATIO(i_I,i_PH2_c) = 100*N_a_H2(i_I,i_PH2_c)/(I(i_I)/(4*F)-N_c_O2(i_I,i_PH2_c));

        if ANODICRATIO(i_I,i_PH2_c) < 2

            % The reversible potential in Volts
            U_P0_rev = 1.4;
            % The Nernst potential in Volts
            U_N(i_I,i_PH2_c) = U_P0_rev + (R*T)/(2*F)*log((PH2_c(i_PH2_c)*sqrt(PO2_a))/(P_0^(3/2)));
            % Conductivity of the membrane in "1/(Ohm.m)"
            Sigma_m = (0.519*lambda_m - 0.326)*exp(1263*(1/303 - 1/T));
            % Protonic resistance of membrane in Ohm.m2
            R_m = t_m(i_tm)/Sigma_m;
            % The Ohmic overvoltage in Volts
            U_ohm(i_I,i_PH2_c) = (R_m+R_el)*I(i_I);
            % cell voltage of electrolyzer in Volts
            U_cell(i_I,i_PH2_c) = U_N(i_I,i_PH2_c)+ U_ohm(i_I,i_PH2_c) +b*log(I(i_I)/I_0);


            % Power consumed by electrolyzer in W/m2
            Power_el(i_I,i_PH2_c) = U_cell(i_I,i_PH2_c).*I(i_I);

            % Electrolyser Efficiency
            Electrolyser_efficiency(i_I,i_PH2_c) = N_prod_H2(i_I,i_PH2_c)*HHV_H2./Power_el(i_I,i_PH2_c);

            % Power consumed by compressor in W/m2
            Power_comp(i_I,i_PH2_c) = ((gamma/(gamma-1))*R*T*((PH2_del./PH2_c(i_PH2_c))^((gamma-1)/gamma)-1))/Efficiency_comp*N_prod_H2(i_I,i_PH2_c);

            % Total Power Consumption in W/m2
            Power_total(i_I,i_PH2_c) = Power_el(i_I,i_PH2_c) + Power_comp(i_I,i_PH2_c);

            % Overall Efficiency of the System
            Efficiency_overall(i_I,i_PH2_c) = N_prod_H2(i_I,i_PH2_c)*HHV_H2./Power_total(i_I,i_PH2_c);

            
                                                        % Cost Analysis

            % Electrolyzer costing in €/m2

            % Electrolyzer Capex 
            Electrolyzer_Capex(i_I,i_PH2_c) = Capex_el*Power_el(i_I,i_PH2_c);
            % Electrolyser OPEX excluding electricity
            Electrolyzer_Opex(i_I,i_PH2_c)=Celopex*Electrolyzer_Capex(i_I,i_PH2_c)*LTel/(365*24);
            % Cost of electricity for the electrolyser
            Power_cost_el(i_I,i_PH2_c)=Power_el(i_I,i_PH2_c)*LTel/1000*Electricity_price(i_Ep);

            % Compressor Costing in €/m2

            % Compressor Capex 
            Compressor_Capex(i_I,i_PH2_c) = Capex_comp*Power_comp(i_I,i_PH2_c)*LTel/LTcomp;
            % Compressor Opex excluding electriciy
            Compressor_Opex(i_I,i_PH2_c)=Ccompopex*Compressor_Capex(i_I,i_PH2_c)*LTel/(365*24);
            % Cost of electricity for the compressor
            Power_cost_comp(i_I,i_PH2_c) =Power_comp(i_I,i_PH2_c)*LTel/1000*Electricity_price(i_Ep);

            % Hydrogen Production

            % Hydrogen Production over the lifetime of electrolyzer in kg/m2
            Total_H2_produced(i_I,i_PH2_c) = N_prod_H2(i_I,i_PH2_c)*M_H2*3600*LTel;


            % Combined Capex in €/kg
            H2_cost_cap(i_I,i_PH2_c)=(Compressor_Capex(i_I,i_PH2_c)+Electrolyzer_Capex(i_I,i_PH2_c))./Total_H2_produced(i_I,i_PH2_c);
            % Combined Opex excluding electiricty in €/kg
            H2_cost_op(i_I,i_PH2_c)=(Compressor_Opex(i_I,i_PH2_c)+Electrolyzer_Opex(i_I,i_PH2_c))./Total_H2_produced(i_I,i_PH2_c);
            % Combined cost of electricity in €/kg
            H2_cost_el(i_I,i_PH2_c)=(Power_cost_comp(i_I,i_PH2_c)+Power_cost_el(i_I,i_PH2_c))./Total_H2_produced(i_I,i_PH2_c);


            % Electricity cost in €/kg
            Power_cost_comp_f(i_I,i_PH2_c) = Power_cost_comp(i_I,i_PH2_c)./Total_H2_produced(i_I,i_PH2_c);
            Power_cost_el_f(i_I,i_PH2_c) = Power_cost_el(i_I,i_PH2_c)./Total_H2_produced(i_I,i_PH2_c);
            % Opex excluding electricity in €/kg
            Compressor_Opex_f(i_I,i_PH2_c)=Compressor_Opex(i_I,i_PH2_c)./Total_H2_produced(i_I,i_PH2_c);
            Electrolyzer_Opex_f(i_I,i_PH2_c)=Electrolyzer_Opex(i_I,i_PH2_c)./Total_H2_produced(i_I,i_PH2_c);
            % Capex in €/kg
            Compressor_Capex_f(i_I,i_PH2_c)=Compressor_Capex(i_I,i_PH2_c)./Total_H2_produced(i_I,i_PH2_c);
            Electrolyzer_Capex_f(i_I,i_PH2_c)=Electrolyzer_Capex(i_I,i_PH2_c)./Total_H2_produced(i_I,i_PH2_c);


            % Levelized cost of Hydrogen in €/kg
            LCOH(i_I,i_PH2_c)=(Electrolyzer_Capex(i_I,i_PH2_c)+Compressor_Capex(i_I,i_PH2_c)+Electrolyzer_Opex(i_I,i_PH2_c)+Compressor_Opex(i_I,i_PH2_c)+Power_cost_el(i_I,i_PH2_c)+Power_cost_comp(i_I,i_PH2_c))./Total_H2_produced(i_I,i_PH2_c);
            
            %LCOH(i_I,i_PH2_c)=(Electrolyzer_Capex(i_I,i_PH2_c)+Electrolyzer_Opex(i_I,i_PH2_c)+Power_cost_el(i_I,i_PH2_c))./Total_H2_produced(i_I,i_PH2_c);



        else
            N_a_H2(i_I,i_PH2_c) = NaN;
            N_c_O2(i_I,i_PH2_c) = NaN;
            N_prod_H2(i_I,i_PH2_c) = NaN;
            ANODICRATIO(i_I,i_PH2_c) = NaN;
            U_N(i_I,i_PH2_c) = NaN;
            U_ohm(i_I,i_PH2_c) = NaN;
            U_cell(i_I,i_PH2_c) = NaN;
            Power_el(i_I,i_PH2_c) = NaN;
            Electrolyser_efficiency(i_I,i_PH2_c) = NaN;
            Power_comp(i_I,i_PH2_c) = NaN;
            Power_total(i_I,i_PH2_c) = NaN;
            Efficiency_overall(i_I,i_PH2_c) = NaN;

            Electrolyzer_Capex(i_I,i_PH2_c) = NaN;
            Electrolyzer_Opex(i_I,i_PH2_c) = NaN;
            Power_cost_el(i_I,i_PH2_c) = NaN;

            Compressor_Capex(i_I,i_PH2_c) = NaN;
            Compressor_Opex(i_I,i_PH2_c) = NaN;
            Power_cost_comp(i_I,i_PH2_c) = NaN;

            Total_H2_produced(i_I,i_PH2_c) = NaN;

            H2_cost_cap(i_I,i_PH2_c) = NaN;
            H2_cost_op(i_I,i_PH2_c) = NaN;
            H2_cost_el(i_I,i_PH2_c) = NaN;

            Power_cost_comp_f(i_I,i_PH2_c) = NaN;
            Power_cost_el_f(i_I,i_PH2_c) = NaN;
            Compressor_Opex_f(i_I,i_PH2_c) = NaN;
            Electrolyzer_Opex_f(i_I,i_PH2_c) = NaN;
            Compressor_Capex_f(i_I,i_PH2_c) = NaN;
            Electrolyzer_Capex_f(i_I,i_PH2_c) = NaN;

            LCOH(i_I,i_PH2_c) = NaN;
       
        end
    end
    % Optimum Parameters Calculation

    % Optimum levelized cost of Hydrogen
    LCOH_opt(i_PH2_c) = min(LCOH(:,i_PH2_c));

    % Cell Address of Optimum Parameters
    [row,col] = find(LCOH(:,i_PH2_c) == LCOH_opt(i_PH2_c));

    % Optimum Current Density
    I_opt(i_PH2_c) = I(1,row);

    % Optimum Hydrogen Production
    Total_H2_produced_opt(i_PH2_c) = Total_H2_produced(row,col);

    % Optimum Costs in €/kg
    Power_cost_comp_f_opt(i_PH2_c) = Power_cost_comp_f(row,i_PH2_c);
    Power_cost_el_f_opt(i_PH2_c) = Power_cost_el_f(row,i_PH2_c);
    Compressor_Opex_f_opt(i_PH2_c) = Compressor_Opex_f(row,i_PH2_c);
    Compressor_Capex_f_opt(i_PH2_c) = Compressor_Capex_f(row,i_PH2_c);
    Electrolyzer_Opex_f_opt(i_PH2_c) = Electrolyzer_Opex_f(row,i_PH2_c);
    Electrolyzer_Capex_f_opt(i_PH2_c) = Electrolyzer_Capex_f(row,i_PH2_c);

    % Optimum Electrolyzer Efficiency
    Electrolyser_efficiency_opt(i_PH2_c) = Electrolyser_efficiency(row,i_PH2_c);

    % Optimum Overall Efficiency
    Efficiency_overall_opt(i_PH2_c) = Efficiency_overall(row, i_PH2_c);



end

% Optimum Parameters for each membrane
LCOH_OPT(i_tm) = min(LCOH_opt);
[row1,col1] = find(LCOH_opt == LCOH_OPT(i_tm));

I_OPT(i_tm) = I_opt(row1, col1);
PH2_c_OPT(i_tm) = PH2_c(row1, col1);

end

%Optimum Cost based on electricity price

LCOH_OPT_EP(i_Ep) = min(LCOH_OPT(i_tm));

[row2,col2] = find(LCOH_OPT == LCOH_OPT_EP(i_Ep));

I_OPT_Ep(i_Ep) = I_OPT(row2, col2);
PH2_c_OPT_Ep(i_Ep) = PH2_c_OPT(row2, col2);

end

figure(1)
plot(PH2_c./1e5, LCOH,'LineWidth', 2)
xlabel("Electrolyzer Pressure [Bar]")
ylabel("LCOH [€/kg]")
legend("I = 0.05 A/cm^2", "I = 0.15 A/cm^2", "I = 0.25 A/cm^2", "I = 0.35 A/cm^2", "I = 0.45 A/cm^2", "I = 0.55 A/cm^2", "I = 0.65 A/cm^2");
%legend("t_m = 60 μm", "t_m = 147 μm", "t_m = 185 μm", "t_m = 203 μm", "t_m = 220 μm", "t_m = 237 μm", "t_m = 265 μm");
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
hold on


% figure(3)
% plot(PH2_del./1e5, LCOH_OPT_Del, 'LineWidth', 2, 'LineStyle','-')
% xlabel("Electrolyzer Pressure [Bar]")
% ylabel("Optimum LCOH [€/kg]")
% legend("Direct Compression", "15 bar + Compression")
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman')
% hold on
% 
% figure(3)
% plot(Electricity_price, LCOH_OPT_EP, 'LineWidth', 2, 'LineStyle','--')
% xlabel("Electricity Price [€/kWh]")
% ylabel("Optimum LCOH [€/kg]")
% legend("30 bar (Direct)","30 bar (Compressor)", "62 bar (Direct)", "62 bar (Compressor)", "200 bar (Direct)", "200 bar (Compressor)", "350 bar (Direct)", "350 bar (Compressor)")
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman')
% hold on

% figure(1)
% plot(PH2_c./1e5,LCOH_opt,'LineWidth', 2)
% hold on
% plot(PH2_c./1e5,Power_cost_el_f_opt,'LineWidth', 2)
% plot(PH2_c./1e5,Power_cost_comp_f_opt,'LineWidth', 2)
% plot(PH2_c./1e5,Compressor_Opex_f_opt,'LineWidth', 2)
% plot(PH2_c./1e5,Compressor_Capex_f_opt,'LineWidth', 2)
% plot(PH2_c./1e5,Electrolyzer_Opex_f_opt,'LineWidth', 2)
% plot(PH2_c./1e5,Electrolyzer_Capex_f_opt,'LineWidth', 2)
% xlabel("Electrolyzer Pressure [Bar]")
% ylabel("Optimum Costs (€/kg)")
% legend("LCOH", "Electrolyzer Electricity", "Compressor Electricity", "Compressor Opex", "Compressor Capex", "Electrolyzer Opex", "Electrolyzer Capex");
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
% grid
% hold on

% figure(2)
% plot(I./1e4,H2_cost_op,'LineWidth', 2)
% xlabel("Current Density [A/cm^2]")
% ylabel("H2 Opex (€/kg)")
% legend("PH2_c = 1 bar", "PH2_c = 100 bar", "PH2_c = 200 bar", "PH2_c = 300 bar", "PH2_c = 400 bar", "PH2_c = 500 bar");
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');


% figure(1)
% plot(I./1e4, ANODICRATIO,'LineWidth', 2)
% xlabel("Current Density [A/cm^2]")
% ylabel("ANODIC RATIO (%)")
% legend("t_m = 60 μm", "t_m = 147 μm", "t_m = 185 μm", "t_m = 203 μm", "t_m = 220 μm", "t_m = 237 μm", "t_m = 265 μm");
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
% hold on


% figure(2)
% subplot(1,2,1)
% hold on
% plot(I./1e4, Efficiency_overall, 'LineWidth', 2)
% xlabel("Current Density [A/cm^2]")
% ylabel("Overall Efficiency (%)")
% legend("t_m = 60 μm", "t_m = 147 μm", "t_m = 185 μm", "t_m = 203 μm", "t_m = 220 μm", "t_m = 237 μm", "t_m = 265 μm");
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
% 
% subplot(1,2,2)
% hold on
% plot(I./1e4, Electrolyser_efficiency, 'LineWidth', 2)
% xlabel("Current Density [A/cm^2]")
% ylabel("Electrolyzer Efficiency (%)")
% legend("t_m = 60 μm", "t_m = 147 μm", "t_m = 185 μm", "t_m = 203 μm", "t_m = 220 μm", "t_m = 237 μm", "t_m = 265 μm");
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
% hold on

% figure(4)
% plot(PH2_c./1e5, U_cell, 'LineWidth', 2)
% xlabel("Electrolyzer Pressure [Bar]")
% ylabel("Cell Voltage (V)")
% legend("t_m = 60 μm", "t_m = 147 μm", "t_m = 185 μm", "t_m = 203 μm", "t_m = 220 μm", "t_m = 237 μm", "t_m = 265 μm");
% hold on

% figure(2)
% plot(I./1e4, N_a_H2)
% xlabel("Current Density [A/cm^2]")
% ylabel("Hydrogen Diffusion [mol/(m^2.s)]")
% 
% figure(3)
% plot(PH2_c./1e5, N_prod_H2)
% xlabel("Electrolyzer Pressure [Bars]")
% ylabel("Hydrogen Production [mol/(m^2.s)]")
% 
% figure(4)
% 
% 
% figure(6)
% plot(PH2_c./1e5, I_opt./1e4);
% xlabel("Electrolyzer Pressure [Bars]")
% ylabel("Current Density [A/cm^2]")
% 
% figure(7)
% plot(PH2_c./1e5, LCOH_opt);
% hold on
% plot(PH2_c./1e5, Power_cost_comp_f_opt);
% plot(PH2_c./1e5, Power_cost_el_f_opt);
% plot(PH2_c./1e5, Compressor_Opex_f_opt);
% plot(PH2_c./1e5, Compressor_Capex_f_opt);
% plot(PH2_c./1e5, Electrolyzer_Opex_f_opt);
% plot(PH2_c./1e5, Electrolyzer_Capex_f_opt);
% xlabel("Electrolyzer Pressure [Bars]")
% ylabel("Cost [E/kg]")
% legend("LCOH", "Compressor Electricity", "Electrolyzer Electricity", "Comp Opex","Comp Capex", "Electrolyzer Opex", "Electrolyzer Capex")
% 
% figure(8)
% subplot(1,2,1)
% hold on
% plot(PH2_c./1e5, Electrolyser_efficiency_opt-0.02);
% xlabel("Pressure")
% ylabel("Electrolyzer Efficiency")
% 
% subplot(1,2,2)
% plot(PH2_c./1e5, Efficiency_overall_opt);
% xlabel("Pressure")
% ylabel("Overall Efficiency")
% 
% yyaxis left
% plot(PH2_c./1e5, Power_el, 'LineWidth', 2);
% ylabel("Electrolyzer Power (W/m^2)")
% yyaxis right
% plot(PH2_c./1e5, Power_comp, 'LineWidth', 2);
% ylabel("Compressor Power (W/m^2)")
% xlabel("Electrolyzer Pressure [Bars]")
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');


% plot(PH2_c./1e5, PH2_del./PH2_c, 'LineWidth', 2);
% ylabel("Compressor Pressure Ratio")
% xlabel("Electrolyzer Pressure [Bars]")
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');


