% Design Parameters

%Load Power in watts

P_Load = 60;

% For PV Cell
AP = 2 ; % Area of PV cells in m2
npv = 0.20; % PV cell Efficiency
Pvoltage = 12; %PV cell nominal voltage

% For Fuel Cell
Nc= 99; %No. of single cells in stack
Vsco=0.82; %Open circuit voltage
asc=0.18;  %Ohmic drop
Vswd=0.7;  %Diode-switch voltage drop
Vfc=12.6; %equal to V_battery @ 0.6 SOC + 0.7V (diode voltage drop)
LHV = 120000; %in KJ/Kg
nfc=0.48; %For PEM
m = 0.002; % molar mass of H2 in kg
v = 22.4; % molar volume of a gas at STP

%Inverter
ninv= 0.95; %Inverter Efficiency

% For Battery
Vbattery= 12; %Nominal Battery voltage in Volts
Vbo=11; %Open circuit voltage
ab=1.5; % No-load voltage drop
nb= 0.95; %Battery efficiency
Battery_Power_Req = P_Load/(nb*ninv);
Qmax=150;

%Initializing Matrices
FinalSolarpowerNancyS2 = zeros(8761,10);
Power_produced_by_panel = zeros(8761,1);
Solar_Current_To_Battery = zeros(8761,1);
Fuel_Cell_Current = zeros(8761,1);
Battery_SOC = zeros(8761,1);
Battery_Q = zeros(8761,1);
Battery_Voltage= zeros(8761,1);
Current_To_Load = zeros(8761,1);
Fuel_Cell_Power = zeros(8761,1);
Hydrogen_Consumption = zeros(8761,1);
% Battery SOC must be between 30% and 70%

%Fetchinhg Excel Sheet Data

Time=xlsread('Solarpower-Nancy.xlsx', 'Data&Results','B2:B8762' );
Solar_Irr=xlsread('Solarpower-Nancy.xlsx', 'Data&Results', 'C2:C8762' );


%Assigning Initial Values to Matrices
% Initially, battery SOC is 0.6 on 1st Jan 0:00 AM
Battery_SOC(1,1)= 0.6;
Battery_Q(1,1)= 150.*Battery_SOC(1,1);

% Inserting Equations 

for k=1:1:8760
n =k
Power_produced_by_panel(n,1) = npv.*AP.*Solar_Irr(n,1)
Solar_Current_To_Battery(n,1) = Power_produced_by_panel(n,1)./Pvoltage
Battery_Voltage(n,1)=Vbo+ab.*Battery_SOC(n,1)
Current_To_Load(n,1) = Battery_Power_Req./Battery_Voltage(n,1)
if Battery_SOC(n,1)<0.6
    Fuel_Cell_Current(n,1)=3.845;
else
    Fuel_Cell_Current(n,1)=0;
end

%Fuel Cell Power in KWatt
Fuel_Cell_Power(n,1)=(Vfc.*Fuel_Cell_Current(n,1))./1000 

if ((Battery_SOC(n,1))+((Power_produced_by_panel(n,1)+Fuel_Cell_Power(n,1).*1000-Battery_Power_Req)./(Battery_Voltage(n,1).*Qmax))) > 0.7
    
    Battery_SOC(n+1,1)=0.7
else
    Battery_SOC(n+1,1)=((Battery_SOC(n,1))+((Power_produced_by_panel(n,1)+Fuel_Cell_Power(n,1).*1000-Battery_Power_Req)./(Battery_Voltage(n,1).*Qmax)))
end 

%Hydrogen Consumption in Litres/hr
Hydrogen_Consumption(n,1)=v.*(3600.*Fuel_Cell_Power(n,1)./(LHV*nfc))./m;

%Battery_SOC(n+1,1)=Battery_Q(n+1,1)./Qmax;
end

FinalSolarpowerNancyS2(:,1)=Time;
FinalSolarpowerNancyS2(:,2)=Solar_Irr;
FinalSolarpowerNancyS2(:,3)=Power_produced_by_panel;
FinalSolarpowerNancyS2(:,4)=Solar_Current_To_Battery;
FinalSolarpowerNancyS2(:,5)=Fuel_Cell_Current;
FinalSolarpowerNancyS2(:,6)=Battery_Voltage;
FinalSolarpowerNancyS2(:,7)=Current_To_Load;
FinalSolarpowerNancyS2(:,8)=Battery_SOC;
FinalSolarpowerNancyS2(:,9)=Fuel_Cell_Power;
FinalSolarpowerNancyS2(:,10)=Hydrogen_Consumption;

