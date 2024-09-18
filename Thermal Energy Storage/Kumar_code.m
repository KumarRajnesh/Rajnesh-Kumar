% General Data

% tank material is stainless steel
rho = 7900;              % kg/m3
Cp = 510;                % J/kg.K
lambda_ss = 17/1000;     % kW/m.K
Th_ss = 6e-3;            % m

% Insulation material: mineral wool
lambda_ins = 0.035/1000; % kW/m.K
Th_ins = 10e-2;          % m

% Tank dimensions
h = 2.5;                % m
d_int = 1.8;            % m            
d_ss = d_int + 2*Th_ss  % m
d_ext = d_ss + 2*Th_ins % m

% Water properties
Cp_w = 4.18;            % K.J/kg.K
rho_w = 1000;           % kg/m3
V_w = pi*(d_int^2)/4*h  % m3

% Ambient properties
T_amb = 20 + 273.15;    % K
alpha_int = 23/1000;    % kW/m2.K
alpha_ext = 12/1000;    % kW/m2.K

% Resistances & Overall Conductivity
% Energy loss of top and bottom are neglected

Res_int = 1./(alpha_int*pi*d_int*h);
Res_ss = (log(d_ss/d_int))./(2*pi*Th_ss*lambda_ss);
Res_ins = (log(d_ext/d_ss))./(2*pi*Th_ins*lambda_ins);
Res_ext = 1./(alpha_ext*pi*d_ext*h);

Res_total = Res_int + Res_ss + Res_ins + Res_ext;

UA = 1./Res_total;            % kW/K

% CASE:1 specific data

% Source data
T_A_source = 90 +273.15;      % K

% Load data
T_A_load = 20 + 273.15;       % K


% Initializing Matrices
step_size = 0.1;              % hour   
Total_days = 5;               % days
Matrix_rows = (Total_days*24/step_size) +1;

Days = zeros(Matrix_rows,1)
Time = zeros(Matrix_rows,1);
Acc_Time = zeros(Matrix_rows,1);
T_tank = zeros(Matrix_rows,1);
T_B_source = zeros(Matrix_rows,1);
T_B_load = zeros(Matrix_rows,1);
deltaE_source = zeros(Matrix_rows,1);
deltaE_load = zeros(Matrix_rows,1);
deltaE_convection = zeros(Matrix_rows,1);
Results = zeros(Matrix_rows,7);


% Standard case

% Initial conditions
n = 1;
Temp_tank = 20 + 273.15;    % K
Time(n,1) = 0*3600;         % s
deltaE_source(n,1) = 0;     % kJ
deltaE_load(n,1) = 0;       % kJ
deltaE_convection(n,1) = 0; % kJ
Acc_Time(n,1) = 0;          % hour
T_tank(n,1) = 20 + 273.15;  % K


% 5 days "for loop"
for d = 1:1:Total_days

% 24h "for loop"
    for t = 0.1:step_size:24
        Time(n+1, 1) = t;   % hour
        Days(n,1) = d;

        if Time(n,1) < 8
            m_source = 0;
            m_load = 0;

        elseif Time(n,1) >=8 && Time(n,1) < 13
            m_source = 0.6;
            m_load = 0;

        elseif Time(n,1) >= 13 && Time(n,1) <= 16
            m_source = 0.6;
            m_load = 0.4;

        else Time(n,1) > 16 && Time(n,1) <= 24
            m_source = 0;
            m_load = 0.4;

        end   
    % Energy Balance: Energy in (source) - Energy out (load and convection) = Change in Energy of tank
    T_B_source(n,1) = T_tank(n,1); %Assumption
    T_B_load(n,1) = T_tank(n,1);   %Assumption

    deltaE_source(n,1) = m_source.*Cp_w.*(T_A_source - T_B_source(n,1)).*step_size.*3600;
    deltaE_load(n,1) = m_load.*Cp_w.*(T_B_load(n,1) - T_A_load).*step_size.*3600;
    deltaE_convection(n,1) = UA.*(T_tank(n,1) - T_amb).*step_size.*3600;

    T_tank(n+1,1) = T_tank(n,1) + (deltaE_source(n,1) - deltaE_load(n,1) - deltaE_convection(n,1))./(rho_w.*V_w.*Cp_w);

    Acc_Time(n+1,1) = Acc_Time(n,1) + step_size;

    n = n+1;    
    end
    

end

Results(:,1) = Days;
Results(:,2) = Time;
Results(:,3) = Acc_Time;
Results(:,4) = T_tank;
Results(:,5) = deltaE_source;
Results(:,6) = deltaE_load;
Results(:,7) = deltaE_convection;

figure(1);
plot(Acc_Time, T_tank-273.15);
xlabel('Time (hrs)');
ylabel('Tank Temperature (C)')
title('Evolution of Tank Temperature with Time')

figure(2);
plot(Acc_Time,deltaE_source, 'b');
hold 
plot(Acc_Time,deltaE_load, 'r');
hold off
xlabel('Time (hrs)');
ylabel('Energy (kJ)');
title('Evolution of Energy with Time')
legend('deltaE_source', 'deltaE,load');


% Parametric Study

% Effect of m_source on Tank Temperature

for m_s = 0.4:0.1:0.7

% Initial conditions
n = 1
Time(n,1) = 0*3600;         % s
deltaE_source(n,1) = 0;     % kJ
deltaE_load(n,1) = 0;       % kJ
deltaE_convection(n,1) = 0; % kJ
Acc_Time(n,1) = 0;          % hour
T_tank(n,1) = 20 + 273.15;  % K


for d = 1:1:Total_days

    for t = 0.1:step_size:24
        Time(n+1, 1) = t;   % hour
        Days(n,1) = d;

        if Time(n,1) < 8
            m_source = 0;
            m_load = 0;

        elseif Time(n,1) >=8 && Time(n,1) < 13
            m_source = m_s;
            m_load = 0;

        elseif Time(n,1) >= 13 && Time(n,1) <= 16
            m_source = m_s;
            m_load = 0.4;

        else Time(n,1) > 16 && Time(n,1) <= 24
            m_source = 0;
            m_load = 0.4;

        end   
    % Energy Balance: Energy in (source) - Energy out (load and convection) = Change in Energy of tank
    T_B_source(n,1) = T_tank(n,1); %Assumption
    T_B_load(n,1) = T_tank(n,1);   %Assumption

    deltaE_source(n,1) = m_source.*Cp_w.*(T_A_source - T_B_source(n,1)).*step_size.*3600;
    deltaE_load(n,1) = m_load.*Cp_w.*(T_B_load(n,1) - T_A_load).*step_size.*3600;
    deltaE_convection(n,1) = UA.*(T_tank(n,1) - T_amb).*step_size.*3600;

    T_tank(n+1,1) = T_tank(n,1) + (deltaE_source(n,1) - deltaE_load(n,1) - deltaE_convection(n,1))./(rho_w.*V_w.*Cp_w);

    Acc_Time(n+1,1) = Acc_Time(n,1) + step_size;

    n = n+1;    
    end
    

end
figure(3);
plot(Acc_Time, T_tank-273.15);
hold on
title("Evolution of Tank Temperature with Time for different flow rates of source");
xlabel('Time (Hrs)');
ylabel('Tank Temperature (C)');
legend('m,source = 0.4', 'm,source = 0.5', 'm,source = 0.6', 'm,source = 0.7');

figure(4);
plot(Acc_Time, deltaE_source);
title("Evolution of Energy extracted from the source");
xlabel('Time (Hrs)');
ylabel('Energy (kJ)');
legend('m,source = 0.4', 'm,source = 0.5', 'm,source = 0.6', 'm,source = 0.7');
hold on

figure(5);
plot(Acc_Time, deltaE_load);
title("Evolution of Energy supplied to the load");
xlabel('Time (Hrs)');
ylabel('Energy (kJ)');
legend('m,source = 0.4', 'm,source = 0.5', 'm,source = 0.6', 'm,source = 0.7');
hold on

end
hold off
figure


% Effect of m_load on Tank Temperature

for m_l = 0.2:0.1:0.5

% Initial conditions
n = 1
Time(n,1) = 0*3600;         % s
deltaE_source(n,1) = 0;     % kJ
deltaE_load(n,1) = 0;       % kJ
deltaE_convection(n,1) = 0; % kJ
Acc_Time(n,1) = 0;          % hour
T_tank(n,1) = 20 + 273.15;  % K


for d = 1:1:Total_days

    for t = 0.1:step_size:24
        Time(n+1, 1) = t;   % hour
        Days(n,1) = d;

        if Time(n,1) < 8
            m_source = 0;
            m_load = 0;

        elseif Time(n,1) >=8 && Time(n,1) < 13
            m_source = 0.6;
            m_load = 0;

        elseif Time(n,1) >= 13 && Time(n,1) <= 16
            m_source = 0.6;
            m_load = m_l;

        else Time(n,1) > 16 && Time(n,1) <= 24
            m_source = 0;
            m_load = m_l;

        end   
    % Energy Balance: Energy in (source) - Energy out (load and convection) = Change in Energy of tank
    T_B_source(n,1) = T_tank(n,1); %Assumption
    T_B_load(n,1) = T_tank(n,1);   %Assumption

    deltaE_source(n,1) = m_source.*Cp_w.*(T_A_source - T_B_source(n,1)).*step_size.*3600;
    deltaE_load(n,1) = m_load.*Cp_w.*(T_B_load(n,1) - T_A_load).*step_size.*3600;
    deltaE_convection(n,1) = UA.*(T_tank(n,1) - T_amb).*step_size.*3600;

    T_tank(n+1,1) = T_tank(n,1) + (deltaE_source(n,1) - deltaE_load(n,1) - deltaE_convection(n,1))./(rho_w.*V_w.*Cp_w);

    Acc_Time(n+1,1) = Acc_Time(n,1) + step_size;

    n = n+1;    
    end
    

end
figure(6);
plot(Acc_Time, T_tank-273.15);
hold on
title("Evolution of Tank Temperature with Time for different flowrates of load");
xlabel('Time (Hrs)');
ylabel('Tank Temperature (C)');
legend('m,load = 0.2', 'm,load = 0.3', 'm,load = 0.4', 'm,load = 0.5');

figure(7);
plot(Acc_Time, deltaE_source);
title("Evolution of Energy extracted from the source");
xlabel('Time (Hrs)');
ylabel('Energy (kJ)');
legend('m,load = 0.2', 'm,load = 0.3', 'm,load = 0.4', 'm,load = 0.5');
hold on

figure(8);
plot(Acc_Time, deltaE_load);
title("Evolution of Energy supplied to the load");
xlabel('Time (Hrs)');
ylabel('Energy (kJ)');
legend('m,load = 0.2', 'm,load = 0.3', 'm,load = 0.4', 'm,load = 0.5');
hold on

end
hold off
figure



% Effect of T_A_source on Tank Temperature

for T_s = 70:10:100
    T_A_source = T_s + 273.15;

% Initial conditions
n = 1
Time(n,1) = 0*3600;         % s
deltaE_source(n,1) = 0;     % kJ
deltaE_load(n,1) = 0;       % kJ
deltaE_convection(n,1) = 0; % kJ
Acc_Time(n,1) = 0;          % hour
T_tank(n,1) = 20 + 273.15;  % K

for d = 1:1:Total_days

    for t = 0.1:step_size:24
        Time(n+1, 1) = t;   % hour
        Days(n,1) = d;

        if Time(n,1) < 8
            m_source = 0;
            m_load = 0;

        elseif Time(n,1) >=8 && Time(n,1) < 13
            m_source = 0.6;
            m_load = 0;

        elseif Time(n,1) >= 13 && Time(n,1) <= 16
            m_source = 0.6;
            m_load = m_l;

        else Time(n,1) > 16 && Time(n,1) <= 24
            m_source = 0;
            m_load = m_l;

        end   
    % Energy Balance: Energy in (source) - Energy out (load and convection) = Change in Energy of tank
    T_B_source(n,1) = T_tank(n,1); %Assumption
    T_B_load(n,1) = T_tank(n,1);   %Assumption

    deltaE_source(n,1) = m_source.*Cp_w.*(T_A_source - T_B_source(n,1)).*step_size.*3600;
    deltaE_load(n,1) = m_load.*Cp_w.*(T_B_load(n,1) - T_A_load).*step_size.*3600;
    deltaE_convection(n,1) = UA.*(T_tank(n,1) - T_amb).*step_size.*3600;

    T_tank(n+1,1) = T_tank(n,1) + (deltaE_source(n,1) - deltaE_load(n,1) - deltaE_convection(n,1))./(rho_w.*V_w.*Cp_w);

    Acc_Time(n+1,1) = Acc_Time(n,1) + step_size;

    n = n+1;    
    end
    

end
figure(9);
plot(Acc_Time, T_tank-273.15);
hold on
title("Evolution of Tank Temperature with Time for different inlet temperatures of source");
xlabel('Time (Hrs)');
ylabel('Tank Temperature (C)');
legend('T,A,source = 70', 'T,A,source  = 80', 'T,A,source  = 90', 'T,A,source  = 100');

figure(10);
plot(Acc_Time, deltaE_source);
title("Evolution of Energy extracted from the source");
xlabel('Time (Hrs)');
ylabel('Energy (kJ)');
legend('T,A,source = 70', 'T,A,source  = 80', 'T,A,source  = 90', 'T,A,source  = 100');
hold on

figure(11);
plot(Acc_Time, deltaE_load);
title("Evolution of Energy supplied to the load");
xlabel('Time (Hrs)');
ylabel('Energy (kJ)');
legend('T,A,source = 70', 'T,A,source  = 80', 'T,A,source  = 90', 'T,A,source  = 100');
hold on

end
hold off
figure


% Effect of T_A_load on Tank Temperature

for T_l = 20:5:35
    T_A_load = T_l + 273.15;

% Initial conditions
n = 1
Time(n,1) = 0*3600;         % s
deltaE_source(n,1) = 0;     % kJ
deltaE_load(n,1) = 0;       % kJ
deltaE_convection(n,1) = 0; % kJ
Acc_Time(n,1) = 0;          % hour
T_tank(n,1) = 20 + 273.15;  % K


for d = 1:1:Total_days

    for t = 0.1:step_size:24
        Time(n+1, 1) = t;   % hour
        Days(n,1) = d;

        if Time(n,1) < 8
            m_source = 0;
            m_load = 0;

        elseif Time(n,1) >=8 && Time(n,1) < 13
            m_source = 0.6;
            m_load = 0;

        elseif Time(n,1) >= 13 && Time(n,1) <= 16
            m_source = 0.6;
            m_load = m_l;

        else Time(n,1) > 16 && Time(n,1) <= 24
            m_source = 0;
            m_load = m_l;

        end   
    % Energy Balance: Energy in (source) - Energy out (load and convection) = Change in Energy of tank
    T_B_source(n,1) = T_tank(n,1); %Assumption
    T_B_load(n,1) = T_tank(n,1);   %Assumption

    deltaE_source(n,1) = m_source.*Cp_w.*(T_A_source - T_B_source(n,1)).*step_size.*3600;
    deltaE_load(n,1) = m_load.*Cp_w.*(T_B_load(n,1) - T_A_load).*step_size.*3600;
    deltaE_convection(n,1) = UA.*(T_tank(n,1) - T_amb).*step_size.*3600;

    T_tank(n+1,1) = T_tank(n,1) + (deltaE_source(n,1) - deltaE_load(n,1) - deltaE_convection(n,1))./(rho_w.*V_w.*Cp_w);

    Acc_Time(n+1,1) = Acc_Time(n,1) + step_size;

    n = n+1;    
    end
    

end
figure(12);
plot(Acc_Time, T_tank-273.15);
hold on
title("Evolution of Tank Temperature with Time for different inlet temperatures of load");
xlabel('Time (Hrs)');
ylabel('Tank Temperature (C)');
legend('T,A,load = 20', 'T,A,load  = 25', 'T,A,load  = 30', 'T,A,load  = 35');

figure(13);
plot(Acc_Time, deltaE_source);
title("Evolution of Energy extracted from the source");
xlabel('Time (Hrs)');
ylabel('Energy (kJ)');
legend('T,A,load = 20', 'T,A,load  = 25', 'T,A,load  = 30', 'T,A,load  = 35');
hold on

figure(14);
plot(Acc_Time, deltaE_load);
title("Evolution of Energy supplied to the load");
xlabel('Time (Hrs)');
ylabel('Energy (kJ)');
legend('T,A,load = 20', 'T,A,load  = 25', 'T,A,load  = 30', 'T,A,load  = 35');
hold on

end
hold off
