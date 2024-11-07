%% Assignment 1: 1D steady conduction heat transfer case in cylindrical walls

clc
clear all
format long
%% Numerical Solution
% Input Data - Physical
% Geometry
R_A = 0.25; %m
R_B_range = 0.5; %m
H = 1; %m

% Material Properties - Stainless Steel
lambda_range = 15; % W/mK

% Internal Heat Source
qv_range = 2000; % W/m3

% Boundary Conditions
alpha_A = 30; % W/m2K
alpha_B = 10; % W/m2K
T_A = 80+273.15; % K
T_B = 40+273.15; % K

% Input Data - Numerical
N_list = 100; % number of control volumes

cas_N = 1; % for different mesh sizes
for N=N_list

cas_qv = 1; % for different cases of qv values
for qv = qv_range

cas_lambda = 1; % for different cases of lambda values
for lambda = lambda_range

cas_R_B = 1; % for different cases of R_B vales
for R_B = R_B_range

max_Res = 1e-12;
max_ite = 1e8;
% Mesh Definitions
Dr(cas_N) = (R_B-R_A)./N;
r = linspace(R_A,R_B,N+2);

% Defining Control Volumes
for i = 1:N+1
    rcv(i) = R_A+ (i-1).*Dr(cas_N); 
end

% Defining nodes
rP(1) = R_A;
for i = 2:N+1
rP(i) = (rcv(i)+rcv(i-1))/2;
end
rP(N+2) = R_B;

% Preallocate Coefficients
ap = zeros(1,N+2);
ae = zeros(1,N+2);
aw = zeros(1,N+2);
bp = zeros(1,N+2);
Se = zeros(1,N+2);
Sw = zeros(1,N+2);
Vp = zeros(1,N+2);

% Areas & Control Volumes

for i = 2:N+1
    Sw(i) = 2*pi*rcv(i-1)*H;
    Se(i) = 2*pi*rcv(i)*H;
    Vp(i) = pi*(rcv(i)^2-rcv(i-1)^2)*H;
end

% Discretization Coefficients
for i = 2:N+1
    ae(i) = lambda*Se(i)/(rP(i+1)-rP(i));
    aw(i) = lambda*Sw(i)/(rP(i)-rP(i-1));
    ap(i) = ae(i) + aw(i);
    bp(i) = qv*Vp(i);
end
% Boundary Conditions
% first node
ae(1) = lambda/(rP(2)-rP(1));
ap(1) = ae(1)+alpha_A;
bp(1) = alpha_A*T_A;

% Last node
aw(end) = lambda/(rP(end)-rP(end-1));
ap(end) = aw(end) + alpha_B;
bp(end) = alpha_B*T_B;
% Initilization of Temperatures
T_Guass = zeros(1,N+2);
T_star_guass = zeros(1,N+2);

ite_guass = 0;
res_guass = max_Res +1;




%% Guass-Seidel Iteration Method
% apTp = aeTp + awTw + bp

while res_guass>max_Res && ite_guass < max_ite
   
    T_Guass(1) = ((ae(1)*T_Guass(2) + bp(1))/ap(1));
    for i = 2:N+1
        T_Guass(i) = ((ae(i)*T_Guass(i+1) + aw(i)*T_Guass(i-1) + bp(i))/ap(i));
    end
    T_Guass(N+2) = ((aw(end)*T_Guass(N+2-1) + bp(end))/ap(end));
for i = 1:N+2
loc_res_guass(i) = T_Guass(i)-T_star_guass(i);
end

res_guass = max(abs(loc_res_guass));
T_star_guass = T_Guass;   
ite_guass = ite_guass+1;
end

%% Code Verification Guass Method
Qin_guass = alpha_A*(T_A-T_Guass(1))*2*pi*R_A*H;
Qout_guass = alpha_B*(T_Guass(end)-T_B)*2*pi*R_B*H;
Qv_guass = qv*pi*(R_B^2-R_A^2)*H;
fprintf('Global Energy Balance Guass Method: %.10f\n', Qin_guass - Qout_guass + Qv_guass);


%% TDMA Method

% Initlialization of vectors
P = zeros(1, N+2);
R = zeros(1, N+2);

P(1) = ae(1)/ap(1);
R(1) = bp(1)/ap(1);

for i = 2:N+1
    P(i) = (ae(i)/(ap(i)-aw(i)*P(i-1)));
    R(i) = (bp(i)-aw(i)*R(i-1))/(ap(i)-aw(i)*P(i-1));
end
P(end) = ae(end)/(ap(end)-aw(end)*P(end-1));
R(end) = (bp(end)+aw(end)*R(end-1))/(ap(end)-aw(end)*P(end-1));

% Temperature
T_TDMA = zeros(1,N+2);
T_star_TDMA = zeros(1,N+2);
ite_TDMA = 0;
res_TDMA = max_Res +1;

%Calculation of P and R
P(1) = ae(1)/ap(1);
R(1) = bp(1)/ap(1);

for i = 2:N+1
    P(i) = ae(i)/(ap(i)-aw(i)*P(i-1));
    R(i) = (bp(i)+aw(i)*R(i-1))/(ap(i)-aw(i)*P(i-1));
end
P(end) = 0;
R(end) = (bp(end)+aw(end)*R(end-1))/(ap(end)-aw(end)*P(end-1));

while res_TDMA>max_Res && ite_TDMA < max_ite    
    T_TDMA(end) = R(end);
    for i = N+1:-1:2
        T_TDMA(i) = P(i)*T_TDMA(i+1)+R(i);
    end
    T_TDMA(1) = P(1)*T_TDMA(2)+R(1);
for i = 1:N+2
loc_res_TDMA(i) = T_TDMA(i)-T_star_TDMA(i);
end

res_TDMA = max(abs(loc_res_TDMA));
T_star_TDMA = T_TDMA;
ite_TDMA = ite_TDMA+1;
end

%% Code Verification TDMA Method
Qin_TDMA = alpha_A*(T_A-T_TDMA(1))*2*pi*R_A*H;
Qout_TDMA = alpha_B*(T_TDMA(end)-T_B)*2*pi*R_B*H;
Qv_TDMA = qv*pi*(R_B^2-R_A^2)*H;
fprintf('Global Energy Balance TDMA Method: %.10f\n', Qin_TDMA - Qout_TDMA + Qv_TDMA);

%% Analytical Solution
A = [log(R_A), 1, -1, 0; log(R_B), 1, 0, -1; (lambda/R_A), 0, - alpha_A, 0; (lambda/R_B), 0, 0, alpha_B];
B = [qv*R_A^2/(4*lambda); qv*R_B^2/(4*lambda); R_A*qv/2 - alpha_A*T_A; R_B*qv/2 + alpha_B*T_B];
X = inv(A)*B;

C1 = X(1);
C2 = X(2);
Tw1 = X(3);
Tw2 = X(4);

% Analytical function
T_ana = @(rP) -qv*rP.^2/(4*lambda)+ C1*log(rP) + C2;


truncation_error_TDMA(cas_N)=max(abs(T_ana(rP)-T_TDMA));
truncation_error_Gauss(cas_N)=max(abs(T_ana(rP)-T_Guass));

% Maximum Temperature for each heat generation value
max_temperautre_RB(cas_R_B) = max(abs(T_TDMA));
cas_R_B = cas_R_B+1;

% Uncomment the results and select increase the range of the variable to see the impact of that variable

% % Impact of Mesh Size
% figure(7)
% plot(rP,T_TDMA-273.15);
% xlabel('Position of Node (m)');
% xlim([min(rP),max(rP)]);
% ylabel('Temperature (°C)');
% title("Temperature Profile");
% legend("N=10","N=35","N=60","N=85","N=110","N=135","N=160","N=185")
% hold on

% % Impact of heat generation
% figure(8)
% plot(rP,T_TDMA-273.15);
% xlabel('Position of Node (m)');
% xlim([min(rP),max(rP)]);
% ylabel('Temperature (°C)');
% title("Temperature Profile");
% legend("qv=2000","qv=2200","qv=2400","qv=2600","qv=2800","qv=3000")
% hold on

% % Impact of lambda
% figure(9)
% plot(rP,T_TDMA-273.15);
% xlabel('Position of Node (m)');
% xlim([min(rP),max(rP)]);
% ylabel('Temperature (°C)');
% title("Temperature Profile");
% legend("λ=15","λ=20","λ=25","λ=30","λ=35","λ=40","λ=45")
% hold on
% 
% % Impact of External Radius
% figure(10)
% plot(rP,T_TDMA-273.15);
% xlabel('Position of Node (m)');
% xlim([min(rP),max(rP)]);
% ylabel('Temperature (°C)');
% title("Temperature Profile");
% legend("R_B=0.5","R_B=0.6","R_B=0.7","R_B=0.8","R_B=0.9","R_B=1")
% hold on


end
% Maximum Temperature for each heat generation value
max_temperautre_lambda(cas_lambda) = max(abs(T_TDMA));
cas_lambda = cas_lambda+1;

end
% Maximum Temperature for each heat generation value
max_temperautre_qv(cas_qv) = max(abs(T_TDMA));
cas_qv = cas_qv+1;

end
cas_N = cas_N +1;
end


%% Results

% Temperature profile with Guass Solver
figure(1)
plot(rP,T_Guass-273.15);
xlabel('Position of Node (m)');
xlim([min(rP),max(rP)]);
ylabel('Temperature (°C)');
ylim([min(T_Guass-273.15)-5, max(T_Guass-273.15)+5]);
title("Temperature Profile - Guass");

% Temperature Profile with TDMA solver
figure(2)
plot(rP,T_TDMA-273.15);
xlabel('Position of Node (m)');
xlim([min(rP),max(rP)]);
ylabel('Temperature (°C)');
ylim([min(T_TDMA-273.15)-5, max(T_TDMA-273.15)+5]);
title("Temperature Profile - TDMA");

% Comparison with Analytical solution - TDMA
figure(3)
hold on;
plot(rP, T_TDMA-273.15);
plot(rP, T_ana(rP)-T_TDMA);
xlabel('Position of Node (m)');

yyaxis left
ylabel('Temperature (°C)');
ylim([min(T_TDMA-273.15)-5, max(T_TDMA-273.15)+5]);

yyaxis right
ylabel('Error');
ylim([min(T_ana(rP)-T_TDMA), max(T_ana(rP)-T_TDMA)]);

legend('Numerical Temperature', 'Error');
title("Comparison of Numerical and Analytical Solutions - TDMA");
hold off

% Comparison with Analytical solution - Guass
figure(4)
hold on;
plot(rP, T_Guass-273.15);
plot(rP, T_ana(rP)-T_Guass);
xlabel('Position of Node (m)');

yyaxis left
ylabel('Temperature (°C)');
ylim([min(T_Guass-273.15)-5, max(T_Guass-273.15)+5]);

yyaxis right
ylabel('Error');
ylim([min(T_ana(rP)-T_Guass), max(T_ana(rP)-T_Guass)]);

legend('Numerical Temperature', 'Error');
title("Comparison of Numerical and Analytical Solutions - Gauss");
hold off

% Error Analysis for differnet mesh sizes
figure(5)
loglog(Dr,truncation_error_TDMA,'-o','LineWidth',1)
xlabel('Mesh size (Δr)')
ylabel('Error')
title('Error Analysis')

% Impact of qv on max temperatures of cylinder
figure(6)
plot(qv_range, max_temperautre_qv-273.15,'-o','LineWidth',1);
xlabel('Heat Generation (W/m3)')
ylabel('Maximum Temperature (°C)')
title('Impact of Heat Generation')

% Impact of lambda on max temperatures of cylinder
figure(7)
plot(lambda_range, max_temperautre_lambda-273.15,'-o','LineWidth',1);
xlabel('Thermal Conductivity (W/m.K)')
ylabel('Maximum Temperature (°C)')
title('Impact of Thermal Conductivity')

% Impact of thermal conductivity on max temperatures of cylinder
figure(8)
plot(R_B_range, max_temperautre_RB-273.15,'-o','LineWidth',1);
xlabel('External Radius (m)')
ylabel('Maximum Temperature (°C)')
title('Impact of External Radius')

% Error Analysis for differnet mesh sizes
figure(15)
loglog(Dr,truncation_error_Gauss,'-o','LineWidth',1)
xlabel('Mesh size (Δr)')
ylabel('Error')
title('Error Analysis')


