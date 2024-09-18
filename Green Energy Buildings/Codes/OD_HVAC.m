clear all
close all
clc
%% Month Selection

months = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
selected_month = menu('Select Month', months);
%% Dimensions
% walls, roof, materials thicknesses;
% dimensions are taken from the drawings and autocad file provided
% Assumption: building height is constant
    
%l for length
%w for width
%h for height
%t for thickness
%a for area

%External Glass Wall
length_wall = 40.5; %m
width_wall = 21.5; %m
height_wall = 3.5; %m, v6e v6i
area_wallA = length_wall * height_wall; %m2
area_wallB = width_wall * height_wall; %m2 the 'door' --> will only use v6e for now
thickness_lam_saf_glass = 0.02; %m, 20mm
thickness_lam_saf_glass_air = 0.81; %m
thickness_low_emiss_glass1 = 0.008; %m
thickness_low_emiss_glass_air = 0.012; %m
thickness_low_emiss_glass2 = 0.012; %m

% Brick wall
thickness_brick_wall = 0.2; %m

%Internal Glass Wall
height_intwall = 3.35; %m
thickness_intwall_v1 = 0.012; %m, facing corridor
thickness_intwall_v2 = 0.016; %m, facing adjacent rooms

%Floor
area_floor = length_wall * width_wall;
thickness_floor = 0.12; %m, 12 cm

%Skylight (laminated glass?)
thickness_skylight = 0.02; %m
length_skylight = 30.4; %m
width_skylight = 3; %m
area_skylight = length_skylight * width_skylight;

%Roof
area_roof = length_wall * width_wall - area_skylight;
thickness_roof_sandwich1 = 0.001; %m, steel sheet
thickness_roof_vaporbar = 0.001; %m
thickness_roof_rockwool = 0.08; %m, 80mm
thickness_roof_sandwich2 = 0.0006; %m, 0.6mm steel sheet
thickness_roofglass = 0.02; %m


%% Material Properties
% density, specific heat capacity and thermal conductivity

% 1. Superior Cover (Roof)
% 1.1. Density [kg/m^3]: rho
rho_sandwich_cover=7800;
rho_vapor_barrier=941;
rho_rock_wool=110;

% 1.2. Specific Heat [J/kg.K]: cp
cp_sandwich_cover=0.46*1000;    % previous 0.9
cp_vapor_barrier=1.55*1000;    % 1820
cp_rock_wool=0.9*1000;

% 1.3. Thermal Conductivity [W/m.K]: lambda
lambda_sandwich_cover=60;
lambda_vapor_barrier=0.4;
lambda_rock_wool=0.033;

% 2. External Glass Facade
% 2.1. Density [kg/m^3]: rho
rho_lam_saf_glass=2480;

% 2.2. Specific Heat [J/kg.K]: cp
cp_lam_saf_glass=1.9*1000;                      %double check

% 2.3. Thermal Conductivity [W/m.K]: lambda
lambda_lam_saf_glass=1;  % 1                     %double check

% 3. Internal Glass Facade
% 3.1. Density [kg/m^3]: rho
rho_low_emiss_glass=2500;
rho_air_layer=1.225;

% 3.2. Specific Heat [J/kg.K]: cp
cp_low_emiss_glass=840;
cp_air_layer=1.005*1000;

% 3.3. Thermal Conductivity [W/m.K]: lambda
lambda_low_emiss_glass=1;   % 1 
lambda_air_layer=0.026;

% 4. Internal Glass Wall (V1)/Same as 2. External Glass Facade: Laminated Safety Glass

% 5. Internal Glass Wall (V2)/Same as 2. External Glass Facade: Laminated Safety Glass

% 6. Floor (Corrugated Sheet Steel and Concrete)
% 6.1. Density [kg/m^3]: rho  (light weight concrete)
rho_floor=1200;                                 %double check

% 6.2. Specific Heat [J/kg.K]: cp
cp_floor=475;                                   %double check

% 6.3. Thermal Conductivity [W/m.K]: lambda
lambda_floor=0.38;                              %double check

%7. Brick wall (light weight concrete)
% 7.1 Thermal Conductivity
lambda_brick = 1; %0.38

% 7.2 Specific heat capacity
cp_brick = 840;

% 7.3 density
rho_brick = 480;


%% Resitances and Overall heat transfer coefficients
% Formula R = thickness/Thermal Conductivity
% Assumption: North/South wall = Laminated + Low emissivity glass;
% East/West wall: Only Laminated

% North/South Wall

if ismember(selected_month,[1,2,11,12])
R_northwall = thickness_lam_saf_glass/lambda_lam_saf_glass + thickness_lam_saf_glass_air/lambda_air_layer + thickness_low_emiss_glass1/lambda_low_emiss_glass + thickness_low_emiss_glass_air/lambda_air_layer + thickness_low_emiss_glass2/lambda_low_emiss_glass;
R_southwall = R_northwall;

else
R_northwall = thickness_lam_saf_glass/lambda_lam_saf_glass + thickness_low_emiss_glass1/lambda_low_emiss_glass + thickness_low_emiss_glass_air/lambda_air_layer + thickness_low_emiss_glass2/lambda_low_emiss_glass;    
R_southwall = R_northwall;

end

% East Wall
R_eastwall = thickness_brick_wall/lambda_brick;

% West Wall

R_westwall = thickness_brick_wall/lambda_brick;

% Roof
R_roof_opaque = thickness_roof_sandwich1/lambda_sandwich_cover + thickness_roof_vaporbar/lambda_vapor_barrier + thickness_roof_rockwool/lambda_rock_wool + thickness_roof_sandwich2/lambda_sandwich_cover;

R_roof_glazed = thickness_roofglass/lambda_lam_saf_glass;

% Floor

R_floor = thickness_floor/lambda_floor;

% External & Internal convection resistances m2.K/W

R_se = 0.04;
R_si_wall = 0.13; 
R_si_roof = 0.17;
R_si_floor = 0.1; 

%% Overall heat transfer coefficients U [W/m2.K]
% U = 1/Total resistance

% North Wall

U_northwall = 1/(R_si_wall + R_northwall + R_se);

% South Wall

U_southwall = 1/(R_si_wall + R_southwall + R_se);

% East Wall

U_eastwall = 1/(R_si_wall + R_eastwall + R_se);

% West Wall

U_westwall = 1/(R_si_wall + R_westwall + R_se);

% Roof

U_roof_opaque = 1/(R_si_roof + R_roof_opaque + R_se);
U_roof_glazed = 1/(R_si_roof + R_roof_glazed + R_se);

% Floor

U_floor = 1/(R_si_floor + R_floor + R_se)*0;


% Sum of UA
UA_opaque_roof = U_roof_opaque*area_roof;

UA_overall = (U_northwall*area_wallA + U_southwall*area_wallA + U_eastwall*area_wallB + U_westwall*area_wallB + U_roof_opaque*area_roof + U_roof_glazed*area_skylight + U_floor*area_floor)*1.1;

%% Capacitances [J/K]
% C = rho*Cp*V

% laminated safety glass (walls/roof)

C_lam_saf_glass = rho_lam_saf_glass*cp_lam_saf_glass*(2*area_wallA*thickness_lam_saf_glass + area_skylight*thickness_lam_saf_glass);

% low emissivity glass (walls)

C_low_emiss_glass = rho_low_emiss_glass*cp_low_emiss_glass*(2*area_wallA*(thickness_low_emiss_glass1+thickness_low_emiss_glass2));

% Air (total volume)

C_air_layer = rho_air_layer*cp_air_layer*(area_floor*height_wall);

% Steel (roof)

C_sandwich_cover = rho_sandwich_cover*cp_sandwich_cover*(area_roof*(thickness_roof_sandwich2 + thickness_roof_sandwich1));

% Rockwool (roof)

C_rock_wool = rho_rock_wool*cp_rock_wool*(area_roof*thickness_roof_rockwool);

% Polyethylene (roof)

C_vapor_barrier =rho_vapor_barrier*cp_vapor_barrier*(area_roof*thickness_roof_vaporbar);

% Concrete (Floor)

C_floor = rho_floor*cp_floor*(area_floor*thickness_floor);

% Brick (wall)
C_brick = rho_brick*cp_brick*(2*area_wallB*thickness_brick_wall);

% Overall Capacitance
C_overall = (C_floor+C_vapor_barrier+C_rock_wool+C_sandwich_cover+C_air_layer+C_low_emiss_glass+C_lam_saf_glass+C_brick);

%% Loading Irradiation, temperature, people & equipment profiles
% Temperature and Irradiance from NASA website (hourly data)

% we are choosing july because its the hottest month
days = 15;

delta_time = 60*15/3600;   % hour (1 second time step)
time_data = 0:24*days;
Time = 0:delta_time:24*days;

Profile = xlsread('All_year.xlsx',selected_month);

G_sun_roof = repmat(Profile(1:24,2),days,1);
G_sun_north = repmat(Profile(1:24,3),days,1);
G_sun_south = repmat(Profile(1:24,4),days,1);
G_sun_east = repmat(Profile(1:24,5),days,1);
G_sun_west = repmat(Profile(1:24,6),days,1);

T_ext = repmat(Profile(1:24,7),days,1);
N_people = repmat(Profile(1:24,8),days,1);
N_equip = repmat(Profile(1:24,9),days,1);

f_facade = repmat(Profile(1:24,10),days,1);

G_sun_roof(24*days+1,1) = 0.5*(G_sun_roof(1,1) + G_sun_roof(24*days,1));
G_sun_north(24*days+1,1) = 0.5*(G_sun_north(1,1) + G_sun_north(24*days,1));
G_sun_south(24*days+1,1) = 0.5*(G_sun_south(1,1) + G_sun_south(24*days,1));
G_sun_east(24*days+1,1) = 0.5*(G_sun_east(1,1) + G_sun_east(24*days,1));
G_sun_west(24*days+1,1) = 0.5*(G_sun_west(1,1) + G_sun_west(24*days,1));

T_ext(24*days+1,1) = 0.5*(T_ext(1,1) + T_ext(24*days,1));
N_people(24*days+1,1) = 0.5*(N_people(1,1) + N_people(24*days,1));
N_equip(24*days+1,1) = 0.5*(N_equip(1,1) + N_equip(24*days,1));

f_facade(24*days+1,1) = 0.5*(f_facade(1,1) + f_facade(24*days,1));

T_ext_interpolated = interp1(time_data,T_ext,Time);

G_sun_north_interpolated = interp1(time_data,G_sun_north,Time);
G_sun_south_interpolated = interp1(time_data,G_sun_south,Time);
G_sun_east_interpolated = interp1(time_data,G_sun_east,Time);
G_sun_west_interpolated = interp1(time_data,G_sun_west,Time);
G_sun_roof_interpolated = interp1(time_data,G_sun_roof,Time);

N_people_interpolated = interp1(time_data,N_people,Time);
N_equip_interpolated=interp1(time_data,N_equip,Time);

f_facade_interpolated = interp1(time_data,f_facade,Time);

Time_24 = mod(Time,24);
%% Heat balance loop

% Parameters
n_vent = 0.35/3600; % number of air shifts
f_sh = 1; % shading factor
tau = 0.5; %transmittance of the material to transmit the solar radiation
tau_garden = 0.2;
alpha_roof = 0.475*1; %absorptance of the roof
alpha_brick = 0.475;
A_opaque_roof = alpha_roof*R_se*UA_opaque_roof*1;
A_opaque_westwall = alpha_brick*R_se*U_westwall*area_wallB*1;
A_opaque_eastwall = alpha_brick*R_se*U_eastwall*area_wallB*1;

Q_people = 100; % Watts
Q_equip = 55; % Watts

% Initial conditions

if ismember(selected_month,[12,1,2])
    T_int(1) = 19.28;
elseif ismember(selected_month,[3,4,5])
    T_int(1) = 20.86;
elseif ismember(selected_month,[6,7,8])
    T_int(1) = 23.66;
else ismember(selected_month,[9,10,11])
    T_int(1) = 21.66;
end
%T_int(1) = 22.5;%T_ext_interpolated(1); % Celsius

for n = 1:length(Time)

    A_glazed_roof = tau_garden*area_skylight;
    A_glazed_northwall(n) = tau*area_wallA*(1-f_facade_interpolated(n));
    A_glazed_southwall(n) = tau*area_wallA*1*(1-f_facade_interpolated(n));

    Q_cond_conv(n) = UA_overall*(T_ext_interpolated(n)-T_int(n));
    Q_vent(n) = cp_air_layer*rho_air_layer*(area_floor*height_wall)*n_vent*(T_ext_interpolated(n)-T_int(n));

    Q_rad_opaque_roof(n) = f_sh*A_opaque_roof*G_sun_roof_interpolated(n);
    Q_rad_opaque_west(n) = f_sh*A_opaque_westwall*G_sun_west_interpolated(n);
    Q_rad_opaque_east(n) = f_sh*A_opaque_eastwall*G_sun_east_interpolated(n);

    Q_rad_opaque(n) = Q_rad_opaque_roof(n)+ Q_rad_opaque_west(n) + Q_rad_opaque_east(n);

    Q_rad_glazed_north(n) = f_sh*A_glazed_northwall(n)*G_sun_north_interpolated(n);
    Q_rad_glazed_south(n) = f_sh*A_glazed_southwall(n)*G_sun_south_interpolated(n);
    Q_rad_glazed_roof(n) = f_sh*A_glazed_roof*G_sun_roof_interpolated(n);

    Q_rad_glazed(n) = Q_rad_glazed_north(n) + Q_rad_glazed_south(n)+Q_rad_glazed_roof(n);

    Q_int_people(n) = N_people_interpolated(n)*Q_people;
    Q_int_equip(n)= N_equip_interpolated(n)*Q_equip;

%% HVAC
    T_comfortable = 22.5;
    T_tolerance = 2.5;
    T_change = 0.5;

    if Time_24(n)>=6 & Time_24(n)<=19
        if T_int(n) > T_comfortable + T_tolerance
            T_int(n+1) = T_int(n) - T_change;
        
        elseif T_int(n) < T_comfortable - T_tolerance
            T_int(n+1) = T_int(n) + T_change;
        
        else
            Q_HVAC(n) = 0;
            T_int(n+1) = T_int(n) + delta_time*3600*(Q_cond_conv(n)+Q_vent(n)+Q_rad_opaque(n)+Q_rad_glazed(n)+Q_int_people(n)+Q_int_equip(n)-Q_HVAC(n))/C_overall;
            
        end
        Q_HVAC(n) = (C_overall*(T_int(n+1)-T_int(n))/(delta_time*3600)-(Q_cond_conv(n)+Q_vent(n)+Q_rad_opaque(n)+Q_rad_glazed(n)+Q_int_people(n)+Q_int_equip(n)));
    else
        Q_HVAC(n) = 0;
        T_int(n+1) = T_int(n) + delta_time*3600*(Q_cond_conv(n)+Q_vent(n)+Q_rad_opaque(n)+Q_rad_glazed(n)+Q_int_people(n)+Q_int_equip(n)-Q_HVAC(n))/C_overall;
        
    end
    
    Energy_HVAC = sum(abs(Q_HVAC))*delta_time/1000;
    %Heat per surfaces
    
    Q_northwall(n)= U_northwall*area_wallA *(T_ext_interpolated(n)-T_int(n)) + Q_rad_glazed_north(n);

    Q_southwall(n)= U_southwall*area_wallA *(T_ext_interpolated(n)-T_int(n)) + Q_rad_glazed_south(n);

    Q_eastwall(n)= U_eastwall*area_wallB *(T_ext_interpolated(n)-T_int(n)) + Q_rad_opaque_east(n);

    Q_westwall(n)= U_westwall*area_wallB *(T_ext_interpolated(n)-T_int(n)) + Q_rad_opaque_west(n);

    Q_roof_opaque(n)= U_roof_opaque*area_roof *(T_ext_interpolated(n)-T_int(n)) + Q_rad_opaque_roof(n);

    Q_roof_glazed(n)= U_roof_glazed*area_skylight *(T_ext_interpolated(n)-T_int(n)) + Q_rad_glazed_roof(n);
    
    Q_floor(n)=  U_floor*area_floor*(T_ext_interpolated(n)-T_int(n));


end

Time=[Time, 24*days+delta_time];
T_ext_interpolated(n+1)=0.5*(T_ext_interpolated(1)+T_ext_interpolated(n));


T_comfortable(n+1) = 22.5;


% %% Figures

% figure(1)
% plot(Time,T_int, LineWidth=2)
% xlabel('Time [h]')
% ylabel('Temperature (°C)')
% axis([0 24*days 0 inf],'auto y')
% grid
% 
% hold on
% 
% plot(Time ,T_ext_interpolated, LineWidth=1.5)
% axis([0 24*days 0 inf],'auto y')
% xlabel('Time [h]')
% ylabel('Temperature [°C]')
% grid
% legend('Inside temperature', 'Outside temperature')
% 
% figure(2)
% plot(Time(1:end-1), Q_cond_conv, LineWidth=1.5)
% hold on
% plot(Time(1:end-1), Q_int_equip,LineWidth=1.5)
% plot(Time(1:end-1), Q_int_people,LineWidth=1.5)
% plot(Time(1:end-1), Q_rad_glazed,LineWidth=1.5)
% plot(Time(1:end-1), Q_rad_opaque,LineWidth=1.5)
% plot(Time(1:end-1), Q_vent,LineWidth=1.5)
% axis([0 24*days 0 inf],'auto y')
% xlabel('Time [h]')
% ylabel('Heat [W]')
% legend("Cond/Conv", "Equipment", "People", "Radiation Glazed", "Radiation Opaque", "Ventilation")
% hold off
% 
% figure(3)
% plot(Time(1:end-1), Q_northwall, LineWidth=1.5)
% hold on
% plot(Time(1:end-1), Q_southwall, LineWidth=1.5)
% plot(Time(1:end-1), Q_eastwall, LineWidth=1.5)
% plot(Time(1:end-1), Q_westwall, LineWidth=1.5)
% plot(Time(1:end-1), Q_roof_glazed, LineWidth=1.5)
% plot(Time(1:end-1), Q_roof_opaque, LineWidth=1.5)
% plot(Time(1:end-1), Q_floor, LineWidth=1.5)
% axis([0 24*days 0 inf],'auto y')
% xlabel('Time [h]')
% ylabel('Heat [W]')
% legend("North Wall", "South wall", "East wall", "West wall", "Roof glazed", "Roof opaque")
% 
% figure(4)
% plot(Time(1:end-1),Q_HVAC/1000, LineWidth=1.5)
% xlabel('Time [h]')
% ylabel('Cooling load [W]')
% axis([0 24*days 0 inf],'auto y')
% 
% figure(5)
% plot(Time(1:end-1), G_sun_north_interpolated, LineWidth=1.5)
% hold on
% plot(Time(1:end-1), G_sun_south_interpolated,LineWidth=1.5)
% plot(Time(1:end-1), G_sun_east_interpolated,LineWidth=1.5)
% plot(Time(1:end-1), G_sun_west_interpolated,LineWidth=1.5)
% plot(Time(1:end-1), G_sun_roof_interpolated,LineWidth=1.5)
% legend("North","South","East","West","Roof")
% hold off

capacity_heatpump = max(abs(Q_HVAC)/1000)
