clear all
close all
clc
%% Month Selection

months = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
selected_month = menu('Select Month', months);

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

%% Dimensions in meter

% Area floor
area_room_floor(1) = 70.38; 
area_room_floor(2) = 34.35;
area_room_floor(3) = 32.89;
area_room_floor(4) = 35.8;
area_room_floor(5) = 34.35;
area_room_floor(6) = 34.35;
area_room_floor(7) = 33.99;
area_room_floor(8) = 223.42;
area_room_floor(9) = 70.38;
area_room_floor(10) = 34.35;
area_room_floor(11) = 32.88;
area_room_floor(12) = 36.85;
area_room_floor(13) = 81.01;
area_room_floor(14) = 34.1;

total_internal_area = sum(area_room_floor);

% Area walls outer North/South

area_wall_outer(1) = 26.703;
area_wall_outer(2) = 13.5;
area_wall_outer(3) = 12.96;
area_wall_outer(4) = 14.04;
area_wall_outer(5) = 13.5;
area_wall_outer(6) = 13.5;
area_wall_outer(7) = 12.83;
area_wall_outer(9) = area_wall_outer(1);
area_wall_outer(10) = area_wall_outer(2);
area_wall_outer(11) = area_wall_outer(3);
area_wall_outer(12) = area_wall_outer(4);
area_wall_outer(13) = 26.81;
area_wall_outer(14) = area_wall_outer(7);

% Area walls inner North/South

area_wall_inner(1) = 34.62;
area_wall_inner(2) = 17.50;
area_wall_inner(3) = 16.8;
area_wall_inner(4) = 18.2;
area_wall_inner(5) = 17.5;
area_wall_inner(6) = 17.5;
area_wall_inner(7) = 16.63;
area_wall_inner(9) = area_wall_inner(1);
area_wall_inner(10) = area_wall_inner(2);
area_wall_inner(11) = area_wall_inner(3);
area_wall_inner(12) = area_wall_inner(4);
area_wall_inner(13) = 34.76;
area_wall_inner(14) = area_wall_inner(7);

% Area walls East/West

area_wall_east_room = 22.54;
area_wall_east_garden = 19.64;

% Previous Dimensions

% Walls thicknesses

% North/South walls
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

average_height_room = (3.5+2.7)/2;

thickness_lam_saf_glass_corridor = 0.012;

% East/West walls

thickness_brick_wall = 0.2; %m
thickness_lam_saf_glass_adjacent = 0.016;

% Floor
area_floor = length_wall * width_wall;
thickness_floor = 0.12;

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

area_roof_garden_opaque = area_room_floor(8)-area_skylight;


%% Resitances and Overall heat transfer coefficients

% Formula R = thickness/Thermal Conductivity
% Assumption: North/South wall = Laminated + Low emissivity glass;
% East/West wall: Only Laminated

% New Resistances

R_eastwall_inner = thickness_lam_saf_glass_adjacent/lambda_lam_saf_glass;

R_northwall_inner = thickness_lam_saf_glass_corridor/lambda_lam_saf_glass;

R_eastwall_garden = thickness_lam_saf_glass_adjacent/lambda_lam_saf_glass;

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

R_roof_glazed = thickness_skylight/lambda_lam_saf_glass;

% Floor

R_floor = thickness_floor/lambda_floor;

% External & Internal convection resistances m2.K/W

R_se = 0.04;
R_si_wall = 0.13; 
R_si_roof = 0.17;
R_si_floor = 0.1; 

%% Overall heat transfer coefficients U [W/m2.K]
% U = 1/Total resistance


% New coefficients

U_eastwall_inner = 1/(R_si_wall + R_eastwall_inner + R_si_wall);

U_eastwall_inner_brick = 1/(R_si_wall + R_eastwall + R_si_wall);

U_eastwall_garden = 1/(R_si_wall + R_eastwall_garden + R_se);

U_northwall_inner = 1/(R_si_wall + R_northwall_inner + R_si_wall);

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
%UA_opaque_roof = U_roof_opaque*area_roof;

%UA_overall = (U_northwall*area_wallA + U_southwall*area_wallA + U_eastwall*area_wallB + U_westwall*area_wallB + U_roof_opaque*area_roof + U_roof_glazed*area_skylight + U_floor*area_floor)*1.1;

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

T_ext_interpolated = interp1(time_data,T_ext,Time)';

G_sun_north_interpolated = interp1(time_data,G_sun_north,Time)';
G_sun_south_interpolated = interp1(time_data,G_sun_south,Time)';
G_sun_east_interpolated = interp1(time_data,G_sun_east,Time)';
G_sun_west_interpolated = interp1(time_data,G_sun_west,Time)';
G_sun_roof_interpolated = interp1(time_data,G_sun_roof,Time)';

N_people_interpolated = interp1(time_data,N_people,Time)';
N_equip_interpolated=interp1(time_data,N_equip,Time)';

f_facade_interpolated = interp1(time_data,f_facade,Time)';

%% Heat balance loop

% Parameters
n_vent = 0.35/3600; % number of air shifts
f_sh = 1; % shading factor
tau = 0.5; %transmittance of the material to transmit the solar radiation
tau_garden = 0.2;
alpha_roof = 0.475*1; %absorptance of the roof
alpha_brick = 0.475*1;

Q_people = 100; % Watts
Q_equip = 55; % Watts

% % Initial Condition
% T_int(1,1) = 22.58;%T_ext_interpolated(1); % Celsius
% T_int(1,2) = 22.58;
% T_int(1,3) = 22.58;
% T_int(1,4) = 22.58;
% T_int(1,5) = 22.58;
% T_int(1,6) = 22.58;
% T_int(1,7) = 22.58;
% T_int(1,8) = 22.58;
% T_int(1,9) = 22.58;
% T_int(1,10) = 22.58;
% T_int(1,11) = 22.58;
% T_int(1,12) = 22.58;
% T_int(1,13) = 22.58;
% T_int(1,14) = 22.58;

% % Initial Condition
% T_int(1,1) = T_ext_interpolated(1); % Celsius
% T_int(1,2) = T_ext_interpolated(1); % Celsius
% T_int(1,3) = T_ext_interpolated(1); % Celsius
% T_int(1,4) = T_ext_interpolated(1); % Celsius
% T_int(1,5) = T_ext_interpolated(1); % Celsius
% T_int(1,6) = T_ext_interpolated(1); % Celsius
% T_int(1,7) = T_ext_interpolated(1); % Celsius
% T_int(1,8) = T_ext_interpolated(1); % Celsius
% T_int(1,9) = T_ext_interpolated(1); % Celsius
% T_int(1,10) = T_ext_interpolated(1); % Celsius
% T_int(1,11) = T_ext_interpolated(1); % Celsius
% T_int(1,12) = T_ext_interpolated(1); % Celsius
% T_int(1,13) = T_ext_interpolated(1); % Celsius
% T_int(1,14) = T_ext_interpolated(1); % Celsius

if ismember(selected_month,[12,1,2])
    T_int(1:14) = 19.28;
elseif ismember(selected_month,[3,4,5])
    T_int(1:14) = 20.86;
elseif ismember(selected_month,[6,7,8])
    T_int(1:14) = 23.66;
else ismember(selected_month,[9,10,11])
    T_int(1:14) = 21.66;
end

for n = 1:length(Time)

    % Room 1
    
    Q_cond_conv(n,1) = U_eastwall_inner*area_wall_east_room*(T_int(n,2)-T_int(n,1))+U_northwall_inner*area_wall_inner(1)*(T_int(n,8)-T_int(n,1))+(U_westwall*area_wall_east_room+U_northwall*area_wall_outer(1)+U_roof_opaque*area_room_floor(1))*(T_ext_interpolated(n,1)-T_int(n,1));
    Q_int_people(n,1) = 0*Q_people;
    Q_int_equip(n,1)= 0*Q_equip;
    Q_vent(n,1) = cp_air_layer*rho_air_layer*(area_room_floor(1)*average_height_room)*n_vent*(T_ext_interpolated(n,1)-T_int(n,1));

    Q_rad_opaque_brick(n,1) = f_sh*(alpha_brick*R_se*U_westwall*area_wall_east_room)*G_sun_west_interpolated(n,1);
    Q_rad_opaque_roof(n,1) = f_sh*(alpha_roof*R_se*U_roof_opaque*area_room_floor(1))*G_sun_roof_interpolated(n,1);

    Q_rad_glazed(n,1) = f_sh*(tau*area_wall_outer(1)*(1-f_facade_interpolated(n,1)))*G_sun_north_interpolated(n,1);

    C_room(1) = C_overall*area_room_floor(1)/total_internal_area;
    T_int(n+1,1) = T_int(n,1) + delta_time*3600/C_room(1)*(Q_vent(n,1)+Q_cond_conv(n,1)+Q_int_people(n,1)+Q_int_equip(n,1)+Q_rad_opaque_brick(n,1)+Q_rad_opaque_roof(n,1)+Q_rad_glazed(n,1));

    % Room 2

    Q_cond_conv(n,2) = U_eastwall_inner*area_wall_east_room*(T_int(n,1)-T_int(n,2))+U_eastwall_inner*area_wall_east_room*(T_int(n,3)-T_int(n,2))+U_northwall_inner*area_wall_inner(2)*(T_int(n,8)-T_int(n,2))+(U_northwall*area_wall_outer(2)+U_roof_opaque*area_room_floor(2))*(T_ext_interpolated(n,1)-T_int(n,2));
    Q_int_people(n,2) = 3*Q_people;
    Q_int_equip(n,2)= 3*Q_equip;
    Q_vent(n,2) = cp_air_layer*rho_air_layer*(area_room_floor(2)*average_height_room)*n_vent*(T_ext_interpolated(n,1)-T_int(n,2));

    Q_rad_opaque_roof(n,2) = f_sh*(alpha_roof*R_se*U_roof_opaque*area_room_floor(2))*G_sun_roof_interpolated(n,1);

    Q_rad_glazed(n,2) = f_sh*(tau*area_wall_outer(2)*(1-f_facade_interpolated(n,1)))*G_sun_north_interpolated(n,1);

    C_room(2) = C_overall*area_room_floor(2)/total_internal_area;
    T_int(n+1,2) = T_int(n,2) + delta_time*3600/C_room(2)*(Q_vent(n,2)+Q_cond_conv(n,2)+Q_int_people(n,2)+Q_int_equip(n,2)+Q_rad_opaque_roof(n,2)+Q_rad_glazed(n,2));


    % Room 3

    Q_cond_conv(n,3) = U_eastwall_inner*area_wall_east_room*(T_int(n,2)-T_int(n,3))+U_eastwall_inner*area_wall_east_room*(T_int(n,4)-T_int(n,3))+U_northwall_inner*area_wall_inner(3)*(T_int(n,8)-T_int(n,3))+(U_northwall*area_wall_outer(3)+U_roof_opaque*area_room_floor(3))*(T_ext_interpolated(n,1)-T_int(n,3));
    Q_int_people(n,3) = 3*Q_people;
    Q_int_equip(n,3)= 3*Q_equip;
    Q_vent(n,3) = cp_air_layer*rho_air_layer*(area_room_floor(3)*average_height_room)*n_vent*(T_ext_interpolated(n,1)-T_int(n,3));

    Q_rad_opaque_roof(n,3) = f_sh*(alpha_roof*R_se*U_roof_opaque*area_room_floor(3))*G_sun_roof_interpolated(n,1);

    Q_rad_glazed(n,3) = f_sh*(tau*area_wall_outer(3)*(1-f_facade_interpolated(n,1)))*G_sun_north_interpolated(n,1);

    C_room(3) = C_overall*area_room_floor(3)/total_internal_area;
    T_int(n+1,3) = T_int(n,3) + delta_time*3600/C_room(3)*(Q_vent(n,3)+Q_cond_conv(n,3)+Q_int_people(n,3)+Q_int_equip(n,3)+Q_rad_opaque_roof(n,3)+Q_rad_glazed(n,3));

    % Room 4

    Q_cond_conv(n,4) = U_eastwall_inner*area_wall_east_room*(T_int(n,3)-T_int(n,4))+U_eastwall_inner*area_wall_east_room*(T_int(n,5)-T_int(n,4))+U_northwall_inner*area_wall_inner(4)*(T_int(n,8)-T_int(n,4))+(U_northwall*area_wall_outer(4)+U_roof_opaque*area_room_floor(4))*(T_ext_interpolated(n,1)-T_int(n,4));
    Q_int_people(n,4) = 3*Q_people;
    Q_int_equip(n,4)= 3*Q_equip;
    Q_vent(n,4) = cp_air_layer*rho_air_layer*(area_room_floor(4)*average_height_room)*n_vent*(T_ext_interpolated(n,1)-T_int(n,4));

    Q_rad_opaque_roof(n,4) = f_sh*(alpha_roof*R_se*U_roof_opaque*area_room_floor(4))*G_sun_roof_interpolated(n,1);

    Q_rad_glazed(n,4) = f_sh*(tau*area_wall_outer(4)*(1-f_facade_interpolated(n,1)))*G_sun_north_interpolated(n,1);

    C_room(4) = C_overall*area_room_floor(4)/total_internal_area;
    T_int(n+1,4) = T_int(n,4) + delta_time*3600/C_room(4)*(Q_vent(n,4)+Q_cond_conv(n,4)+Q_int_people(n,4)+Q_int_equip(n,4)+Q_rad_opaque_roof(n,4)+Q_rad_glazed(n,4));
    
    % Room 5

    Q_cond_conv(n,5) = U_eastwall_inner*area_wall_east_room*(T_int(n,4)-T_int(n,5))+U_eastwall_inner*area_wall_east_room*(T_int(n,6)-T_int(n,5))+U_northwall_inner*area_wall_inner(5)*(T_int(n,8)-T_int(n,5))+(U_northwall*area_wall_outer(5)+U_roof_opaque*area_room_floor(5))*(T_ext_interpolated(n,1)-T_int(n,5));
    Q_int_people(n,5) = 3*Q_people;
    Q_int_equip(n,5)= 3*Q_equip;
    Q_vent(n,5) = cp_air_layer*rho_air_layer*(area_room_floor(5)*average_height_room)*n_vent*(T_ext_interpolated(n,1)-T_int(n,5));

    Q_rad_opaque_roof(n,5) = f_sh*(alpha_roof*R_se*U_roof_opaque*area_room_floor(5))*G_sun_roof_interpolated(n,1);

    Q_rad_glazed(n,5) = f_sh*(tau*area_wall_outer(5)*(1-f_facade_interpolated(n,1)))*G_sun_north_interpolated(n,1);

    C_room(5) = C_overall*area_room_floor(5)/total_internal_area;
    T_int(n+1,5) = T_int(n,5) + delta_time*3600/C_room(5)*(Q_vent(n,5)+Q_cond_conv(n,5)+Q_int_people(n,5)+Q_int_equip(n,5)+Q_rad_opaque_roof(n,5)+Q_rad_glazed(n,5));

    % Room 6

    Q_cond_conv(n,6) = U_eastwall_inner*area_wall_east_room*(T_int(n,5)-T_int(n,6))+U_eastwall_inner*area_wall_east_room*(T_int(n,7)-T_int(n,6))+U_northwall_inner*area_wall_inner(6)*(T_int(n,8)-T_int(n,6))+(U_northwall*area_wall_outer(6)+U_roof_opaque*area_room_floor(6))*(T_ext_interpolated(n,1)-T_int(n,6));
    Q_int_people(n,6) = 3*Q_people;
    Q_int_equip(n,6)= 3*Q_equip;
    Q_vent(n,6) = cp_air_layer*rho_air_layer*(area_room_floor(6)*average_height_room)*n_vent*(T_ext_interpolated(n,1)-T_int(n,6));

    Q_rad_opaque_roof(n,6) = f_sh*(alpha_roof*R_se*U_roof_opaque*area_room_floor(6))*G_sun_roof_interpolated(n,1);

    Q_rad_glazed(n,6) = f_sh*(tau*area_wall_outer(6)*(1-f_facade_interpolated(n,1)))*G_sun_north_interpolated(n,1);

    C_room(6) = C_overall*area_room_floor(6)/total_internal_area;
    T_int(n+1,6) = T_int(n,6) + delta_time*3600/C_room(6)*(Q_vent(n,6)+Q_cond_conv(n,6)+Q_int_people(n,6)+Q_int_equip(n,6)+Q_rad_opaque_roof(n,6)+Q_rad_glazed(n,6));
    
    % Room 7

    Q_cond_conv(n,7) = U_eastwall_inner*area_wall_east_room*(T_int(n,6)-T_int(n,7))+U_northwall_inner*area_wall_inner(7)*(T_int(n,8)-T_int(n,7))+(U_westwall*area_wall_east_room+U_northwall*area_wall_outer(7)+U_roof_opaque*area_room_floor(7))*(T_ext_interpolated(n,1)-T_int(n,7));
    Q_int_people(n,7) = 0*Q_people;
    Q_int_equip(n,7)= 0*Q_equip;
    Q_vent(n,7) = cp_air_layer*rho_air_layer*(area_room_floor(7)*average_height_room)*n_vent*(T_ext_interpolated(n,1)-T_int(n,7));

    Q_rad_opaque_brick(n,7) = f_sh*(alpha_brick*R_se*U_westwall*area_wall_east_room)*G_sun_east_interpolated(n,1);
    Q_rad_opaque_roof(n,7) = f_sh*(alpha_roof*R_se*U_roof_opaque*area_room_floor(7))*G_sun_roof_interpolated(n,1);

    Q_rad_glazed(n,7) = f_sh*(tau*area_wall_outer(7)*(1-f_facade_interpolated(n,1)))*G_sun_north_interpolated(n,1);

    C_room(7) = C_overall*area_room_floor(7)/total_internal_area;
    T_int(n+1,7) = T_int(n,7) + delta_time*3600/C_room(7)*(Q_vent(n,7)+Q_cond_conv(n,7)+Q_int_people(n,7)+Q_int_equip(n,7)+Q_rad_opaque_brick(n,7)+Q_rad_opaque_roof(n,7)+Q_rad_glazed(n,7));
    
    % Room 8
    Q_cond_conv(n,8) = U_northwall_inner*(area_wall_inner(1)*(T_int(n,1)-T_int(n,8))+area_wall_inner(2)*(T_int(n,2)-T_int(n,8))+area_wall_inner(3)*(T_int(n,3)-T_int(n,8)) ...
        +area_wall_inner(4)*(T_int(n,4)-T_int(n,8))+area_wall_inner(5)*(T_int(n,5)-T_int(n,8))+area_wall_inner(6)*(T_int(n,6)-T_int(n,8))+area_wall_inner(7)*(T_int(n,7)-T_int(n,8)) ...
        +area_wall_inner(9)*(T_int(n,9)-T_int(n,8))+area_wall_inner(10)*(T_int(n,10)-T_int(n,8))+area_wall_inner(11)*(T_int(n,11)-T_int(n,8))+area_wall_inner(12)*(T_int(n,12)-T_int(n,8)) ...
        +area_wall_inner(13)*(T_int(n,13)-T_int(n,8))+area_wall_inner(14)*(T_int(n,14)-T_int(n,8)))+(2*U_eastwall_garden*area_wall_east_garden+U_roof_opaque*area_roof_garden_opaque+U_roof_glazed*area_skylight)*(T_ext_interpolated(n,1)-T_int(n,8));
    
    Q_vent(n,8) = cp_air_layer*rho_air_layer*(area_room_floor(8)*height_wall)*n_vent*(T_ext_interpolated(n,1)-T_int(n,8));
    Q_rad_opaque_roof(n,8) = f_sh*(alpha_roof*R_se*U_roof_opaque*area_roof_garden_opaque)*G_sun_roof_interpolated(n,1);
    Q_rad_glazed(n,8) = f_sh*(tau_garden*area_skylight)*G_sun_roof_interpolated(n,1)+f_sh*(tau*area_wall_east_garden)*(G_sun_west_interpolated(n,1)+G_sun_east_interpolated(n,1));

    C_room(8) = 1.5*C_overall*area_room_floor(8)/total_internal_area;
    T_int(n+1,8) = T_int(n,8) + delta_time*3600/C_room(8)*(Q_vent(n,8)+Q_cond_conv(n,8)+Q_rad_opaque_roof(n,8)+Q_rad_glazed(n,8));

    % Room 9
    
    Q_cond_conv(n,9) = U_eastwall_inner*area_wall_east_room*(T_int(n,10)-T_int(n,9))+U_northwall_inner*area_wall_inner(9)*(T_int(n,8)-T_int(n,9))+(U_westwall*area_wall_east_room+U_northwall*area_wall_outer(9)+U_roof_opaque*area_room_floor(9))*(T_ext_interpolated(n,1)-T_int(n,9));
    Q_int_people(n,9) = 0*Q_people;
    Q_int_equip(n,9)= 0*Q_equip;
    Q_vent(n,9) = cp_air_layer*rho_air_layer*(area_room_floor(9)*average_height_room)*n_vent*(T_ext_interpolated(n,1)-T_int(n,9));

    Q_rad_opaque_brick(n,9) = f_sh*(alpha_brick*R_se*U_westwall*area_wall_east_room)*G_sun_west_interpolated(n,1);
    Q_rad_opaque_roof(n,9) = f_sh*(alpha_roof*R_se*U_roof_opaque*area_room_floor(9))*G_sun_roof_interpolated(n,1);

    Q_rad_glazed(n,9) = f_sh*(tau*area_wall_outer(9)*(1-f_facade_interpolated(n,1)))*G_sun_south_interpolated(n,1);

    C_room(9) = C_overall*area_room_floor(9)/total_internal_area;
    T_int(n+1,9) = T_int(n,9) + delta_time*3600/C_room(9)*(Q_vent(n,9)+Q_cond_conv(n,9)+Q_int_people(n,9)+Q_int_equip(n,9)+Q_rad_opaque_brick(n,9)+Q_rad_opaque_roof(n,9)+Q_rad_glazed(n,9));
    
    % Room 10

    Q_cond_conv(n,10) = U_eastwall_inner*area_wall_east_room*(T_int(n,9)-T_int(n,10))+U_eastwall_inner*area_wall_east_room*(T_int(n,11)-T_int(n,10))+U_northwall_inner*area_wall_inner(10)*(T_int(n,8)-T_int(n,10))+(U_northwall*area_wall_outer(10)+U_roof_opaque*area_room_floor(10))*(T_ext_interpolated(n,1)-T_int(n,10));
    Q_int_people(n,10) = 3*Q_people;
    Q_int_equip(n,10)= 3*Q_equip;
    Q_vent(n,10) = cp_air_layer*rho_air_layer*(area_room_floor(10)*average_height_room)*n_vent*(T_ext_interpolated(n,1)-T_int(n,10));

    Q_rad_opaque_roof(n,10) = f_sh*(alpha_roof*R_se*U_roof_opaque*area_room_floor(10))*G_sun_roof_interpolated(n,1);

    Q_rad_glazed(n,10) = f_sh*(tau*area_wall_outer(10)*(1-f_facade_interpolated(n,1)))*G_sun_south_interpolated(n,1);

    C_room(10) = C_overall*area_room_floor(10)/total_internal_area;
    T_int(n+1,10) = T_int(n,10) + delta_time*3600/C_room(10)*(Q_vent(n,10)+Q_cond_conv(n,10)+Q_int_people(n,10)+Q_int_equip(n,10)+Q_rad_opaque_roof(n,10)+Q_rad_glazed(n,10));

    % Room 11

    Q_cond_conv(n,11) = U_eastwall_inner*area_wall_east_room*(T_int(n,10)-T_int(n,11))+U_eastwall_inner_brick*area_wall_east_room*(T_int(n,12)-T_int(n,11))+U_northwall_inner*area_wall_inner(11)*(T_int(n,8)-T_int(n,11))+(U_northwall*area_wall_outer(11)+U_roof_opaque*area_room_floor(11))*(T_ext_interpolated(n,1)-T_int(n,11));
    Q_int_people(n,11) = 3*Q_people;
    Q_int_equip(n,11)= 3*Q_equip;
    Q_vent(n,11) = cp_air_layer*rho_air_layer*(area_room_floor(11)*average_height_room)*n_vent*(T_ext_interpolated(n,1)-T_int(n,11));

    Q_rad_opaque_roof(n,11) = f_sh*(alpha_roof*R_se*U_roof_opaque*area_room_floor(11))*G_sun_roof_interpolated(n,1);

    Q_rad_glazed(n,11) = f_sh*(tau*area_wall_outer(11)*(1-f_facade_interpolated(n,1)))*G_sun_south_interpolated(n,1);

    C_room(11) = C_overall*area_room_floor(11)/total_internal_area;
    T_int(n+1,11) = T_int(n,11) + delta_time*3600/C_room(11)*(Q_vent(n,11)+Q_cond_conv(n,11)+Q_int_people(n,11)+Q_int_equip(n,11)+Q_rad_opaque_roof(n,11)+Q_rad_glazed(n,11));

    
    % Room 12

    Q_cond_conv(n,12) = U_eastwall_inner_brick*area_wall_east_room*(T_int(n,11)-T_int(n,12))+U_eastwall_inner_brick*area_wall_east_room*(T_int(n,13)-T_int(n,12))+U_northwall_inner*area_wall_inner(12)*(T_int(n,8)-T_int(n,12))+(U_eastwall_inner_brick*area_wall_outer(12)+U_roof_opaque*area_room_floor(12))*(T_ext_interpolated(n,1)-T_int(n,12));
    Q_int_people(n,12) = 0*Q_people;
    Q_int_equip(n,12)= 0*Q_equip;
    Q_vent(n,12) = cp_air_layer*rho_air_layer*(area_room_floor(12)*average_height_room)*n_vent*(T_ext_interpolated(n,1)-T_int(n,12));

    Q_rad_opaque_roof(n,12) = f_sh*(alpha_roof*R_se*U_roof_opaque*area_room_floor(12))*G_sun_roof_interpolated(n,1);
    Q_rad_opaque_brick(n,12) = f_sh*(alpha_brick*R_se*U_westwall*area_wall_outer(12))*G_sun_south_interpolated(n,1);

    C_room(12) = C_overall*area_room_floor(12)/total_internal_area;
    T_int(n+1,12) = T_int(n,12) + delta_time*3600/C_room(12)*(Q_vent(n,12)+Q_cond_conv(n,12)+Q_int_people(n,12)+Q_int_equip(n,12)+Q_rad_opaque_roof(n,12)+Q_rad_opaque_brick(n,12));

    % Room 13

    Q_cond_conv(n,13) = U_eastwall_inner*area_wall_east_room*(T_int(n,14)-T_int(n,13))+U_eastwall_inner_brick*area_wall_east_room*(T_int(n,12)-T_int(n,13))+U_northwall_inner*area_wall_inner(13)*(T_int(n,8)-T_int(n,13))+(U_northwall*area_wall_outer(13)+U_roof_opaque*area_room_floor(13))*(T_ext_interpolated(n,1)-T_int(n,13));
    Q_int_people(n,13) = 0*Q_people;
    Q_int_equip(n,13)= 0*Q_equip;
    Q_vent(n,13) = cp_air_layer*rho_air_layer*(area_room_floor(13)*average_height_room)*n_vent*(T_ext_interpolated(n,1)-T_int(n,13));

    Q_rad_opaque_roof(n,13) = f_sh*(alpha_roof*R_se*U_roof_opaque*area_room_floor(13))*G_sun_roof_interpolated(n,1);

    Q_rad_glazed(n,13) = f_sh*(tau*area_wall_outer(13)*(1-f_facade_interpolated(n,1)))*G_sun_south_interpolated(n,1);

    C_room(13) = C_overall*area_room_floor(13)/total_internal_area;
    T_int(n+1,13) = T_int(n,13) + delta_time*3600/C_room(13)*(Q_vent(n,13)+Q_cond_conv(n,13)+Q_int_people(n,13)+Q_int_equip(n,13)+Q_rad_opaque_roof(n,13)+Q_rad_glazed(n,13));


    % Room 14

    Q_cond_conv(n,14) = U_eastwall_inner*area_wall_east_room*(T_int(n,13)-T_int(n,14))+U_northwall_inner*area_wall_inner(14)*(T_int(n,8)-T_int(n,14))+(U_westwall*area_wall_east_room+U_northwall*area_wall_outer(14)+U_roof_opaque*area_room_floor(14))*(T_ext_interpolated(n,1)-T_int(n,14));
    Q_int_people(n,14) = 3*Q_people;
    Q_int_equip(n,14)= 3*Q_equip;
    Q_vent(n,14) = cp_air_layer*rho_air_layer*(area_room_floor(14)*average_height_room)*n_vent*(T_ext_interpolated(n,1)-T_int(n,14));

    Q_rad_opaque_brick(n,14) = f_sh*(alpha_brick*R_se*U_westwall*area_wall_east_room)*G_sun_east_interpolated(n,1);
    Q_rad_opaque_roof(n,14) = f_sh*(alpha_roof*R_se*U_roof_opaque*area_room_floor(14))*G_sun_roof_interpolated(n,1);

    Q_rad_glazed(n,14) = f_sh*(tau*area_wall_outer(14)*(1-f_facade_interpolated(n,1)))*G_sun_south_interpolated(n,1);

    C_room(14) = C_overall*area_room_floor(14)/total_internal_area;
    T_int(n+1,14) = T_int(n,14) + delta_time*3600/C_room(14)*(Q_vent(n,14)+Q_cond_conv(n,14)+Q_int_people(n,14)+Q_int_equip(n,14)+Q_rad_opaque_brick(n,14)+Q_rad_opaque_roof(n,14)+Q_rad_glazed(n,14));

    
    
    %T_avg(n,1) = sum(T_int(n,:))/14; %(T_int(n,1)+T_int(n,2)+T_int(n,3)+T_int(n,4)+T_int(n,5)+T_int(n,6)+T_int(n,7)+T_int(n,8)+T_int(n,9)+T_int(n,10)+T_int(n,11)+T_int(n,12)+T_int(n,13)+T_int(n,14));
    T_avg(n,1) = (T_int(n,1)+T_int(n,2)+T_int(n,3)+T_int(n,4)+T_int(n,5)+T_int(n,6)+T_int(n,7)+T_int(n,9)+T_int(n,10)+T_int(n,11)+T_int(n,13)+T_int(n,14))/12;


end
Time=[Time, 24*days+delta_time];
T_ext_interpolated(n+1)=0.5*(T_ext_interpolated(1)+T_ext_interpolated(n));

% figure(1)
% plot(Time,T_int, LineWidth=2)
% hold on
% plot(Time,T_ext_interpolated, LineWidth=2)
% xlabel('Time [h]')
% ylabel('Temperature (°C)')
% axis([0 24*days 0 inf],'auto y')
% legend("room 1", "room 2", "room 3", "room 4", "room 5", "room 6", "room 7", "room 8", "room 9", "room 10", "room 11", "room 12", "room 13", "room 14")
% grid
% 
% figure(2)
% plot(Time(1:end-1),Q_vent, LineWidth=2)
% legend("room 1", "room 2", "room 3", "room 4", "room 5", "room 6", "room 7", "room 8", "room 9", "room 10", "room 11", "room 12", "room 13", "room 14")
% grid
% 
% figure(4)
% plot(Time(1:end-1),Q_cond_conv, LineWidth=2)
% legend("room 1", "room 2", "room 3", "room 4", "room 5", "room 6", "room 7", "room 8", "room 9", "room 10", "room 11", "room 12", "room 13", "room 14")
% grid
% 
% figure(5)
% plot(Time(1:end-1),Q_rad_glazed, LineWidth=2)
% legend("room 1", "room 2", "room 3", "room 4", "room 5", "room 6", "room 7", "room 8", "room 9", "room 10", "room 11", "room 12", "room 13", "room 14")
% grid
% 
% figure(3)
% plot(Time(1:end-1),T_avg, LineWidth=2)
% 
% figure(6)
% plot(Time(1:end-1), Q_cond_conv(:,8), LineWidth=1.5)
% hold on
% plot(Time(1:end-1), Q_int_equip(:,8),LineWidth=1.5)
% plot(Time(1:end-1), Q_int_people(:,8),LineWidth=1.5)
% plot(Time(1:end-1), Q_rad_glazed(:,8),LineWidth=1.5)
% plot(Time(1:end-1), Q_rad_opaque_roof(:,8),LineWidth=1.5)
% plot(Time(1:end-1), Q_vent(:,8),LineWidth=1.5)
% axis([0 24*days 0 inf],'auto y')
% xlabel('Time [h]')
% ylabel('Heat [W]')
% legend("Cond/Conv", "Equipment", "People", "Radiation Glazed", "Radiation Opaque", "Ventilation")
% hold off

