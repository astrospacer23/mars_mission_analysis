clc; clear; close all;

% Input Parameters %

R = 10; % Radius of Spherical Airship in m
r = 0.05; % Radius of cylindrical airbeam in m
n = 5; % number of airbeams
rho = 0.018; % martian atmospheric density in kg/m^3
lambda_sphere = 0.025; % membrane material density in kg/m^3
lambda_airbeam = 0.15; % airbeam density in kg/m^3
rho_gas = 0.045; % hydrogen gas density in kg/m^3
g_mars = 3.71; % gravitational acceleration on mars in m/s^2

% Output Parameters %

% 1.Sphere Geometry %

V = (4/3) * pi * (R^3); % volume of sphere in m^3
S = 4 * pi *(R^2); % surface area of sphere in m^2


% 2.Airbeam Geometry %

l = pi * R; % length of airbeam in m
S_airbeam = 2 * (pi^2) * r * R; % surface area per airbeam in m^2
S_total_airbeam = n * S_airbeam; % TSA of airbeam in m^2
V_airbeam = ((pi*r)^2) * R; % Volume per airbeam
V_total_airbeam = n * V_airbeam; % total airbeam volume

% 3.Mass Calculations %

M_sphere = lambda_sphere * S; % mass of sphere membrane in kg
M_airbeam = lambda_airbeam * S_total_airbeam; % mass of airbeams in kg
M_gas = rho_gas * V_total_airbeam; % inflated gas mass in kg
M_total = M_sphere + M_airbeam + M_gas; % total mass in kg

% 4.Buoyancy and Lift %

L = (rho - rho_gas) * g_mars * V; % lift 
Liftable_mass = (rho - rho_gas) * (V) - (M_total); % liftable mass in kg

% 5.Display Results %

fprintf('=== MARTIAN AIRSHIP DESIGN CALCULATIONS ===\n\n');

fprintf('INPUT PARAMETERS:\n');
fprintf('Sphere Radius (R): %.2f m\n', R);
fprintf('Airbeam Radius (r): %.3f m\n', r);
fprintf('Number of Airbeams (n): %d\n', n);
fprintf('Martian Atmospheric Density: %.3f kg/m³\n', rho);
fprintf('Membrane Material Density: %.3f kg/m²\n', lambda_sphere);
fprintf('Airbeam Density: %.1f kg/m²\n', lambda_airbeam);
fprintf('Hydrogen Gas Density: %.3f kg/m³\n', rho_gas);
fprintf('Mars Gravity: %.2f m/s²\n\n', g_mars);

fprintf('SPHERE GEOMETRY:\n');
fprintf('Volume (V): %.2f m³\n', V);
fprintf('Surface Area (S): %.2f m²\n\n', S);

fprintf('AIRBEAM GEOMETRY:\n');
fprintf('Airbeam Length (l): %.2f m\n', l);
fprintf('Surface Area per Airbeam: %.2f m²\n', S_airbeam);
fprintf('Total Airbeam Surface Area: %.2f m²\n', S_total_airbeam);
fprintf('Volume per Airbeam: %.4f m³\n', V_airbeam);
fprintf('Total Airbeam Volume: %.4f m³\n\n', V_total_airbeam);

fprintf('MASS CALCULATIONS:\n');
fprintf('Sphere Membrane Mass: %.2f kg\n', M_sphere);
fprintf('Airbeam Mass: %.2f kg\n', M_airbeam);
fprintf('Gas Mass: %.4f kg\n', M_gas);
fprintf('Total Mass: %.2f kg\n\n', M_total);

fprintf('BUOYANCY AND LIFT:\n');
fprintf('Total Lift Force: %.2f N\n', L);
fprintf('Liftable Mass: %.2f kg\n\n', Liftable_mass);

