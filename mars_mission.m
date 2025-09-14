clc; clear; close all;

%% === 1. Constants ===
mu_earth = 3.986004418e5; % km^3/s^2
mu_mars = 4.282837e4;     % km^3/s^2
R_earth = 6378.1;         % km
R_mars = 3389.5;          % km

% J2
J2_earth = 1.08263e-3;
J2_mars = 1.96045e-3;

% Solar
AU = 1.495978707e8; % km
mu_sun = 1.32712440018e11;

% Vehicle
mass = 1000; % kg
Cd = 1.8;
A = 10; % m^2

%% === 2. Earth Departure ===
r_leo = R_earth + 200; % km
v_circ = sqrt(mu_earth / r_leo);
C3 = 12; % km^2/s^2
v_inf = sqrt(C3);
v_hyper = sqrt(v_inf^2 + 2 * mu_earth / r_leo);
delta_v_tmi = v_hyper - v_circ;
fprintf('TMI Δv = %.2f km/s\n', delta_v_tmi);

%% === 3. Hohmann Transfer ===
r1 = AU;
r2 = 1.524 * AU;
a_trans = (r1 + r2) / 2;
T_trans = pi * sqrt(a_trans^3 / mu_sun);
fprintf('Transfer Time: %.1f days\n', T_trans / 86400);

%% === 4. Mars Arrival & MOI ===
r_p = R_mars + 250;
r_a = R_mars + 10000;
a_m = (r_p + r_a) / 2;
v_circ_mars = sqrt(mu_mars / r_p);
v_hyper_mars = sqrt(v_inf^2 + 2 * mu_mars / r_p);
v_elliptic = sqrt(2 * (mu_mars / r_p - mu_mars / (2 * a_m)));
delta_v_moi = v_hyper_mars - v_elliptic;
fprintf('MOI Δv = %.2f km/s\n', delta_v_moi);

%% === 5. Orbit Propagation with Mars J2 ===
oe0 = [a_m; (r_a - r_p)/(r_a + r_p); deg2rad(30); 0; 0; 0];
[r0, v0] = oe2rv(mu_mars, oe0);
tspan = [0 10*24*3600]; % 10 days
y0 = [r0; v0];

[t_orb, y_orb] = ode45(@(t,y) twobody_J2_mars(t, y, mu_mars, R_mars, J2_mars), tspan, y0);
figure; plot3(y_orb(:,1), y_orb(:,2), y_orb(:,3));
xlabel('X'); ylabel('Y'); zlabel('Z'); title('Mars Orbit with J2');

%% === 6. Atmospheric Entry with Drag ===
h_entry = 125;      % km
v_entry = 5.8;      % km/s
gamma = -12 * pi/180; % entry angle
alt_hover = 10;     % km target
tspan_entry = [0 300]; % sec

% Entry conditions
r_entry = R_mars + h_entry;
vx_entry = v_entry * cos(gamma);
vz_entry = -v_entry * sin(gamma);
y0_entry = [0; r_entry; vx_entry; vz_entry; h_entry]; % x, z, vx, vz, h

[t_entry, Y_entry] = ode45(@(t,y) entry_with_drag(t,y,mass,Cd,A,R_mars,mu_mars), tspan_entry, y0_entry);

% Plot trajectory
figure;
plot(Y_entry(:,1)/1e3, Y_entry(:,2)/1e3);
xlabel('Downrange [km]'); ylabel('Altitude [km]');
title('Atmospheric Entry Trajectory');

%% === 7. Stable Hovering at Target Altitude with PID ===
alt_target = alt_hover;   % km
hover_time = 120;         % seconds
g_mars = mu_mars / (R_mars + alt_target)^2; % km/s²
mass_kg = mass;           % confirm consistent use

% PID Gains (tuned)
kP = 0.5;
kI = 0.02;
kD = 0.4;

% Time loop
dt = 0.1;
t_hover = 0:dt:hover_time;
N = length(t_hover);

% Initialization
alt = alt_target + 0.05; % slight offset to simulate descent
vel = -0.02;             % small descent rate
err_int = 0;
err_prev = alt - alt_target;
alt_log = zeros(1, N);
vel_log = zeros(1, N);
thrust_log = zeros(1, N);

for i = 1:N
    err = alt - alt_target;
    derr = (err - err_prev) / dt;
    err_int = err_int + err * dt;

    % PID output (acceleration command)
    a_cmd = - (kP*err + kI*err_int + kD*derr);

    % Compute required thrust
    acc_total = g_mars + a_cmd;
    thrust = mass_kg * acc_total; % in km·kg/s²

    % Safety clamp (no negative thrust, limit to 2g)
    thrust = min(max(thrust, 0), 2 * mass_kg * g_mars);

    % Update dynamics
    acc = (thrust / mass_kg) - g_mars; % net acceleration (km/s²)
    vel = vel + acc * dt;
    alt = alt + vel * dt;

    % Store logs
    err_prev = err;
    alt_log(i) = alt;
    vel_log(i) = vel;
    thrust_log(i) = thrust;
end

% Plot Hover Performance
figure;
plot(t_hover, alt_log, 'b-', 'LineWidth', 1.5); hold on;
yline(alt_target, 'r--', 'Target Altitude');
xlabel('Time [s]'); ylabel('Altitude [km]');
title('PID Hovering Altitude Control');
legend('Altitude','Target');
grid on;


%% === Helper Functions ===

function dydt = twobody_J2_mars(~, y, mu, R, J2)
    r = y(1:3); v = y(4:6); r_norm = norm(r);
    a_tb = -mu * r / r_norm^3;
    z2 = r(3)^2; r2 = r_norm^2;
    tx = r(1)/r_norm*(5*z2/r2 - 1);
    ty = r(2)/r_norm*(5*z2/r2 - 1);
    tz = r(3)/r_norm*(5*z2/r2 - 3);
    a_j2 = 1.5 * J2 * mu * R^2 / r_norm^4 * [tx; ty; tz];
    dydt = [v; a_tb + a_j2];
end

function [r,v] = oe2rv(mu, oe)
    a = oe(1); e = oe(2); i = oe(3); RAAN = oe(4); omega = oe(5); theta = oe(6);
    p = a * (1 - e^2);
    r_pf = [p*cos(theta)/(1+e*cos(theta)); p*sin(theta)/(1+e*cos(theta)); 0];
    v_pf = [-sqrt(mu/p)*sin(theta); sqrt(mu/p)*(e + cos(theta)); 0];
    R3_W = [cos(RAAN) sin(RAAN) 0; -sin(RAAN) cos(RAAN) 0; 0 0 1];
    R1_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
    R3_w = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];
    Q_pX = (R3_W') * (R1_i') * (R3_w');
    r = Q_pX * r_pf; v = Q_pX * v_pf;
end

function dydt = entry_with_drag(~, y, ~, Cd, A, R, mu)
    x = y(1); z = y(2); vx = y(3); vz = y(4); h = y(5);
    r = sqrt(x^2 + z^2); v = sqrt(vx^2 + vz^2);
    alt = r - R;
    rho = mars_density(alt);
    drag = 0.5 * rho * v^2 * Cd * A / 1000; % km/s²
    ax = -drag * vx/v; az = -mu*z/r^3 - drag * vz/v;
    dydt = [vx; vz; ax; az; alt];
end

function rho = mars_density(alt)
    % GRAM-like exponential model
    if alt > 125
        rho = 0;
    elseif alt > 80
        rho = 1.2e-9 * exp(-(alt - 80)/6);
    elseif alt > 50
        rho = 5e-8 * exp(-(alt - 50)/8);
    elseif alt > 0
        rho = 1e-6 * exp(-alt/10);
    else
        rho = 1e-5; % ground-level
    end
end
%% === 6B. 3D Animation of Atmospheric Entry ===
% Extract position components
x_km = Y_entry(:,1);          % Downrange x in km
z_km = Y_entry(:,2);          % Altitude + Mars radius (z) in km
y_km = zeros(size(x_km));     % Assuming planar motion for now

% Convert altitude (z) from radius to height above Mars surface
z_km_alt = z_km - R_mars;

figure;
hold on;
comet3(x_km, y_km, z_km_alt);
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Altitude [km]');
title('3D Animation of Atmospheric Entry');
grid on;
view(45, 25);

% Plot Mars surface sphere for reference
[xx, yy, zz] = sphere(50);
surf(xx*R_mars, yy*R_mars, zz*R_mars, ...
     'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [1 0.4 0.1]);
%% === 8. Δv Budget and Mission Summary ===
delta_v_total = delta_v_tmi + delta_v_moi; % km/s
fprintf('\n===== MISSION SUMMARY =====\n');
fprintf('TMI Delta-V:        %.2f km/s\n', delta_v_tmi);
fprintf('Transfer Time:      %.1f days\n', T_trans/86400);
fprintf('MOI Delta-V:        %.2f km/s\n', delta_v_moi);
fprintf('Total Delta-V:      %.2f km/s\n', delta_v_total);
fprintf('Hover Altitude:     %.1f km\n', alt_hover);
fprintf('Final Hover Error:  %.3f km\n', abs(alt_log(end) - alt_hover));

if abs(alt_log(end) - alt_hover) < 0.01
    fprintf('✅ Mission Success: Stable hover achieved at target altitude!\n');
else
    fprintf('⚠️ Mission Warning: Hover deviation exceeds tolerance.\n');
end

%% === 9. Exporting Figures ===
fprintf('\nExporting figures...\n');
figs = findall(groot, 'Type', 'figure');
labels = {'Mars_J2_Orbit', 'Entry_Profile', 'Entry_3D', 'Hover_PID'};

for k = 1:min(length(figs), length(labels))
    saveas(figs(k), ['MarsMission_', labels{k}, '.png']);
end


