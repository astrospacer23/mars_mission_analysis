% Plotting script for Mars mission simulation variables
% Plots: alt_log, Y_entry, t_entry, t_hover, t_orb, thrust_log, vel_log, 
% x_km, xx, y_km, yy, y_orb, z_km, z_km_alt, zz

clc; % Clear command window
close all; % Close existing figures

% Note: This script assumes all variables (alt_log, Y_entry, etc.) are 
% available in the workspace from running the original simulation code.
% Run the original code first to generate these variables.

%% 1. Hover Phase Plots: alt_log, t_hover, thrust_log, vel_log
figure('Name', 'Hover Phase Dynamics');
% Subplot 1: Altitude vs. Time
subplot(2, 1, 1);
plot(t_hover, alt_log, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Altitude');
hold on;
yline(alt_hover, 'r--', 'Target Altitude', 'LineWidth', 1.5, 'DisplayName', 'Target');
xlabel('Time [s]');
ylabel('Altitude [km]');
title('Hover Altitude Control');
legend('show');
grid on;

% Subplot 2: Thrust and Velocity vs. Time
subplot(2, 1, 2);
yyaxis left;
plot(t_hover, thrust_log, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Thrust');
ylabel('Thrust [km·kg/s²]');
yyaxis right;
plot(t_hover, vel_log, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Velocity');
ylabel('Velocity [km/s]');
xlabel('Time [s]');
title('Thrust and Velocity During Hover');
legend('show');
grid on;

%% 2. Atmospheric Entry Plots: Y_entry, t_entry, x_km, y_km, z_km, z_km_alt
% Y_entry contains [x, z, vx, vz, h]
% x_km, y_km, z_km, z_km_alt are derived from Y_entry
figure('Name', 'Atmospheric Entry Trajectories');
% Subplot 1: Downrange (x_km) vs. Altitude (z_km_alt)
subplot(2, 1, 1);
plot(x_km/1e3, z_km_alt/1e3, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Trajectory');
xlabel('Downrange [km]');
ylabel('Altitude [km]');
title('Entry Trajectory: Downrange vs. Altitude');
legend('show');
grid on;

% Subplot 2: Altitude (z_km_alt) vs. Time (t_entry)
subplot(2, 1, 2);
plot(t_entry, z_km_alt/1e3, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Altitude');
xlabel('Time [s]');
ylabel('Altitude [km]');
title('Altitude vs. Time During Entry');
legend('show');
grid on;

% 3D Plot: x_km, y_km, z_km_alt
figure('Name', '3D Atmospheric Entry Trajectory');
plot3(x_km/1e3, y_km/1e3, z_km_alt/1e3, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Entry Path');
hold on;
xlabel('Downrange X [km]');
ylabel('Crossrange Y [km]');
zlabel('Altitude [km]');
title('3D Atmospheric Entry Trajectory');
legend('show');
grid on;

%% 3. Mars Orbit Propagation: t_orb, y_orb
% y_orb contains [x, y, z, vx, vy, vz] in Mars inertial frame
figure('Name', 'Mars Orbit with J2');
plot3(y_orb(:,1)/1e3, y_orb(:,2)/1e3, y_orb(:,3)/1e3, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Orbit');
hold on;
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('Mars Orbit with J2 Perturbation (10 Days)');
legend('show');
grid on;
% Plot initial point
plot3(y_orb(1,1)/1e3, y_orb(1,2)/1e3, y_orb(1,3)/1e3, 'ro', 'MarkerSize', 5, 'DisplayName', 'Start');

%% 4. Mars Surface Sphere: xx, yy, zz
% xx, yy, zz are from sphere(50) for Mars surface visualization
figure('Name', 'Mars Surface and Entry Trajectory');
% Plot Mars surface
surf(xx*R_mars, yy*R_mars, zz*R_mars, ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [1 0.4 0.1], ...
    'DisplayName', 'Mars Surface');
hold on;
% Overlay entry trajectory for context
plot3(x_km/1e3, y_km/1e3, z_km_alt/1e3, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Entry Path');
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Altitude [km]');
title('Mars Surface with Entry Trajectory');
legend('show');
grid on;
view(45, 25);

%% 5. Additional Entry Data from Y_entry
% Y_entry = [x, z, vx, vz, h]
figure('Name', 'Entry Dynamics from Y_entry');
% Subplot 1: Position Components (x, z) vs. Time
subplot(2, 1, 1);
plot(t_entry, Y_entry(:,1)/1e3, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Downrange X');
hold on;
plot(t_entry, Y_entry(:,2)/1e3, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Radial Distance Z');
xlabel('Time [s]');
ylabel('Distance [km]');
title('Entry Position Components');
legend('show');
grid on;

% Subplot 2: Velocity Components (vx, vz) vs. Time
subplot(2, 1, 2);
plot(t_entry, Y_entry(:,3), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Velocity X');
hold on;
plot(t_entry, Y_entry(:,4), 'm-', 'LineWidth', 1.5, 'DisplayName', 'Velocity Z');
xlabel('Time [s]');
ylabel('Velocity [km/s]');
title('Entry Velocity Components');
legend('show');
grid on;

% Save all figures
fprintf('Saving figures...\n');
figs = findall(groot, 'Type', 'figure');
labels = {'Hover_Dynamics', 'Entry_Trajectories', 'Entry_3D', 'Mars_Orbit', ...
          'Mars_Surface', 'Entry_Dynamics'};
for k = 1:min(length(figs), length(labels))
    saveas(figs(length(figs)-k+1), ['MarsMission_', labels{k}, '.png']);
end
fprintf('Figures saved successfully.\n');