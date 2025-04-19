%% Isodense Contours
clear all;close all;clc;

% Define ranges
u = linspace(0.1, 5, 1000);
xi = linspace(0.1, 90, 1000);

[U, Xi] = meshgrid( u, xi );

% Parker equation
C = U.^2 - log(U.^2) - 4*log(Xi) - 4./Xi;

C_values = -6:1:4;

% Plot the contours
figure;
contour_handle = contour(Xi, U, C, C_values, 'ShowText','On', 'LineWidth', 1.5);
hold on;

plot(1, 1, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % Critical Point

text(0.2, 1.3, 'I','FontSize',16, 'FontWeight', 'bold');
text(3, 1.3, 'II','FontSize',16, 'FontWeight', 'bold');
text(1, 1.9, 'III','FontSize',16, 'FontWeight', 'bold');
text(1, 1.3, 'IV','FontSize',16, 'FontWeight', 'bold');
text(1, 0.6, 'V','FontSize',16, 'FontWeight', 'bold');

xlim([0 6]);
ylim([0 3]);
xlabel('Height \xi (normalized)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Velocity \it{u} (normalized)', 'FontSize', 16, 'FontWeight', 'bold');
title('Isodense Contours of Parker Equation Solutions', 'FontSize', 18, 'FontWeight', 'bold');
colorbar;
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.5);
grid on;
hold off;

%% Solution Choice

close all; clc;

C_target = -3;

contour_data = contour_handle;

% Find the section corresponding to C = -3
xi_vals = 0;
u_vals = 0;

i = 1;
while i < size(contour_data, 2)
    C_level = contour_data(1, i);
    n_points = contour_data(2, i);
    
    if C_level == C_target
        for j = 1:n_points
            xi_temp = contour_data(1, i + j);
            u_temp = contour_data(2, i + j);

            if (xi_temp <= 1 && u_temp <= 1) || (xi_temp > 1.01 && u_temp > 1.01)
                xi_vals = [xi_vals, xi_temp]; % Append valid values
                u_vals = [u_vals, u_temp];
            end
        end
    end
    i = i + n_points + 1; % Next contour level
end

% Select 20 points from the correct solution
num_points = min(20, length(xi_vals));
indices = round(linspace(1, length(xi_vals), num_points));
xi_selected = xi_vals(indices);
u_selected = u_vals(indices);

% Plot extracted points
figure;
plot(xi_selected, u_selected, 'ko', 'LineWidth', 1, 'MarkerSize', 8, 'MarkerFaceColor', 'r');

xlabel('Height \xi (normalized)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Velocity \it{u} (normalized)', 'FontSize', 16, 'FontWeight', 'bold');
title('Filtered Points for C = -3', 'FontSize', 18, 'FontWeight', 'bold');
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.5);
grid on;

%% Interpolation

close all; clc;

xi_interp = linspace(min(xi_selected), max(xi_selected), 500);

% Linear interpolation
u_interp = interp1(xi_selected, u_selected, xi_interp, 'linear');

figure;
plot(xi_selected, u_selected, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Original points
hold on;
plot(xi_interp, u_interp, 'k-', 'LineWidth', 2); % Interpolated curve

ylim([0 5]);
xlabel('Height \xi (normalized)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Velocity u (normalized)', 'FontSize', 16, 'FontWeight', 'bold');
title('Interpolated Parker Wind Solution', 'FontSize', 18, 'FontWeight', 'bold');
grid on;
legend('Data points', 'Linear Interpolation', 'Location', 'best');
hold off;

%% Velocity (km/s) and Distance (R_sol) for different Temperatures of the solar corona

close all; clc;

G = 6.67428E-11;
M_sol = 1.9884E30;
R_sol = 6.955E8; % [m]
k_b = 1.38E-23; % Boltzmann constant
m_p = 1.67E-27; % Proton mass

solar = char(9737);

e_orb = 149.6E9 / R_sol; % Average distance of Earth's orbit from the Sun in solar radii

T_vals = [0.5E6, 1E6, 2E6, 3E6];
colors = ['r', 'g', 'b', 'm'];

figure;
hold on;

legend_handles = [];
legend_labels = {};

for j = 1:length(T_vals)
    T = T_vals(j);
    
    % Calculate velocity
    cs = sqrt (2 * k_b * T / m_p); % sound's velocity [m/s]
    u_final = u_selected * cs / 1000; % [km/s]
    
    % Calculate distance
    rc = 0.5 * (G * M_sol / cs^2) * 1.27; % critical point where u = cs [m]
    xi_final  = xi_selected * rc / R_sol; % [solar radii]
    
    xi_interp = linspace(min(xi_final), max(xi_final), 100);
    u_interp = interp1(xi_final, u_final, xi_interp, 'spline');
    
    plot(xi_final, u_final, 'o', 'Color', colors(j), 'LineWidth', 1.5);
    
    % Interpolated curve
    h = plot(xi_interp, u_interp, '-', 'Color', colors(j), 'LineWidth', 2);
    legend_handles = [legend_handles, h]; % Store handle
    legend_labels{end+1} = ['T = ' num2str(T/1e6) 'MK'];
end

xline(e_orb, '--k', 'Orbit of Earth', 'FontSize', 14, 'FontWeight', 'bold', ...
    'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'middle');

xlim([0 250]);
xlabel(['Distance r in solar radii (R_',solar, ')'], 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Velocity \it{u} (km/s)', 'FontSize', 16, 'FontWeight', 'bold');
title('Parker Wind Solution for Different Temperatures', 'FontSize', 18, 'FontWeight', 'bold');
legend(legend_handles, legend_labels, 'Location', 'best');
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.5);
grid on;
hold off;

%% Density vs distance

close all; clc;

n_0 = 1E8; % Plasma Density in the solar corona [cm^(-3)]
r_0 = 1; % [solar radii]

T_vals = [0.5E6, 1E6, 2E6, 3E6];
colors = ['r', 'g', 'b', 'm'];

figure;
hold on;

legend_handles = [];
legend_labels = {};

for j = 1:length(T_vals)
    T = T_vals(j);
    
    % Calculate critical point
    cs = sqrt (2 * k_b * T / m_p); % [m/s]
    rc = 0.5 * (G * M_sol / cs^2) * 1.27; % [m]

    % Starting velocity
    k = 4 * log(R_sol / rc) + 4 * (rc / R_sol) - 3;
    u_0_norm = exp(-k / 2); % normalised
    
    fprintf('T = %.1f MK -> u_0 = %.4f\n', T / 1e6, u_0_norm);  
    
    % Calculate distance
    xi_final = xi_selected * rc; % [m]   
    
    % Calculate vleocity
    u_final = u_selected * cs; % [m/s]
    u_0 = u_0_norm * cs; % [m/s]
   
    % Density
    n_final = n_0 * (R_sol ./ xi_final).^2 .* (u_0 ./ u_final);
    
    xi_interp = linspace(min(xi_final), max(xi_final), 100);
    n_interp = interp1(xi_final, n_final, xi_interp, 'spline', 'extrap');
    
    plot(xi_final(1:3) / R_sol, n_final(1:3), 'o', 'Color', colors(j), 'LineWidth', 2);
    plot(xi_final(3:end) / R_sol, n_final(3:end), 'o-', 'Color', colors(j), 'LineWidth', 2);
    
    % Interpolated curve
    h = plot(xi_interp, n_interp, '-', 'Color', colors(j), 'LineWidth', 2);
    legend_handles = [legend_handles, h];
    legend_labels{end+1} = ['T = ' num2str(T/1e6) 'MK'];
end

xline(e_orb, '--k', 'Orbit of Earth', 'FontSize', 14, 'FontWeight', 'bold', ...
    'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');

set(gca, 'YScale', 'log');
xlim([0 250]);
xlabel(['Distance r in solar radii (R_',solar, ')'], 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Density  \it{n} (cm^{-3})', 'FontSize', 16, 'FontWeight', 'bold');
title('Density of Solar Wind for Different Temperatures', 'FontSize', 18, 'FontWeight', 'bold');
legend(legend_handles, legend_labels, 'Location', 'best');
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.5);
grid on;
hold off;

%% Save density and frequency values up to any distance in solar radii

close all; clc;

r_desired = 1:9:250;

for j = 1:length(T_vals)
    T = T_vals(j);
    
    cs = sqrt (2 * k_b * T / m_p); % [m/s]
    
    % Calculate critical point
    rc = 0.5 * (G * M_sol / cs^2) * 1.27; % [m]
    xi_final  = xi_selected * rc; % [m]
    
    % Starting velocity
    k = 4. * log(R_sol / rc) + 4. * rc / R_sol - 3;
    u_0 = exp(-k / 2);

    % Calculate density
    n_final = n_0 * (R_sol ./ xi_final).^2 .* (u_0 ./ u_selected);
        
    r_filtered = xi_final / R_sol; % Distance in solar radii
        
    n_at_r = interp1(r_filtered, n_final, r_desired, 'linear','extrap');
    
    % Calculate plasma frequency
    fp = 8.98E3 * sqrt(n_at_r); % [Hz]

    filename = sprintf('density_freq_profile_in_250_rsol_T%.1fMK.txt', T / 1e6);
    fileID = fopen(filename, 'w');
    fprintf(fileID, 'Distance (R_solar)\t Density (cm^-3) \t Plasma Frequency (Hz)\n');
    
    for i = 1:length(r_desired)
        fprintf(fileID, '%.2f \t\t %.2e \t\t %.2e\n', r_desired(i), n_at_r(i), fp(i));
    end

    fclose(fileID);
    fprintf('Data saved to %s\n\n', filename);    
end

%% Atmospheric models (Newkirk, Saito, Leblanc and Vrsnak)

close all; clc;

r_desired = 1:10:250; % [solar redii]
    
    % Density in cm^(-3)
    n_newk = 4.2 * 1E4 * exp(9.95 ./ r_desired); % Newkirk

    n_saito = 1E8 * (3.09 * r_desired.^(-16) + 1.58 * r_desired.^(-6) + 0.0251 * r_desired.^(-2.5)); % Saito for ö = 0

    n_leb = 1E5 * (3.3 * r_desired.^(-2) + 41 * r_desired.^(-4) + 800 * r_desired.^(-6)); % Leblanc

    n_vr = 1E8 * (0.0033 * r_desired.^(-2) + r_desired.^(-4) + 3.16 * r_desired.^(-6) + 15.45 * r_desired.^(-16)); % Vrsnak

    % Calculate plasma frequency
    fp_newk = 8.98E3 * sqrt(n_newk);
    
    fp_saito = 8.98E3 * sqrt(n_saito);
    
    fp_leb = 8.98E3 * sqrt(n_leb);
    
    fp_vr = 8.98E3 * sqrt(n_vr);
    
    filename = sprintf('models.txt');
    fileID = fopen(filename, 'w');
    fprintf(fileID, 'Distance (R_solar)\t Newkirk (cm^-3) \t Saito (cm^-3) \t Leblanc (cm^-3) \t Vrsnak (cm^-3) \t\n');
    
    for i = 1:length(r_desired)
        fprintf(fileID, '%.2f \t\t %.2e \t\t %.2e \t\t %.2e \t\t %.2e\n', r_desired(i), n_newk(i), n_saito(i), n_leb(i), n_vr(i));
    end
    
    fprintf(fileID, '\n Distance (R_solar)\t Newkirk (Hz) \t Saito (Hz) \t Leblanc (Hz) \t Vrsnak (Hz) \t\n');
    
    for i = 1:length(r_desired)
        fprintf(fileID, '%.2f \t\t %.2e \t\t %.2e \t\t %.2e \t\t %.2e\n', r_desired(i), fp_newk(i), fp_saito(i), fp_leb(i), fp_vr(i));
    end

    fclose(fileID);
    fprintf('Data saved to %s\n\n', filename);
    
%% Frequency vs distance for each temperature

close all; clc;

T_vals = [0.5E6, 1E6, 2E6, 3E6];
colors = ['r', 'g', 'b', 'm'];

figure;
hold on;
legend_labels = {};

for j = 1:length(T_vals)
    T = T_vals(j);
    
    filename = sprintf('density_freq_profile_in_250_rsol_T%.1fMK.txt', T / 1e6);
    data = dlmread(filename, '', 1, 0);  % Skip 1 row

    r_vals = data(:,1);       % Distance (R_solar)
    freq_vals = data(:,3);    % Frequency (Hz)
    
    plot(r_vals, freq_vals, 'o-', 'Color', colors(j), 'LineWidth', 2);
    legend_labels{end+1} = sprintf('T = %.1f MK', T / 1e6);
end

set(gca, 'YScale', 'log');
xlabel(['Distance r (R_',solar, ')'], 'FontSize', 14, 'FontWeight', 'bold');
ylabel('log_{10} f_p (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
title('Plasma Frequency vs Distance for Different Temperatures', 'FontSize', 16, 'FontWeight', 'bold');
legend(legend_labels, 'Location', 'best');
grid on;
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'LineWidth', 1.2);
xlim([0 250]);
hold off;

%% Plot density (+ theoritical models) vs distance for T = 2 MK

close all; clc;

sim_filename = 'density_freq_profile_in_250_rsol_T2.0MK.txt';
sim_data = dlmread(sim_filename, '', 1, 0);

r_sim = sim_data(:,1);      % Distance (R_solar)
n_sim = sim_data(:,2);      % Denisty (cm^(-3))

% Load theoretical model data
theo_filename = 'models_d.txt';
theo_data = dlmread(theo_filename, '', 1, 0);
% data_limited = theo_data(1:24, :);

r_theo = theo_data(:,1);
n_newk = theo_data(:,2);
n_saito = theo_data(:,3);
n_leblanc = theo_data(:,4);
n_vrsnak = theo_data(:,5);

figure;
hold on;

plot(r_sim, n_sim, 'k-o', 'LineWidth', 3, 'DisplayName', 'T = 2 MK (Simulated)');

plot(r_theo, n_newk, '--r', 'LineWidth', 2, 'DisplayName', 'Newkirk');
plot(r_theo, n_saito, '--g', 'LineWidth', 2, 'DisplayName', 'Saito');
plot(r_theo, n_leblanc, '--b', 'LineWidth', 2, 'DisplayName', 'Leblanc');
plot(r_theo, n_vrsnak, '--m', 'LineWidth', 2, 'DisplayName', 'Vrsnak');

set(gca, 'YScale', 'log');
xlabel(['Distance r (R_',solar, ')'], 'FontSize', 14, 'FontWeight', 'bold');
ylabel('log_{10} n (cm^{-3})', 'FontSize', 14, 'FontWeight', 'bold');
title('Simulated vs Theoretical Solar Wind Densities at T = 2 MK', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'northeast');
grid on;
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'LineWidth', 1.2);
xlim([0 250]);

hold off;

%% Plot frequency (+ theoretical models) vs distance for T = 2 MK

close all; clc;

sim_filename = 'density_freq_profile_in_250_rsol_T2.0MK.txt';
sim_data = dlmread(sim_filename, '', 1, 0);

r_sim = sim_data(:,1);      % Distance (R_solar)
f_sim = sim_data(:,3);      % Frequency (Hz)

% Load theoretical model data
theo_filename = 'models.txt';
theo_data = dlmread(theo_filename, '', 28, 0);

r_theo = theo_data(:,1);
f_newk = theo_data(:,2);
f_saito = theo_data(:,3);
f_leblanc = theo_data(:,4);
f_vrsnak = theo_data(:,5);

figure;
hold on;

plot(r_sim, f_sim, 'k-o', 'LineWidth', 3, 'DisplayName', 'T = 2 MK (Simulated)');

plot(r_theo, f_newk, '--r', 'LineWidth', 2, 'DisplayName', 'Newkirk');
plot(r_theo, f_saito, '--g', 'LineWidth', 2, 'DisplayName', 'Saito');
plot(r_theo, f_leblanc, '--b', 'LineWidth', 2, 'DisplayName', 'Leblanc');
plot(r_theo, f_vrsnak, '--m', 'LineWidth', 2, 'DisplayName', 'Vrsnak');

set(gca, 'YScale', 'log');
xlabel(['Distance r (R_',solar, ')'], 'FontSize', 14, 'FontWeight', 'bold');
ylabel('log_{10} f_p (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
title('Simulated vs Theoretical Plasma Frequencies at T = 2 MK', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'northeast');
grid on;
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'LineWidth', 1.2);
xlim([0 250]);

hold off;