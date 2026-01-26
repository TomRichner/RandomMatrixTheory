
close all
clear all
clc

rng(2)

set(groot, 'DefaultAxesFontSize', 16);
set(groot, 'DefaultTextFontSize', 16);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultAxesLineWidth', 1.5);

%% Struct to hold RMT objects
G = struct();

%% 1. G(1): n=600, density=50%
g = 1;
n = 600;
G(g).rmt = RMT2(n);
G(g).rmt.set_mean_parameters(1/sqrt(600)); % mu_E = 1. mu_I calculated automatically.
G(g).rmt.set_f(0.5);
b_E = 0.5/sqrt(600);
b_eff = 0.5/sqrt(600);
b_I = G(g).rmt.compute_sigma_I_from_eff(b_eff, b_E);
G(g).rmt.set_stdev_parameters(b_E, b_I);
G(g).rmt.set_density(1);
G(g).rmt.description = 'density 100%';
G(g).rmt.display_parameters();

%% 2.
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
% G(g).rmt.set_mu(0);
G(g).rmt.set_density(.75);
G(g).rmt.description = 'density 75%';
G(g).rmt.display_parameters();

%% 3
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
% G(g).rmt.set_mu(0);
G(g).rmt.set_density(.5);
G(g).rmt.description = 'density 50%';
G(g).rmt.display_parameters();

%% 4
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
% G(g).rmt.set_mu(0);
G(g).rmt.set_density(.25);
G(g).rmt.description = 'density 25%';
G(g).rmt.display_parameters();

%% 1. G(1): n=600, density=50%
g = g + 1;
n = 300;
G(g).rmt = RMT2(n);
G(g).rmt.set_mean_parameters(1/sqrt(600)); % mu_E = 1. mu_I calculated automatically.
G(g).rmt.set_f(0.5);
b_E = 0.5/sqrt(600);
b_eff = 0.5/sqrt(600);
b_I = G(g).rmt.compute_sigma_I_from_eff(b_eff, b_E);
G(g).rmt.set_stdev_parameters(b_E, b_I);
G(g).rmt.set_density(0.50);
G(g).rmt.description = '300, 50%';
G(g).rmt.display_parameters();

%% 1. G(1): n=600, density=50%
g = g + 1;
n = 600;
G(g).rmt = RMT2(n);
G(g).rmt.set_mean_parameters(1/sqrt(600)); % mu_E = 1. mu_I calculated automatically.
G(g).rmt.set_f(0.5);
b_E = 0.5/sqrt(600);
b_eff = 0.5/sqrt(600);
b_I = G(g).rmt.compute_sigma_I_from_eff(b_eff, b_E);
G(g).rmt.set_stdev_parameters(b_E, b_I);
G(g).rmt.set_density(0.25);
G(g).rmt.description = '600, 25%';
G(g).rmt.display_parameters();

%% 1. G(1): n=600, density=50%
g = g + 1;
n = 1200;
G(g).rmt = RMT2(n);
G(g).rmt.set_mean_parameters(1/sqrt(600)); % mu_E = 1. mu_I calculated automatically.
G(g).rmt.set_f(0.5);
b_E = 0.5/sqrt(600);
b_eff = 0.5/sqrt(600);
b_I = G(g).rmt.compute_sigma_I_from_eff(b_eff, b_E);
G(g).rmt.set_stdev_parameters(b_E, b_I);
G(g).rmt.set_density(0.125);
G(g).rmt.description = '1200 12.5%';
G(g).rmt.display_parameters();

%% 1. G(1): n=600, density=50%
g = g + 1;
n = 2400;
G(g).rmt = RMT2(n);
G(g).rmt.set_mean_parameters(1/sqrt(600)); % mu_E = 1. mu_I calculated automatically.
G(g).rmt.set_f(0.5);
b_E = 0.5/sqrt(600);
b_eff = 0.5/sqrt(600);
b_I = G(g).rmt.compute_sigma_I_from_eff(b_eff, b_E);
G(g).rmt.set_stdev_parameters(b_E, b_I);
G(g).rmt.set_density(0.125/2);
G(g).rmt.description = '2400 12.5/2%';
G(g).rmt.display_parameters();

% %% 2. G(2): n=1200, density=25%
% g = 2;
% n = 600;
% G(g).rmt = RMT2(n);
% G(g).rmt.set_mean_parameters(1/sqrt(600)); % mu_E = 1. mu_I calculated automatically.
% G(g).rmt.set_f(0.5);
% G(g).rmt.set_stdev_parameters(b_E, b_I);
% G(g).rmt.set_density(0.5);
% G(g).rmt.description = sprintf('n=%d, p=0.25', n);
% G(g).rmt.display_parameters();

%% Compute and Plot
f1 = figure(1);
set(f1, 'Position', [100   200   990   600], 'Color', 'white');
t = tiledlayout(2, ceil(g/2), 'TileSpacing', 'compact', 'Padding', 'compact');

ax = gobjects(length(G), 1);

fprintf('\nExpected Radius vs Computed Spectral Radius:\n');

for i = 1:length(G)
    % Compute eigenvalues
    G(i).rmt.compute_eigenvalues();

    % Calculate expected radius
    r = G(i).rmt.compute_expected_radius();

    % Get max abs eigenvalue
    max_lambda = max(abs(G(i).rmt.eigenvalues));
    fprintf('%s: Expected r = %.4f, Max |lambda| = %.4f\n', G(i).rmt.description, r, max_lambda);

    % Set plot circle (center is determined by shift)
    % Note: shift is 0 unless set otherwise. The example code had "shift = -r" commented out or set.
    % The example code L27-39 does NOT set shift, but later L46 does. We stick to L27-39.
    x_center = G(i).rmt.shift;
    G(i).rmt.set_plot_circle(r, x_center);

    % Enable outlier coloring
    G(i).rmt.set_color_outliers(true);
    G(i).rmt.set_outlier_factor(1.04);

    % Plot
    ax(i) = nexttile;
    G(i).rmt.plot_circle(ax(i));
    hold on;
    G(i).rmt.plot_eigenvalue_distribution(ax(i));
    hold off;

    title(G(i).rmt.description, 'Interpreter', 'none', 'FontWeight', 'bold');
end

% Determine global scale
max_radius_x = 0;
max_radius_y = 0;

for i = 1:length(G)
    axis(ax(i), 'normal');
    axis(ax(i), 'tight');
    current_xlim = xlim(ax(i));
    current_ylim = ylim(ax(i));
    dist_x = max(abs(current_xlim - G(i).rmt.shift));
    dist_y = max(abs(current_ylim));
    max_radius_x = max(max_radius_x, dist_x);
    max_radius_y = max(max_radius_y, dist_y);
end

margin = 1.1;
common_radius_x = max_radius_x * margin;
common_radius_y = max_radius_y * margin;

for i = 1:length(G)
    center_x = G(i).rmt.shift;
    xlim(ax(i), [center_x - common_radius_x, center_x + common_radius_x]);
    ylim(ax(i), [-common_radius_y, common_radius_y]);
    daspect(ax(i), [1 1 1]);

    axes(ax(i));
    x_limit = xlim;
    y_limit = ylim;
    axis off;
    hold on;
    plot(x_limit, [0,0], 'k', 'LineWidth', 1.5);
    plot([0,0], y_limit, 'k', 'LineWidth', 1.5);
    hold off;
end
