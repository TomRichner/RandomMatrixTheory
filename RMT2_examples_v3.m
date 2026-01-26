close all
clear all

clc

% if ~exist('seed', 'var')
seed = 1
% else
%     seed = seed + 1
% end
% clearvars -except seed

rng(seed)

set(groot, 'DefaultAxesFontSize', 16);
set(groot, 'DefaultTextFontSize', 16);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultAxesLineWidth', 1.5);

save_figs = false;

n = 600;

%% Struct to hold RMT objects
G = struct();

%% 1. Standard Circular Law with outlier and shift
g = 1;
G(g).rmt = RMT2(n);
b_E = 1/sqrt(n);
b_eff = 1/sqrt(n);
% G(g).rmt.set_mu(1.2*b_eff);
% Calculate b_I using the method from RMT2
b_I = G(g).rmt.compute_sigma_I_from_eff(b_eff, b_E);
G(g).rmt.set_stdev_parameters(b_E, b_I);
% r = G(g).rmt.compute_expected_radius();
% G(g).rmt.shift = -r;
G(g).rmt.description = 'Standard Circular Law with outlier and shift';
G(g).rmt.display_parameters();

%% 2. Shifted
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
% G(g).rmt.set_mu(0);
r = G(g).rmt.compute_expected_radius();
G(g).rmt.shift = -r;
G(g).rmt.description = 'Shifted';
G(g).rmt.display_parameters();

%% 3. Dale's law
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_mean_parameters(1/sqrt(n)); % mu_E = 1. mu_I calculated automatically.
G(g).rmt.set_f(0.5);
G(g).rmt.description = 'Dale''s Law';
G(g).rmt.display_parameters();

% %% 4. Zero Row Sum
% g = g + 1;
% G(g).rmt = G(g-1).rmt.copy();
% G(g).rmt.set_row_sum_zero(true);
% G(g).rmt.description = 'Zero Row Sum';
% G(g).rmt.display_parameters();

%% 5. Different standard deviations
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
% b_E = 1.1;
% b_eff = 1;
b_E = 0.25/sqrt(n);
b_eff = 0.5/sqrt(n);
% Calculate b_I using the method from RMT2
b_I = G(g).rmt.compute_sigma_I_from_eff(b_eff, b_E);
G(g).rmt.set_stdev_parameters(b_E, b_I);
G(g).rmt.description = 'Different std devs';
G(g).rmt.display_parameters();

%% 6. 50% Sparse
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_row_sum_zero(false);
G(g).rmt.set_density(0.5);
G(g).rmt.description = '50% sparse';
G(g).rmt.display_parameters();

% %% 7. 50% Sparse with Post-Sparsification ZRS
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_post_sparsification_zrs(true);
G(g).rmt.description = 'Sparse w/ Post-ZRS';
G(g).rmt.display_parameters();

%% Compute and Plot
f1 = figure(1);
set(f1, 'Position', [100   200   990   600], 'Color', 'white');
t = tiledlayout(2, ceil(length(G)/2), 'TileSpacing', 'compact', 'Padding', 'compact');

ax = gobjects(length(G), 1);

for i = 1:length(G)
    % Compute eigenvalues
    G(i).rmt.compute_eigenvalues();

    % Calculate expected radius
    r = G(i).rmt.compute_expected_radius();

    % Set plot circle (center is determined by shift)
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

    % Add description as title
    % title(G(i).rmt.description, 'Interpreter', 'none', 'FontWeight', 'bold');
end

% Determine global scale centered on each plot's shift
max_radius_x = 0;
max_radius_y = 0;

for i = 1:length(G)
    % Temporarily relax aspect ratio to measure actual data extent
    axis(ax(i), 'normal');
    axis(ax(i), 'tight');

    % Get current limits
    current_xlim = xlim(ax(i));
    current_ylim = ylim(ax(i));

    % Center for this plot
    center_x = G(i).rmt.shift;

    % Calculate max distance from center
    dist_x = max(abs(current_xlim - center_x));
    dist_y = max(abs(current_ylim)); % y is centered at 0

    max_radius_x = max(max_radius_x, dist_x);
    max_radius_y = max(max_radius_y, dist_y);
end

% Apply common scale with margin, but centered individually
margin = 1.1;
common_radius_x = max_radius_x * margin;
common_radius_y = max_radius_y * margin;

for i = 1:length(G)
    center_x = G(i).rmt.shift;
    xlim(ax(i), [center_x - common_radius_x, center_x + common_radius_x]);
    ylim(ax(i), [-common_radius_y, common_radius_y]);

    % Ensure x and y units are visually proportional (circles look circular)
    daspect(ax(i), [1 1 1]);
end

% Format axes
for i = 1:length(G)
    axes(ax(i));
    x_lim = xlim;
    y_lim = ylim;
    axis off;

    hold on;
    h_x = plot(x_lim, [0,0], 'k', 'LineWidth', 1.5);
    h_y = plot([0,0], y_lim, 'k', 'LineWidth', 1.5);
    uistack([h_x, h_y], 'bottom');

    % text(x_lim(2), 0, ' Re($\lambda$)', 'Interpreter', 'latex', 'VerticalAlignment', 'middle');
    % text(0, y_lim(2), 'Im($\lambda$)', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    text(x_lim(2), 0, ' Re', 'Interpreter', 'latex', 'VerticalAlignment', 'middle','FontSize',16);
    text(0, y_lim(2), 'Im', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center','FontSize',16);

    text(x_lim(1), y_lim(1), G(i).rmt.description, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontWeight', 'normal', 'Rotation', 90);

    xlim(x_lim);
    ylim(y_lim);

    hold off;
end

% Add letters (a), (b), etc. to subplots
letters = arrayfun(@(c) sprintf('(%s)', c), 'a':'z', 'UniformOutput', false);
% HShift/VShift: positive moves Right/Down (into plot), negative moves Left/Up (outside)
AddLetters2Plots(num2cell(ax), letters, 'FontSize', 18, 'FontWeight', 'normal', 'HShift', -0.01, 'VShift', -0.03);

%% Figure 2: Weight Histograms
f2 = figure(2);
set(f2, 'Position', [150   150   990   600], 'Color', 'white');
t2 = tiledlayout(2, ceil(length(G)/2), 'TileSpacing', 'compact', 'Padding', 'compact');

ax2 = gobjects(length(G), 1);

for i = 1:length(G)
    ax2(i) = nexttile;
    G(i).rmt.plot_weight_histogram(ax2(i));
    title(ax2(i), G(i).rmt.description, 'Interpreter', 'none', 'FontWeight', 'normal');
end

%% Figure 3: Discrete-Time Eigenvalue Distributions
f3 = figure(3);
set(f3, 'Position', [200   100   990   600], 'Color', 'white');
t3 = tiledlayout(2, ceil(length(G)/2), 'TileSpacing', 'compact', 'Padding', 'compact');

ax3 = gobjects(length(G), 1);
dt = 0.01;

for i = 1:length(G)
    ax3(i) = nexttile;
    G(i).rmt.plot_discrete_eigenvalue_distribution(ax3(i), dt);
    title(ax3(i), G(i).rmt.description, 'Interpreter', 'none', 'FontWeight', 'normal');
end

% Determine global scale centered on (1, 0) for discrete systems
max_radius = 0;

for i = 1:length(G)
    % Temporarily relax aspect ratio to measure actual data extent
    axis(ax3(i), 'normal');
    axis(ax3(i), 'tight');

    % Get current limits
    current_xlim = xlim(ax3(i));
    current_ylim = ylim(ax3(i));

    % Center at (1, 0) since exp(0) = 1 is the neutral point
    dist_x = max(abs(current_xlim - 1));
    dist_y = max(abs(current_ylim));

    max_radius = max(max_radius, max(dist_x, dist_y));
end

% Apply common scale with margin, centered at (1, 0)
margin = 1.1;
common_radius = max_radius * margin;

for i = 1:length(G)
    xlim(ax3(i), [1 - common_radius, 1 + common_radius]);
    ylim(ax3(i), [-common_radius, common_radius]);
    daspect(ax3(i), [1 1 1]);
end

% Add title to the tiled layout
title(t3, sprintf('Discrete System Eigenvalues (dt = %.4f)', dt), 'FontWeight', 'bold');

if save_figs
    save_folder = 'test_RMT_figs';
    save_name = 'test_RMT2_example';
    fig_vec = [];
    fig_type = {'fig', 'png', 'svg'};
    save_some_figs_to_folder_2(save_folder, save_name, fig_vec, fig_type)
end
