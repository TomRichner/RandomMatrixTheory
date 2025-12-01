close all
clear all
clc

set(groot, 'DefaultAxesFontSize', 14);
set(groot, 'DefaultTextFontSize', 14);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultAxesLineWidth', 1.5);

n = 500;

%% Struct to hold RMT objects
G = struct();

%% 1. Standard Circular Law with outlier
g = 1;
G(g).rmt = RMT2(n);
G(g).rmt.set_mu(30/n);
G(g).rmt.description = 'Standard Circular Law with outlier';

%% 2. Shifted
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_mu(0);
r = G(g).rmt.compute_expected_radius();
G(g).rmt.shift = -r;
G(g).rmt.description = 'Shifted';

%% 3. Dale's law
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_f(0.5);
G(g).rmt.set_mean_parameters(1); % mu_E = 1. mu_I calculated automatically.
G(g).rmt.description = 'Dale''s Law';

%% 4. Zero Row Sum
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_row_sum_zero(true);
G(g).rmt.description = 'Zero Row Sum';

%% 5. Different standard deviations
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
% b_E = 1.25;
% b_eff = 1.0;
% % Calculate b_I using the method from RMT2
% b_I = G(g).rmt.compute_sigma_I_from_eff(b_eff, b_E);
% G(g).rmt.set_stdev_parameters(b_E, b_I);
G(g).rmt.description = 'Different std devs';

%% 6. 50% Sparse
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_sparsity(0.5);
G(g).rmt.description = '50% sparse';

%% 7. 50% Sparse with Post-Sparsification ZRS
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_post_sparsification_zrs(true);
G(g).rmt.description = 'Sparse w/ Post-ZRS';

%% Compute and Plot
f1 = figure(1);
set(f1, 'Position', [100 100 1200 700], 'Color', 'white');
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
    
    % Plot
    ax(i) = nexttile;
    G(i).rmt.plot_circle(ax(i));
    hold on;
    G(i).rmt.plot_eigenvalue_distribution(ax(i));
    hold off;
    
    % Add description as title
    % title(G(i).rmt.description, 'Interpreter', 'none', 'FontWeight', 'bold');
end

% Link axes to ensure same scale
linkaxes(ax, 'xy');

% Adjust limits to be symmetric and contain all data
xl = xlim(ax(1));
yl = ylim(ax(1));
max_val = max([abs(xl), abs(yl)]);
new_lim = [-max_val, max_val] * 1.1; % Add 10% margin
xlim(ax, new_lim);
ylim(ax, new_lim);

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
    
    text(x_lim(2), 0, ' Re($\lambda$)', 'Interpreter', 'latex', 'VerticalAlignment', 'middle');
    text(0, y_lim(2), 'Im($\lambda$)', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    
    text(x_lim(1), y_lim(1), G(i).rmt.description, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontWeight', 'normal', ...
        'Rotation', 90);
        
    xlim(x_lim);
    ylim(y_lim);
    
    hold off;
end

