close all
% if ~exist('seed', 'var')
    seed = 78
% else
%     seed = seed + 1
% end
clearvars -except seed
clc

rng(seed)

set(groot, 'DefaultAxesFontSize', 16);
set(groot, 'DefaultTextFontSize', 16);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultAxesLineWidth', 1.5);

save_figs = true;

n = 600;

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
b_E = 1.375;
b_eff = 1.0;
% Calculate b_I using the method from RMT2
b_I = G(g).rmt.compute_sigma_I_from_eff(b_eff, b_E);
G(g).rmt.set_stdev_parameters(b_E, b_I);
G(g).rmt.description = 'Different std devs';

%% 6. 50% Sparse
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_sparsity(0.5);
G(g).rmt.description = '50% sparse';

% %% 7. 50% Sparse with Post-Sparsification ZRS
% g = g + 1;
% G(g).rmt = G(g-1).rmt.copy();
% G(g).rmt.set_post_sparsification_zrs(true);
% G(g).rmt.description = 'Sparse w/ Post-ZRS';

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
    G(i).rmt.set_outlier_factor(1.035);

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
    
    % text(x_lim(1), y_lim(1), G(i).rmt.description, ...
    %     'HorizontalAlignment', 'left', ...
    %     'VerticalAlignment', 'top', ...
    %     'FontWeight', 'normal', ...
    %     'Rotation', 90);
        
    xlim(x_lim);
    ylim(y_lim);
    
    hold off;
end

% Add letters (a), (b), ... to subplots
addpath('old');
letters = arrayfun(@(c) sprintf('(%s)', c), 'a':'z', 'UniformOutput', false);
% HShift/VShift: positive moves Right/Down (into plot), negative moves Left/Up (outside)
AddLetters2Plots(num2cell(ax), letters, 'FontSize', 18, 'FontWeight', 'normal', 'HShift', -0.01, 'VShift', -0.03);

if save_figs
    save_folder = 'new_RMT_figs';
    save_name = 'RMT2_example';
    fig_vec = [];
    fig_type = {'fig', 'png', 'svg'};
    save_some_figs_to_folder_2(save_folder, save_name, fig_vec, fig_type)
end
