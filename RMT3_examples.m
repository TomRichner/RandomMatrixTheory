% RMT3_examples.m
% Illustrative examples for RMT3 class based on Harris et al. (2023)

clear; close all; clc;

N = 500;
G = struct();
g = 0;

%% Dense unbalanced Network
% Expectation: All eigenvalues within circle. No outliers.
g = g + 1;
G(g).rmt = RMT3(N);
% set_params(mu_e, mu_i, sigma_e, sigma_i, f, alpha)
G(g).rmt.set_params(0, -0.3/sqrt(N), 0.5/sqrt(N), 0.5/sqrt(N), 0.5, 1.0);
G(g).rmt.set_zrs_mode('none');
G(g).rmt.description = 'Dense unbalanced';
G(g).rmt.display_parameters();

%% Dense Unbalanced Dales
% Expectation: Outlier at -2 (since N*(-2/N + 0)) if we set means that way?
% Let's use mu_e=1/N, mu_i=-3/N. (Not strictly balanced).
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
% Scale means by 1/N for outlier O(1)
G(g).rmt.set_params(1/sqrt(N), -1.3/sqrt(N), 0.5/sqrt(N), 0.5/sqrt(N), 0.5);
G(g).rmt.set_zrs_mode('none');
G(g).rmt.description = 'Dense Unbalanced Dales';
G(g).rmt.display_parameters();

%% Dense Unbalanced Dales ZRS
% Expectation: Outlier at -2 (since N*(-2/N + 0)) if we set means that way?
% Let's use mu_e=1/N, mu_i=-3/N. (Not strictly balanced).
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
% Scale means by 1/N for outlier O(1)
G(g).rmt.set_zrs_mode('ZRS');
G(g).rmt.description = 'Dense Unbalanced Dales ZRS';
G(g).rmt.display_parameters();


%% Sparse Unbalanced
% Expectation: Sparse matrix (alpha=0.1), balanced means. SZRS removes local outliers.
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_alpha(0.5);
G(g).rmt.set_zrs_mode('none');
G(g).rmt.description = 'Sparse Unbalanced';
G(g).rmt.display_parameters();


%% Sparse Unbalanced SZRS
% Expectation: Imbalance (Outlier) preserved. Local outliers removed (Radius clean).
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
% Use Partial SZRS to keep the outlier from mu_e/mu_i but clean up bulk
G(g).rmt.set_zrs_mode('SZRS');
G(g).rmt.description = 'Sparse Unbalanced SZRS';
G(g).rmt.display_parameters();

%% Sparse Unbalanced Partial SZRS
% Expectation: Imbalance (Outlier) preserved. Local outliers removed (Radius clean).
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
% Use Partial SZRS to keep the outlier from mu_e/mu_i but clean up bulk
G(g).rmt.set_zrs_mode('Partial_SZRS');
G(g).rmt.description = 'Sparse Unbalanced Partial SZRS';
G(g).rmt.display_parameters();


%% Compute and Plot
f1 = figure(1);
set(f1, 'Position', [100   200   990   600], 'Color', 'white');
t = tiledlayout(2, ceil(length(G)/2), 'TileSpacing', 'compact', 'Padding', 'compact');

ax = gobjects(length(G), 1);

for i = 1:length(G)
    % Compute eigenvalues
    G(i).rmt.compute_eigenvalues();

    % Plot
    ax(i) = nexttile;
    G(i).rmt.plot_spectrum(ax(i));
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

fprintf('All examples verified.\n');
