% RMT4_examples.m
% Illustrative examples for RMT4 class based on Harris et al. (2023)
% Testing equations 15-18, 24-25, 30-31 with ZRS, SZRS, and Partial SZRS

clear; close all; clc;

N = 500;
G = struct();
g = 0;

%% 1. Dense unbalanced Network (no ZRS)
% Expectation: Outlier visible, bulk within radius
g = g + 1;
G(g).rmt = RMT4(N);
% set_params(mu_tilde_e, mu_tilde_i, sigma_tilde_e, sigma_tilde_i, f, alpha)
G(g).rmt.set_params(0, -0/sqrt(N), 0.5/sqrt(N), 0.5/sqrt(N), 0.5, 1.0);
G(g).rmt.set_zrs_mode('none');
G(g).rmt.description = 'Dense unbalanced (none)';
G(g).rmt.display_parameters();

%% 2. Dense Unbalanced with means
% Expectation: Visible outlier from E-I imbalance
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_params(1/sqrt(N), -1/sqrt(N), 0.5/sqrt(N), 0.5/sqrt(N), 0.5);
G(g).rmt.set_zrs_mode('none');
G(g).rmt.description = 'Dense Unbalanced Dales';
G(g).rmt.display_parameters();

%% 3. Dense Unbalanced with ZRS (Eq 24-25)
% Expectation: Projection removes random outliers, preserves mean structure
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_zrs_mode('ZRS');
G(g).rmt.description = 'Dense Unbalanced ZRS';
G(g).rmt.display_parameters();

%% 4. Sparse Unbalanced (no correction)
% Expectation: Local outliers may escape radius
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_alpha(0.5);
G(g).rmt.set_zrs_mode('none');
G(g).rmt.description = 'Sparse Unbalanced (none)';
G(g).rmt.display_parameters();

%% 5. Sparse Unbalanced with SZRS (Eq 30-31)
% Expectation: Full SZRS forces row sums to zero, removes E-I outlier
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_zrs_mode('SZRS');
G(g).rmt.description = 'Sparse Unbalanced SZRS';
G(g).rmt.display_parameters();

%% 6. Sparse Unbalanced with Partial SZRS (Eq 32)
% Expectation: Outlier preserved, local random outliers removed
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_zrs_mode('Partial_SZRS');
G(g).rmt.description = 'Sparse Partial SZRS';
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
    axis(ax(i), 'normal');
    axis(ax(i), 'tight');
    current_xlim = xlim(ax(i));
    current_ylim = ylim(ax(i));
    center_x = G(i).rmt.shift;
    dist_x = max(abs(current_xlim - center_x));
    dist_y = max(abs(current_ylim));
    max_radius_x = max(max_radius_x, dist_x);
    max_radius_y = max(max_radius_y, dist_y);
end

% Apply common scale with margin
margin = 1.1;
common_radius_x = max_radius_x * margin;
common_radius_y = max_radius_y * margin;

for i = 1:length(G)
    center_x = G(i).rmt.shift;
    xlim(ax(i), [center_x - common_radius_x, center_x + common_radius_x]);
    ylim(ax(i), [-common_radius_y, common_radius_y]);
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

    text(x_lim(2), 0, ' Re', 'Interpreter', 'latex', 'VerticalAlignment', 'middle', 'FontSize', 16);
    text(0, y_lim(2), 'Im', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 16);
    text(x_lim(1), y_lim(1), G(i).rmt.description, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontWeight', 'normal', 'Rotation', 90);

    xlim(x_lim);
    ylim(y_lim);
    hold off;
end

% Add letters (a), (b), etc. to subplots
letters = arrayfun(@(c) sprintf('(%s)', c), 'a':'z', 'UniformOutput', false);
AddLetters2Plots(num2cell(ax), letters, 'FontSize', 18, 'FontWeight', 'normal', 'HShift', -0.01, 'VShift', -0.03);

fprintf('All RMT4 examples completed.\n');
