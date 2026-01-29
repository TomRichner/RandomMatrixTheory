% Illustrative examples for RMT4 class based on Harris et al. 2023

clear; close all; clc;

save_figs = true;

rng(100)

N = 600;
G = struct();
g = 0;
E_W = +0.05/sqrt(N);

%% Dense unbalanced Network
% Expectation: Outlier visible, bulk within radius
g = g + 1;
G(g).rmt = RMT4(N);
% set_params(mu_tilde_e, mu_tilde_i, sigma_tilde_e, sigma_tilde_i, f, alpha)
G(g).rmt.set_params(0+E_W, 0+E_W, 1/sqrt(N), 1/sqrt(N), 1, 1.0);
G(g).rmt.set_zrs_mode('none');
G(g).rmt.description = 'Dense unbalanced (none)';
G(g).rmt.display_parameters();

%% Dense unbalanced shifted Network
% Expectation: Outlier visible, bulk within radius
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_params(0,0); % remove imbalance
R = G(g).rmt.R;
G(g).rmt.shift = -R;  % Shift by spectral radius
G(g).rmt.description = 'Dense unbalanced shifted';
G(g).rmt.display_parameters();

%% Dense Unbalanced Dales
% Expectation: Visible outlier from E-I imbalance
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_params(1/sqrt(N), -1/sqrt(N), 1/sqrt(N), 1/sqrt(N), 0.5);
G(g).rmt.description = 'Dense Dales';
G(g).rmt.display_parameters();

%% Dense Unbalanced Dalse ZRS
% Expectation: Visible outlier from E-I imbalance
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_zrs_mode('ZRS');
G(g).rmt.description = 'Dense Dales ZRS';
G(g).rmt.display_parameters();

%% Dense Unbalanced Dales different sigmas
% Expectation: Non-uniform eigenvalue density when sigma_e != sigma_i
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
% Set a different sigma_tilde_e, then compute sigma_tilde_i to maintain target variance
G(g).rmt.sigma_tilde_e = 0.35/sqrt(N);
target_variance = 1/N;  % Same total variance as before
G(g).rmt.sigma_tilde_i = G(g).rmt.compute_sigma_tilde_i_for_target_variance(target_variance);
G(g).rmt.description = 'Dense Dales different sigmas';
G(g).rmt.display_parameters();

%% Sparse Unbalanced with Partial SZRS (Eq 32)
% Expectation: Outlier preserved, local random outliers removed
g = g + 1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.set_alpha(0.5);
G(g).rmt.set_zrs_mode('Partial_SZRS');
G(g).rmt.description = 'Sparse Unbalanced Dales Partial SZRS';
G(g).rmt.display_parameters();

%% Compute and Plot
f1 = figure(1);
set(f1, 'Position', [100   300   760   500], 'Color', 'white');
t = tiledlayout(2, ceil(length(G)/2), 'TileSpacing', 'tight', 'Padding', 'compact');

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
max_R = 0;

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
    max_R = max(max_R, G(i).rmt.R);
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
    y_lim_axis = min(0.75*y_lim,1.1*max_R);%y_lim)
    axis off;

    hold on;
    h_x = plot(x_lim, [0,0], 'k', 'LineWidth', 1.5);
    h_y = plot([0,0], y_lim_axis, 'k', 'LineWidth', 1.5);
    uistack([h_x, h_y], 'bottom');

    text(x_lim(2), 0, ' Re', 'Interpreter', 'latex', 'VerticalAlignment', 'middle', 'FontSize', 16);
    text(0, y_lim_axis(2), 'Im', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 16);
    % text(x_lim(1), y_lim(1), G(i).rmt.description, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontWeight', 'normal', 'Rotation', 90);

    xlim(x_lim);
    ylim(y_lim);
    hold off;
end

% Add letters (a), (b), etc. to subplots
letters = arrayfun(@(c) sprintf('(%s)', c), 'a':'z', 'UniformOutput', false);
AddLetters2Plots(num2cell(ax), letters, 'FontSize', 14, 'FontWeight', 'normal', 'HShift', +0.005, 'VShift', +0.005);

drawnow;
set(f1, 'Position', [100   300   760   500]);

if save_figs
    save_some_figs_to_folder_2('figs', 'RMT4_examples', [], {'fig', 'svg', 'png', 'jp2'});
end