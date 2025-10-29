close all
clear all
clc

set(groot, 'DefaultAxesFontSize', 20);
set(groot, 'DefaultTextFontSize', 18);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultAxesLineWidth', 1.5);

%% default parameters
n = 500;
b = 1;
mu = 0;
f = 0.5;
mu_E = 1; % relative to b
mean_indegree = n; % desired mean number of connections in and out if sparse

%% Struct to hold RMT objects
G = struct();

%% 1 Base Ginibre ensemble
g = 1;
G(g).rmt = RMT(n, b, mu);
G(g).rmt.description = ['n=' num2str(n) ', dense, \mu = ' num2str(mu)];

%% 2 Bigger Ginibre ensemble
g = g+1;
G(g).rmt = RMT(2*n, b, mu);
G(g).rmt.description = ['n=' num2str(2*n) ', dense, \mu = ' num2str(mu)];

%% 3 Sparse
g = g+1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.apply_sparsity(mean_indegree);
G(g).rmt.description = ['n=' num2str(2*n) ', k_{in} = ' num2str(mean_indegree) ', \mu = ' num2str(mu)];

%% 4 EI imbalanced
g = g+1;
G(g).rmt = G(g-1).rmt.copy();
mu2 = 0.025;
G(g).rmt.add_constant(mu2);
G(g).rmt.description = ['n=' num2str(2*n) ', k_{in} = ' num2str(mean_indegree) ', \mu = ' num2str(mu2)];

%% 5 Rajan
g = g+1;
G(g).rmt = G(1).rmt.copy();
G(g).rmt.set_rajan_means(f, mu_E);
G(g).rmt.description = {['n=' num2str(n) ', dense'], ['\mu_E = 1, \mu_I = -1']};

%% 6 Rajan zero sum
g = g+1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.row_sum_to_zero();
G(g).rmt.description = {['n=' num2str(n) ', dense'], ['\mu_E = 1, \mu_I = -1, rows zeroed']};

%% Ginibre with shift to edge
% g = g+1;
% G(g).rmt = G(1).rmt.copy();
% shift = -G(g).rmt.b*sqrt(G(g).rmt.n)*sqrt(G(g).rmt.density);
% G(g).rmt.shift_diagonal(shift);
% G(g).rmt.description = 'Ginibre with shift';



%% Compute eigenvalues and set circles for all
for i_g=1:length(G)
    G(i_g).rmt.compute_eigenvalues();
    r = G(i_g).rmt.b*sqrt(G(i_g).rmt.n)*sqrt(G(i_g).rmt.density);
    x_center = 0;
    G(i_g).rmt.set_plot_circle(r, x_center);
end

%% Make figures of different scenarios
f1 = figure(1);
set(f1, 'Position', [-1715 -114 640 1060])
tiledlayout(ceil(length(G)/2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = gobjects(length(G), 1);
for i_g = 1:length(G)
    ax(i_g) = nexttile;
    G(i_g).rmt.plot_circle(ax(i_g));
    hold on
    G(i_g).rmt.plot_eigenvalue_distribution(ax(i_g));
    hold off
end

linkaxes(ax,'y')

%% Make figure of A matrices
% Collect all A matrix values to determine clims
all_A_values = [];
for i_g = 1:length(G)
    all_A_values = [all_A_values; G(i_g).rmt.A(:)];
end

% Calculate symmetric clims based on 5% percentile
max_abs_vals = abs(all_A_values);
clim_val = prctile(max_abs_vals, 95);
clims = [-clim_val, clim_val];

% Create custom colormap: blue (negative) -> black (zero) -> red (positive)
n_colors = 256;
half = n_colors / 2;
% Blue to black for negative values
blues = [linspace(0, 0, half)', linspace(0, 0, half)', linspace(1, 0, half)'];
% Black to red for positive values
reds = [linspace(0, 1, half)', linspace(0, 0, half)', linspace(0, 0, half)'];
custom_cmap = [blues; reds];

f2 = figure(2);
set(f2, 'Position', [-1715 -114 640 1060])
tiledlayout(ceil(length(G)/2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax2 = gobjects(length(G), 1);
for i_g = 1:length(G)
    ax2(i_g) = nexttile;
    
    % Create A matrix with NaN for sparse (zero) elements
    A_plot = G(i_g).rmt.A;
    A_plot(~G(i_g).rmt.dense_mask) = NaN;
    
    % Plot with transparency for NaN values
    h = imagesc(ax2(i_g), A_plot);
    box off
    set(h, 'AlphaData', ~isnan(A_plot));
    
    % Set background color to bright green for sparse elements
    set(ax2(i_g), 'Color', [1 1 1]);
    
    colormap(ax2(i_g), custom_cmap);
    caxis(ax2(i_g), clims);
    axis(ax2(i_g), 'equal', 'tight');
    title(ax2(i_g), G(i_g).rmt.description, 'FontWeight', 'normal');
    colorbar(ax2(i_g));
    
    % Hide the axes completely
    set(ax2(i_g), 'Visible', 'off');
    
end
