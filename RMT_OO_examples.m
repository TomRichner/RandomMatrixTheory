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

%% Create concatenated matrix with constant pixel size
% Parameters for layout
n_cols = 2;
n_rows = ceil(length(G) / n_cols);
row_spacing = 50;  % NaN pixels between rows
col_spacing = 50;  % NaN pixels between columns

% Find maximum matrix size
max_size = 0;
for i_g = 1:length(G)
    max_size = max(max_size, size(G(i_g).rmt.A, 1));
end

% Create cell array to hold padded matrices
padded_matrices = cell(length(G), 1);
for i_g = 1:length(G)
    A_plot = G(i_g).rmt.A;
    A_plot(~G(i_g).rmt.dense_mask) = NaN;
    
    % Pad matrix to max_size, centered
    current_size = size(A_plot, 1);
    pad_amount = max_size - current_size;
    pad_before = floor(pad_amount / 2);
    pad_after = pad_amount - pad_before;
    
    padded_matrices{i_g} = NaN(max_size, max_size);
    padded_matrices{i_g}(pad_before+1:pad_before+current_size, ...
                         pad_before+1:pad_before+current_size) = A_plot;
end

% Build the concatenated matrix row by row
concat_matrix = [];
for i_row = 1:n_rows
    row_matrices = [];
    for i_col = 1:n_cols
        idx = (i_row - 1) * n_cols + i_col;
        if idx <= length(G)
            row_matrices = [row_matrices, padded_matrices{idx}];
        else
            % Fill empty space with NaN
            row_matrices = [row_matrices, NaN(max_size, max_size)];
        end
        
        % Add column spacing (except after last column)
        if i_col < n_cols
            row_matrices = [row_matrices, NaN(max_size, col_spacing)];
        end
    end
    
    concat_matrix = [concat_matrix; row_matrices];
    
    % Add row spacing (except after last row)
    if i_row < n_rows
        concat_matrix = [concat_matrix; NaN(row_spacing, size(row_matrices, 2))];
    end
end

% Plot the concatenated matrix
f2 = figure(2);
set(f2, 'Position', [-1715 -114 900 1060])
ax2 = axes('Parent', f2);

h = imagesc(ax2, concat_matrix);
set(h, 'AlphaData', ~isnan(concat_matrix));
set(ax2, 'Color', [1 1 1]);

colormap(ax2, custom_cmap);
caxis(ax2, clims);
axis(ax2, 'equal', 'tight');
box off

% Add single colorbar
cb = colorbar(ax2);
ylabel(cb, 'Connection Strength');

% Add titles at appropriate positions
% Calculate center positions for each matrix
title_offset = max_size * 0.05;  % Position titles above matrices
for i_g = 1:length(G)
    i_row = ceil(i_g / n_cols);
    i_col = mod(i_g - 1, n_cols) + 1;
    
    % Calculate position
    x_pos = (i_col - 1) * (max_size + col_spacing) + max_size / 2;
    y_pos = (i_row - 1) * (max_size + row_spacing) - title_offset;
    
    % Add text annotation
    desc = G(i_g).rmt.description;
    if iscell(desc)
        desc = strjoin(desc, ', ');
    end
    text(ax2, x_pos, y_pos, desc, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 14, ...
        'FontWeight', 'normal');
end

set(ax2, 'Visible', 'off');
