close all
clear all
clc

set(groot, 'DefaultAxesFontSize', 20);
set(groot, 'DefaultTextFontSize', 20);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultAxesLineWidth', 1.5);

%% default parameters
n = 500;
scale_n = 1.5; % scale factor for larger matrices
assert(mod(scale_n * n, 1) == 0, 'scale_n * n must be an integer');
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
n_large = scale_n * n;
G(g).rmt = RMT(n_large, b, mu);
G(g).rmt.description = ['n=' num2str(n_large) ', dense, \mu = ' num2str(mu)];

%% 3 Sparse
g = g+1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.apply_sparsity(mean_indegree);
G(g).rmt.description = ['n=' num2str(n_large) ', k_{in} = ' num2str(mean_indegree) ', \mu = ' num2str(mu)];

%% 4 EI imbalanced
g = g+1;
G(g).rmt = G(g-1).rmt.copy();
mu2 = 0.025;
G(g).rmt.add_constant(mu2);
G(g).rmt.description = ['n=' num2str(n_large) ', k_{in} = ' num2str(mean_indegree) ', \mu = ' num2str(mu2)];

%% 5 Rajan
g = g+1;
G(g).rmt = G(1).rmt.copy();
G(g).rmt.set_rajan_means(f, mu_E);
G(g).rmt.description = {['n=' num2str(n) ', dense'], '\mu_E = 1, \mu_I = -1'};

%% 6 Rajan zero sum
g = g+1;
G(g).rmt = G(g-1).rmt.copy();
G(g).rmt.row_sum_to_zero();
G(g).rmt.description = {['n=' num2str(n) ', dense'], '\mu_E = 1, \mu_I = -1', 'rows zeroed'};

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
set(f1, 'Position', [-1715 -114 640 1060], 'Color', 'white')
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
row_spacing = 250;  % NaN pixels between rows
col_spacing = 250;  % NaN pixels between columns

% Find maximum matrix size
max_size = 0;
for i_g = 1:length(G)
    max_size = max(max_size, size(G(i_g).rmt.A, 1));
end

% Create cell array to hold padded matrices and track padding amounts
padded_matrices = cell(length(G), 1);
pad_before_vert_array = zeros(length(G), 1);
pad_before_horiz_array = zeros(length(G), 1);
for i_g = 1:length(G)
    A_plot = G(i_g).rmt.A;
    A_plot(~G(i_g).rmt.dense_mask) = NaN;
    
    % Pad matrix to max_size
    % Horizontal: centered (left-right)
    % Vertical: aligned to top (no padding above, all padding below)
    current_size = size(A_plot, 1);
    pad_amount = max_size - current_size;
    
    % Horizontal padding (centered)
    pad_before_horiz = floor(pad_amount / 2);
    pad_after_horiz = pad_amount - pad_before_horiz;
    
    % Vertical padding (aligned to top, all padding below)
    pad_before_vert = 0;
    pad_after_vert = pad_amount;
    
    % Store padding amounts for label positioning
    pad_before_vert_array(i_g) = pad_before_vert;  % Will be 0 for top alignment
    pad_before_horiz_array(i_g) = pad_before_horiz;  % Track horizontal padding
    
    padded_matrices{i_g} = NaN(max_size, max_size);
    padded_matrices{i_g}(pad_before_vert+1:pad_before_vert+current_size, ...
                         pad_before_horiz+1:pad_before_horiz+current_size) = A_plot;
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
set(f2, 'Position', [-1715 -114 900 1060], 'Color', 'white')
ax2 = axes('Parent', f2);

h = imagesc(ax2, concat_matrix);
set(h, 'AlphaData', ~isnan(concat_matrix));
set(ax2, 'Color', [1 1 1]);

colormap(ax2, custom_cmap);
caxis(ax2, clims);
axis(ax2, 'equal', 'tight');
box off

% Add ylabel-style labels at appropriate positions
% Calculate positions for each matrix
ylabel_offset = -15;  % Position labels this many pixels from the left edge of actual data
for i_g = 1:length(G)
    i_row = ceil(i_g / n_cols);
    i_col = mod(i_g - 1, n_cols) + 1;
    
    % Calculate position
    % x position: at the left edge of the non-NaN submatrix
    x_pos = (i_col - 1) * (max_size + col_spacing) + pad_before_horiz_array(i_g) + ylabel_offset;
    % y position: at the top of the matrix (no vertical padding above)
    y_pos = (i_row - 1) * (max_size + row_spacing) + 1;
    
    % Add text annotation as ylabel
    desc = G(i_g).rmt.description;
    text(ax2, x_pos, y_pos, desc, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom', ...
        'FontWeight', 'normal', ...
        'Rotation', 90);
end

set(ax2, 'Visible', 'off');

%% Figure 3: colorbar
% Create separate colorbar figure
% Create modified colormap with thin white line at center for colorbar only
custom_cmap_white = custom_cmap;
% Replace center 2-3 pixels with white
center_idx = round(size(custom_cmap, 1) / 2);
white_line_width = 2; % number of pixels for white line
white_indices = center_idx + (-floor(white_line_width/2):floor(white_line_width/2));
custom_cmap_white(white_indices, :) = repmat([1 1 1], length(white_indices), 1);

show_colorbar_ticks = false; % set to false to hide all ticks

f3 = figure(3);
set(f3, 'Position', [-794   608   291   213], 'Color', 'white')
ax3 = axes('Parent', f3);
colormap(ax3, custom_cmap_white);
cb = colorbar(ax3, 'Location', 'west');
caxis(ax3, clims);
set(ax3, 'Visible', 'off');
if show_colorbar_ticks
    set(cb, 'Box', 'off', 'TickLength', 0, 'Ticks', [-1, 0, 1]);
else
    set(cb, 'Box', 'off', 'TickLength', 0, 'Ticks', []);
end
ylabel(cb, 'Connection Strength');

%% Figure 4: Histograms of A matrices (ignoring NaNs)
f4 = figure(4);
set(f4, 'Position', [-1715 -114 640 1060], 'Color', 'white')
tiledlayout(ceil(length(G)/2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax4 = gobjects(length(G), 1);
for i_g = 1:length(G)
    ax4(i_g) = nexttile;
    A_plot = G(i_g).rmt.A;
    A_plot(~G(i_g).rmt.dense_mask) = NaN;  % Set non-dense elements to NaN
    
    % Check if this RMT has E/I structure defined
    if ~isempty(G(i_g).rmt.E) && ~isempty(G(i_g).rmt.I)
        % Separate E and I column values
        A_E_cols = A_plot(:, G(i_g).rmt.E);  % Excitatory columns
        A_I_cols = A_plot(:, G(i_g).rmt.I);  % Inhibitory columns
        
        A_E_vals = A_E_cols(:);
        A_E_vals = A_E_vals(~isnan(A_E_vals));  % Remove NaNs
        
        A_I_vals = A_I_cols(:);
        A_I_vals = A_I_vals(~isnan(A_I_vals));  % Remove NaNs
        
        % Plot overlapping histograms with transparency
        hold on;
        histogram(A_E_vals, 50, 'FaceColor', [0 0 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);  % Blue for E
        histogram(A_I_vals, 50, 'FaceColor', [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);  % Red for I
        hold off;
    else
        % Original gray histogram when no E/I structure
        A_vals = A_plot(:);
        A_vals = A_vals(~isnan(A_vals));  % Remove NaNs
        histogram(A_vals, 50, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
    end
    
    desc = G(i_g).rmt.description;
    % title(desc);
    grid off;
    axis off;
end

% Link x and y axes across all subplots
linkaxes(ax4, 'xy');

%% Figure 5: Row sums as vertical plots
f5 = figure(5);
set(f5, 'Position', [-1715 -114 640 1060], 'Color', 'white')
tiledlayout(ceil(length(G)/2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax5 = gobjects(length(G), 1);
for i_g = 1:length(G)
    ax5(i_g) = nexttile;
    A_curr = G(i_g).rmt.A;
    A_curr(~G(i_g).rmt.dense_mask) = NaN;  % Set non-dense elements to NaN
    n_curr = size(A_curr, 1);
    row_sums = sum(A_curr, 2, 'omitnan');  % Sum across columns, ignoring NaNs
    plot(row_sums, 1:n_curr, 'k');  % 'k' for black color
    xlabel('Row Sum');
    ylabel('Row Index');
    desc = G(i_g).rmt.description;
    % title(desc);
    grid off;
    set(gca, 'YDir', 'reverse');  % Make row 1 at top
    axis off;
end

% Link x and y axes across all subplots
linkaxes(ax5, 'xy');
% Set the xlims of ax5(1) to 3x its current xlim
xlim(ax5(1), 3*xlim(ax5(1)));
