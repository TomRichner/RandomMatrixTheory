close all
clear all
clc

set(groot, 'DefaultAxesFontSize', 20);
set(groot, 'DefaultTextFontSize', 20);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultAxesLineWidth', 1.5);

%% default parameters
make_imagesc_nans_transparent = false; % if false, NaNs are treated as zeros
n = 1000;
scale_n = 1; % scale factor for larger matrices
assert(mod(scale_n * n, 1) == 0, 'scale_n * n must be an integer');
b = 1;
mu = 0;
f = 0.66667;
mu_E = 1; % relative to b
mean_indegree = 500; % desired mean number of connections in and out if sparse


G = struct();
%% 1 Base Ginibre ensemble
g = 1;
G(g).rmt = RMT(n, b, mu);
G(g).rmt.description = ['G(0,1)'];

%% Rajan
g = 2;
G(g).rmt = G(1).rmt.copy();
G(g).rmt.set_rajan_means(f, mu_E);
G(g).rmt.description = {['Rajan, f = 2/3, mu_E = 1, mu_I = -2']};

%% make network sparse
g = 3;
G(g).rmt = G(2).rmt.copy();
G(g).rmt.apply_sparsity(mean_indegree);
G(g).rmt.description = ['Sparse connections'];

activation_one = rand(n,1)<0.5;
activation_two= rand(n,1)<0.5;

%% apply sparse activation one
g = 4;
G(g).rmt = G(3).rmt.copy();
G(g).rmt.zero_columns(activation_one);
G(g).rmt.description = ['Sparse Activation One'];

%% homeostasis row zero
g = 5;
G(g).rmt = G(4).rmt.copy();
G(g).rmt.row_sum_to_zero();
row_sum_Z_diff = G(g).rmt.row_sum_Z_diff; % these were the modified weights 
abscissa5 = max(real(eig(G(g).rmt.A))); % abscissa of this matrix
G(g).rmt.description = ['Rows summed to zero'];

%% shifted to EOS, This will be subplot(1,2,1)
g = 6;
G(g).rmt = G(5).rmt.copy();
G(g).rmt.shift_diagonal(-abscissa5)
G(g).rmt.compute_eigenvalues();
G(g).rmt.set_plot_circle(abscissa5, -abscissa5);
G(g).rmt.description = ['Activation One, row-zeroed, and shifted to EOS'];

%% apply sparse activation two, This will be subplot(1,2,2)
g = 7;
G(g).rmt = G(3).rmt.copy(); % reset to 3, before activation one was applied
G(g).rmt.A = G(g).rmt.A + row_sum_Z_diff; % apply the previous homeostatic row-zero matrix
G(g).rmt.zero_columns(activation_two); % zero according to the new activation
G(g).rmt.shift_diagonal(-abscissa5) % shift accordin previous homeostatic shift
G(g).rmt.compute_eigenvalues();
G(g).rmt.set_plot_circle(abscissa5, -abscissa5);
G(g).rmt.description = ['Activation Two'];

%% figure 1, eigenvalue distribution
f1 = figure(1);
% set(f1, 'Position', [100 100 710 1060], 'Color', 'white')
set(f1, 'Color', 'white')

% tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
tiledlayout(1, 2, 'TileSpacing', 'compact');

g_to_plot = [6, 7]
ax = gobjects(2, 1);
for i_g = 1:length(g_to_plot)
    g = g_to_plot(i_g)
    ax(i_g) = nexttile;
    G(g).rmt.plot_circle(ax(i_g));
    hold on
    G(g).rmt.plot_eigenvalue_distribution(ax(i_g));
    hold off
end

linkaxes(ax,'xy')

% Customize axes after linking
for i_g = 1:length(ax)
    axes(ax(i_g));
    x_lim = xlim;
    y_lim = ylim;
    axis off;
    
    hold on;
    h_x = plot(x_lim, [0,0], 'k', 'LineWidth', 1.5);
    h_y = plot([0,0], y_lim, 'k', 'LineWidth', 1.5);
    uistack([h_x, h_y], 'bottom');

    text(x_lim(2), 0, ' Re($\lambda$)', 'Interpreter', 'latex', 'VerticalAlignment', 'middle');
    text(0, y_lim(2), 'Im($\lambda$)', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

    % text(x_lim(1), y_lim(1), G(i_g).rmt.description, ...
    %     'HorizontalAlignment', 'left', ...
    %     'VerticalAlignment', 'top', ...
    %     'FontWeight', 'normal', ...
    %     'Rotation', 90);

    xlim(x_lim);
    ylim(y_lim);

    hold off
end
axis equal

%% figure 2, modified imagesc plot of A
% Collect all A matrix values to determine clims
all_A_values = [];
for i_g = 1:length(g_to_plot)
    g = g_to_plot(i_g);
    all_A_values = [all_A_values; G(g).rmt.A(:)];
end

% Calculate symmetric clims based on 95% percentile
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

%% Create concatenated matrix with constant pixel size (1x2 layout)
% Parameters for layout
n_cols = 2;
n_rows = 1;
col_spacing = 250;  % NaN pixels between columns

% Find maximum matrix size
max_size = 0;
for i_g = 1:length(g_to_plot)
    g = g_to_plot(i_g);
    max_size = max(max_size, size(G(g).rmt.A, 1));
end

% Create cell array to hold padded matrices and track padding amounts
padded_matrices = cell(length(g_to_plot), 1);
pad_before_horiz_array = zeros(length(g_to_plot), 1);
for i_g = 1:length(g_to_plot)
    g = g_to_plot(i_g);
    A_plot = G(g).rmt.A;
    if make_imagesc_nans_transparent
        A_plot(~G(g).rmt.dense_mask) = NaN;
    else
        A_plot(~G(g).rmt.dense_mask) = 0;
    end
    
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
    pad_before_horiz_array(i_g) = pad_before_horiz;
    
    padded_matrices{i_g} = NaN(max_size, max_size);
    padded_matrices{i_g}(pad_before_vert+1:pad_before_vert+current_size, ...
                         pad_before_horiz+1:pad_before_horiz+current_size) = A_plot;
end

% Build the concatenated matrix (1 row with 2 columns)
concat_matrix = [];
for i_col = 1:n_cols
    concat_matrix = [concat_matrix, padded_matrices{i_col}];
    
    % Add column spacing (except after last column)
    if i_col < n_cols
        concat_matrix = [concat_matrix, NaN(max_size, col_spacing)];
    end
end

% Plot the concatenated matrix
f2 = figure(2);
set(f2, 'Position', [100 100 1200 600], 'Color', 'white')
ax2 = axes('Parent', f2);

h = imagesc(ax2, concat_matrix);
set(h, 'AlphaData', ~isnan(concat_matrix));
set(ax2, 'Color', [1 1 1]);

colormap(ax2, custom_cmap);
caxis(ax2, clims);
axis(ax2, 'equal', 'tight');
box off

% Add ylabel-style labels at appropriate positions
ylabel_offset = -15;  % Position labels this many pixels from the left edge of actual data
for i_g = 1:length(g_to_plot)
    g = g_to_plot(i_g);
    i_col = i_g;
    
    % Calculate position
    % x position: at the left edge of the non-NaN submatrix
    x_pos = (i_col - 1) * (max_size + col_spacing) + pad_before_horiz_array(i_g) + ylabel_offset;
    % y position: at the top of the matrix
    y_pos = 1;
    
    % Add text annotation as ylabel
    desc = G(g).rmt.description;
    text(ax2, x_pos, y_pos, desc, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom', ...
        'FontWeight', 'normal', ...
        'Rotation', 90);
end

set(ax2, 'Visible', 'off');
