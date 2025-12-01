close all
clear all
clc

set(groot, 'DefaultAxesFontSize', 20);
set(groot, 'DefaultTextFontSize', 20);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultAxesLineWidth', 1.5);

%% default parameters
make_imagesc_nans_transparent = false; % if false, NaNs are treated as zeros
n = 500;
scale_n = 1.5; % scale factor for larger matrices
assert(mod(scale_n * n, 1) == 0, 'scale_n * n must be an integer');
b = 1;
mu = 0;
f = 0.5;
mu_E = 1; % relative to b
mean_indegree = n; % desired mean number of connections in and out if sparse

% % % Figures 1-5 examples of changing n, changing density, changing mu, changing mu_E/mu_I and applying row zero
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
mu2 = 0.04;
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
set(f1, 'Position', [100 200 710 1060], 'Color', 'white')
tiledlayout(ceil(length(G)/2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = gobjects(length(G), 1);
for i_g = 1:length(G)
    ax(i_g) = nexttile;
    G(i_g).rmt.plot_circle(ax(i_g));
    hold on
    G(i_g).rmt.plot_eigenvalue_distribution(ax(i_g));
    hold off
end

linkaxes(ax,'xy')

% ensure left and right xlims are the same and +/- max(abs(xlim))
xl = xlim(ax(1));
new_lim = max(abs(xl));
xlim(ax, [-new_lim, new_lim]);

% Customize axes after linking
for i_g = 1:length(G)
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

    text(x_lim(1), y_lim(1), G(i_g).rmt.description, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontWeight', 'normal', ...
        'Rotation', 90);

    xlim(x_lim);
    ylim(y_lim);

    hold off
end

%% Figure 2 Make figure of A matrices
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
    if make_imagesc_nans_transparent
        A_plot(~G(i_g).rmt.dense_mask) = NaN;
    else
        A_plot(~G(i_g).rmt.dense_mask) = 0;
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
set(f2, 'Position', [100 200 900 1060], 'Color', 'white')
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
% % Replace center 2-3 pixels with white
if make_imagesc_nans_transparent % otherwise don't add the white line
    center_idx = round(size(custom_cmap, 1) / 2);
    white_line_width = 2; % number of pixels for white line
    white_indices = center_idx + (-floor(white_line_width/2):floor(white_line_width/2));
    custom_cmap_white(white_indices, :) = repmat([1 1 1], length(white_indices), 1);
end

show_colorbar_ticks = false; % set to false to hide all ticks

f3 = figure(3);
set(f3, 'Position', [100 200   291   213], 'Color', 'white')
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
set(f4, 'Position', [100 200 640 1060], 'Color', 'white')
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
        hold on
        histogram(A_vals, 50, 'FaceColor', [0 0 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % overlapping identical histograms to match the colors above
        histogram(A_vals, 50, 'FaceColor', [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        hold off
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
set(f5, 'Renderer', 'painters');
set(f5, 'Position', [100 200 640 1060], 'Color', 'white')
tiledlayout(ceil(length(G)/2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax5 = gobjects(length(G), 1);
for i_g = 1:length(G)
    ax5(i_g) = nexttile;
    A_curr = G(i_g).rmt.A;
    A_curr(~G(i_g).rmt.dense_mask) = NaN;  % Set non-dense elements to NaN
    n_curr = size(A_curr, 1);
    row_sums = sum(A_curr, 2, 'omitnan');  % Sum across columns, ignoring NaNs
    
    % Plot faint red-grey line at zero
    if i_g == 6
        % For ax5(6), plot red line first so row_sum appears on top
        % plot([0 0], [1 n_curr], 'Color', [1 0 0], 'LineWidth', 3);
        % hold on
        plot(row_sums, 1:n_curr, 'Color',[0 0 0]', 'LineWidth', 1.5);  
        % hold off;
    else
        plot(row_sums, 1:n_curr, 'Color',[0 0 0], 'LineWidth', 1.5);
        hold on
        plot([0 0], [1 n_curr], 'Color', [1 0 0], 'LineWidth', 3);
        hold off;
    end
    
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
xlim(ax5(1), 5*xlim(ax5(1)));

% % % Figures 11-15 new set of figures, similar to figures 1-5 above in content.  Except this figure will have 1 x 4 subplots.  
% subplot 1 will be H(2), subplot 2 will be H(3), subplot 3 will be H(4), and subplot 4 will be H(5).  Note that H(1) will NOT be in the figures.
% first, we need a RMT with mu_E and mu_I similar to G(6) above, call it H(1)
% In this new set of figure, we want to simulate the sparse activation of a network
% Lets assume the neurons have an internal bias vector bias_b = randn(n,1) which is due to their own intrinsic excitability and average input from other neurons
% Neurons can't have negative firing rates, so if the bias is negative, we must zero out the corresponding columns of the A matrix to form H(2)
% for H(3) we compute the actual spectral radius of H(2) and shift H(3) such that it is on the edge of stability using RMT.shift_diagonal()
% this is to simulate homeostatic mechanisms which put it on the EOC
% now we add a random, sparse stimulus u to bias_b, which simulates a change in environmental stimulus. This instantaneously changes the activation of the 
% neurons.  Copy H(1) to make H(4), then zero out its columns for negative u+bias_b values.  Then shift it the same amount as we shifted H(2) to get H(3)
% this is because it takes a while for homeostatic mechanisms to change, so the shift can't change in time.
% due to the imbalance caused by u, we should now have outlier eigenvalues in H(4)
% for H(5) copy H(4), reduce the magnitude of all non-diagonal elements enough such that the new spectral radius is on the edge of stability.
% figure 11 are the distribution of eigenvalues
% figure 12 are the A matrices imagesc plots
% figure 13 is the colorbar
% figure 14 are the A matrix histograms
% figure 15 are the row sums.

%% Struct to hold RMT objects for figures 11-15
H = struct();

%% H(1) - Base Rajan with row zero sum (similar to G(6))
h = 1;
H(h).rmt = G(1).rmt.copy();
H(h).rmt.set_rajan_means(f, mu_E);
H(h).rmt.row_sum_to_zero();
H(h).rmt.description = {['n=' num2str(n) ', dense'], '\mu_E = 1, \mu_I = -1', 'rows zeroed', 'baseline'};

%% Create bias vector
bias_b = randn(n, 1);

%% H(2) - Zero columns where bias is negative
h = h+1;
H(h).rmt = H(1).rmt.copy();
negative_bias = bias_b < 0;
H(h).rmt.zero_columns(negative_bias);
H(h).rmt.description = {['n=' num2str(n) ', dense'], '\mu_E = 1, \mu_I = -1', 'rows zeroed', 'cols zeroed (neg bias)'};
% Get actual spectral radius of H(2) to use for all circles
H(2).rmt.compute_eigenvalues();
r_H2 = max(real(H(2).rmt.eigenvalues));
H(h).rmt.set_plot_circle(r_H2, 0);

%% H(3) - zero rows and Shift H(2) to edge of stability
h = h+1;
H(h).rmt = H(2).rmt.copy();
H(h).rmt.row_sum_to_zero();
H(h).rmt.compute_eigenvalues();
r_H3 = max(real(H(3).rmt.eigenvalues));
H(h).rmt.shift_diagonal(-r_H3);
H(h).rmt.description = {['n=' num2str(n) ', dense'], '\mu_E = 1, \mu_I = -1', 'rows zeroed', 'cols zeroed, shifted to EOC'};
H(h).rmt.set_plot_circle(r_H3, -r_H3);

%% Create sparse stimulus u
u_density = 0.2;  % 10% of neurons receive stimulus
u = zeros(n, 1);
stimulated_neurons = rand(n, 1) < u_density;
u(stimulated_neurons) = randn(sum(stimulated_neurons), 1);

%% H(4) - Copy H(1), zero columns based on u+bias_b, apply same shift
h = h+1;
H(h).rmt = H(1).rmt.copy();
negative_activation = (u + bias_b) < 0;
H(h).rmt.zero_columns(negative_activation);
H(h).rmt.shift_diagonal(-r_H3);  % Same shift as H(3)
H(h).rmt.description = {['n=' num2str(n) ', dense'], '\mu_E = 1, \mu_I = -1', 'rows zeroed', 'cols zeroed (stimulus), shifted'};
H(h).rmt.set_plot_circle(r_H3, -r_H3);

%% H(5) - Scale non-diagonal elements to edge of stability
h = h+1;
H(h).rmt = H(4).rmt.copy();
H(h).rmt.scale_nondiagonal_to_spectral_radius(r_H2);
H(h).rmt.description = {['n=' num2str(n) ', dense'], '\mu_E = 1, \mu_I = -1', 'rows zeroed', 'cols zeroed, scaled to EOC'};
H(h).rmt.set_plot_circle(r_H3, -r_H3);



%% Figure 11: Eigenvalue distributions
f11 = figure(11);
set(f11, 'Position', [100 200 1420 530], 'Color', 'white')
tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

ax11 = gobjects(4, 1);
for i_h = 2:length(H)
    H(i_h).rmt.compute_eigenvalues();
    ax11(i_h-1) = nexttile;
    H(i_h).rmt.plot_circle(ax11(i_h-1));
    hold on
    H(i_h).rmt.plot_eigenvalue_distribution(ax11(i_h-1));
    hold off
end

linkaxes(ax11,'xy')

% ensure left and right xlims are the same and +/- max(abs(xlim))
xl = xlim(ax11(1));
new_lim = max(abs(xl));
xlim(ax11, [-new_lim, new_lim]);

% Customize axes after linking
for i_h = 2:length(H)
    axes(ax11(i_h-1));
    x_lim = xlim;
    y_lim = ylim;
    axis off;
    
    hold on;
    h_x = plot(x_lim, [0,0], 'k', 'LineWidth', 1.5);
    h_y = plot([0,0], y_lim, 'k', 'LineWidth', 1.5);
    uistack([h_x, h_y], 'bottom');

    text(x_lim(2), 0, ' Re($\lambda$)', 'Interpreter', 'latex', 'VerticalAlignment', 'middle');
    text(0, y_lim(2), 'Im($\lambda$)', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

    text(x_lim(1), y_lim(1), H(i_h).rmt.description, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontWeight', 'normal', ...
        'Rotation', 90);

    xlim(x_lim);
    ylim(y_lim);

    hold off
end

%% Figure 12: A matrices
% Collect all A matrix values to determine clims (reuse from figures 1-5)
all_H_values = [];
for i_h = 2:length(H)
    all_H_values = [all_H_values; H(i_h).rmt.A(:)];
end

% Use same clims as before for consistency
max_abs_vals_H = abs(all_H_values);
clim_val_H = prctile(max_abs_vals_H, 95);
clims_H = [-clim_val_H, clim_val_H];

%% Create concatenated matrix with constant pixel size
% Parameters for layout
n_cols_H = 4;
n_rows_H = 1;
row_spacing_H = 250;
col_spacing_H = 250;

% All H matrices have same size n
max_size_H = n;

% Create cell array to hold padded matrices
padded_matrices_H = cell(4, 1);
pad_before_horiz_array_H = zeros(4, 1);
for i_h = 2:length(H)
    A_plot = H(i_h).rmt.A;
    if make_imagesc_nans_transparent
        A_plot(~H(i_h).rmt.dense_mask) = NaN;
    else
        A_plot(~H(i_h).rmt.dense_mask) = 0;
    end
    padded_matrices_H{i_h-1} = A_plot;
    pad_before_horiz_array_H(i_h-1) = 0;
end

% Build the concatenated matrix
concat_matrix_H = [];
for i_col = 1:n_cols_H
    if i_col == 1
        concat_matrix_H = padded_matrices_H{i_col};
    else
        concat_matrix_H = [concat_matrix_H, NaN(max_size_H, col_spacing_H), padded_matrices_H{i_col}];
    end
end

% Plot the concatenated matrix
f12 = figure(12);
set(f12, 'Position', [100 200 1800 450], 'Color', 'white')
ax12 = axes('Parent', f12);

h12 = imagesc(ax12, concat_matrix_H);
set(h12, 'AlphaData', ~isnan(concat_matrix_H));
set(ax12, 'Color', [1 1 1]);

colormap(ax12, custom_cmap);
caxis(ax12, clims_H);
axis(ax12, 'equal', 'tight');
box off

% Add ylabel-style labels
ylabel_offset_H = -15;
for i_h = 2:length(H)
    i_col = i_h - 1;
    x_pos = (i_col - 1) * (max_size_H + col_spacing_H) + ylabel_offset_H;
    y_pos = 1;
    
    desc = H(i_h).rmt.description;
    text(ax12, x_pos, y_pos, desc, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom', ...
        'FontWeight', 'normal', ...
        'Rotation', 90);
end

set(ax12, 'Visible', 'off');

%% Figure 13: colorbar (reuse same colorbar setup as Figure 3)
% Uses same colormap and clims as figures 1-5

%% Figure 14: Histograms of A matrices
f14 = figure(14);
set(f14, 'Position', [100 200 1420 530], 'Color', 'white')
tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

ax14 = gobjects(4, 1);
for i_h = 2:length(H)
    ax14(i_h-1) = nexttile;
    A_plot = H(i_h).rmt.A;
    A_plot(~H(i_h).rmt.dense_mask) = NaN;
    
    % Separate E and I column values
    A_E_cols = A_plot(:, H(i_h).rmt.E);
    A_I_cols = A_plot(:, H(i_h).rmt.I);
    
    A_E_vals = A_E_cols(:);
    A_E_vals = A_E_vals(~isnan(A_E_vals));
    
    A_I_vals = A_I_cols(:);
    A_I_vals = A_I_vals(~isnan(A_I_vals));
    
    % Plot overlapping histograms with transparency
    hold on;
    histogram(A_E_vals, 50, 'FaceColor', [0 0 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    histogram(A_I_vals, 50, 'FaceColor', [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold off;
    
    grid off;
    axis off;
end

% Link x and y axes across all subplots
linkaxes(ax14, 'xy');

%% Figure 15: Row sums as vertical plots
f15 = figure(15);
set(f15, 'Renderer', 'painters');
set(f15, 'Position', [100 200 1420 530], 'Color', 'white')
tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

ax15 = gobjects(4, 1);
for i_h = 2:length(H)
    ax15(i_h-1) = nexttile;
    A_curr = H(i_h).rmt.A;
    A_curr(~H(i_h).rmt.dense_mask) = NaN;
    n_curr = size(A_curr, 1);
    row_sums = sum(A_curr, 2, 'omitnan');
    
    plot(row_sums, 1:n_curr, 'Color',[0 0 0], 'LineWidth', 1.5);
    hold on
    plot([0 0], [1 n_curr], 'Color', [1 0 0], 'LineWidth', 3);
    hold off;
    
    xlabel('Row Sum');
    ylabel('Row Index');
    grid off;
    set(gca, 'YDir', 'reverse');
    axis off;
end

% Link x and y axes across all subplots
linkaxes(ax15, 'xy');
% Set the xlims to reasonable values
xlim(ax15(1), 5*xlim(ax15(1)));

